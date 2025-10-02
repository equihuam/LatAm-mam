import os, io, tarfile, re, zipfile, tempfile, requests, json
import geopandas as gpd
import pandas as pd
from shapely.geometry import box
from xml.etree import ElementTree as ET

# --- OUTPUT PATH ---
OUT_GPKG = "latin_america_biomes.gpkg"

# --- 1) Download IUCN GET Level-3 maps (GeoJSON tarball) + map-details.xml ---
ZENODO_README = "https://zenodo.org/records/10081251/preview/README.md?include_deleted=0"
ZENODO_TARBALL = "https://zenodo.org/records/10081251/files/all-maps-vector-geojson.tar.bz2?download=1"
ZENODO_XML = "https://zenodo.org/records/10081251/files/map-details.xml?download=1"

r = requests.get(ZENODO_TARBALL)
r.raise_for_status()
tar_bytes = io.BytesIO(r.content)

rxml = requests.get(ZENODO_XML)
rxml.raise_for_status()
xml_root = ET.fromstring(rxml.content)

# Parse XML to map EFG code -> EFG name
efg_code_to_name = {}
for m in xml_root.findall(".//Map"):
    code = m.attrib.get("efg_code")
    fname = m.findtext("Functional_group")
    if code and fname:
        # e.g., "F1.1 Permanent upland streams" -> ("F1.1", "Permanent upland streams")
        efg_code_to_name[code] = re.sub(r"^[A-Z]+\d+\.\d+\s+", "", fname)

# --- 2) Define Latin America mask (countries) using Natural Earth ---
NE_ADMIN = "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip"
z = zipfile.ZipFile(io.BytesIO(requests.get(NE_ADMIN).content))
tmpdir = tempfile.mkdtemp()
z.extractall(tmpdir)
admin_path = [os.path.join(tmpdir, f) for f in z.namelist() if f.endswith(".shp")][0]

admin = gpd.read_file(admin_path)
admin = admin.to_crs(3857)

# ISO_A3 list for Latin America & Caribbean (UN M49 region 419) + Mexico
LAC_ISO3 = {
    'MEX','BLZ','GTM','SLV','HND','NIC','CRI','PAN',
    'BHS','CUB','JAM','HTI','DOM','TTO','BRB','DMA','GRD','KNA','LCA','VCT','ATG','BRB','BRB',
    'COL','VEN','GUY','SUR','GUF','BRA','ECU','PER','BOL','PRY','URY','ARG','CHL',
    'ABW','BES','CUW','SXM','AIA','VGB','TCA','PRI','MTQ','GLP','MSR','CYM','BMU',
    'BOL','FLK'
}
# Normalize: keep sovereign Latin America + Caribbean; fall back by SUBREGION labels too
lac = admin[(admin['ISO_A3'].isin(LAC_ISO3)) |
            (admin['SUBREGION'].isin(['South America','Central America','Caribbean'])) |
            (admin['NAME_EN'] == 'Mexico')].copy()

# Dissolve to a single polygon for land clip
lac_land = lac.dissolve().geometry.iloc[0]

# --- 3) EEZ for marine clip (Marine Regions) ---
EEZ_URL = "https://www.marineregions.org/downloads/World_EEZ_v12_20231025/eez_v12.zip"
eez_zip = zipfile.ZipFile(io.BytesIO(requests.get(EEZ_URL).content))
eez_zip.extractall(tmpdir)
eez_shp = [os.path.join(tmpdir, f) for f in eez_zip.namelist() if f.endswith(".shp")][0]
eez = gpd.read_file(eez_shp).to_crs(3857)

# Filter EEZ to Latin America sovereigns (match by TERRITORY1 or ISO_Ter1 where possible)
# We'll use a broad spatial filter: intersect with lac_land buffered offshore
eez_latam = eez[eez.intersects(lac_land.buffer(500000))].copy()
eez_latam = eez_latam.dissolve().geometry.iloc[0]

# --- 4) Load all IUCN GET EFG GeoJSONs from tarball, tag realm/biome, clip ---
gdfs = []
with tarfile.open(fileobj=tar_bytes, mode='r:bz2') as tf:
    members = [m for m in tf.getmembers() if m.name.endswith(".json")]
    for m in members:
        f = tf.extractfile(m); 
        if f is None: 
            continue
        data = json.load(io.TextIOWrapper(f, encoding='utf-8'))
        g = gpd.GeoDataFrame.from_features(data['features'], crs="EPSG:4326").to_crs(3857)
        # derive fields from filename and attributes
        # filenames like: T4.1.web.mix_v1.0.json
        base = os.path.basename(m.name)
        code_match = re.match(r"([A-Z]\d+\.\d+)", base)
        if not code_match:
            continue
        efg_code = code_match.group(1)
        realm = efg_code[0]  # T/F/M/S/...
        biome_code = re.match(r"([A-Z]\d+)", efg_code).group(1)  # e.g., "T4"
        efg_name = efg_code_to_name.get(efg_code, None)

        # occurrence: major=1, minor=2 (per README)
        # some layers store it in "occurrence" or "value" fields; harmonize:
        occ = None
        for cand in ["occurrence","value","class","VAL","val"]:
            if cand in g.columns:
                occ = g[cand]
                break
        if occ is None:
            occ = 1  # default if missing
        g["occurrence"] = g.get("occurrence", occ).replace({1:"major", 2:"minor"})

        g["realm"] = realm
        g["biome_code"] = biome_code
        g["efg_code"] = efg_code
        g["efg_name"] = efg_name if efg_name else efg_code

        # Clip by realm: land realms (T,F,SF,FM,MT?) to land; marine (M, SM) to EEZ
        if realm in ("M","SM"):   # marine realms (M=Marine; SM=Submarine/Seamount/Sea-ice grids)
            g = g[g.geometry.notnull()].copy()
            g["geometry"] = g.geometry.intersection(eez_latam)
        else:
            g = g[g.geometry.notnull()].copy()
            g["geometry"] = g.geometry.intersection(lac_land)

        g = g[g.geometry.is_valid & ~g.geometry.is_empty]
        if len(g):
            gdfs.append(g[["realm","biome_code","efg_code","efg_name","occurrence","geometry"]])

efg_latam = pd.concat(gdfs, ignore_index=True) if gdfs else gpd.GeoDataFrame(geometry=[], crs=3857)

# Map biome_code -> friendly L2 biome names (concise; based on IUCN GET v2.0 profiles)
biome_name_map = {
    # Terrestrial (T*)
    "T1":"Tropical-subtropical forests","T2":"Temperate-boreal forests & woodlands",
    "T3":"Shrublands & shrubby woodlands","T4":"Savannas & dry woodlands",
    "T5":"Temperate grasslands/savannas/shrublands","T6":"Montane grasslands & shrublands",
    "T7":"Polar & alpine","T8":"Flooded grasslands & savannas","T9":"Deserts & xeric shrublands",
    "T10":"Forested wetlands & peatlands","T11":"Non-forest wetlands","T12":"Ice/rock/bare",
    "T13":"Inland saline & ephemeral","T14":"Coastal & tidal systems",
    "T15":"Mosaic/complex land mosaics","T16":"Intensive land-use systems",
    # Freshwater (F*)
    "F1":"Rivers & streams","F2":"Lakes","F3":"Palustrine/wetlands",
    # Marine (M*)
    "M1":"Coastal & shelf pelagic/benthic","M2":"Oceanic pelagic","M3":"Deep-sea benthic",
    "M4":"Sea-ice & polar marine",
    # Subterranean / Submarine / Misc (S*, SF*, SM*, MT*, FM*)
    "S1":"Subterranean karst & caves","S2":"Subterranean groundwater",
    "SF1":"Saline/estuarine wetlands","SF2":"Inland saline wetlands",
    "SM1":"Seamounts/undersea features","MT1":"Tidal/fluvial transitional","MT2":"Delta/estuary",
    "FM1":"Floodplains & riparian"
}
efg_latam["biome_name"] = efg_latam["biome_code"].map(biome_name_map).fillna(efg_latam["biome_code"])

# --- 5) Dissolve to Level-2 biomes for quick maps ---
biomes_latam = efg_latam.dissolve(by=["realm","biome_code","biome_name"])

# --- 6) Write the GeoPackage ---
if os.path.exists(OUT_GPKG):
    os.remove(OUT_GPKG)

efg_latam.to_crs(4326).to_file(OUT_GPKG, layer="get_efg_latam", driver="GPKG")
biomes_latam.to_crs(4326).to_file(OUT_GPKG, layer="get_biomes_latam", driver="GPKG")
lac.to_crs(4326).to_file(OUT_GPKG, layer="countries_latam", driver="GPKG")

# Write EEZ reference (singlepart dissolve already)
gpd.GeoDataFrame(geometry=[eez_latam], crs=3857).to_crs(4326).to_file(OUT_GPKG, layer="eez_latam", driver="GPKG")

print(f"✅ Wrote {OUT_GPKG} with layers: get_efg_latam, get_biomes_latam, countries_latam, eez_latam")
