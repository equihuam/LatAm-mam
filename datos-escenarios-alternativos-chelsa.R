# ---------- Paquetes ----------
# install.packages(c("pastclim","terra","sf","rnaturalearth","rnaturalearthdata","tmap"))
pacman::p_load(pastclim, terra, sf, rnaturalearth,
               rnaturalearthdata, tmap)

# ---------- Parámetros ----------
ssp  <- "ssp370"  # elige: "ssp126" | "ssp370" | "ssp585"
gcms <- c("GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-HR","MRI-ESM2-0","UKESM1-0-LL")

# Carpeta local para datasets (cámbiala si quieres)
pastclim::set_data_path("chelsa_store")
options(timeout = max(3600, getOption("timeout", 60)))

# ---------- Máscara de América Latina ----------
# Filtramos por región del Banco Mundial "Latin America & Caribbean"
latam_sf <- ne_countries(scale = "medium", returnclass = "sf")
latam_sf <- subset(latam_sf, region_wb == "Latin America & Caribbean")
latam_union <- st_union(latam_sf)
latam_v <- vect(latam_union)  # a SpatVector (terra)

# ---------- Variables mensuales necesarias ----------
present_ds <- "CHELSA_2.1_0.5m_vsi"
# Lista de futuros para cada GCM con VSI (descarga virtual via /vsicurl)
future_ds  <- sprintf("CHELSA_2.1_%s_%s_0.5m_vsi", gcms, ssp)

# Nombres de variables mensuales (evitamos min/max: solo promedio mensual)
temp_months <- sprintf("temperature_%02d", 1:12)
prec_months <- sprintf("precipitation_%02d", 1:12)
vars_needed <- c(temp_months, prec_months)

# ---------- Descargas (configura VRTs/archivos virtuales) ----------
# Presente
download_dataset(dataset = present_ds, bio_variables = vars_needed)
# Futuro (se crea VRT por cada combinación GCM-SSP)
for (ds in future_ds) download_dataset(dataset = ds, bio_variables = vars_needed)

# ---------- Slice (raster) de Presente (1990 = 1981–2000) ----------
present_slice <- region_slice(
  time_ce      = 1990,
  dataset      = present_ds,
  bio_variables = vars_needed,
  region       = latam_union   # recorta el área para traer solo lo necesario
)

# Separamos temperatura y precipitación (12 capas cada uno)
present_T <- present_slice[[grep("^temperature_\\d{2}$", names(present_slice))]]
present_P <- present_slice[[grep("^precipitation_\\d{2}$", names(present_slice))]]

# Anuales (°C y mm/año)
T_present_ann <- mean(present_T)     # promedio de las 12 capas
P_present_ann <- sum(present_P)      # suma de las 12 capas

# Aseguramos recorte/máscara por si quedó algo fuera
T_present_ann <- mask(crop(T_present_ann, latam_v), latam_v)
P_present_ann <- mask(crop(P_present_ann, latam_v), latam_v)

# ---------- Slices de Futuro (2025 = 2011–2040) y deltas por GCM ----------
dT_list <- list()
dP_list <- list()

for (i in seq_along(future_ds)) {
  ds <- future_ds[i]
  
  fut_slice <- region_slice(
    time_ce       = 2025,
    dataset       = ds,
    bio_variables = vars_needed,
    region        = latam_union
  )
  
  fut_T <- fut_slice[[grep("^temperature_\\d{2}$", names(fut_slice))]]
  fut_P <- fut_slice[[grep("^precipitation_\\d{2}$", names(fut_slice))]]
  
  T_fut_ann <- mean(fut_T)
  P_fut_ann <- sum(fut_P)
  
  # Alineamos a la malla del presente para evitar desplazamientos
  T_fut_ann <- project(T_fut_ann, T_present_ann, method = "bilinear")
  P_fut_ann <- project(P_fut_ann, P_present_ann, method = "bilinear")
  
  # Deltas
  dT  <- T_fut_ann - T_present_ann                                   # °C
  dP  <- ((P_fut_ann / P_present_ann) - 1) * 100                      # %
  
  # Guardamos
  names(dT) <- paste0("dT_", gcms[i])
  names(dP) <- paste0("dP_", gcms[i])
  dT_list[[gcms[i]]] <- dT
  dP_list[[gcms[i]]] <- dP
}

# ---------- Ensamble multi-modelo (promedio) ----------
dT_stack <- rast(dT_list)  # capas: una por GCM
dP_stack <- rast(dP_list)

dT_ens <- mean(dT_stack, na.rm = TRUE)  # °C
dP_ens <- mean(dP_stack, na.rm = TRUE)  # %

names(dT_ens) <- "Delta_T_C_2011_2040_vs_1981_2000"
names(dP_ens) <- "Delta_Perc_P_2011_2040_vs_1981_2000"

# ---------- Salidas GeoTIFF ----------
dir.create("salidas_maps", showWarnings = FALSE)
writeRaster(dT_ens, "salidas_maps/CHELSA_CMIP6_Ens_DeltaT_2025_vs_1990_LATAM.tif", overwrite = TRUE)
writeRaster(dP_ens, "salidas_maps/CHELSA_CMIP6_Ens_DeltaPpct_2025_vs_1990_LATAM.tif", overwrite = TRUE)

# ---------- Mapas rápidos con tmap ----------
tmap_mode("plot")

# Paletas sugeridas: divergente para dT (azul-rojo invertido), y secuencial/centrada para dP%
m1 <- tm_shape(dT_ens) +
  tm_raster(title = "ΔT (°C)\n2011–2040 vs 1981–2000",
            palette = "-RdBu", style = "cont", n = 9) +
  tm_layout(title = paste0("CHELSA v2.1 CMIP6 (", ssp, ") – Ensamble 5 GCMs"),
            legend.outside = TRUE)
tmap_save(m1, "salidas_maps/map_DeltaT_LATAM.png", width = 2000, height = 1400, dpi = 200)

m2 <- tm_shape(dP_ens) +
  tm_raster(title = "ΔP (%)\n2011–2040 vs 1981–2000",
            palette = "PuOr", style = "cont", n = 9) +
  tm_layout(title = paste0("CHELSA v2.1 CMIP6 (", ssp, ") – Ensamble 5 GCMs"),
            legend.outside = TRUE)
tmap_save(m2, "salidas_maps/map_DeltaPpct_LATAM.png", width = 2000, height = 1400, dpi = 200)

# ---------- Extras: guarda deltas por GCM (opcional) ----------
# writeRaster(dT_stack, "salidas_maps/_byGCM_DeltaT_stack.tif", overwrite = TRUE)
# writeRaster(dP_stack, "salidas_maps/_byGCM_DeltaPpct_stack.tif", overwrite = TRUE)
