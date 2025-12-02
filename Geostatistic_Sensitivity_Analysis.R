# ==============================================================
# COMPLETE, CORRECTED PIPELINE: Kriging, Sensitivity, SD maps, Reports
# Processes sheets: df1, df2, df3, df4, df5
# ==============================================================
# 0. PACKAGES --------------------------------------------------
library(sf)
library(sp)
library(gstat)
library(dplyr)
library(raster)
library(ggplot2)
library(viridis)
library(readxl)
library(cowplot)
library(ggspatial)
library(gridExtra)
library(rmarkdown)

# 1. CLEAN AND SET WORKDIR ------------------------------------
rm(list = ls())
#setwd("")  

# Create output folders (silently)
dirs <- c("PNG","SVG","TIF","CSV","REPORT")
for (d in dirs) if (!dir.exists(d)) dir.create(d)

# 2. GLOBAL SETTINGS -------------------------------------------
SCALE_C13 <- c(-28.5, -14.5)   # fixed color range for predicted maps
SD_GLOBAL_RANGE <- NULL        # computed after first pass
GRID_RES <- 10                 # grid resolution in meters
LAG_DIST <- 165                # variogram lag width
EXPAND_RATIO <- 0.00           # expand bbox by 0% (set >0 to add margin)
MODELS_TO_TEST <- c("Sph","Exp","Gau")

# 3. Helper: extract legend (cowplot)
get_legend <- function(mygg) {
  gg <- ggplotGrob(mygg)
  guides <- which(sapply(gg$grobs, function(x) x$name) == "guide-box")
  if (length(guides) == 0) return(NULL)
  gg$grobs[[guides]]
}

# 4. Processing function (single sheet) - robust and self-contained
process_soil_sheet <- function(sheet_name) {
  cat("\n===================================================================\n")
  cat("Processing sheet:", sheet_name, "\n")
  cat("===================================================================\n")
  
  # 4.1 Load data
  df <- read_excel("Soil_DB.xlsx", sheet = sheet_name)
  required_cols <- c("Coord_X","Coord_Y","C13")
  if (!all(required_cols %in% names(df))) stop("Sheet ", sheet_name, " must contain: ", paste(required_cols, collapse = ", "))
  
  # 4.2 Create sf and ensure validity
  sf_pts <- st_as_sf(df, coords = c("Coord_X","Coord_Y"), crs = 31982, remove = FALSE)
  sf_pts <- sf_pts[!st_is_empty(sf_pts), ]
  sf_pts <- st_make_valid(sf_pts)
  
  # Force attribute transfer and create sp object robustly
  sp_pts <- as(sf_pts, "Spatial")
  # Ensure @data matches df rows and columns
  sp_pts@data <- as.data.frame(sf_pts)[, setdiff(names(as.data.frame(sf_pts)), "geometry"), drop = FALSE]
  
  # 4.3 Coordinates and trend
  coords <- coordinates(sp_pts)
  dados_sp <- sp_pts@data
  dados_sp$X <- coords[,1]; dados_sp$Y <- coords[,2]
  
  lm_trend <- lm(C13 ~ X + Y, data = dados_sp)
  dados_sp$resid <- resid(lm_trend)
  sp_pts$resid <- dados_sp$resid  # attach to sp object
  
  # 4.4 Empirical variogram (use sp_pts)
  dist_mat <- sp::spDists(coords)
  cutoff <- 0.5 * max(dist_mat, na.rm = TRUE)
  vg_emp <- variogram(resid ~ 1, sp_pts, cutoff = cutoff, width = LAG_DIST)
  
  # Convert vg_emp to dataframe for ggplot
  vg_df <- as.data.frame(vg_emp)
  
  # If empty, warn explicitly
  if (nrow(vg_df) == 0) warning("Empirical variogram is empty for sheet: ", sheet_name)
  
  # 4.5 Save empirical variogram (ggplot) - points only for now
  p_vg_emp <- ggplot(vg_df, aes(x = dist, y = gamma)) +
    geom_point(size = 2) +
    geom_line() +
    theme_minimal(base_size = 14) +
    labs(title = paste("Empirical semivariogram -", sheet_name),
         x = "Distance (m)", y = "Semivariance")
  ggsave(filename = file.path("PNG", paste0("Variogram_Empirical_", sheet_name, ".png")),
         plot = p_vg_emp, width = 8, height = 6, dpi = 300)
  ggsave(filename = file.path("SVG", paste0("Variogram_Empirical_", sheet_name, ".svg")),
         plot = p_vg_emp, width = 8, height = 6, device = "svg")
  
  # 4.6 Test models and LOOCV metrics, capture vgm params; store raw numbers
  results_raw <- data.frame(Model=character(), RMSE=numeric(), Bias=numeric(), MAE=numeric(), R2=numeric(),
                            Nugget=numeric(), Range=numeric(), Sill=numeric(), stringsAsFactors = FALSE)
  
  for (m in MODELS_TO_TEST) {
    init <- vgm(psill = var(dados_sp$resid, na.rm = TRUE), model = m, range = cutoff/3, nugget = 0)
    fit_try <- try(fit.variogram(vg_emp, model = init), silent = TRUE)
    if (inherits(fit_try, "try-error")) {
      results_raw <- rbind(results_raw, data.frame(Model=m, RMSE=NA, Bias=NA, MAE=NA, R2=NA, Nugget=NA, Range=NA, Sill=NA))
      next
    }
    fit <- fit_try
    cv <- krige.cv(resid ~ 1, sp_pts, model = fit, nfold = nrow(sp_pts))
    rmse <- sqrt(mean(cv$residual^2, na.rm = TRUE))
    bias <- mean(cv$residual, na.rm = TRUE)
    mae <- mean(abs(cv$residual), na.rm = TRUE)
    r2  <- cor(cv$observed, cv$observed - cv$residual, use = "complete.obs")^2
    # variogram params
    nugget <- fit$psill[1]
    sill   <- sum(fit$psill)
    # pick last positive range
    rngs <- fit$range[which(!is.na(fit$range) & fit$range > 0)]
    rng    <- ifelse(length(rngs) > 0, tail(rngs, 1), NA)
    
    results_raw <- rbind(results_raw,
                         data.frame(Model=m, RMSE=rmse, Bias=bias, MAE=mae, R2=r2, Nugget=nugget, Range=rng, Sill=sill))
  }
  
  # 4.7 Round sensitivity numbers to 2 decimal places (as requested)
  results <- results_raw %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  # Save per-sheet sensitivity (rounded)
  write.csv(results, file.path("CSV", paste0("Sensitivity_", sheet_name, ".csv")), row.names = FALSE)
  
  # 4.8 Select best model by RMSE (using raw values to avoid rounding tie issues)
  best_idx <- which.min(results_raw$RMSE)
  if (length(best_idx) == 0 || is.na(results_raw$RMSE[best_idx])) {
    stop("No successful variogram fit for sheet: ", sheet_name)
  }
  best_model_name <- results_raw$Model[best_idx]
  cat("Best model selected (by RMSE):", best_model_name, "\n")
  
  # Fit best model for plotting and kriging
  best_init <- vgm(psill = var(dados_sp$resid, na.rm = TRUE), model = best_model_name, range = cutoff/3, nugget = 0)
  fitted_vgm <- fit.variogram(vg_emp, model = best_init)
  
  # 4.9 Build theoretical variogram line for overlay (fine resolution)
  vg_line <- variogramLine(fitted_vgm, maxdist = cutoff, n = 200)
  vg_line_df <- as.data.frame(vg_line)
  
  # 4.10 Plot empirical + fitted variogram (ggplot) and save
  p_vg_both <- ggplot() +
    geom_point(data = vg_df, aes(x = dist, y = gamma), size = 2) +
    geom_line(data = vg_line_df, aes(x = dist, y = gamma), color = "red", size = 1) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Empirical semivariogram and fitted model -", sheet_name),
         subtitle = paste0("Model: ", best_model_name),
         x = "Distance (m)", y = "Semivariance")
  ggsave(file.path("PNG", paste0("Variogram_Empirical_Fitted_", sheet_name, ".png")), p_vg_both, width = 8, height = 6, dpi = 300)
  ggsave(file.path("SVG", paste0("Variogram_Empirical_Fitted_", sheet_name, ".svg")), p_vg_both, width = 8, height = 6, device = "svg")
  
  # 4.11 Final LOOCV with fitted_vgm (raw metrics)
  cv_final <- krige.cv(resid ~ 1, sp_pts, model = fitted_vgm, nfold = nrow(sp_pts))
  rmse_final <- sqrt(mean(cv_final$residual^2, na.rm = TRUE))
  bias_final <- mean(cv_final$residual, na.rm = TRUE)
  mae_final <- mean(abs(cv_final$residual), na.rm = TRUE)
  r2_final <- cor(cv_final$observed, cv_final$observed - cv_final$residual, use = "complete.obs")^2
  
  # Save validation summary (rounded)
  validation_summary <- data.frame(Sheet = sheet_name,
                                   Model = best_model_name,
                                   RMSE  = round(rmse_final, 2),
                                   Bias  = round(bias_final, 2),
                                   MAE   = round(mae_final, 2),
                                   R2    = round(r2_final, 2),
                                   stringsAsFactors = FALSE)
  write.csv(validation_summary, file.path("CSV", paste0("ValidationSummary_", sheet_name, ".csv")), row.names = FALSE)
  
  # Save observed vs predicted (raw)
  obs_pred_df <- data.frame(Observed = cv_final$observed,
                            Predicted = cv_final$observed - cv_final$residual)
  write.csv(obs_pred_df, file.path("CSV", paste0("ObsPred_", sheet_name, ".csv")), row.names = FALSE)
  
  # Save fitted variogram parameters (rounded)
  vgm_params <- data.frame(Sheet = sheet_name,
                           Model = best_model_name,
                           Nugget = round(fitted_vgm$psill[1], 2),
                           PartialSill = round(fitted_vgm$psill[2], 2),
                           Range = round(fitted_vgm$range[2], 2),
                           Sill = round(sum(fitted_vgm$psill), 2),
                           stringsAsFactors = FALSE)
  write.csv(vgm_params, file.path("CSV", paste0("VGM_Params_", sheet_name, ".csv")), row.names = FALSE)
  
  # 4.12 Kriging interpolation (grid based on shapefile)
  area <- st_read("GLEBA_2.shp", quiet = TRUE)
  area_utm <- sf::st_transform(area, crs = 31982)
  bbox <- st_bbox(area_utm)
  # optional expand
  xspan <- as.numeric(bbox["xmax"] - bbox["xmin"]); yspan <- as.numeric(bbox["ymax"] - bbox["ymin"])
  bbox["xmin"] <- bbox["xmin"] - xspan * EXPAND_RATIO
  bbox["xmax"] <- bbox["xmax"] + xspan * EXPAND_RATIO
  bbox["ymin"] <- bbox["ymin"] - yspan * EXPAND_RATIO
  bbox["ymax"] <- bbox["ymax"] + yspan * EXPAND_RATIO
  
  grd <- expand.grid(x = seq(bbox["xmin"], bbox["xmax"], by = GRID_RES),
                     y = seq(bbox["ymin"], bbox["ymax"], by = GRID_RES))
  grd_sf <- st_as_sf(grd, coords = c("x","y"), crs = st_crs(sf_pts))
  grd_sp <- as(grd_sf, "Spatial")
  
  krig_resid <- krige(resid ~ 1, sp_pts, grd_sp, model = fitted_vgm)
  pred_trend <- predict(lm_trend, newdata = data.frame(X = coordinates(grd_sp)[,1], Y = coordinates(grd_sp)[,2]))
  krig_resid$C13_pred <- krig_resid$var1.pred + pred_trend
  pred_df <- as.data.frame(krig_resid)
  pred_df$SD <- sqrt(pred_df$var1.var)
  
  # Rasters
  r_pred <- rasterFromXYZ(pred_df[, c("coords.x1","coords.x2","C13_pred")])
  crs(r_pred) <- CRS(SRS_string = "EPSG:31982")
  r_sd   <- rasterFromXYZ(pred_df[, c("coords.x1","coords.x2","SD")])
  crs(r_sd)   <- CRS(SRS_string = "EPSG:31982")
  
  area_sp <- as(area_utm, "Spatial")
  r_pred_mask <- raster::mask(raster::crop(r_pred, area_sp), area_sp)
  r_sd_mask   <- raster::mask(raster::crop(r_sd, area_sp), area_sp)
  
  # Save rasters
  writeRaster(r_pred_mask, file.path("TIF", paste0("C13_Kriging_", sheet_name, ".tif")), overwrite = TRUE)
  writeRaster(r_sd_mask,   file.path("TIF", paste0("C13_SD_", sheet_name, ".tif")), overwrite = TRUE)
  
  # Return objects and max SD for global scaling
  list(sheet = sheet_name,
       validation = validation_summary,
       sensitivity = results,
       vgm_params = vgm_params,
       obs_pred = obs_pred_df,
       maxSD = max(pred_df$SD, na.rm = TRUE))
}

# 5. FIRST PASS: process all sheets and collect SD maxima
sheets <- c("df1","df2","df3","df4","df5")
first_pass <- list(); sd_max <- numeric(length(sheets))
for (i in seq_along(sheets)) {
  s <- sheets[i]
  res <- process_soil_sheet(s)
  first_pass[[s]] <- res
  sd_max[i] <- res$maxSD
}

# Set global SD range (0 .. max)
SD_GLOBAL_RANGE <- c(0, max(sd_max, na.rm = TRUE))
cat("\nGlobal SD range (fixed):", SD_GLOBAL_RANGE, "\n")

# 6. POST-PROCESS: generate maps, obs×pred plots, and combine tables
map_list <- list(); sd_list <- list(); obs_list <- list()
all_sensitivity <- data.frame(); all_validation <- data.frame()

for (s in sheets) {
  cat("\nPost-processing:", s, "\n")
  # read rasters and CSVs
  r_pred <- raster(file.path("TIF", paste0("C13_Kriging_", s, ".tif")))
  r_sd   <- raster(file.path("TIF", paste0("C13_SD_", s, ".tif")))
  sens   <- read.csv(file.path("CSV", paste0("Sensitivity_", s, ".csv")), stringsAsFactors = FALSE)
  val    <- read.csv(file.path("CSV", paste0("ValidationSummary_", s, ".csv")), stringsAsFactors = FALSE)
  obs_df <- read.csv(file.path("CSV", paste0("ObsPred_", s, ".csv")), stringsAsFactors = FALSE)
  vgm_p  <- read.csv(file.path("CSV", paste0("VGM_Params_", s, ".csv")), stringsAsFactors = FALSE)
  
  sens$Sheet <- s; val$Sheet <- s
  all_sensitivity <- rbind(all_sensitivity, sens)
  all_validation  <- rbind(all_validation, val)
  
  # Prepare dataframes for ggplot
  r_pred_df <- as.data.frame(rasterToPoints(r_pred)); colnames(r_pred_df) <- c("x","y","C13")
  r_sd_df   <- as.data.frame(rasterToPoints(r_sd));   colnames(r_sd_df)   <- c("x","y","SD")
  area <- st_read("GLEBA_2.shp", quiet = TRUE); area_utm <- st_transform(area, crs = 31982)
  
  # C13 map with fixed scale and map elements
  p_map <- ggplot() +
    geom_raster(data = r_pred_df, aes(x=x, y=y, fill=C13)) +
    scale_fill_viridis(name = "δ13C", limits = SCALE_C13, option = "D", na.value = "transparent") +
    geom_sf(data = st_geometry(area_utm), fill = NA, color = "black", linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = 0.4) +
    annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
    coord_sf() + theme_minimal(base_size = 14) +
    labs(title = paste("δ13C (Kriging -", s, ")"), x = "Easting (m)", y = "Northing (m)") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # SD map with global SD range fixed
  p_sd <- ggplot() +
    geom_raster(data = r_sd_df, aes(x=x, y=y, fill=SD)) +
    scale_fill_viridis(name = "SD (δ13C)", limits = SD_GLOBAL_RANGE, option = "C", na.value = "transparent") +
    geom_sf(data = st_geometry(area_utm), fill = NA, color = "black", linewidth = 0.5) +
    annotation_scale(location = "bl", width_hint = 0.4) +
    annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
    coord_sf() + theme_minimal(base_size = 14) +
    labs(title = paste("Prediction SD -", s), x = "Easting (m)", y = "Northing (m)") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # Observed vs Predicted with 95% CI and metrics in subtitle
  rmse_v <- val$RMSE; bias_v <- val$Bias; mae_v <- val$MAE; r2_v <- val$R2
  p_obs <- ggplot(obs_df, aes(x = Observed, y = Predicted)) +
    geom_point(aes(color = abs(Observed - Predicted)), size = 3, alpha = 0.85) +
    scale_color_viridis_c(name = "|Residual|") +
    geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "black", fill = "gray70") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_equal() +
    theme_minimal(base_size = 14) +
    labs(title = paste("Observed vs Predicted δ13C -", s),
         subtitle = paste0("R²=", round(r2_v, 2), " | RMSE=", round(rmse_v, 2), " | MAE=", round(mae_v, 2), " | Bias=", round(bias_v,2)),
         x = "Observed δ13C", y = "Predicted δ13C")
  
  # Save individual plots
  ggsave(file.path("PNG", paste0("C13_Kriging_Map_", s, ".png")), p_map, width = 6, height = 5, dpi = 300)
  ggsave(file.path("PNG", paste0("C13_SD_Map_", s, ".png")), p_sd, width = 6, height = 5, dpi = 300)
  ggsave(file.path("PNG", paste0("C13_ObsPred_", s, ".png")), p_obs, width = 6, height = 5, dpi = 300)
  
  # Save SVGs too
  ggsave(file.path("SVG", paste0("C13_Kriging_Map_", s, ".svg")), p_map, width = 6, height = 5, device = "svg")
  ggsave(file.path("SVG", paste0("C13_SD_Map_", s, ".svg")), p_sd, width = 6, height = 5, device = "svg")
  ggsave(file.path("SVG", paste0("C13_ObsPred_", s, ".svg")), p_obs, width = 6, height = 5, device = "svg")
  
  # Store objects for combination
  map_list[[s]] <- p_map
  sd_list[[s]] <- p_sd
  obs_list[[s]] <- p_obs
}

# 7. Merge sensitivity & validation into master CSVs (round numeric to 2 decimals)
all_sensitivity <- all_sensitivity %>% mutate()  # placeholder
# Read back per-sheet sensitivity files (ensures consistent ordering and rounding)
sens_merge <- do.call(rbind, lapply(sheets, function(s) {
  tmp <- read.csv(file.path("CSV", paste0("Sensitivity_", s, ".csv")), stringsAsFactors = FALSE)
  tmp$Sheet <- s
  # enforce numeric rounding (already rounded on save but ensure)
  tmp <- tmp %>% mutate(across(where(is.numeric), ~round(., 2)))
  tmp
}))
write.csv(sens_merge, file.path("CSV", "Sensitivity_Results_All_Depths.csv"), row.names = FALSE)

val_merge <- do.call(rbind, lapply(sheets, function(s) {
  tmp <- read.csv(file.path("CSV", paste0("ValidationSummary_", s, ".csv")), stringsAsFactors = FALSE)
  tmp$Sheet <- s
  tmp <- tmp %>% mutate(across(where(is.numeric), ~round(., 2)))
  tmp
}))
write.csv(val_merge, file.path("CSV", "Validation_Results_All_Depths.csv"), row.names = FALSE)

cat("\nSaved master sensitivity and validation CSVs.\n")

# 8. Create combined panels with shared legends
# Combined predicted maps
shared_leg_map <- get_legend(map_list[[1]] + theme(legend.position = "right"))
map_grid <- plot_grid(plotlist = lapply(map_list, function(x) x + theme(legend.position = "none")), ncol = 2)
final_map_panel <- plot_grid(map_grid, shared_leg_map, rel_widths = c(1, 0.12))
ggsave(filename = file.path("PNG", "C13_Kriging_All_Combined.png"), final_map_panel, width = 12, height = 10, dpi = 300)

# Combined SD maps
shared_leg_sd <- get_legend(sd_list[[1]] + theme(legend.position = "right"))
sd_grid <- plot_grid(plotlist = lapply(sd_list, function(x) x + theme(legend.position = "none")), ncol = 2)
final_sd_panel <- plot_grid(sd_grid, shared_leg_sd, rel_widths = c(1, 0.12))
ggsave(filename = file.path("PNG", "C13_SD_All_Combined.png"), final_sd_panel, width = 12, height = 10, dpi = 300)

# Combined Obs×Pred
obs_grid <- plot_grid(plotlist = obs_list, ncol = 2)
ggsave(filename = file.path("PNG", "C13_ObsPred_All_Combined.png"), obs_grid, width = 12, height = 10, dpi = 300)

cat("\nCombined PNGs created.\n")

# 9. Optional PDF report via rmarkdown (produces REPORT/C13_Report.pdf if possible)
rmd_file <- file.path("REPORT", "C13_Report.Rmd")
pdf_file <- file.path("REPORT", "C13_Report.pdf")
rmd_text <- '
