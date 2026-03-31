# =============================================================================
# SOIL ANALYSIS — FINAL SCRIPT
# Covers: Cstock (total + by layer), Texture (Sand/Silt/Clay),
#         Bulk Density (BD), Carbon concentration (C)
# Structure:
#   SECTION A — Between-treatment comparisons per depth layer (df1–df5)
#   SECTION B — Within-treatment comparisons across depth (df6, per vegetation)
#   SECTION C — Supplementary Table S1 (assumption audit)
#   SECTION D — Plots (Cstock boxplots + δ¹³C profile)
# =============================================================================

library(tidyverse)
library(rstatix)
library(dplyr)
library(ggplot2)
library(FSA)
library(multcompView)
library(agricolae)
library(readxl)
library(gridExtra)

set.seed(1234)
rm(list = ls())
setwd("")


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Shapiro-Wilk — guards against non-finite values (e.g. log(0) = -Inf)
test_normality <- function(x, label = "") {
  x_clean <- x[is.finite(x) & !is.na(x)]
  if (length(x_clean) < 3) {
    cat(sprintf("[%s] SKIPPED — fewer than 3 finite values.\n", label))
    return(invisible(list(statistic = NA, p.value = 1)))
  }
  res    <- shapiro.test(x_clean)
  status <- ifelse(res$p.value < 0.05, "does NOT FOLLOW", "FOLLOWS")
  cat(sprintf("[%s] %s a normal distribution (W = %.4f, p = %.6f)\n",
              label, status, res$statistic, res$p.value))
  invisible(res)
}

#' Levene's test
test_levene <- function(data, formula, label = "") {
  res    <- levene_test(data, formula)
  status <- ifelse(res$p < 0.05, "NOT homogeneous", "homogeneous")
  cat(sprintf("[%s] Levene: variances are %s (p = %.4f)\n",
              label, status, res$p))
  invisible(res)
}

#' ANOVA + Tukey HSD; back-transforms means when log_transformed = TRUE
run_anova_tukey <- function(data, response, group,
                            alpha = 0.05, log_transformed = FALSE) {
  formula <- as.formula(paste(response, "~", group))
  model   <- aov(formula, data = data)
  hsd     <- HSD.test(model, group, alpha = alpha, group = TRUE, console = FALSE)

  if (log_transformed) {
    cat("\n--- ANOVA + Tukey HSD (log scale) ---\n")
    cat("Significance letters:\n")
    print(hsd$groups[, "groups", drop = FALSE])
    gm <- data %>%
      group_by(.data[[group]]) %>%
      summarise(n = n(),
                geometric_mean = round(exp(mean(.data[[response]], na.rm = TRUE)), 4),
                .groups = "drop") %>%
      left_join(tibble(!!group := rownames(hsd$groups),
                       letter  = hsd$groups$groups), by = group)
    cat("Geometric means (back-transformed):\n"); print(gm)
  } else {
    cat("\n--- ANOVA + Tukey HSD ---\n"); print(hsd)
  }
  invisible(list(model = model, tukey = TukeyHSD(model), hsd = hsd))
}

#' Welch ANOVA + Games-Howell (normal but heterogeneous variances)
run_welch_gh <- function(data, response, group, log_transformed = FALSE) {
  formula <- as.formula(paste(response, "~", group))
  welch   <- oneway.test(formula, data = data, var.equal = FALSE)
  cat(sprintf("\n--- Welch ANOVA: F(%.2f, %.2f) = %.4f, p = %.6f ---\n",
              welch$parameter[1], welch$parameter[2],
              welch$statistic, welch$p.value))
  gh <- games_howell_test(data, formula)
  cat("Games-Howell:\n"); print(gh)
  if (log_transformed) {
    gm <- data %>%
      group_by(.data[[group]]) %>%
      summarise(n = n(),
                geometric_mean = round(exp(mean(.data[[response]], na.rm = TRUE)), 4),
                .groups = "drop")
    cat("Geometric means (back-transformed):\n"); print(gm)
  }
  invisible(list(welch = welch, games_howell = gh))
}

#' Kruskal-Wallis + Dunn (Bonferroni); reports median + IQR
run_kruskal_dunn <- function(data, response, group,
                             p.adjust = "bonferroni", label = "") {
  desc <- data %>%
    group_by(.data[[group]]) %>%
    summarise(n      = sum(!is.na(.data[[response]])),
              Median = round(median(.data[[response]], na.rm = TRUE), 4),
              Q25    = round(quantile(.data[[response]], .25, na.rm = TRUE), 4),
              Q75    = round(quantile(.data[[response]], .75, na.rm = TRUE), 4),
              IQR    = round(IQR(.data[[response]], na.rm = TRUE), 4),
              .groups = "drop")
  cat("\nDescriptive (median + IQR):\n"); print(desc)

  formula <- as.formula(paste(response, "~", group))
  kr      <- kruskal_test(data, formula)
  cat(sprintf("[%s] Kruskal-Wallis: chi2(%d) = %.4f, p = %.6f\n",
              label, kr$df, kr$statistic, kr$p))
  if (kr$p < 0.05) {
    ph <- dunn_test(data, formula, p.adjust.method = p.adjust)
    cat(sprintf("Dunn post-hoc (p.adjust = '%s'):\n", p.adjust)); print(ph)
  } else {
    cat("Not significant — Dunn post-hoc not necessary.\n")
  }
  invisible(kr)
}

# -----------------------------------------------------------------------------
#' Core pipeline:
#'   1. Coerce to numeric
#'   2. Check can_log (skip if zeros/negatives)
#'   3. Shapiro-Wilk original
#'   4. Log-transform + re-test if needed and possible
#'   5. Levene on chosen column
#'   6. Route → ANOVA / Welch / Kruskal
analyze_variable <- function(data, response, group = "Vegetation",
                             label = "", p.adjust = "bonferroni") {
  cat(sprintf("\n%s\n  Variable: %s | Group: %s | Label: %s\n%s\n",
              strrep("=", 60), response, group, label, strrep("=", 60)))

  data[[response]] <- suppressWarnings(as.numeric(data[[response]]))
  if (all(is.na(data[[response]]))) {
    cat(sprintf("  SKIPPED: '%s' is entirely NA or non-numeric.\n", response))
    return(invisible(NULL))
  }

  log_col         <- paste0(response, "_log")
  log_transformed <- FALSE
  can_log         <- all(data[[response]] > 0, na.rm = TRUE)

  if (!can_log)
    cat(sprintf("  Log-transform skipped: '%s' contains zero or negative values.\n",
                response))

  sw     <- test_normality(data[[response]], label)
  normal <- !is.na(sw$p.value) && sw$p.value >= 0.05

  if (!normal && can_log) {
    data[[log_col]] <- log(data[[response]])
    sw_log  <- test_normality(data[[log_col]], paste(label, "[log]"))
    normal  <- !is.na(sw_log$p.value) && sw_log$p.value >= 0.05
    if (normal) { log_transformed <- TRUE; use_col <- log_col
    } else        use_col <- response
  } else {
    use_col <- response
  }

  lev         <- test_levene(data, as.formula(paste(use_col, "~", group)), label)
  homogeneous <- lev$p >= 0.05

  if (!normal)
    run_kruskal_dunn(data, response, group, p.adjust = p.adjust, label = label)
  else if (homogeneous)
    run_anova_tukey(data, use_col, group, log_transformed = log_transformed)
  else {
    cat(sprintf("  ⚠ Heterogeneous variances (Levene p = %.4f) → Welch + Games-Howell\n",
                lev$p))
    run_welch_gh(data, use_col, group, log_transformed = log_transformed)
  }
}

#' Run analyze_variable for multiple variables in one layer
analyze_layer <- function(data, variables, group = "Vegetation",
                          layer_name = "", p.adjust = "bonferroni") {
  data[[group]] <- as.factor(data[[group]])
  for (var in variables)
    if (var %in% names(data))
      analyze_variable(data, var, group = group,
                       label    = paste(layer_name, var),
                       p.adjust = p.adjust)
}

#' Within-treatment depth comparison for one vegetation type
analyze_depth_within_veg <- function(df6, veg_name, response,
                                     depth_levels = c("P20","P40","P60","P80","P100"),
                                     depth_labels = c("20cm","40cm","60cm","80cm","100cm")) {
  cat(sprintf("\n%s\n  Variable: %s | Vegetation: %s\n%s\n",
              strrep("=", 60), response, veg_name, strrep("=", 60)))

  if (!response %in% names(df6)) {
    cat(sprintf("  SKIPPED: '%s' not found in sheet 7.\n", response))
    return(invisible(NULL))
  }

  df_veg <- df6 %>%
    filter(Vegetation == veg_name) %>%
    mutate(Depth = factor(Depth, levels = depth_levels, labels = depth_labels))
  df_veg[[response]] <- suppressWarnings(as.numeric(df_veg[[response]]))

  analyze_variable(df_veg, response, group = "Depth",
                   label    = paste(veg_name, response),
                   p.adjust = "bonferroni")
}

# -----------------------------------------------------------------------------
#' Assumption audit — one row per variable × layer (for Supp Table S1)
audit_variable <- function(data, response, group = "Vegetation",
                           layer = "", variable = "") {
  data[[response]] <- suppressWarnings(as.numeric(data[[response]]))
  x       <- data[[response]][!is.na(data[[response]])]
  if (length(x) < 3)
    return(data.frame(Layer = layer, Variable = variable,
                      W_original = NA, p_original = NA, Normal_orig = NA,
                      W_log = NA, p_log = NA, Normal_log = NA,
                      Transform = "Skipped — insufficient data",
                      Levene_p = NA, Homogeneous = NA,
                      Test_Used = "Skipped", stringsAsFactors = FALSE))

  can_log     <- all(x > 0)
  x_finite    <- x[is.finite(x)]
  sw_orig     <- if (length(x_finite) >= 3) shapiro.test(x_finite) else list(statistic = NA, p.value = NA)
  normal_orig <- !is.na(sw_orig$p.value) && sw_orig$p.value >= 0.05

  if (can_log) {
    x_log_fin  <- log(x)[is.finite(log(x))]
    sw_log     <- if (length(x_log_fin) >= 3) shapiro.test(x_log_fin) else list(statistic = NA, p.value = NA)
    normal_log <- !is.na(sw_log$p.value) && sw_log$p.value >= 0.05
  } else {
    sw_log     <- list(statistic = NA, p.value = NA)
    normal_log <- FALSE
  }

  use_col         <- if (normal_orig) response else paste0(response, "_log_tmp")
  data[[use_col]] <- if (normal_orig) data[[response]] else {
    if (can_log) log(data[[response]]) else data[[response]]
  }

  lev <- tryCatch(
    levene_test(data, as.formula(paste(use_col, "~", group))),
    error = function(e) list(p = NA))
  homogeneous <- !is.na(lev$p) && lev$p >= 0.05

  transform <- if (normal_orig) "None" else if (!can_log) "Not applicable (zero/negative)" else if (normal_log) "log(x)" else "log(x) — insufficient"

  test_used <-
    if      (normal_orig && homogeneous)  "ANOVA + Tukey HSD"
    else if (normal_orig && !homogeneous) "Welch ANOVA + Games-Howell"
    else if (normal_log  && homogeneous)  "ANOVA + Tukey HSD (on log)"
    else if (normal_log  && !homogeneous) "Welch ANOVA + Games-Howell (on log)"
    else                                  "Kruskal-Wallis + Dunn (Bonferroni)"

  data.frame(
    Layer       = layer,    Variable    = variable,
    W_original  = round(sw_orig$statistic, 4),
    p_original  = round(sw_orig$p.value,   4),
    Normal_orig = normal_orig,
    W_log       = ifelse(is.na(sw_log$statistic), NA, round(sw_log$statistic, 4)),
    p_log       = ifelse(is.na(sw_log$p.value),   NA, round(sw_log$p.value,   4)),
    Normal_log  = normal_log,
    Transform   = transform,
    Levene_p    = ifelse(is.na(lev$p), NA, round(lev$p, 4)),
    Homogeneous = ifelse(is.na(lev$p), NA, homogeneous),
    Test_Used   = test_used,
    stringsAsFactors = FALSE
  )
}

build_audit_table <- function(layers_list, vars_per_layer, group = "Vegetation") {
  rows <- list()
  for (lyr in layers_list) {
    lyr$data[[group]] <- as.factor(lyr$data[[group]])
    for (var in vars_per_layer)
      if (var %in% names(lyr$data))
        rows <- append(rows, list(
          audit_variable(lyr$data, var, group, lyr$name, var)))
  }
  do.call(rbind, rows)
}


# =============================================================================
# DATA LOADING
# =============================================================================

df  <- read_excel("Soil_DB.xlsx", sheet = 1)   # Total Cstock (0–100 cm)
df1 <- read_excel("Soil_DB.xlsx", sheet = 2)   # Layer 1 (0–20 cm)
df2 <- read_excel("Soil_DB.xlsx", sheet = 3)   # Layer 2 (20–40 cm)
df3 <- read_excel("Soil_DB.xlsx", sheet = 4)   # Layer 3 (40–60 cm)
df4 <- read_excel("Soil_DB.xlsx", sheet = 5)   # Layer 4 (60–80 cm)
df5 <- read_excel("Soil_DB.xlsx", sheet = 6)   # Layer 5 (80–100 cm)
df6 <- read_excel("Soil_DB.xlsx", sheet = 7)   # All layers combined

DEPTH_LEVELS <- c("P20", "P40", "P60", "P80", "P100")
DEPTH_LABELS <- c("20cm", "40cm", "60cm", "80cm", "100cm")
VEG_COLORS   <- c("SSF" = "skyblue", "DWS" = "lightgreen", "DA" = "coral")

# Variables analysed between treatments (Section A)
VARS_BETWEEN <- c("Cstock", "Sand", "Silt", "Clay", "BD", "C", "N", "N15", "C/N")

# Variables analysed within treatments across depth (Section B)
VARS_WITHIN  <- c("Cstock", "C", "Clay","Silt", "Sand", "BD", "N", "N15", "C/N")   # only vars present in sheet 7

layers <- list(
  list(data = df1, name = "Layer 1 (0-20 cm)"),
  list(data = df2, name = "Layer 2 (20-40 cm)"),
  list(data = df3, name = "Layer 3 (40-60 cm)"),
  list(data = df4, name = "Layer 4 (60-80 cm)"),
  list(data = df5, name = "Layer 5 (80-100 cm)")
)


# =============================================================================
# SECTION D — PLOTS (built first so df6 Depth labels are set for plots only)
# =============================================================================

## -- Plot 1: Total Cstock boxplot (0–100 cm) ---------------------------------
df$Vegetation <- as.factor(df$Vegetation)
anova_tot  <- aov(Cstock_Total ~ Vegetation, data = df)
tukey_tot  <- TukeyHSD(anova_tot)
tukey_p    <- tukey_tot$Vegetation[, "p adj"]
names(tukey_p) <- rownames(tukey_tot$Vegetation)
cld        <- multcompLetters(tukey_p, threshold = 0.05)

letters_df <- data.frame(Vegetation = names(cld$Letters),
                         Letters    = toupper(cld$Letters))
label_pos  <- df %>%
  group_by(Vegetation) %>%
  summarise(y_pos = quantile(Cstock_Total, 0.75, na.rm = TRUE) + 1) %>%
  left_join(letters_df, by = "Vegetation")

p1 <- ggplot(df, aes(x = Vegetation, y = Cstock_Total, fill = Vegetation)) +
  geom_boxplot() +
  geom_text(data = label_pos,
            aes(x = Vegetation, y = y_pos, label = Letters),
            position = position_nudge(x = 0.2, y = 0.7),
            vjust = 0, size = 7) +
  labs(x = "Land Use and Land Cover", y = "SOC Stocks (Mg C ha\u207b\u00b9)") +
  scale_fill_manual(values = VEG_COLORS) +
  theme_bw() +
  theme(axis.title  = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
p2 <- p1 + labs(tag = "A")

## -- Plot 2: Cstock by depth interval ----------------------------------------
df6_plot <- df6 %>%
  mutate(Depth = factor(Depth,
                        levels = DEPTH_LEVELS,
                        labels = c("0-20 cm","20-40 cm","40-60 cm","60-80 cm","80-100 cm")))

p3 <- df6_plot %>%
  filter(Cstock >= 1, Cstock <= 50) %>%
  ggplot(aes(x = Cstock, y = Vegetation, fill = Depth)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(y = "Land Use and Land Cover",
       x = "SOC Stocks (Mg C ha\u207b\u00b9)", fill = "Depth") +
  coord_flip() + theme_bw() +
  theme(axis.title   = element_text(size = 20),
        legend.text  = element_text(size = 18),
        legend.title = element_text(size = 18))
p4 <- p3 + labs(tag = "B")
grid.arrange(p2, p4, widths = 1)

## -- Plot 3: delta13C profile -------------------------------------------------
df6_iso <- df6 %>%
  mutate(Depth     = factor(Depth, levels = DEPTH_LEVELS, labels = DEPTH_LABELS),
         Vegetation = as.factor(Vegetation)) %>%
  group_by(Vegetation, Depth) %>%
  summarise(mean = mean(C13, na.rm = TRUE),
            lci  = t.test(C13, conf.level = 0.95)$conf.int[1],
            uci  = t.test(C13, conf.level = 0.95)$conf.int[2],
            .groups = "drop")

ggplot(df6_iso, aes(x = Depth, y = mean, group = Vegetation, color = Vegetation)) +
  geom_line(position = position_dodge(0.3), linewidth = 1) +
  geom_point(position = position_dodge(0.3), size = 3.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci),
                width = 0.2, linewidth = 0.5,
                position = position_dodge(0.3)) +
  coord_flip() +
  scale_x_discrete(limits = rev(DEPTH_LABELS)) +
  scale_color_manual(values = VEG_COLORS) +
  labs(x = "Depth (cm)",
       y = expression(paste("Mean Carbon Isotope (", delta^{13}, "C \u2030)")),
       color = "Vegetation Type") +
  theme_classic(base_size = 20)


# =============================================================================
# SECTION A — BETWEEN-TREATMENT COMPARISONS PER LAYER
# =============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION A — Between-treatment comparisons (DA vs DWS vs SSF)\n")
cat(strrep("#", 70), "\n\n")

# Total Cstock (0-100 cm) — uses df, not layer sheets
cat(strrep("=", 60), "\n")
cat("  Total Cstock (0-100 cm)\n")
cat(strrep("=", 60), "\n")
analyze_variable(df, "Cstock_Total", group = "Vegetation", label = "Total Cstock")

# Per-layer analysis
for (lyr in layers)
  analyze_layer(lyr$data, VARS_BETWEEN, group = "Vegetation", layer_name = lyr$name)


# =============================================================================
# SECTION B — WITHIN-TREATMENT DEPTH COMPARISONS
# =============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION B — Within-treatment depth comparisons\n")
cat(strrep("#", 70), "\n")

df6_clean <- read_excel("Soil_DB.xlsx", sheet = 7)   # reload without factor relabels

for (veg in c("SSF", "DWS", "DA"))
  for (var in VARS_WITHIN)
    analyze_depth_within_veg(df6_clean, veg, var,
                             depth_levels = DEPTH_LEVELS,
                             depth_labels = DEPTH_LABELS)


# =============================================================================
# SECTION C — SUPPLEMENTARY TABLE S1: ASSUMPTION AUDIT
# =============================================================================

cat("\n", strrep("#", 70), "\n")
cat("  SECTION C — Generating Supplementary Table S1\n")
cat(strrep("#", 70), "\n")

# Between-treatment audit (df1–df5)
supp_between <- build_audit_table(layers, VARS_BETWEEN)

# Add total Cstock row
df$Vegetation <- as.factor(df$Vegetation)
row_total <- audit_variable(df, "Cstock_Total", group = "Vegetation",
                            layer = "All layers (0-100 cm)", variable = "Cstock_Total")

# Within-treatment audit (df6 filtered by vegetation)
df6_clean2 <- read_excel("Soil_DB.xlsx", sheet = 7)
within_rows <- list()
for (veg in c("SSF", "DWS", "DA")) {
  df_veg <- df6_clean2 %>%
    filter(Vegetation == veg) %>%
    mutate(Depth = factor(Depth, levels = DEPTH_LEVELS, labels = DEPTH_LABELS))
  for (var in VARS_WITHIN)
    if (var %in% names(df_veg))
      within_rows <- append(within_rows, list(
        audit_variable(df_veg, var, group = "Depth",
                       layer    = paste("Within", veg),
                       variable = var)
      ))
}
supp_within <- do.call(rbind, within_rows)

# Combine all sections
supp_table <- rbind(row_total, supp_between, supp_within)

# Add human-readable comparison type column
supp_table <- supp_table %>%
  mutate(Comparison = case_when(
    grepl("Within", Layer) ~ "Within treatment (across depth)",
    grepl("0-100",  Layer) ~ "Between treatments (total)",
    TRUE                   ~ "Between treatments (by layer)"
  )) %>%
  select(Comparison, Layer, Variable, everything())

write.csv(supp_table,
          "Supp_Table_S1_Assumption_Checks.csv",
          row.names = FALSE)

cat("\n✔ Supplementary Table S1 saved: Supp_Table_S1_Assumption_Checks.csv\n")
cat(sprintf("  Total rows: %d\n", nrow(supp_table)))
print(supp_table)
