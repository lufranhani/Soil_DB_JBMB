# =============================================================================
# SOIL CARBON STOCK ANALYSIS - REFACTORED VERSION
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
setwd("C:/Users/lufra/OneDrive/buffet/Documentos/Análise Solos")


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Test normality with Shapiro-Wilk and print result
test_normality <- function(x, label = "") {
  result <- shapiro.test(x)
  status <- ifelse(result$p.value < 0.05, "does NOT FOLLOW", "FOLLOWS")
  cat(sprintf("[%s] Variable %s a normal distribution (p = %.8f)\n",
              label, status, result$p.value))
  invisible(result)
}

#' Test variance homogeneity with Levene's Test
test_levene <- function(data, formula, label = "") {
  result <- levene_test(data, formula)
  status <- ifelse(result$p < 0.05, "NOT homogeneous", "homogeneous")
  cat(sprintf("[%s] Levene's Test: variances are %s (p = %.3f)\n",
              label, status, result$p))
  invisible(result)
}

#' Run ANOVA + Tukey HSD and return compact letter display
run_anova_tukey <- function(data, response, group, alpha = 0.05) {
  formula  <- as.formula(paste(response, "~", group))
  model    <- aov(formula, data = data)
  tukey    <- TukeyHSD(model)
  tukey_p  <- tukey[[group]][, "p adj"]
  names(tukey_p) <- rownames(tukey[[group]])
  hsd      <- HSD.test(model, group, alpha = alpha, group = TRUE, console = FALSE)
  print(hsd)
  invisible(list(model = model, tukey = tukey, hsd = hsd))
}

#' Run Kruskal-Wallis + optional Dunn post-hoc test
run_kruskal_dunn <- function(data, response, group,
                             p.adjust = "none", label = "") {
  formula <- as.formula(paste(response, "~", group))
  kruskal <- kruskal_test(data, formula)
  cat(sprintf("\n[%s] Kruskal-Wallis: p = %.4f\n", label, kruskal$p))
  print(kruskal)

  if (kruskal$p < 0.05) {
    posthoc <- dunn_test(data, formula, p.adjust.method = p.adjust)
    print(posthoc)
  } else {
    cat("Kruskal-Wallis not significant. Dunn post-hoc not necessary.\n")
  }
  invisible(kruskal)
}

#' Full normality-aware analysis pipeline for one variable
#' Tries log-transform if original data is non-normal.
#' Chooses ANOVA or Kruskal based on normality of (transformed) data.
analyze_variable <- function(data, response, group = "Vegetation",
                             label = "", p.adjust = "bonferroni") {
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("  Variable: %s | Layer: %s\n", response, label))
  cat(strrep("=", 60), "\n")

  # 1. Test normality on original data
  sw <- test_normality(data[[response]], label)
  normal <- sw$p.value >= 0.05

  # 2. If not normal, try log-transform
  log_col <- paste0(response, "_log")
  if (!normal) {
    data[[log_col]] <- log(data[[response]])
    sw_log <- test_normality(data[[log_col]], paste(label, "log"))
    normal <- sw_log$p.value >= 0.05
    use_col <- if (normal) log_col else response
  } else {
    use_col <- response
  }

  # 3. Test homogeneity
  formula_lev <- as.formula(paste(use_col, "~", group))
  test_levene(data, formula_lev, label)

  # 4. Choose test
  if (normal) {
    run_anova_tukey(data, use_col, group)
  } else {
    run_kruskal_dunn(data, response, group,
                     p.adjust = p.adjust, label = label)
  }
}

#' Bulk analysis of multiple variables in one layer's dataframe
analyze_layer <- function(data, variables, group = "Vegetation",
                          layer_name = "", p.adjust = "bonferroni") {
  data[[group]] <- as.factor(data[[group]])
  for (var in variables) {
    analyze_variable(data, var, group = group,
                     label = paste(layer_name, var),
                     p.adjust = p.adjust)
  }
}


# =============================================================================
# DATA LOADING
# =============================================================================

df  <- read_excel("Soil_DB.xlsx", sheet = 1)   # Total C stocks
df1 <- read_excel("Soil_DB.xlsx", sheet = 2)   # Layer 1 (0–20 cm)
df2 <- read_excel("Soil_DB.xlsx", sheet = 3)   # Layer 2 (20–40 cm)
df3 <- read_excel("Soil_DB.xlsx", sheet = 4)   # Layer 3 (40–60 cm)
df4 <- read_excel("Soil_DB.xlsx", sheet = 5)   # Layer 4 (60–80 cm)
df5 <- read_excel("Soil_DB.xlsx", sheet = 6)   # Layer 5 (80–100 cm)
df6 <- read_excel("Soil_DB.xlsx", sheet = 7)   # All layers combined

DEPTH_LEVELS <- c("P20", "P40", "P60", "P80", "P100")
DEPTH_LABELS <- c("0–20 cm", "20–40 cm", "40–60 cm", "60–80 cm", "80–100 cm")
DEPTH_LABELS_SHORT <- c("20cm", "40cm", "60cm", "80cm", "100cm")
VEG_COLORS <- c("SSF" = "skyblue", "DWS" = "lightgreen", "DA" = "coral")

# Relabel depth in df6
df6$Depth <- factor(df6$Depth, levels = DEPTH_LEVELS, labels = DEPTH_LABELS)


# =============================================================================
# SECTION 1 – TOTAL C STOCKS (0–100 cm)
# =============================================================================

df$Vegetation <- as.factor(df$Vegetation)
analyze_variable(df, "Cstock_Total", group = "Vegetation", label = "Total")

# --- Build CLD for boxplot labels ---
anova_total  <- aov(Cstock_Total ~ Vegetation, data = df)
tukey_total  <- TukeyHSD(anova_total)
tukey_p_tot  <- tukey_total$Vegetation[, "p adj"]
names(tukey_p_tot) <- rownames(tukey_total$Vegetation)
cld          <- multcompLetters(tukey_p_tot, threshold = 0.05)

letters_df <- data.frame(
  Vegetation = names(cld$Letters),
  Letters    = toupper(cld$Letters)
)

label_pos <- df %>%
  group_by(Vegetation) %>%
  summarise(y_pos = quantile(Cstock_Total, 0.75, na.rm = TRUE) + 1) %>%
  left_join(letters_df, by = "Vegetation")

p1 <- ggplot(df, aes(x = Vegetation, y = Cstock_Total, fill = Vegetation)) +
  geom_boxplot() +
  geom_text(data = label_pos,
            aes(x = Vegetation, y = y_pos, label = Letters),
            position = position_nudge(x = 0.2, y = 0.7),
            vjust = 0, size = 7) +
  labs(x = "Land Use and Land Cover", y = "SOC Stocks (Mg C ha⁻¹)") +
  scale_fill_manual(values = VEG_COLORS) +
  theme_bw() +
  theme(axis.title  = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))

p2 <- p1 + labs(tag = "A")


# =============================================================================
# SECTION 2 – C STOCKS BY DEPTH (boxplot panel B)
# =============================================================================

p3 <- df6 %>%
  filter(Cstock >= 1, Cstock <= 50) %>%
  ggplot(aes(x = Cstock, y = Vegetation, fill = Depth)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  labs(y = "Land Use and Land Cover",
       x = "SOC Stocks (Mg C ha⁻¹)", fill = "Depth") +
  coord_flip() +
  theme_bw() +
  theme(axis.title   = element_text(size = 20),
        legend.text  = element_text(size = 18),
        legend.title = element_text(size = 18))

p4 <- p3 + labs(tag = "B")

grid.arrange(p2, p4, widths = 1)


# =============================================================================
# SECTION 3 – CARBON ISOTOPE PROFILE (δ¹³C)
# =============================================================================

if (!"Depth_raw" %in% names(df6))
  df6$Depth_raw <- rep(DEPTH_LEVELS, length.out = nrow(df6))

df6$Depth_short <- factor(df6$Depth_raw,
                          levels = DEPTH_LEVELS,
                          labels = DEPTH_LABELS_SHORT)
df6$Vegetation  <- as.factor(df6$Vegetation)

dt_iso <- df6 %>%
  group_by(Vegetation, Depth_short) %>%
  summarise(mean = mean(C13, na.rm = TRUE),
            lci  = t.test(C13, conf.level = 0.95)$conf.int[1],
            uci  = t.test(C13, conf.level = 0.95)$conf.int[2],
            .groups = "drop")

ggplot(dt_iso, aes(x = Depth_short, y = mean,
                   group = Vegetation, color = Vegetation)) +
  geom_line(position = position_dodge(0.3), linewidth = 1) +
  geom_point(position = position_dodge(0.3), size = 3.5) +
  geom_errorbar(aes(ymin = lci, ymax = uci),
                width = 0.2, linewidth = 0.5,
                position = position_dodge(0.3)) +
  coord_flip() +
  scale_x_discrete(limits = rev(DEPTH_LABELS_SHORT)) +
  scale_color_manual(values = VEG_COLORS) +
  labs(x = "Depth (cm)",
       y = expression(paste("Mean Carbon Isotope (", delta^{13}, "C \u2030)")),
       color = "Vegetation Type") +
  theme_classic(base_size = 20)


# =============================================================================
# SECTION 3.5 – SUPPLEMENTARY TABLE: NORMALITY & HOMOGENEITY AUDIT
# =============================================================================

#' Build one row of the assumption-checking audit table
audit_variable <- function(data, response, group = "Vegetation", 
                           layer = "", variable = "") {
  x <- data[[response]]
  x <- x[!is.na(x)]
  
  # Shapiro-Wilk on original
  sw_orig  <- shapiro.test(x)
  normal_orig <- sw_orig$p.value >= 0.05
  
  # Log transform attempt
  x_log       <- log(x)
  sw_log      <- shapiro.test(x_log)
  normal_log  <- sw_log$p.value >= 0.05
  
  # Levene on whichever column passed (or log if neither passed)
  use_col <- if (normal_orig) response else paste0(response, "_log_tmp")
  data[[use_col]] <- if (normal_orig) data[[response]] else log(data[[response]])
  
  formula_lev <- as.formula(paste(use_col, "~", group))
  lev <- tryCatch(
    levene_test(data, formula_lev),
    error = function(e) list(p = NA)
  )
  homogeneous <- !is.na(lev$p) && lev$p >= 0.05
  
  # Decision logic
  if (normal_orig) {
    transform  <- "None"
    sw_w_used  <- round(sw_orig$statistic, 4)
    sw_p_used  <- round(sw_orig$p.value, 4)
    test_used  <- ifelse(homogeneous, "ANOVA + Tukey HSD", "ANOVA + Tukey HSD*")
    rationale  <- ifelse(homogeneous,
                         "Normal + homogeneous → ANOVA",
                         "Normal but heterogeneous → ANOVA (Welch recommended)")
  } else if (normal_log) {
    transform  <- "log(x)"
    sw_w_used  <- round(sw_log$statistic, 4)
    sw_p_used  <- round(sw_log$p.value, 4)
    test_used  <- ifelse(homogeneous, "ANOVA + Tukey HSD", "ANOVA + Tukey HSD*")
    rationale  <- ifelse(homogeneous,
                         "Log normalized + homogeneous → ANOVA on log",
                         "Log normalized but heterogeneous → ANOVA on log (Welch recommended)")
  } else {
    transform  <- "log(x) — insufficient"
    sw_w_used  <- round(sw_log$statistic, 4)
    sw_p_used  <- round(sw_log$p.value, 4)
    test_used  <- "Kruskal-Wallis + Dunn"
    rationale  <- "Neither original nor log met normality → non-parametric"
  }
  
  data.frame(
    Layer           = layer,
    Variable        = variable,
    W_original      = round(sw_orig$statistic, 4),
    p_original      = round(sw_orig$p.value,   4),
    W_log           = round(sw_log$statistic,  4),
    p_log           = round(sw_log$p.value,    4),
    Transform_Used  = transform,
    Levene_p        = ifelse(is.na(lev$p), "NA", round(lev$p, 4)),
    Homogeneous     = ifelse(is.na(lev$p), "NA", ifelse(homogeneous, "Yes", "No")),
    Test_Used       = test_used,
    Rationale       = rationale,
    stringsAsFactors = FALSE
  )
}

#' Run audit across all layers and variables, return combined table
build_audit_table <- function(layers_list, vars_per_layer, group = "Vegetation") {
  rows <- list()
  for (lyr in layers_list) {
    lyr$data[[group]] <- as.factor(lyr$data[[group]])
    for (var in vars_per_layer) {
      if (var %in% names(lyr$data)) {
        row <- audit_variable(lyr$data, var, group = group,
                              layer    = lyr$name,
                              variable = var)
        rows <- append(rows, list(row))
      }
    }
  }
  do.call(rbind, rows)
}

layers <- list(
  list(data = df1, name = "Layer 1 (0–20 cm)"),
  list(data = df2, name = "Layer 2 (20–40 cm)"),
  list(data = df3, name = "Layer 3 (40–60 cm)"),
  list(data = df4, name = "Layer 4 (60–80 cm)"),
  list(data = df5, name = "Layer 5 (80–100 cm)")
)

audit_vars  <- c("Cstock", "C", "Sand", "Silt", "Clay", "BD")
audit_table <- build_audit_table(layers, audit_vars)

# Total C stock (df, all depths)
df$Vegetation <- as.factor(df$Vegetation)
row_total <- audit_variable(df, "Cstock_Total", group = "Vegetation",
                            layer = "All layers (0–100 cm)", variable = "Cstock_Total")
audit_table <- rbind(row_total, audit_table)

# Export to CSV
write.csv(audit_table,
          "Supplementary_Table_Assumption_Checks.csv",
          row.names = FALSE)

cat("\n✔ Supplementary table saved: Supplementary_Table_Assumption_Checks.csv\n")
print(audit_table)


# =============================================================================
# SECTION 4 – SOIL TEXTURE & BULK DENSITY PER LAYER
# =============================================================================

TEXTURE_VARS <- c("Sand", "Silt", "Clay")
BD_VAR       <- "BD"

layers <- list(
  list(data = df1, name = "Layer1"),
  list(data = df2, name = "Layer2"),
  list(data = df3, name = "Layer3"),
  list(data = df4, name = "Layer4"),
  list(data = df5, name = "Layer5")
)

for (lyr in layers) {
  analyze_layer(lyr$data, TEXTURE_VARS, layer_name = lyr$name)
  analyze_layer(lyr$data, BD_VAR,       layer_name = lyr$name)
}


# =============================================================================
# SECTION 5 – CARBON CONCENTRATION PER LAYER
# =============================================================================

C_VAR <- "C"

for (lyr in layers) {
  analyze_layer(lyr$data, C_VAR, layer_name = lyr$name)
}


# =============================================================================
# SECTION 6 – C STOCKS BY DEPTH WITHIN EACH VEGETATION TYPE
# =============================================================================

# Reload df6 cleanly (rm() removed — utility functions must stay in memory)
df6 <- read_excel("Soil_DB.xlsx", sheet = 7)

DEPTH_LEVELS_2 <- c("P20", "P40", "P60", "P80", "P100")
DEPTH_LABELS_2 <- c("20cm", "40cm", "60cm", "80cm", "100cm")

#' Analyze C stocks across depth intervals for a single vegetation type
analyze_depth_within_veg <- function(df6, veg_name) {
  cat(sprintf("\n%s\n  Vegetation: %s\n%s\n",
              strrep("=", 60), veg_name, strrep("=", 60)))

  df_veg <- df6 %>%
    filter(Vegetation == veg_name) %>%
    mutate(Depth = factor(Depth,
                          levels = DEPTH_LEVELS_2,
                          labels = DEPTH_LABELS_2))

  sw <- test_normality(df_veg$Cstock, veg_name)

  df_veg$Cstock_log <- log(df_veg$Cstock)
  sw_log <- test_normality(df_veg$Cstock_log, paste(veg_name, "log"))
  normal <- sw_log$p.value >= 0.05

  if (normal) {
    lev <- levene_test(df_veg, Cstock_log ~ Depth)
    cat(sprintf("Levene's p = %.3f\n", lev$p))
    run_anova_tukey(df_veg, "Cstock_log", "Depth")
  } else {
    run_kruskal_dunn(df_veg, "Cstock", "Depth",
                     p.adjust = "bonferroni", label = veg_name)
  }
}

for (veg in c("SSF", "DWS", "DA")) {
  analyze_depth_within_veg(df6, veg)
}
