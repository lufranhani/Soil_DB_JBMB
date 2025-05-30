library(tidyverse)    
library(rstatix)      
library(dplyr)        
library(ggplot2)      
library(FSA)          
library(multcompView) 
library (agricolae)
library (readxl)
library(gridExtra)
set.seed(1234)

### ----- Clean R's workspaace
rm(list=ls())

### ----- Set working directory
#setwd("")

#Read table in Directory
df <- read_excel("Soil_DB.xlsx")
df1 <- read_excel("Soil_DB.xlsx", sheet = 2)
df2 <- read_excel("Soil_DB.xlsx", sheet = 3)
df3 <- read_excel("Soil_DB.xlsx", sheet = 4)
df4 <- read_excel("Soil_DB.xlsx", sheet = 5)
df5 <- read_excel("Soil_DB.xlsx", sheet = 6)
df6 <- read_excel("Soil_DB.xlsx", sheet = 7)

#Carbon Stocks plots and 0 - 100 cm statistical comparison #######################
#Perform the Shapiro-Wilk test on the variable
shapiro_test_result <- shapiro.test(df$Cstock_Total)

# Print the result
if (shapiro_test_result$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_test_result$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_test_result$p.value, 8), ")."))}

qqnorm(df$Cstock_Total); qqline(df$Cstock_Total)
hist(df$Cstock_Total, main = "Histogram of Cstock", xlab = "Cstock_Total")

#Convert Vegetation to a factor
df$Vegetation <- as.factor(df$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_result <- df %>% 
  levene_test(Cstock_Total ~ Vegetation)

# Extract p-value
p_val <- levene_result$p

# Create interpretation phrase
if (p_val < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val, 3), ").")
}

# Print result
print(message)

# Run ANOVA
anova_result <- aov(Cstock_Total ~ Vegetation, data = df)

# Tukey's HSD
tukey_result <- TukeyHSD(anova_result)

# Extract adjusted p-values
tukey_p <- tukey_result$Vegetation[, "p adj"]
names(tukey_p) <- rownames(tukey_result$Vegetation)

# Compact Letter Display (CLD) - Uppercase
cld <- multcompLetters(tukey_p, threshold = 0.05)
letters_df <- data.frame(
  Vegetation = names(cld$Letters),
  Letters = toupper(cld$Letters))

# Calculate 75th percentile (upper quartile) for label positioning
label_positions <- df %>%
  group_by(Vegetation) %>%
  summarise(y_pos = quantile(Cstock_Total, 0.75, na.rm = TRUE) + 1) %>%
  left_join(letters_df, by = "Vegetation")

# Create boxplot with uppercase letters above
p1 <- ggplot(df, aes(x = Vegetation, y = Cstock_Total, fill = Vegetation)) +
  geom_boxplot() +
  geom_text(
    data = label_positions,
    aes(x = Vegetation, y = y_pos, label = Letters),
    position = position_nudge(x = 0.2),  
    vjust = 0,
    size = 5,
    fontface = "plain") +
  labs(x = "Treatments", y = "SOC Stocks (Mg C ha⁻¹)") +
  theme_bw() +
  theme (axis.title.x = element_text(size = 20), 
         axis.title.y = element_text(size = 20), 
         legend.text = element_text(size = 18), 
         legend.title = element_text(size = 18)) +
  scale_fill_manual(values = c("SSF" = "skyblue", "DWS" = "lightgreen", "DA" = "coral"))
p1 
p2 = p1 + labs(tag = "a")

# Defining new name for the variables
new_labels  <- c("0 - 20 cm", "20 - 40 cm", "40 - 60 cm", "60 - 80 cm", "80 - 100 cm")
# Defining new order to the levels 
order_levels <- c("P20", "P40", "P60", "P80", "P100")
# Converting Depth to a factor with new names 
df6$Depth <- factor(df6$Depth, levels = order_levels, labels = new_labels)

# PLOT Carbon Stocks by Depth Intervals 
p3 <- df6 %>%
  filter(Cstock >= 1L & Cstock <= 50L) %>%
  ggplot() +
  aes(x = Cstock, y = Vegetation, fill = Depth) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdYlGn", 
                    direction = -1) +
  labs(y = "Treatments", x = 'SOC Stocks (Mg C ha⁻¹)', fill = "Depth") +
  coord_flip() +
  theme_bw()  +
  theme (axis.title.x = element_text(size = 20), 
         axis.title.y = element_text(size = 20),
         legend.text = element_text(size = 18), 
         legend.title = element_text(size = 18))
p3 
p4 = p3 + labs(tag = 'b')
p4 

# Arrange the plots one above the other
grid.arrange(p2, p4, widths = 1)

# 1st Soil Layer - Cstocks #######################
# Perform the Shapiro-Wilk test on the variable
shapiro_test_result1 <- shapiro.test(df1$Cstock)

# Print the result
if (shapiro_test_result1$p.value < 0.05) { print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
round(shapiro_test_result1$p.value, 8), ")."))} else { print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_test_result1$p.value, 8), ")."))}

# Plotting distribution
qqnorm(df1$Cstock); qqline(df1$Cstock)
hist(df1$Cstock, main = "Histogram of Cstock by Depth Intervals", xlab = "Cstock")

# Create a new column with log from Cstock 
df1$Cstock_log <- log(df1$Cstock)

# Perform the Shapiro-Wilk test on the variable
shapiro_test_result1 <- shapiro.test(df1$Cstock_log)

# Print the result
if (shapiro_test_result1$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_test_result1$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_test_result1$p.value, 8), ")."))}

# Plotting distribution
qqnorm(df1$Cstock_log); qqline(df1$Cstock_log)
hist(df1$Cstock_log, main = "Histogram of Cstock by Depth Intervals - 1 layer", xlab = "Cstock")

# Convert Vegetation to a factor
df1$Vegetation <- as.factor(df1$Vegetation)

# Main test
kruskal_result <- df1 %>% kruskal_test(Cstock ~ Vegetation)

# Result
print(kruskal_result)

# Verifique se o p-valor é significativo antes de continuar
if(kruskal_result$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df1 %>% dunn_test(Cstock ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# 2nd Soil Layer - Cstocks#######################
# Perform the Shapiro-Wilk test on the variable
shapiro_test_result2 <- shapiro.test(df2$Cstock)

# Print the result
if (shapiro_test_result2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", round(shapiro_test_result2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =",round(shapiro_test_result2$p.value, 8), ")."))}

#Plotting distribution
qqnorm(df2$Cstock); qqline(df2$Cstock)
hist(df2$Cstock, main = "Histogram of Cstock by Depth Intervals - 2 layer", xlab = "Cstock")

# Convert Cstock to log to minimize effect and meet assumption
df2$Cstock_log <- log(df2$Cstock)

# Shapiro 
shapiro_test_result2 <- shapiro.test(df2$Cstock_log)

# Print the result
if (shapiro_test_result2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", round(shapiro_test_result2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =",round(shapiro_test_result2$p.value, 8), ")."))}

#Convert Vegetation to a factor
df2$Vegetation <- as.factor(df2$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_result2 <- df2%>%levene_test(Cstock_log ~ Vegetation)

# Extract p-value
p_val2 <- levene_result2$p

# Create interpretation phrase
if (p_val2 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val2, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val2, 3), ").")
}

# Print result
print(message)

#Perfom ANOVA
anova2 <- aov(Cstock_log ~ Vegetation, data = df2)

# Tukey's HSD
tukey_result2 <- TukeyHSD(anova2)

# Extract adjusted p-values
tukey_p2 <- tukey_result2$Vegetation[, "p adj"]
names(tukey_p2) <- rownames(tukey_result2$Vegetation)
print (tukey_p2, 3)
hsd_agricolae2 <- HSD.test(anova2, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_agricolae2)

# 3rd Soil Layer - Cstocks#######################
# Perform the Shapiro-Wilk test on the variable
shapiro_test_result3 <- shapiro.test(df3$Cstock)

# Print the result
if (shapiro_test_result3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", round(shapiro_test_result3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =",round(shapiro_test_result3$p.value, 8), ")."))}

#Plotting distribution
qqnorm(df3$Cstock); qqline(df3$Cstock)
hist(df3$Cstock, main = "Histogram of Cstock by Depth Intervals - 3 layer", xlab = "Cstock")

#Convert Vegetation to a factor
df3$Vegetation <- as.factor(df3$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_result3 <- df3%>%levene_test(Cstock ~ Vegetation)

# Extract p-value
p_val3 <- levene_result3$p

# Create interpretation phrase
if (p_val3 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val3, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val3, 3), ").")
}

# Print result
print(message)

kruskal3 <- kruskal(df3$Cstock, df3$Vegetation, alpha = 0.05, p.adj=c("none"), group = TRUE, main = NULL, console = FALSE)
print (kruskal3)

# 4th Soil Layer - Cstocks#######################
# Perform the Shapiro-Wilk test on the variable
shapiro_test_result4 <- shapiro.test(df4$Cstock)

# Print the result
if (shapiro_test_result4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", round(shapiro_test_result4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =",round(shapiro_test_result4$p.value, 8), ")."))}

#Plotting distribution
qqnorm(df4$Cstock); qqline(df4$Cstock)
hist(df4$Cstock, main = "Histogram of Cstock by Depth Intervals - 4 Layer", xlab = "Cstock")

#Convert Vegetation to a factor
df4$Vegetation <- as.factor(df4$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_result4 <- df4%>%levene_test(Cstock ~ Vegetation)

# Extract p-value
p_val4 <- levene_result4$p

# Create interpretation phrase
if (p_val4 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val4, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val4, 3), ").")
}

# Print result
print(message)

#Perfom ANOVA
anova4 <- aov(Cstock ~ Vegetation, data = df4)

# Tukey's HSD
tukey_result4 <- TukeyHSD(anova4)

# Extract adjusted p-values
tukey_p4 <- tukey_result4$Vegetation[, "p adj"]
names(tukey_p4) <- rownames(tukey_result4$Vegetation)
print (tukey_p4, 3)
hsd_agricolae4 <- HSD.test(anova4, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_agricolae4)

# 5th Soil Layer - Cstocks#######################
# Perform the Shapiro-Wilk test on the variable
shapiro_test_result5 <- shapiro.test(df5$Cstock)

# Print the result
if (shapiro_test_result5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", round(shapiro_test_result5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =",round(shapiro_test_result5$p.value, 8), ")."))}

#Plotting distribution
qqnorm(df5$Cstock); qqline(df5$Cstock)
hist(df5$Cstock, main = "Histogram of Cstock by Depth Intervals - 5 Layer", xlab = "Cstock")

#Convert Vegetation to a factor
df5$Vegetation <- as.factor(df5$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_result5 <- df5%>%levene_test(Cstock ~ Vegetation)

# Extract p-value
p_val5 <- levene_result5$p

# Create interpretation phrase
if (p_val5 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val5, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val5, 3), ").")
}

# Print result
print(message)

#Perfom ANOVA
anova5 <- aov(Cstock ~ Vegetation, data = df5)

# Tukey's HSD
tukey_result5 <- TukeyHSD(anova5)

# Extract adjusted p-values
tukey_p5 <- tukey_result5$Vegetation[, "p adj"]
names(tukey_p5) <- rownames(tukey_result5$Vegetation)
print (tukey_p5, 3)
hsd_agricolae5 <- HSD.test(anova5, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_agricolae5)

# Carbon Isotope Plot #######################
#Carbon Isotope 
#Calculating 
if (!"Depth_raw" %in% names(df6)) {
  df6$Depth_raw <- rep(c("P20", "P40", "P60", "P80", "P100"), length.out = nrow(df6))
}

new_labels_1 <- c("20cm", "40cm", "60cm", "80cm", "100cm")
order_levels_1 <- c("P20", "P40", "P60", "P80", "P100")
df6$Depth <- factor(df6$Depth_raw, levels = order_levels_1, labels = new_labels_1)
df6$Vegetation <- as.factor(df6$Vegetation)

dt_grouped <- df6 %>%
  dplyr::group_by(Vegetation, Depth) %>%
  dplyr::summarise(
    mean = mean(C13, na.rm = TRUE),
    lci = t.test(C13, conf.level = 0.95)$conf.int[1],
    uci = t.test(C13, conf.level = 0.95)$conf.int[2],
    .groups = 'drop'
  )


collors <- c("DA" = "coral", 
            "DWS" = "lightgreen",
            "SSF" = "skyblue")


pl2_grouped <- ggplot(data = dt_grouped,
                                aes(x = Depth,
                                    y = mean,
                                    group = Vegetation,
                                    color = Vegetation)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.3), size = 3.5) +
  geom_errorbar(
    aes(ymin = lci, ymax = uci),
    width = 0.2,
    linewidth = 0.5,
    position = position_dodge(width = 0.3)
  ) +
  theme_classic(base_size = 20) +
  labs(
    title = "",
    x = "Depth (cm)",
    y = expression(paste("Mean Carbon Isotope (", delta^{13}, "C \u2030)")),
    color = "Vegetation Type"
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev(new_labels_1)) +
  scale_color_manual(values = collors)

print(pl2_grouped)

# Comparison of soil textures between treatments 1st Layer #######################
# Sand
# Perform the Shapiro-Wilk test on the variable
shapiro_sand <- shapiro.test(df1$Sand)

# Print the result
if (shapiro_sand$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand$p.value, 8), ")."))}

qqnorm(df1$Sand); qqline(df1$Sand)
hist(df1$Sand, main = "Histogram of Sand", xlab = "Sand")

# Create a new column with log from Sand 
df1$Sand_log <- log(df1$Sand)

# Perform the Shapiro-Wilk test on the variable
shapiro_sand_log <- shapiro.test(df1$Sand_log)

# Print the result
if (shapiro_sand_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand_log$p.value, 8), ")."))}

qqnorm(df1$Sand_log); qqline(df1$Sand_log)
hist(df1$Sand_log, main = "Histogram of Sand_log", xlab = "Sand_log")

# Convert Vegetation to a factor
df1$Vegetation <- as.factor(df1$Vegetation)

# Main test
kruskal_Sand <- df1 %>% kruskal_test(Sand ~ Vegetation)

# Result
print(kruskal_Sand)

# Verifique se o p-valor é significativo antes de continuar
if(kruskal_Sand$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df1 %>% dunn_test(Sand ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# Silt concentration
# Perform the Shapiro-Wilk test on the variable
shapiro_Silt <- shapiro.test(df1$Silt)

# Print the result
if (shapiro_Silt$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Silt$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Silt$p.value, 8), ")."))}

#Convert Vegetation to a factor
df1$Vegetation <- as.factor(df1$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_Silt <- df1 %>% 
  levene_test(Silt ~ Vegetation)

# Extract p-value
p_val_Silt <- levene_Silt$p

# Create interpretation phrase
if (p_val_Silt < 0.05) {
  message_Silt <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_Silt, 3), ").")
} else {
  message_Silt <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_Silt, 3), ").")
}

# Print result
print(message_Silt)

#Perfom ANOVA
anova_Silt <- aov(Silt ~ Vegetation, data = df1)

# Tukey's HSD
tukey_Silt <- TukeyHSD(anova_Silt)

# Extract adjusted p-values
tukey_Silt_R <- tukey_Silt$Vegetation[, "p adj"]
names(tukey_Silt_R) <- rownames(tukey_Silt$Vegetation)
print (tukey_Silt_R, 3)
hsd_Silt <- HSD.test(anova_Silt, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_Silt)

# Clay concentration
# Perform the Shapiro-Wilk test on the variable
shapiro_Clay <- shapiro.test(df1$Clay)

# Print the result
if (shapiro_Clay$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Clay$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Clay$p.value, 8), ")."))}

qqnorm(df1$Clay); qqline(df1$Clay)
hist(df1$Clay, main = "Histogram of Clay", xlab = "Clay")

# Create a new column with log from Clay 
df1$Clay_log <- log(df1$Clay)

# Perform the Shapiro-Wilk test on the variable
shapiro_Clay_log <- shapiro.test(df1$Clay_log)

# Print the result
if (shapiro_Clay_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Clay_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Clay_log$p.value, 8), ")."))}

qqnorm(df1$Clay_log); qqline(df1$Clay_log)
hist(df1$Clay_log, main = "Histogram of Clay_log", xlab = "Clay")

#Convert Vegetation to a factor
df1$Vegetation <- as.factor(df1$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_Clay <- df1 %>% 
  levene_test(Clay_log ~ Vegetation)

# Extract p-value
p_val_Clay <- levene_Clay$p

# Create interpretation phrase
if (p_val_Clay < 0.05) {
  message_Clay <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_Clay, 3), ").")
} else {
  message_Clay <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_Clay, 3), ").")
}

# Print result
print(message_Clay)

anova_Clay <- aov(Clay_log ~ Vegetation, data = df1)

# Tukey's HSD
tukey_Clay <- TukeyHSD(anova_Clay)

# Extract adjusted p-values
tukey_Clay_R <- tukey_Clay$Vegetation[, "p adj"]
names(tukey_Clay_R) <- rownames(tukey_Clay$Vegetation)
print (tukey_Clay_R, 3)
hsd_Clay <- HSD.test(anova_Clay, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_Clay)

# Comparison of soil textures between treatments 2nd Layer #######################
#Sand 
shapiro_sand2 <- shapiro.test(df2$Sand)

if (shapiro_sand2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand2$p.value, 8), ")."))}

qqnorm(df2$Sand); qqline(df2$Sand)
hist(df2$Sand, main = "Histogram of Sand", xlab = "Sand")

df2$Vegetation <- as.factor(df2$Vegetation)

levene_Sand2 <- df2 %>% 
  levene_test(Sand ~ Vegetation)

p_val_Sand2 <- levene_Sand2$p

if (p_val_Sand2 < 0.05) {
  message_Sand2 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_Sand2, 3), ").")
} else {
  message_Sand2 <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_Sand2, 3), ").")
}

print(message_Sand2)

anova_Sand2 <- aov (Sand ~ Vegetation, data = df2)
tukey_Sand2 <- TukeyHSD(anova_Sand2)
tukey_Sand2_R <- tukey_Sand2$Vegetation[, "p adj"]
names(tukey_Sand2_R) <- rownames(tukey_Sand2$Vegetation)
print (tukey_Sand2_R, 3)
hsd_Sand2 <- HSD.test(anova_Sand2, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_Sand2)


#Silt 
shapiro_Silt2 <- shapiro.test(df2$Silt)

if (shapiro_Silt2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Silt2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Silt2$p.value, 8), ")."))}


qqnorm(df2$Silt); qqline(df2$Silt)
hist(df2$Silt, main = "Histogram of Silt", xlab = "Silt")

df2$Silt_log2 <- log(df2$Silt)

shapiro_Silt_log2 <- shapiro.test(df2$Silt_log2)

if (shapiro_Silt_log2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Silt_log2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Silt_log2$p.value, 8), ")."))}

df2$Vegetation <- as.factor(df2$Vegetation)

levene_Silt_log2 <- df2 %>% 
  levene_test(Silt_log2 ~ Vegetation)

p_val_Silt2 <- levene_Silt_log2$p

if (p_val_Silt2 < 0.05) {
  message_Silt2 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_Silt2, 3), ").")
} else {
  message_Silt2 <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_Silt2, 3), ").")
}

print(message_Silt2)

anova_Silt_log2 <- aov (Silt_log2 ~ Vegetation, data = df2)
tukey_Silt_log2 <- TukeyHSD(anova_Silt_log2)
tukey_Silt_log2_R <- tukey_Silt_log2$Vegetation[, "p adj"]
names(tukey_Silt_log2_R) <- rownames(tukey_Silt_log2$Vegetation)
print (tukey_Silt_log2_R, 3)
hsd_Silt_log2 <- HSD.test(anova_Silt_log2, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_Silt_log2)

#Clay 
Shapiro_Clay2 <- shapiro.test(df2$Clay)
if (Shapiro_Clay2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Clay2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Clay2$p.value, 8), ")."))}

df2$Vegetation <- as.factor(df2$Vegetation)

levene_Clay2 <- df2 %>% 
  levene_test(Clay ~ Vegetation)

p_val_Clay2 <- levene_Clay2$p

if (p_val_Clay2 < 0.05) {
  message_Clay2 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_Clay2, 3), ").")
} else {
  message_Clay2 <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_Clay2, 3), ").")
}
print(message_Clay2)

anova_Clay2 <- aov (Clay ~ Vegetation, data = df2)
tukey_Clay2 <- TukeyHSD(anova_Clay2)
tukey_Clay2_R <- tukey_Clay2$Vegetation[, "p adj"]
names(tukey_Clay2_R) <- rownames(tukey_Clay2$Vegetation)
print (tukey_Clay2_R, 3)
HSD_Clay2 <- HSD.test(anova_Clay2, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print(HSD_Clay2)

# Comparison of soil textures between treatments 3rd Layer #######################
# Sand 
shapiro_sand3 <- shapiro.test(df3$Sand)

if (shapiro_sand3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand3$p.value, 8), ")."))}

qqnorm(df3$Sand); qqline(df3$Sand)
hist(df3$Sand, main = "Histogram of Sand", xlab = "Sand")

df3$Vegetation <- as.factor(df3$Vegetation)

levene_Sand3 <- df3 %>% levene_test(Sand ~ Vegetation)

p_val_Sand3 <- levene_Sand3$p

if (p_val_Sand3 < 0.05) {
  message_Sand3 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =",
                         round(p_val_Sand3, 3), ").")
} else {
  message_Sand3 <- paste("Levene’s Test indicates that variances are homogeneous (p =", 
                         round(p_val_Sand3, 3), ").")
}

print(message_Sand3)

anova_Sand3 <- aov (Sand ~ Vegetation, data = df3)
tukey_Sand3 <- TukeyHSD(anova_Sand3)
tukey_Sand3_R <- tukey_Sand3$Vegetation[, "p adj"]
names(tukey_Sand3_R) <- rownames(tukey_Sand3$Vegetation)
print (tukey_Sand3_R, 3)
hsd_Sand3 <- HSD.test(anova_Sand3, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_Sand3)

#Silt 
shapiro_Silt3 <- shapiro.test(df3$Silt)

if (shapiro_Silt3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Silt3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Silt3$p.value, 8), ")."))}


qqnorm(df3$Silt); qqline(df3$Silt)
hist(df3$Silt, main = "Histogram of Silt", xlab = "Silt")

df3$Vegetation <- as.factor(df3$Vegetation)

levene_Silt3 <- df3 %>% levene_test(Silt ~ Vegetation)

p_val_Silt3 <- levene_Silt3$p

if (p_val_Silt3 < 0.05) {
  message_Silt3 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", 
                         round(p_val_Silt3, 3), ").")
} else {
  message_Silt3 <- paste("Levene’s Test indicates that variances are homogeneous (p =", 
                         round(p_val_Silt3, 3), ").")
}

print(message_Silt3)

anova_Silt3 <- aov (Silt ~ Vegetation, data = df3)
tukey_Silt3 <- TukeyHSD(anova_Silt3)
tukey_Silt3_R <- tukey_Silt3$Vegetation[, "p adj"]
names(tukey_Silt3_R) <- rownames(tukey_Silt3$Vegetation)
print (tukey_Silt3_R, 3)
HSD_Silt3 <- HSD.test(anova_Silt3, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (HSD_Silt3)

#Clay 
Shapiro_Clay3 <- shapiro.test(df3$Clay)
if (Shapiro_Clay3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Clay3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Clay3$p.value, 8), ")."))}

df3$Vegetation <- as.factor(df3$Vegetation)

levene_Clay3 <- df3 %>%  levene_test(Clay ~ Vegetation)

p_val_Clay3 <- levene_Clay3$p

if (p_val_Clay3 < 0.05) {
  message_Clay3 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", 
                         round(p_val_Clay3, 3), ").")
} else {
  message_Clay3 <- paste("Levene’s Test indicates that variances are homogeneous (p =",
                         round(p_val_Clay3, 3), ").")
}
print(message_Clay3)

anova_Clay3 <- aov (Clay ~ Vegetation, data = df3)
tukey_Clay3 <- TukeyHSD(anova_Clay3)
tukey_Clay3_R <- tukey_Clay3$Vegetation[, "p adj"]
names(tukey_Clay3_R) <- rownames(tukey_Clay3$Vegetation)
print (tukey_Clay3_R, 3)
HSD_Clay3 <- HSD.test(anova_Clay3, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print(HSD_Clay3)


# Comparison of soil textures between treatments 4th Layer #######################
#Sand
shapiro_sand4 <- shapiro.test(df4$Sand)

if (shapiro_sand4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand4$p.value, 8), ")."))}

qqnorm(df4$Sand); qqline(df4$Sand)
hist(df4$Sand, main = "Histogram of Sand", xlab = "Sand")

df4$Sand_log4 <- log(df4$Sand)

shapiro_sand_log4 <- shapiro.test(df4$Sand_log4)

if (shapiro_sand_log4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand_log4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand_log4$p.value, 8), ")."))}

qqnorm(df4$Sand_log4); qqline(df4$Sand_log4)
hist(df4$Sand_log4, main = "Histogram of Sand", xlab = "Sand")

kruskal_Sand4 <- df4 %>% kruskal_test(Sand ~ Vegetation)

print(kruskal_Sand4)

# Confirm p-value before continue
if(kruskal_Sand4$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df4 %>% dunn_test(Sand ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

#Silt 
shapiro_Silt4 <- shapiro.test(df4$Silt)

if (shapiro_Silt4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Silt4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Silt4$p.value, 8), ")."))}


qqnorm(df4$Silt); qqline(df4$Silt)
hist(df4$Silt, main = "Histogram of Silt", xlab = "Silt")

df4$Vegetation <- as.factor(df4$Vegetation)

levene_Silt4 <- df4 %>% levene_test(Silt ~ Vegetation)

p_val_Silt4 <- levene_Silt4$p

if (p_val_Silt4 < 0.05) {
  message_Silt4 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", 
                         round(p_val_Silt4, 3), ").")
} else {
  message_Silt4 <- paste("Levene’s Test indicates that variances are homogeneous (p =", 
                         round(p_val_Silt4, 3), ").")
}

print(message_Silt4)

anova_Silt4 <- aov (Silt ~ Vegetation, data = df4)
tukey_Silt4 <- TukeyHSD(anova_Silt4)
tukey_Silt4_R <- tukey_Silt4$Vegetation[, "p adj"]
names(tukey_Silt4_R) <- rownames(tukey_Silt4$Vegetation)
print (tukey_Silt4_R, 3)
HSD_Silt4 <- HSD.test(anova_Silt4, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (HSD_Silt4)

#Clay 
Shapiro_Clay4 <- shapiro.test(df4$Clay)
if (Shapiro_Clay4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Clay4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Clay4$p.value, 8), ")."))}

df4$Clay_log4 <- log(df4$Clay)

Shapiro_Clay_log4 <- shapiro.test(df4$Clay)
if (Shapiro_Clay_log4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Clay_log4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Clay_log4$p.value, 8), ")."))}

qqnorm(df4$Clay_log4); qqline(df4$Clay_log4)
hist(df4$Clay_log4, main = "Histogram of Clay", xlab = "Clay")

kruskal_Clay4 <- df4 %>% kruskal_test(Clay ~ Vegetation)

print(kruskal_Clay4)

# Confirm p-value before continue
if(kruskal_Clay4$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_Clay4 <- df4 %>% dunn_test(Clay ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result_Clay4)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}


# Comparison of soil textures between treatments 5th Layer #######################
#Sand 
shapiro_sand5 <- shapiro.test(df5$Sand)

if (shapiro_sand5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand5$p.value, 8), ")."))}

qqnorm(df5$Sand); qqline(df5$Sand)
hist(df5$Sand, main = "Histogram of Sand", xlab = "Sand")

df5$Sand_log5 <- log(df5$Sand)

shapiro_sand_log5 <- shapiro.test(df5$Sand_log5)

if (shapiro_sand_log5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_sand_log5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_sand_log5$p.value, 8), ")."))}

qqnorm(df5$Sand_log5); qqline(df5$Sand_log5)
hist(df5$Sand_log5, main = "Histogram of Sand", xlab = "Sand")

kruskal_Sand5 <- df5 %>% kruskal_test(Sand ~ Vegetation)

print(kruskal_Sand5)

# Confirm p-value before continue
if(kruskal_Sand5$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result5 <- df5 %>% dunn_test(Sand ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result5)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

#Silt 
shapiro_Silt5 <- shapiro.test(df5$Silt)

if (shapiro_Silt5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_Silt5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_Silt5$p.value, 8), ")."))}


qqnorm(df5$Silt); qqline(df5$Silt)
hist(df5$Silt, main = "Histogram of Silt", xlab = "Silt")

df5$Vegetation <- as.factor(df5$Vegetation)

levene_Silt5 <- df5 %>% levene_test(Silt ~ Vegetation)

p_val_Silt5 <- levene_Silt5$p

if (p_val_Silt5 < 0.05) {
  message_Silt5 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", 
                         round(p_val_Silt5, 3), ").")
} else {
  message_Silt5 <- paste("Levene’s Test indicates that variances are homogeneous (p =", 
                         round(p_val_Silt5, 3), ").")
}

print(message_Silt5)

anova_Silt5 <- aov (Silt ~ Vegetation, data = df5)
tukey_Silt5 <- TukeyHSD(anova_Silt5)
tukey_Silt5_R <- tukey_Silt5$Vegetation[, "p adj"]
names(tukey_Silt5_R) <- rownames(tukey_Silt5$Vegetation)
print (tukey_Silt5_R, 3)
HSD_Silt5 <- HSD.test(anova_Silt5, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (HSD_Silt5)

#Clay 
Shapiro_Clay5 <- shapiro.test(df5$Clay)
if (Shapiro_Clay5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Clay5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Clay5$p.value, 8), ")."))}

df5$Clay_log5 <- log(df5$Clay)

Shapiro_Clay_log5 <- shapiro.test(df5$Clay)
if (Shapiro_Clay_log5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Clay_log5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Clay_log5$p.value, 8), ")."))}

qqnorm(df5$Clay_log5); qqline(df5$Clay_log5)
hist(df5$Clay_log5, main = "Histogram of Clay", xlab = "Clay")

kruskal_Clay5 <- df5 %>% kruskal_test(Clay ~ Vegetation)

print(kruskal_Clay5)

# Confirm p-value before continue
if(kruskal_Clay5$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_Clay5 <- df5 %>% dunn_test(Clay ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result_Clay5)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}


# Statistical Analysis for Bulk Density (BD) #### 
#Estimation of BD 
df6 <- df6 %>%
  mutate(BD = round(1.56 - 0.0005 * Clay - 0.001 * C + 0.0075 * SB, digits = 2))

# BD for 1st Layer ####################################################
shapiro_BD <- shapiro.test(df1$BD)

if (shapiro_BD$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_BD$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_BD$p.value, 8), ")."))}


qqnorm(df1$BD); qqline(df1$BD)
hist(df1$BD, main = "Histogram of Bulk Density", xlab = "BD")

df1$BD_log <- log(df1$BD)

Shapiro_BD_log <- shapiro.test(df1$BD_log)
if (Shapiro_BD_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_BD_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_BD_log$p.value, 8), ")."))}

qqnorm(df1$BD_log); qqline(df1$BD_log)
hist(df1$BD_log, main = "Histogram of Clay", xlab = "Clay")

Kruskal_BD <- df1 %>% kruskal_test(BD ~ Vegetation)

print(Kruskal_BD)

# Confirm p-value before continue
if(Kruskal_BD$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_BD1 <- df1 %>% dunn_test(BD ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_BD1)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# BD for 2nd Layer ####################################################
shapiro_BD2 <- shapiro.test(df2$BD)

if (shapiro_BD2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_BD2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_BD2$p.value, 8), ")."))}

qqnorm(df2$BD); qqline(df2$BD)
hist(df2$BD, main = "Histogram of Bulk Density", xlab = "BD")

df2$BD_log2 <- log(df2$BD)

Shapiro_BD_log2 <- shapiro.test(df2$BD_log2)
if (Shapiro_BD_log2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_BD_log2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_BD_log2$p.value, 8), ")."))}

qqnorm(df2$BD_log2); qqline(df2$BD_log2)
hist(df2$BD_log2, main = "Histogram of Bulk Density", xlab = "BD")

Kruskal_BD2 <- df2 %>% kruskal_test(BD ~ Vegetation)

print(Kruskal_BD2)

# Confirm p-value before continue
if(Kruskal_BD2$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_BD2 <- df2 %>% dunn_test(BD ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_BD2)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# BD for 3rd layer ####################################################
shapiro_BD3 <- shapiro.test(df3$BD)

if (shapiro_BD3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_BD3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_BD3$p.value, 8), ")."))}

qqnorm(df3$BD); qqline(df3$BD)
hist(df3$BD, main = "Histogram of Bulk Density", xlab = "BD")

df3$BD_log3 <- log(df3$BD)

Shapiro_BD_log3 <- shapiro.test(df3$BD_log3)
if (Shapiro_BD_log3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_BD_log3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_BD_log3$p.value, 8), ")."))}

qqnorm(df3$BD_log3); qqline(df3$BD_log3)
hist(df3$BD_log3, main = "Histogram of Bulk Density", xlab = "BD")

Kruskal_BD3 <- df3 %>% kruskal_test(BD ~ Vegetation)

print(Kruskal_BD3)

# Confirm p-value before continue
if(Kruskal_BD3$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_BD3 <- df3 %>% dunn_test(BD ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_BD3)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# BD for 4th layer ####################################################
shapiro_BD4 <- shapiro.test(df4$BD)

if (shapiro_BD4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_BD4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_BD4$p.value, 8), ")."))}

qqnorm(df4$BD); qqline(df4$BD)
hist(df4$BD, main = "Histogram of Bulk Density", xlab = "BD")

df4$BD_log4 <- log(df4$BD)

Shapiro_BD_log4 <- shapiro.test(df4$BD_log4)
if (Shapiro_BD_log4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_BD_log4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_BD_log4$p.value, 8), ")."))}

qqnorm(df4$BD_log4); qqline(df4$BD_log4)
hist(df4$BD_log4, main = "Histogram of Bulk Density", xlab = "BD")

Kruskal_BD4 <- df4 %>% kruskal_test(BD ~ Vegetation)

print(Kruskal_BD4)

# Confirm p-value before continue
if(Kruskal_BD4$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_BD4 <- df4 %>% dunn_test(BD ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_BD4)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# BD for 5th layer  ####################################################
shapiro_BD5 <- shapiro.test(df5$BD)

if (shapiro_BD5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_BD5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_BD5$p.value, 8), ")."))}

qqnorm(df5$BD); qqline(df5$BD)
hist(df5$BD, main = "Histogram of Bulk Density", xlab = "BD")

df5$BD_log5 <- log(df5$BD)

Shapiro_BD_log5 <- shapiro.test(df5$BD_log5)
if (Shapiro_BD_log5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_BD_log5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_BD_log5$p.value, 8), ")."))}

qqnorm(df5$BD_log5); qqline(df5$BD_log5)
hist(df5$BD_log5, main = "Histogram of Bulk Density", xlab = "BD")

Kruskal_BD5 <- df5 %>% kruskal_test(BD ~ Vegetation)

print(Kruskal_BD5)

# Confirm p-value before continue
if(Kruskal_BD5$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_BD5 <- df5 %>% dunn_test(BD ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_BD5)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}


# Statistical analysis of Carbon concentration ####
#1st Layer 
shapiro_C <- shapiro.test(df1$C)

if (shapiro_C$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C$p.value, 8), ")."))}

df1$C_log1 <- log(df1$C)

Shapiro_C_log1 <- shapiro.test(df1$C_log1)
if (Shapiro_C_log1$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_C_log1$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_C_log1$p.value, 8), ")."))}

Kruskal_C1 <- df1 %>% kruskal_test(C ~ Vegetation)

print(Kruskal_C1)

# Confirm p-value before continue
if(Kruskal_C1$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_C1 <- df1 %>% dunn_test(C ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_C1)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

#2nd Layer
shapiro_C2 <- shapiro.test(df2$C)

if (shapiro_C2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C2$p.value, 8), ")."))}

df2$C_log2 <- log(df2$C)

Shapiro_C_log2 <- shapiro.test(df2$C_log2)
if (Shapiro_C_log2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_C_log2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_C_log2$p.value, 8), ")."))}

Kruskal_C2 <- df2 %>% kruskal_test(C_log2 ~ Vegetation)

print(Kruskal_C2)

# Confirm p-value before continue
if(Kruskal_C2$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_C2 <- df2 %>% dunn_test(C ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_C2)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

#3rd Layer
shapiro_C3 <- shapiro.test(df3$C)

if (shapiro_C3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C3$p.value, 8), ")."))}

#Convert Vegetation to a factor
df3$Vegetation <- as.factor(df3$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_C3 <- df3%>%levene_test(C ~ Vegetation)

# Extract p-value
p_val_C3 <- levene_C3$p

# Create interpretation phrase
if (p_val_C3 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_C3, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_C3, 3), ").")
}

# Print result
print(message)

Kruskal_C3 <- df3 %>% kruskal_test(C ~ Vegetation)

print(Kruskal_C3)

# Confirm p-value before continue
if(Kruskal_C3$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result_C3 <- df3 %>% dunn_test(C ~ Vegetation, p.adjust.method = "bonferroni")
  # Result of post-hoc test
  print(posthoc_result_C3)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

#4th Layer
shapiro_C4 <- shapiro.test(df4$C)

if (shapiro_C4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C4$p.value, 8), ")."))}

#Convert Vegetation to a factor
df4$Vegetation <- as.factor(df4$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_C4 <- df4%>%levene_test(C ~ Vegetation)

# Extract p-value
p_val_C4 <- levene_C4$p

# Create interpretation phrase
if (p_val_C4 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_C4, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_C4, 3), ").")
}

# Print result
print(message)

anova_C4 <- aov(C ~ Vegetation, data = df4)
tukey_C4 <- TukeyHSD(anova_C4)
tukey_C4_R <- tukey_C4$Vegetation[, "p adj"]
names(tukey_C4_R) <- rownames(tukey_C4$Vegetation)
print (tukey_C4_R, 3)
hsd_C4 <- HSD.test(anova_C4, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_C4)

#5th Layer 
shapiro_C5 <- shapiro.test(df5$C)

if (shapiro_C5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C5$p.value, 8), ")."))}

df5$Vegetation <- as.factor(df5$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_C5 <- df5%>%levene_test(C ~ Vegetation)

# Extract p-value
p_val_C5 <- levene_C5$p

# Create interpretation phrase
if (p_val_C5 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_C5, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_C5, 3), ").")
}

print(message)

anova_C5 <- aov(C ~ Vegetation, data = df5)
tukey_C5 <- TukeyHSD(anova_C5)
tukey_C5_R <- tukey_C5$Vegetation[, "p adj"]
names(tukey_C5_R) <- rownames(tukey_C5$Vegetation)
print (tukey_C5_R, 3)
hsd_C5 <- HSD.test(anova_C5, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_C5)

# Comparison of SOC Stocks inside treatments, considering Depth Intervals.  #########
### ----- Clean R's workspaace
rm(list=ls())

df6 <- read_excel("Soil_DB.xlsx", sheet = 7)

#Semideciduous Seasonal Forest (SSF)
DF_SSF <- df6 %>% filter(Vegetation == "SSF")

#Transform Depth Interval in a Factor with 5 levels
DF_SSF$Depth <- factor(DF_SSF$Depth)
# Defining new name for variables 
new_labels_1 <- c("20cm", "40cm", "60cm", "80cm", "100cm")
# Defining order for the levels 
order_levels_1 <- c("P20", "P40", "P60", "P80", "P100")
# Converting 
DF_SSF$Depth <- factor(DF_SSF$Depth, levels = order_levels_1, labels = new_labels_1)

shapiro_SSF <- shapiro.test(DF_SSF$Cstock)

# Print the result
if (shapiro_SSF$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_SSF$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_SSF$p.value, 8), ")."))}

DF_SSF$Cstock_log <- log(DF_SSF$Cstock)

Shapiro_Cstock_SSF_log <- shapiro.test(DF_SSF$Cstock_log)
if (Shapiro_Cstock_SSF_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_Cstock_SSF_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_Cstock_SSF_log$p.value, 8), ")."))}

levene_SSF_log <- DF_SSF%>%levene_test(Cstock_log ~ Depth)

# Extract p-value
p_val_SSF_log <- levene_SSF_log$p

# Create interpretation phrase
if (p_val_SSF_log < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_SSF_log, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_SSF_log, 3), ").")
}

anova_Depth_SSF <- aov(Cstock_log ~ Depth, data = DF_SSF)
tukey_SSF <- TukeyHSD(anova_Depth_SSF)
tukey_SSF_R <- tukey_SSF$Depth[, "p adj"]
names(tukey_SSF_R) <- rownames(tukey_SSF$Depth)
print (tukey_SSF_R, 3)
HSD_SSF_DEPTH <- HSD.test(anova_Depth_SSF, "Depth", alpha = 0.05, group=TRUE, console=FALSE)
print (HSD_SSF_DEPTH)


#Densely Wooded Savanna 
DF_DWS <- df6 %>% filter(Vegetation == "DWS")

DF_DWS$Depth <- factor(DF_DWS$Depth)
new_labels_1 <- c("20cm", "40cm", "60cm", "80cm", "100cm")
order_levels_1 <- c("P20", "P40", "P60", "P80", "P100")
DF_DWS$Depth <- factor(DF_DWS$Depth, levels = order_levels_1, labels = new_labels_1)

shapiro_DWS <- shapiro.test(DF_DWS$Cstock)

if (shapiro_DWS$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DWS$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DWS$p.value, 8), ")."))}

DF_DWS$Cstock_log_DWS <- log(DF_DWS$Cstock)

shapiro_DWS_log <- shapiro.test(DF_DWS$Cstock_log_DWS)

if (shapiro_DWS_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DWS_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DWS_log$p.value, 8), ")."))}

kruskal_DWS <- kruskal(DF_DWS$Cstock, DF_DWS$Depth, alpha = 0.05, p.adj=c("bonferroni"), group = TRUE, main = NULL, console = FALSE)
print (kruskal_DWS)

#Deforested Area (DA)
DF_DA <- df6 %>% filter(Vegetation == "DA")
DF_DA$Depth <- factor(DF_DA$Depth)
new_labels_1 <- c("20cm", "40cm", "60cm", "80cm", "100cm")
order_levels_1 <- c("P20", "P40", "P60", "P80", "P100")
DF_DA$Depth <- factor(DF_DA$Depth, levels = order_levels_1, labels = new_labels_1)

shapiro_DA <- shapiro.test(DF_DA$Cstock)

if (shapiro_DA$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =",
              round(shapiro_DA$p.value, 8), ").")) 
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DA$p.value, 8), ")."))}

DF_DA$Cstock_log_DA <- log(DF_DA$Cstock)

shapiro_DA_log <- shapiro.test(DF_DA$Cstock_log_DA)

if (shapiro_DA_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DA_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DA_log$p.value, 8), ")."))}

# Levene's Test for Homogeneity of Variance
levene_DA_log <- DF_DA%>%levene_test(Cstock_log_DA ~ Depth)

# Extract p-value
p_val_DA_log <- levene_DA_log$p

# Create interpretation phrase
if (p_val_DA_log < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_DA_log, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_DA_log, 3), ").")
}

print(message)

anova_DA_log <- aov(Cstock_log_DA ~ Depth, data = DF_DA)
tukey_DA <- TukeyHSD(anova_DA_log)
tukey_DA_log_R <- tukey_DA$Depth[, "p adj"]
names(tukey_DA_log_R) <- rownames(tukey_DA$Depth)
print (tukey_DA_log_R, 3)
HSD_DA_LOG <- HSD.test(anova_DA_log, "Depth", alpha = 0.05, group=TRUE, console=FALSE)
print (HSD_DA_LOG)
