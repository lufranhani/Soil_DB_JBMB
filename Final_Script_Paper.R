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

# Perform the Shapiro-Wilk test on the variable
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
         axis.title.y = element_text(size = 20)) +
  scale_fill_manual(values = c("SSF" = "skyblue", "DWS" = "lightgreen", "DA" = "coral"))
p1 
p2 = p1 + labs(tag = "a")

# Defining new name for the variables
new_labels  <- c("0 - 20 cm", "20 - 40 cm", "40 - 60 cm", "60 - 80 cm", "80 - 100 cm")
# Defining new order to the levels 
order_levels <- c("P20", "P40", "P60", "P80", "P100")
# Converting Depth to a factor with new names 
df6$Depth <- factor(df6$Depth, levels = order_levels, labels = new_labels)

# PLOT CARBONO POR CAMADAS!!!! INTERVALOS DE 20 EM 20CM 
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
         axis.title.y = element_text(size = 20)) 

p3 
p4 = p3 + labs(tag = 'b')
p4 

# Arrange the plots one above the other
grid.arrange(p2, p4, widths = 1)

# 1st Soil Layer #######################
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

# 2nd Soil Layer #######################
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

# 3rd Soil Layer #######################
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

# 4th Soil Layer #######################
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

# 5th Soil Layer #######################
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

#
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
