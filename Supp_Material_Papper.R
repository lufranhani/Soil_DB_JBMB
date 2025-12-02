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

# pH for 1st Layer ####################################################
shapiro_pH <- shapiro.test(df1$pH)

if (shapiro_pH$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_pH$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_pH$p.value, 8), ")."))}


qqnorm(df1$pH); qqline(df1$pH)
hist(df1$pH, main = "Histogram of pH", xlab = "pH")

df1$pH_log <- log(df1$pH)

Shapiro_pH_log <- shapiro.test(df1$pH_log)

if (Shapiro_pH_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_pH_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_pH_log$p.value, 8), ")."))}


#Convert Vegetation to a factor
df1$Vegetation <- as.factor(df1$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_pH1 <- df1 %>% 
  levene_test(pH_log ~ Vegetation)

# Extract p-value
pval_pH_log <- levene_pH1$p

# Create interpretation phrase
if (pval_pH_log < 0.05) {
  message_pH1 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(pval_pH_log, 3), ").")
} else {
  message_pH1 <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(pval_pH_log, 3), ").")
}

# Print result
print(message_pH1)

anova_pH1 <- aov(pH_log ~ Vegetation, data = df1)

# Tukey's HSD
tukey_pH <- TukeyHSD(anova_pH1)

# Extract adjusted p-values
tukey_pH_R <- tukey_pH$Vegetation[, "p adj"]
names(tukey_pH_R) <- rownames(tukey_pH$Vegetation)
print (tukey_pH_R, 3)
hsd_pH <- HSD.test(anova_pH1, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_pH)


# pH for 2nd Layer ####################################################
shapiro_pH2 <- shapiro.test(df2$pH)

if (shapiro_pH2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_pH2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_pH2$p.value, 8), ")."))}


qqnorm(df2$pH); qqline(df2$pH)
hist(df2$pH, main = "Histogram of pH2", xlab = "pH2")


df2$pH_log2 <- log(df2$pH)

Shapiro_pH_log2 <- shapiro.test(df2$pH_log2)

if (Shapiro_pH_log2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_pH_log2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_pH_log2$p.value, 8), ")."))}


#Convert Vegetation to a factor
df2$Vegetation <- as.factor(df2$Vegetation)

kruskal_pH2 <- df2 %>% kruskal_test(pH ~ Vegetation)

print(kruskal_pH2)

# Confirm p-value before continue
if(kruskal_pH2$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df2 %>% dunn_test(pH ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# pH for 3rd Layer ####################################################
shapiro_pH3 <- shapiro.test(df3$pH)

if (shapiro_pH3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_pH3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_pH3$p.value, 8), ")."))}


qqnorm(df3$pH); qqline(df3$pH)
hist(df3$pH, main = "Histogram of pH3", xlab = "pH3")


df3$pH_log3 <- log(df3$pH)

Shapiro_pH_log3 <- shapiro.test(df3$pH_log3)

if (Shapiro_pH_log3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_pH_log3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_pH_log3$p.value, 8), ")."))}


#Convert Vegetation to a factor
df3$Vegetation <- as.factor(df3$Vegetation)

kruskal_pH3 <- df3 %>% kruskal_test(pH ~ Vegetation)

print(kruskal_pH3)

# Confirm p-value before continue
if(kruskal_pH3$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df3 %>% dunn_test(pH ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# pH for 4th Layer ####################################################
shapiro_pH4 <- shapiro.test(df4$pH)

if (shapiro_pH4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_pH4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_pH4$p.value, 8), ")."))}


qqnorm(df4$pH); qqline(df4$pH)
hist(df4$pH, main = "Histogram of pH4", xlab = "pH4")


df4$pH_log4 <- log(df4$pH)

Shapiro_pH_log4 <- shapiro.test(df4$pH_log4)

if (Shapiro_pH_log4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_pH_log4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_pH_log4$p.value, 8), ")."))}


#Convert Vegetation to a factor
df4$Vegetation <- as.factor(df4$Vegetation)

kruskal_pH4 <- df4 %>% kruskal_test(pH ~ Vegetation)

print(kruskal_pH4)

# Confirm p-value before continue
if(kruskal_pH4$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df4 %>% dunn_test(pH ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# pH for 5th Layer ####################################################
shapiro_pH5 <- shapiro.test(df4$pH)

if (shapiro_pH5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_pH5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_pH5$p.value, 8), ")."))}


qqnorm(df5$pH); qqline(df5$pH)
hist(df5$pH, main = "Histogram of pH5", xlab = "pH5")


df5$pH_log5 <- log(df5$pH)

Shapiro_pH_log5 <- shapiro.test(df5$pH_log5)

if (Shapiro_pH_log5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_pH_log5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_pH_log5$p.value, 8), ")."))}

qqnorm(df5$pH_log5); qqline(df5$pH_log5)
hist(df5$pH_log5, main = "Histogram of pH5", xlab = "pH5")

#Convert Vegetation to a factor
df5$Vegetation <- as.factor(df5$Vegetation)

kruskal_pH5 <- df5 %>% kruskal_test(pH ~ Vegetation)

print(kruskal_pH5)

# Confirm p-value before continue
if(kruskal_pH5$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df5 %>% dunn_test(pH ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}


# Comparison of pH  inside treatments, considering Depth Intervals.  #########
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

shapiro_SSF <- shapiro.test(DF_SSF$pH)

# Print the result
if (shapiro_SSF$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_SSF$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_SSF$p.value, 8), ")."))}

DF_SSF$pH_log <- log(DF_SSF$pH)

Shapiro_pH_SSF_log <- shapiro.test(DF_SSF$pH_log)
if (Shapiro_pH_SSF_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(Shapiro_pH_SSF_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(Shapiro_pH_SSF_log$p.value, 8), ")."))}

levene_SSF_log <- DF_SSF%>%levene_test(pH_log ~ Depth)

# Extract p-value
p_val_SSF_log <- levene_SSF_log$p

# Create interpretation phrase
if (p_val_SSF_log < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_SSF_log, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_SSF_log, 3), ").")
}

print (message)


anova_Depth_SSF <- aov(pH_log ~ Depth, data = DF_SSF)
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

shapiro_DWS <- shapiro.test(DF_DWS$pH)

if (shapiro_DWS$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DWS$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DWS$p.value, 8), ")."))}

DF_DWS$pH_log_DWS <- log(DF_DWS$pH)

shapiro_DWS_log <- shapiro.test(DF_DWS$pH_log_DWS)

if (shapiro_DWS_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DWS_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DWS_log$p.value, 8), ")."))}

kruskal_DWS <- kruskal(DF_DWS$pH, DF_DWS$Depth, alpha = 0.05, p.adj=c("bonferroni"), group = TRUE, main = NULL, console = FALSE)
print (kruskal_DWS)



#Deforested Area (DA)
DF_DA <- df6 %>% filter(Vegetation == "DA")
DF_DA$Depth <- factor(DF_DA$Depth)
new_labels_1 <- c("20cm", "40cm", "60cm", "80cm", "100cm")
order_levels_1 <- c("P20", "P40", "P60", "P80", "P100")
DF_DA$Depth <- factor(DF_DA$Depth, levels = order_levels_1, labels = new_labels_1)

shapiro_DA <- shapiro.test(DF_DA$pH)

if (shapiro_DA$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =",
              round(shapiro_DA$p.value, 8), ").")) 
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DA$p.value, 8), ")."))}

DF_DA$pH_log_DA <- log(DF_DA$pH)

shapiro_DA_log <- shapiro.test(DF_DA$pH_log_DA)

if (shapiro_DA_log$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DA_log$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DA_log$p.value, 8), ")."))}

kruskal_DA <- kruskal(DF_DA$pH, DF_DA$Depth, alpha = 0.05, p.adj=c("bonferroni"), group = TRUE, main = NULL, console = FALSE)
print (kruskal_DA)

# 13C for 1st Layer ####################################################
shapiro_C13 <- shapiro.test(df1$C13)

if (shapiro_C13$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C13$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C13$p.value, 8), ")."))}


qqnorm(df1$C13); qqline(df1$C13)
hist(df1$C13, main = "Histogram of C13", xlab = "C13")


#Convert Vegetation to a factor
df1$Vegetation <- as.factor(df1$Vegetation)

kruskal_C13 <- df1 %>% kruskal_test(C13 ~ Vegetation)

print(kruskal_C13)

# Confirm p-value before continue
if(kruskal_C13$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df2 %>% dunn_test(C13 ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# 13C for 2ND Layer ####################################################
shapiro_C13_2 <- shapiro.test(df2$C13)

if (shapiro_C13_2$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C13_2$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C13_2$p.value, 8), ")."))}


qqnorm(df2$C13); qqline(df2$C13)
hist(df2$C13, main = "Histogram of C13", xlab = "C13")


#Convert Vegetation to a factor
df2$Vegetation <- as.factor(df2$Vegetation)

kruskal_C13_2 <- df2 %>% kruskal_test(C13 ~ Vegetation)

print(kruskal_C13_2)

# Confirm p-value before continue
if(kruskal_C13_2$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df2 %>% dunn_test(C13 ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}

# 13C for 3rd Layer ####################################################
shapiro_C13_3 <- shapiro.test(df3$C13)

if (shapiro_C13_3$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C13_3$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C13_3$p.value, 8), ")."))}


qqnorm(df3$C13); qqline(df3$C13)
hist(df3$C13, main = "Histogram of C13", xlab = "C13")


#Convert Vegetation to a factor
df3$Vegetation <- as.factor(df3$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_C13 <- df3 %>% 
  levene_test(C13 ~ Vegetation)


# Extract p-value
pval_C13 <- levene_C13$p

# Create interpretation phrase
if (pval_C13 < 0.05) {
  message_C13 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(pval_C13, 3), ").")
} else {
  message_C13 <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(pval_C13, 3), ").")
}

# Print result
print(message_pH1)

anova_C13 <- aov(C13 ~ Vegetation, data = df3)

# Tukey's HSD
tukey_C13 <- TukeyHSD(anova_C13)

# Extract adjusted p-values
tukey_C13_R <- tukey_C13$Vegetation[, "p adj"]
names(tukey_C13_R) <- rownames(tukey_C13$Vegetation)
print (tukey_C13_R, 3)
hsd_C13 <- HSD.test(anova_C13, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_C13)

# 13C for 4TH Layer ####################################################
shapiro_C13_4 <- shapiro.test(df4$C13)

if (shapiro_C13_4$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C13_4$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C13_4$p.value, 8), ")."))}


qqnorm(df4$C13); qqline(df4$C13)
hist(df4$C13, main = "Histogram of C13", xlab = "C13")

#Convert Vegetation to a factor
df4$Vegetation <- as.factor(df4$Vegetation)

# Levene's Test for Homogeneity of Variance
levene_C13_4 <- df4 %>% 
  levene_test(C13 ~ Vegetation)


# Extract p-value
pval_C13_4 <- levene_C13$p

# Create interpretation phrase
if (pval_C13_4 < 0.05) {
  message_C13_4 <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(pval_C13_4, 3), ").")
} else {
  message_C13_4 <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(pval_C13_4, 3), ").")
}

# Print result
print(message_C13_4)

anova_C13_4 <- aov(C13 ~ Vegetation, data = df4)

# Tukey's HSD
tukey_C13_4 <- TukeyHSD(anova_C13_4)

# Extract adjusted p-values
tukey_C13_R_4 <- tukey_C13_4$Vegetation[, "p adj"]
names(tukey_C13_R_4) <- rownames(tukey_C13_4$Vegetation)
print (tukey_C13_R_4, 3)
hsd_C13_4 <- HSD.test(anova_C13_4, "Vegetation", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_C13_4)


# 13C for 5TH Layer ####################################################
shapiro_C13_5 <- shapiro.test(df5$C13)

if (shapiro_C13_5$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_C13_5$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_C13_5$p.value, 8), ")."))}


qqnorm(df5$C13); qqline(df5$C13)
hist(df5$C13, main = "Histogram of C13", xlab = "C13")

#Convert Vegetation to a factor
df5$Vegetation <- as.factor(df5$Vegetation)

kruskal_C13_5 <- df5 %>% kruskal_test(C13 ~ Vegetation)

print(kruskal_C13_5)

# Confirm p-value before continue
if(kruskal_C13_5$p < 0.05) {
  # If necessary, run Dunn Test
  posthoc_result <- df5 %>% dunn_test(C13 ~ Vegetation, p.adjust.method = "none")
  # Result of post-hoc test
  print(posthoc_result)
} else {
  print("Kruskal-Wallis test wasnt significative. Dunn Post-Hoc test isnt necessary.")
}


# Comparison of C13 inside treatments, considering Depth Intervals.  #########
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

shapiro_SSF <- shapiro.test(DF_SSF$C13)

# Print the result
if (shapiro_SSF$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_SSF$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_SSF$p.value, 8), ")."))}

levene_SSF_C13 <- DF_SSF%>%levene_test(C13 ~ Depth)

# Extract p-value
p_val_SSF_C13 <- levene_SSF_C13$p

# Create interpretation phrase
if (p_val_SSF_C13 < 0.05) {
  message <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(p_val_SSF_C13, 3), ").")
} else {
  message <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(p_val_SSF_C13, 3), ").")
}

print (message)


anova_Depth_SSF <- aov(C13 ~ Depth, data = DF_SSF)
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

shapiro_DWS <- shapiro.test(DF_DWS$C13)

if (shapiro_DWS$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =", 
              round(shapiro_DWS$p.value, 8), ")."))
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DWS$p.value, 8), ")."))}


kruskal_DWS <- kruskal(DF_DWS$C13, DF_DWS$Depth, alpha = 0.05, p.adj=c("none"), group = TRUE, main = NULL, console = FALSE)
print (kruskal_DWS)

#Deforested Area (DA)
DF_DA <- df6 %>% filter(Vegetation == "DA")
DF_DA$Depth <- factor(DF_DA$Depth)
new_labels_1 <- c("20cm", "40cm", "60cm", "80cm", "100cm")
order_levels_1 <- c("P20", "P40", "P60", "P80", "P100")
DF_DA$Depth <- factor(DF_DA$Depth, levels = order_levels_1, labels = new_labels_1)

shapiro_DA <- shapiro.test(DF_DA$C13)

if (shapiro_DA$p.value < 0.05) {
  print(paste("The variable does NOT FOLLOW a normal distribution (p =",
              round(shapiro_DA$p.value, 8), ").")) 
} else {
  print(paste("The variable FOLLOWS a normal distribution (p =", 
              round(shapiro_DA$p.value, 8), ")."))}

#Convert Vegetation to a factor
df6$Depth <- as.factor(df6$Depth)

# Levene's Test for Homogeneity of Variance
levene_C13 <- DF_DA %>% 
  levene_test(C13 ~ Depth)

# Extract p-value
pval_C13_DA<- levene_C13$p

# Create interpretation phrase
if (pval_C13_DA < 0.05) {
  message_C13_DA <- paste("Levene’s Test indicates that variances are NOT homogeneous (p =", round(pval_C13_DA, 3), ").")
} else {
  message_C13_DA <- paste("Levene’s Test indicates that variances are homogeneous (p =", round(pval_C13_DA, 3), ").")
}

# Print result
print(message_C13_DA)

anova_C13_DA <- aov(C13 ~ Depth, data = DF_DA)

# Tukey's HSD
tukey_C13_DA <- TukeyHSD(anova_C13_DA)

# Extract adjusted p-values
tukey_C13_DA_R <- tukey_C13_DA$Depth[, "p adj"]
names(tukey_C13_DA_R) <- rownames(tukey_C13_DA$Depth)
print (tukey_C13_DA_R, 3)
hsd_C13_DA <- HSD.test(anova_C13_DA, "Depth", alpha = 0.05, group=TRUE, console=FALSE)
print (hsd_C13_DA)
