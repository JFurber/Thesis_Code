### This script is used to compute the univariate statistical tests

library(dplyr)
library(PMCMRplus)

# Normality Test
shapiro.test(Noise_GLMM_data$L)

Noise_GLMM_data$Site <- as.factor(Noise_GLMM_data$Site)
Noise_GLMM_data$Month <- as.factor(Noise_GLMM_data$Month)
Noise_GLMM_data$Capture_Year <- as.factor(Noise_GLMM_data$Capture_Year)

# Filter the data based on site
Woodchester <- Noise_GLMM_data %>% filter(Site == "Woodchester")
C2 <- Noise_GLMM_data %>% filter(Site == "C2")
C4 <- Noise_GLMM_data %>% filter(Site == "C4")
F1 <- Noise_GLMM_data %>% filter(Site == "F1")
F2 <- Noise_GLMM_data %>% filter(Site == "F2")
Ireland <- Noise_GLMM_data %>% filter(Site == "Ireland")

# Mann-whitney U test - sex differences
wilcox.test(L ~ Sex, data = Woodchester)
wilcox.test(L ~ Sex, data = C2)
wilcox.test(L ~ Sex, data = C4)
wilcox.test(L ~ Sex, data = F1)
wilcox.test(L ~ Sex, data = F2)
wilcox.test(L ~ Sex, data = Ireland)

# Kruskal Wallis test - month differences
kruskal.test(L~Month, data = Woodchester)
dscfAllPairsTest(L~ Month, data = Woodchester)

kruskal.test(L~Month, data = C2)
# dscfAllPairsTest(L~ Month, data = C2) #- not needed

kruskal.test(L~Month, data = C4)
# dscfAllPairsTest(L~ Month, data = C4) #- not needed

kruskal.test(L~Month, data = F1)
dscfAllPairsTest(L~ Month, data = F1)

kruskal.test(L~Month, data = F2)
dscfAllPairsTest(L~ Month, data = F2)

kruskal.test(L~Month, data = Ireland)
dscfAllPairsTest(L~ Month, data = Ireland)

# Kruskal Wallis test - year differences
kruskal.test(L~Capture_Year, data = Woodchester)
dscfAllPairsTest(L~ Capture_Year, data = Woodchester)

kruskal.test(L~Capture_Year, data = C2)
dscfAllPairsTest(L~ Capture_Year, data = C2)

kruskal.test(L~Capture_Year, data = C4)
dscfAllPairsTest(L~ Capture_Year, data = C4)

kruskal.test(L~Capture_Year, data = F1)
dscfAllPairsTest(L~ Capture_Year, data = F1)

kruskal.test(L~Capture_Year, data = F2)
# dscfAllPairsTest(L~ Capture_Year, data = F2) #- not needed

kruskal.test(L~Capture_Year, data = Ireland)
# dscfAllPairsTest(L~ Capture_Year, data = Ireland) # - not needed


#Kruskal Wallis test - site differences
kruskal.test(L~Site, data = Noise_GLMM_data)
dscfAllPairsTest(L~ Site, data = Noise_GLMM_data)
