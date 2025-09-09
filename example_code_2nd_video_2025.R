#######################################################################################
# Code example for the YouTube video Unmixing science 2. The FingerPro package UPDATED!
#######################################################################################

# GitHub repository: https://github.com/eead-csic-eesa/fingerPro
# Citing FingerPro and its tools: https://github.com/eead-csic-eesa/fingerPro#citing-fingerpro-and-its-tools
# ResearchGate: https://www.researchgate.net/profile/Ivan-Lizaga
#
# Latest development repository (may contain experimental features): https://github.com/ivanlizaga/fingerPro_devel
# GitHub personal repository: https://github.com/ivanlizaga
#
# The creation of this code was supported by the Research Foundation - Flanders, mandate 12V8622N at 
# Ghent University/Department of Green Chemistry and Technology/ISOFYS group
################################################################################

# To ensure you're working with the last version of the package install it either from GitHub 
# or download the .zip file from GitHub repository and install it in your computer

# From CRAN (version 2.0 || 01/09/2025)
install.packages("fingerPro")

# Download the fingerPro_1.3.zip file from GitHub on you computer
setwd("C:/your/file/directory")
install.packages('fingerPro_2.0.zip', repos = NULL)

# From GitHub (version 2.0)
devtools::install_github("ivanlizaga/fingerPro_devel", ref = "master", force = T)


library(fingerPro)
?fingerPro

# Example stored in the help file
################################################################################
# Created by Ivan Lizaga 09/09/2025

# If you want to use your own data,
# I recommend following the format of the datasets included in the package
# setwd("the directory that contains your dataset")
# data <- read.table('your dataset.csv', header = T, sep = ',')

# Example of the data included in the fingerPro package
# Load the dataset called "catchment" 
# "catchment": this dataset has been selected from a Mediterranean catchment
# to test the multiple functions inside the package. 
# It contains high-quality radionuclides and geochemistry data, four sources and one mixture.
# AG (cropland), PI and PI1 (Two types of Pine forest), and SS (subsoil)


# Load the 'catchment' dataset to test the graphical functions
data <- catchment_2025

# Display boxplots and a correlation graph of the loaded dataset
box_plot(data, tracers = 1:6, ncol = 3, colors = c("#663300", "#339900", "darkgreen", "#666666", "purple"))
correlation_plot(data, columns = 1:6, mixtures = TRUE, colors = c("#663300", "#339900", "darkgreen", "#666666", "purple"))

# Check how the selected tracers discriminate between the potential sources
LDA_plot(data[, c(1:10)], P3D = FALSE, text = FALSE, colors = c("#663300", "#339900", "darkgreen", "#666666"))
LDA_plot(data[, c(1:10)], P3D = FALSE, text = TRUE, colors = c("#663300", "#339900", "darkgreen", "#666666")) #adding point information
LDA_plot(data[, c(1:10)], P3D = TRUE, text = FALSE, interactive = TRUE, colors = c("#663300", "#339900", "darkgreen", "#666666")) #3D LDA

# Displaying a PCA graph 
PCA_plot(data, colors = c("#663300", "#339900", "darkgreen", "#666666", "purple"))

################################################################################
rm(list = ls())
dev.off()

################################################################################
# Compute the CI method and plot the triangles
################################################################################
# The following database was extracted from Raigani et al. (2019): https://doi.org/10.1016/j.ejrh.2019.100613
sources.file <- system.file("extdata", "Raigani_2025.csv", package = "fingerPro")
data <- read.csv(sources.file)
rm(sources.file)
data <- read.csv("G:/My Drive/Congress and workshops/0_Fingerprinting_Courses/4_Plymouth_2025/Raigani_2025.csv")

# Run the CI method
  results_CI <- individual_tracer_analysis(data, iter = 4000, seed = 1234567L)
  
  output_CI <- CI(results_CI)

# Save Ternary diagram to PDF
  ternary_diagram(results_CI, tracers = c(1:8), 
                  rows = 2, cols = 4, solution = c(0.05, 0.9, 0.05))

################################################################################
# Compute the CR method
################################################################################

crgeo <- CR(data, debates = 1000, seed = 1234567L)# 2000 it's more recommended
print(crgeo)

################################################################################
# Compute the CTS method
################################################################################
library(dplyr)
# compute pairs/triplets or quartets depending on your number of sources
# Be aware that for four sources, instead of pairs of tracers you'll have triplets
pgeo <- CTS_seeds(data, iter = 1000, seed=1234567)# 2000 it's more recommended
head(pgeo)

#Explore those pairs
#########
# Pair 1
#########
sol <- print(pgeo[pgeo$tracers == "P Sr",])
output_CTS <- CTS_error(data, solution = sol)
CI_CR_CTS_summary <- merge(output_CTS, crgeo, by = "tracer")
CI_CR_CTS_summary <- merge(CI_CR_CTS_summary, output_CI, by = "tracer")
head(CI_CR_CTS_summary)
# Filter the summary data to select only the most robust tracers.
# The criteria are a low CTS error (< 0.025) and a high CR score (> 80).
sel_tracers_pair_1 <- 
  CI_CR_CTS_summary[CI_CR_CTS_summary$CTS_err < 0.025 & CI_CR_CTS_summary$CR_score > 80, ]
print(sel_tracers_pair_1)

#########
# Pair 2
#########
sol <- print(pgeo[pgeo$tracers == "Ba Sr",])
output_CTS <- CTS_error(data, solution = sol)
CI_CR_CTS_summary <- merge(output_CTS, crgeo, by = "tracer")
CI_CR_CTS_summary <- merge(CI_CR_CTS_summary, output_CI, by = "tracer")
head(CI_CR_CTS_summary)
# Filter the summary data to select only the most robust tracers.
# The criteria are a low CTS error (< 0.025) and a high CR score (> 80).
sel_tracers_pair_2 <- 
  CI_CR_CTS_summary[CI_CR_CTS_summary$CTS_err < 0.025 & CI_CR_CTS_summary$CR_score > 80, ]
print(sel_tracers_pair_2)

sol <- pgeo[pgeo$id=="Ba Sr",]
ctsgeo <- cts_3s(source=sgeo, mixture=mgeo, sol=c(sol$w1, sol$w2, sol$w3))
ctsgeo <- ctsgeo %>% right_join(crgeo, by=c("tracer"))
ctsgeo <- ctsgeo[ctsgeo$err<0.025 & ctsgeo$score>80,]
head(ctsgeo)


################################################################################
# Unmix using the Consistent solutions
################################################################################
# First we select the tracers suggested by the CTS and CR together with their SD (D*) for the two first pairs
data_sel_pair_1 <- data %>% dplyr::select(ID, samples, all_of(paste0("mean_", sel_tracers_pair_1$tracer)),
                                               all_of(paste0("sd_", sel_tracers_pair_1$tracer)), n)
# Transform all columns (except the second) to numeric just in case
data_sel_pair_1[, -2] <- lapply(data_sel_pair_1[, -2], as.numeric)

rs_pair_1 <- unmix(data_sel_pair_1, iter = "short", seed = 1234567L)

p_pair_1 <- plot_results(rs_pair_1, violin = F, bounds = c(-0.2, 1.2), colors = c("#996300", "#8e44ad", "#3498db"))



data_sel_pair_2 <- data %>% dplyr::select(ID, samples, all_of(paste0("mean_", sel_tracers_pair_2$tracer)),
                                          all_of(paste0("sd_", sel_tracers_pair_2$tracer)), n)
# Transform all columns (except the second) to numeric just in case
data_sel_pair_2[, -2] <- lapply(data_sel_pair_2[, -2], as.numeric)

rs_pair_2 <- unmix(data_sel_pair_2, iter = "short", seed = 1234567L)

p_pair_2 <- plot_results(rs_pair_2, violin = F, bounds = c(-0.2, 1.2), colors = c("#996300", "#8e44ad", "#3498db"))

library(gridExtra)
grid.arrange(plot_results(rs_pair_1, y_high = 1, violin = F, colors = c("#996300", "#8e44ad", "#3498db")), plot_results(rs_pair_2, y_high = 1, violin = F, colors = c("#996300", "#8e44ad", "#3498db")),
             plot_results(rs_pair_1, violin = T, colors = c("#996300", "#8e44ad", "#3498db")), plot_results(rs_pair_2, violin = T, colors = c("#996300", "#8e44ad", "#3498db")), ncol=2, nrow =2)

# In the previous example we saw two conservative, consensual and consistent solutions that despite being different had
# similar results. This is not always the case as can be seen in the CTS paper. Thus, you can find totally different
# results for your mixture, mostly depending on the nature of your tracers, your sources, transport processes..... 
# and now that you know of their existence you have the opportunity to discuss about why!!!
# How to choose between them?
# 1) Not everything is black or white and there is not unique or define method for that, 
#     display them in your research and highlight this characteristic, let's move towards a less black box unmixing. 
# 2) Use all previous tools combined with the CI, CR and the triangles together with your expert knowledge, 
#     the unmixing is only the icing on the cake after you understand your data.
# 3) New methods and improvements are being constantly uploaded to the repository, so in case of doubts do not hesitate
#     to contact the developers || lizaga.ivan10@gmail.com