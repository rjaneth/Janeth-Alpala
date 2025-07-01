# Configurar CRAN
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# install and load packages only if they are missing
install_and_load <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  }
  invisible(lapply(packages, library, character.only = TRUE))
}
#
#  packages
required_packages <- c(
  "ggplot2", "reshape2", "knitr", "pandoc", "gridExtra", 
  "gtools", "stats4", "rmutil", "scales", "tidyr", "invgamma", 
  "tidyverse", "RColorBrewer", "ggsci", "carData", "ggpubr",  "patchwork", "dplyr", 
  "kableExtra", "ggthemes", "latex2exp", "e1071", "viridis", "nortest", "bookdown","terra", "sf", "pROC", "purrr"
)

# Install and load only missing packages
install_and_load(required_packages)

# 
# source("./Code/entropy_gamma_sar.R")
# source("./Code/entropy_gI0.R")
# source("./Code/gamma_sar_sample.R")
# source("./Code/gi0_sample.R")
# source("./Code/al_omari_1_estimator.R")
# source("./Code/bootstrap_al_omari_1_estimator.R")
# source("./Code/renyi_entropy_estimator_v1.R")
# source("./Code/bootstrap_renyi_entropy_estimator_v1.R")
# source("./Code/entropy_renyi_gamma_sar.R")
# 
# 
# source("./Code/vasicek_estimator.R")
# source("./Code/van_es_estimator.R")
# source("./Code/noughabi_arghami_estimator.R")
# source("./Code/correa_estimator.R")
# source("./Code/ebrahimi_estimator.R")
# source("./Code/bootstrap_correa_estimator.R")
# 
# source("./Code/tsallis_entropy_gamma_sar.R")
# source("./Code/tsallis_estimator.R")
# source("./Code/bootstrap_tsallis_entropy.R")
# source("./Code/tsallis_estimator_optimized.R")
# source("./Code/bootstrap_tsallis_entropy_optimized.R")
# 
# source("./Code/renyi_estimator.R")
# source("./Code/bootstrap_renyi_estimator.R")
# source("./Code/renyi_entropy_optimized.R")
# source("./Code/bootstrap_renyi_estimator_Op.R")
# 
# source("./Code/shannon_alomari_estimator.R")
# source("./Code/bootstrap_shannon_alomari.R")
# source("./Code/shannon_entropy_gamma_sar.R")
# 
# source("./Code/functions_sample_bias_mse.R")
 source("./Code/functions_sample_bias_mse_1.R")
# 
# theme_set(theme_minimal() +
#             theme(text = element_text(family = "serif"),
#                   legend.position = "bottom"))


#source("./Code/read_ENVI_images.R")
