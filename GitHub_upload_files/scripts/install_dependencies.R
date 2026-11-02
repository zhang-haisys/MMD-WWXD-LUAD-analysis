# Optional, explicit installation into the current R environment.
# This is NOT a historical version lock. Review packages.csv first.
if (!file.exists("metadata/packages.csv")) stop("Start in repository root.")
if (!"--install" %in% commandArgs(trailingOnly = TRUE))
  stop("Review metadata/packages.csv, then use: Rscript scripts/install_dependencies.R --install")
p <- read.csv("metadata/packages.csv", colClasses = "character")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
cran <- p$package[p$source == "CRAN_candidate"]
bioc <- p$package[p$source == "Bioconductor"]
if (length(cran)) install.packages(cran, repos = "https://cloud.r-project.org")
if (length(bioc)) BiocManager::install(bioc, update = FALSE, ask = FALSE)
message("Review KEGG.db/oncoPredict separately; they are not automatically installed. Helper-internal dependencies may be missing.")
