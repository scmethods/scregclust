# Load data to be used inside package
human_tfs <- read.csv("datasets/humanTFs.txt", header = FALSE)[, 1]
human_tfs_v2 <- read.csv("datasets/humanTFs_v2.txt", header = FALSE)[, 1]
human_tfs_v3 <- read.csv("datasets/humanTFs_v3.txt", header = FALSE)[, 1]

human_kinases <- read.csv("datasets/humanKinases.txt", header = FALSE)[, 1]

human_regulators <- read.csv("datasets/humanRegulators.txt", header = FALSE)[, 1]

# Create R/sysdata.rda with those datasets
usethis::use_data(
  human_tfs,
  human_tfs_v2,
  human_tfs_v3,
  human_kinases,
  human_regulators,
  internal = TRUE,
  overwrite = TRUE
)

# Use them with e.g. scregclust:::human_tfs
# Note that there are three `:`