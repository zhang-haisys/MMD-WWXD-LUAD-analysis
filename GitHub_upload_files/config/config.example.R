# Copy to config/config.local.R and configure paths. Never commit private paths.
config <- list(
  project_root = file.path(repo_root, "workspace"),
  helper_dir = file.path(repo_root, "R", "vendor"),
  acknowledge_unresolved_legacy_methods = FALSE,
  resources = list(
    resource_032493440b = NA_character_, # oncoPredict_Training_Data
    resource_057bf6ea30 = NA_character_, # TCGA_list.rda
    resource_4ae41aa186 = NA_character_, # TF-Target-information.txt
    resource_4d1dbc6d3e = NA_character_, # starbase_miRNA_lncRNA_High.csv
    resource_d0c0f93a27 = NA_character_ # miRNA_mRNA.csv
  )
)
