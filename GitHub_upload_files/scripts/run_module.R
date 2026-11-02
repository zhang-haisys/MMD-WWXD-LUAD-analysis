# Launch one reviewed historical module; no automatic end-to-end execution.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/run_module.R <module_id>")
repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "metadata", "modules.csv")))
  stop("Start in the repository root.")
source(file.path(repo_root, "config", "config.example.R"), encoding = "UTF-8")
local_config <- file.path(repo_root, "config", "config.local.R")
if (file.exists(local_config)) source(local_config, encoding = "UTF-8")
source(file.path(repo_root, "R", "project.R"), encoding = "UTF-8")
modules <- read.csv(file.path(repo_root, "metadata", "modules.csv"),
                    colClasses = "character", check.names = FALSE)
module <- modules[modules$id == sprintf("%02d", as.integer(args[1])), , drop = FALSE]
if (nrow(module) != 1L) stop("Unknown module.")
if (module$status == "template_not_study_validation")
  stop("This is a placeholder external-validation template, not an executable study model.")
if (!isTRUE(config$acknowledge_unresolved_legacy_methods))
  stop("Read docs/METHODS_AUDIT.md and docs/PATH_REVIEW.md; resolve relevant items before enabling execution in config/config.local.R.")
script <- file.path(repo_root, module$script)
parse(file = script, encoding = "UTF-8")
deps <- strsplit(module$required_helpers, ";", fixed = TRUE)[[1]]
for (dep in deps[nzchar(deps)]) helper_file(dep)
pkgs <- strsplit(module$packages, ";", fixed = TRUE)[[1]]
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
dir.create(config$project_root, recursive = TRUE, showWarnings = FALSE)
stages <- unique(c(modules$stage, "3.Venn", "4.machine_learning_1",
                  "6.亚型分析_1", "8.Phenotype_1", "9.tmb/Mutation", "9.tmb/Mutation_1"))
for (stage in stages) dir.create(project_dir(stage), recursive = TRUE, showWarnings = FALSE)
logs <- file.path(repo_root, "results", "logs")
dir.create(logs, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
logbase <- file.path(logs, paste0(module$id, "_", stamp))
capture.output(sessionInfo(), file = paste0(logbase, "_session_before.txt"))
writeLines(c(paste("module:", module$id), paste("script:", module$script),
             paste("original_sha256:", module$sha256),
             paste("adapted_md5:", unname(tools::md5sum(script))),
             "This records the CURRENT run, not the historical environment."),
           paste0(logbase, "_provenance.txt"))
old <- getwd()
tryCatch({
  setwd(project_dir(module$stage))
  source(script, local = .GlobalEnv, encoding = "UTF-8")
}, finally = {
  setwd(old)
  capture.output(sessionInfo(), file = paste0(logbase, "_session_after.txt"))
})
