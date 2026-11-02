# Shared path configuration. No expression transformations are performed here.
project_dir <- function(...) file.path(config$project_root, ...)
helper_file <- function(name) {
  path <- file.path(config$helper_dir, name)
  if (!file.exists(path)) stop("Missing original helper: ", name,
    ". See metadata/helper_dependencies.csv. Do not substitute an unverified helper.")
  path
}
resource_path <- function(key) {
  value <- config$resources[[key]]
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value))
    stop("Configure external resource: ", key, "; see metadata/external_resources.csv")
  if (!file.exists(value)) stop("Resource does not exist: ", value)
  value
}
