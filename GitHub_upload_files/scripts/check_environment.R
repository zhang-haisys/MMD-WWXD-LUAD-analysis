# Read-only inventory and R syntax check. No analysis or downloads.
root <- normalizePath(getwd(), winslash = "/")
if (!file.exists("metadata/modules.csv")) stop("Start in the repository root.")
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
p <- read.csv("metadata/packages.csv", colClasses = "character")
p$installed <- vapply(p$package, requireNamespace, logical(1), quietly = TRUE)
p$current_version <- vapply(p$package, function(x)
  if (requireNamespace(x, quietly = TRUE)) as.character(packageVersion(x)) else NA_character_, "")
write.csv(p, "results/logs/current_packages.csv", row.names = FALSE)
m <- read.csv("metadata/modules.csv", colClasses = "character")
m$parse_status <- vapply(m$script, function(x) tryCatch({
  parse(file = x, encoding = "UTF-8"); "OK"
}, error = function(e) conditionMessage(e)), "")
write.csv(m[c("id", "script", "parse_status")], "results/logs/parse_status.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "results/logs/current_sessionInfo.txt")
print(m[c("id", "parse_status")], row.names = FALSE)
message("Current environment only; this does not recover historical versions or validate scientific results.")
