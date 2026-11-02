# Record the environment actually used for a new run.
dir.create("results/logs", recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = "results/logs/sessionInfo.txt")
ip <- installed.packages()
write.csv(ip[, intersect(c("Package", "Version", "Built"), colnames(ip)), drop = FALSE],
          "results/logs/installed_packages.csv", row.names = FALSE)
# Create a lock ONLY after a successful, verified run, using its actual environment:
# renv::init(bare = TRUE)
# renv::snapshot(type = "all", prompt = FALSE)
