# Path and object review

The adapted sources replace former server roots with `project_dir()` and absolute
resource references with configuration lookups. Relative paths are deliberately
preserved because aliases may represent different historical runs.

- `3.Venn` versus `3.Venn_1`, `4.machine_learning` versus `_1`,
  `6.亚型分析` versus `_1`, `8.Phenotype` versus `_1`, and Mutation versus Mutation_1
  are inconsistent across files. Select the intended upstream result explicitly;
  do not copy different-run files indiscriminately to satisfy a path.
- Module 04 reads `../0.data` before its later `setwd('4.machine_learning/')`.
  Its early working directory assumptions are inconsistent with later ones.
- Module 01 runs three DEG analyses consecutively. The original helper determines
  whether it changes the working directory and how separate outputs are named.
  Without that helper, expected DEG_PC9/DEG_H1975/DEG_HCC827 outputs are uncertain.
- Module 09 changes directories more than once and references both Mutation and
  Mutation_1. Review each resolved path before enabling it.
- Several downstream modules assume that `colors`, loaded packages, or expression
  objects already exist. Fresh-process execution can expose these dependencies.
- Some functions write files or change the current directory as side effects.
  Missing helpers must be reviewed before relying on launcher working directories.
- Module 11 `library()` does not load analysis packages; explicit requirements
  should be recovered from the original successful run.
- Legacy preprocessing refers to intermediate files not created by all preceding
  lines, and no explicit complete 90-LUAD/5-normal sample selection is shown.

`metadata/input_references.csv` is a static inventory of literal read/load paths,
not a verified file dependency graph. Computed paths, downloads, source helper
dependencies and manually created outputs require additional inspection.
