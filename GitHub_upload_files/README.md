# MMD + WXD / LUAD and osimertinib-resistance analysis code

This repository organizes the 16 supplied R source files for the MMD + WXD study.
It is a **source-code release candidate, not a verified end-to-end reproduction**.
Original files are preserved byte-for-byte in `archive/original/`; path-adapted
versions are under `scripts/`. No analyses were rerun and no study results changed.

## Repository structure

```text
archive/original/       Original supplied scripts (unchanged)
scripts/               Numbered modules and execution/environment utilities
R/project.R            Configurable project, helper and external-resource paths
R/vendor/              Original helper functions must be supplied here
config/                Example config; local config is ignored by Git
metadata/              Module mapping, packages, missing helpers, input references
data/                  Input schemas and data-placement instructions
docs/                  Methods audit, path review and upload instructions
results/logs/          Environment/run records produced on the user's R machine
workspace/             Local input/output layout; large files excluded from Git
```

## Status and prerequisites

Read `docs/METHODS_AUDIT.md`, `docs/PATH_REVIEW.md` and
`metadata/helper_dependencies.csv` before execution. The supplied modules call
custom functions that were not supplied. No replacement functions, fabricated
package versions or synthetic historical `sessionInfo()` are included.
The preprocessing script retains its historical transformations for audit;
its scale/normalization questions have not been resolved by packaging.

R is required. Inspect `metadata/packages.csv` for packages directly referenced
in active code; dependencies inside missing helpers are unknown. No minimum R
version or historical package versions have been verified. Optional installation:

```sh
Rscript scripts/install_dependencies.R --install
Rscript scripts/check_environment.R
```

The installation helper uses current repositories, not historical versions.
`KEGG.db` and `oncoPredict` require manual dependency/source review.

## Configure and execute a reviewed module

1. Copy `config/config.example.R` to `config/config.local.R`.
2. Place study input files under `workspace/0.data/` or point `project_root` to a
   separate copy of the original project. Use a copy: legacy modules write files.
3. Supply original helpers in `R/vendor/` and configure external resources.
4. Resolve the module's path and method issues, then enable the explicit review
   acknowledgement in local config. This acknowledgement does not verify results.
5. Run from the repository root in a fresh R process:

```sh
Rscript scripts/run_module.R 02
```

Do not run every file in one session. Standalone modules may still depend on
unrecorded objects such as `colors`; these are flagged for review. Module 14 is a
placeholder template and the launcher intentionally refuses to execute it.
No new global seed was imposed: adding seeds now would not recover the old run.

## Suggested dependency order (not an automatic pipeline)

00 preprocessing and 15 compound targets -> 01 resistance genes + 02 disease DEGs
-> 03 intersection/enrichment -> 04 feature selection -> 05 diagnostic model and
06 subtypes -> 07 correlation, 08 drug sensitivity, 09 TMB, 10 checkpoints,
11 TIDE, 12 networks, 13 pan-cancer. Module 11 requires externally generated TIDE
output. Module 14 is not verified GSE31210 tissue-classification code.
The filename “single-gene GSEA” was misleading: supplied module 07 performs
correlation analysis; it is named accordingly in the adapted copy.

## Reproducibility and results

Current run utilities save `sessionInfo()` and package versions. Historical
versions remain unknown. A lockfile should be created only from a verified
analysis environment; no invented `renv.lock` is supplied. Package installation,
function defaults, preprocessing scales, contrast direction and source datasets
must be checked before using this release to support manuscript claims.

No expression matrices, clinical tables, mutation files, API tokens or result
datasets are bundled. Public-data provenance and redistribution conditions must
be recorded when inputs are added. Code availability does not by itself establish
analytical validity or reproducibility. See `README_zh-CN.md` for Chinese guidance.

## License and citation

No license is assigned on the authors' behalf. See `LICENSE_NOTE.md`. Add an
author-approved license and the verified publication/repository citation before
public release. No GitHub repository has been created or uploaded by this task.
