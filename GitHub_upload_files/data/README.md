# Data placement and provenance

No study datasets are included. Put local inputs under workspace/0.data, and
configure external databases in config/config.local.R. Large data are excluded.

Known core inputs:

| File | Expected structure |
|---|---|
| exp_combat.csv | Gene rows, sample columns, first column gene identifiers |
| groupinfo2.csv | sample, dataset, group; group values tumor/control |
| PC9.csv, H1975.csv, HCC827.csv | Sample column actually used as gene identifier; remaining columns expression samples; verify source scale |
| PC9_group.csv, H1975_group.csv, HCC827_group.csv | sample and group (Parental/OR) |
| total_clinical.csv | Sample-linked clinical/survival fields; confirm units |
| 药物靶点.txt | Target gene identifiers; verify exact header with source |

The separately inspected processed matrix contains 10,543 genes and 605 samples;
the supplied metadata has TCGA-LUAD 479 tumors/31 controls and GSE11969 90/5.
These counts establish consistency of those files, not correctness of normalization
or completeness of cohort selection. They are not enforced as software constants.

See metadata/input_references.csv and external_resources.csv for further inputs.
Record accession, platform, retrieval date, download query/manifest, processing
scale, sample exclusions and source licensing. TCGA-LUAD/GSE11969/GSE202859/
GSE31210 are named in the manuscript; supplied code also includes GSE54236 and
unverified placeholders. Reconcile these before reproducibility claims.
Do not upload clinical files, database snapshots or training datasets merely to
satisfy missing paths without checking provenance and redistribution conditions.
