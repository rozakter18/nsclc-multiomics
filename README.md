# NSCLC multi-omics integration and knowledge mining

Research pipeline for TCGA lung adenocarcinoma (LUAD) and lung squamous cell
carcinoma (LUSC). It integrates RNA expression, DNA methylation, and
copy-number data with MOFA+, evaluates stage associations, mines biological
evidence, and transfers RNA-weight scores to an independent GEO cohort.

The analysis is exploratory and is not a clinical prediction tool. The
research question, cohort definition, and label rules are described in
[docs/objective.md](docs/objective.md); data retrieval and reproducibility
details are in [docs/provenance.md](docs/provenance.md).

## Pipeline

Run commands from the repository root and in order:

| Step | Command | Purpose |
|---|---|---|
| 1 | `make inventory` | Summarize available TCGA samples and stage by omics layer |
| 2 | `make download` | Download selected Xena inputs and update SHA-256 manifest |
| 3 | `make preprocess` | Define the cohort, select features, and run QC |
| 4 | `make integrate` | Fit MOFA+ models and annotate factors |
| 5 | `make model` | Evaluate stage prediction with patient- and site-based CV |
| 6 | `make knowledge` | Rank factor genes and mine pathway, literature, and curated evidence |
| 7 | `make validate` | Validate factor RNA scores in GEO GSE37745 |

`make all` runs the full sequence. `make verify` checks downloaded files
against `data/raw_manifest.tsv`. These analyses include network downloads and
can take substantial time; the MOFA+, modeling, knowledge-mining, and GEO
steps also depend on outputs from earlier steps.

To run the pipeline without Make, execute the scripts from the repository
root in the same order:

```powershell
python scripts/01_sample_inventory.py
python scripts/02_download.py
python scripts/02_download.py --verify
python scripts/03_preprocess.py
python scripts/04_integrate_mofa.py
python scripts/05_stage_model.py
python scripts/06_knowledge_mining.py
python scripts/07_validate_geo.py --gse GSE37745
```

For available options, run each script with `--help` (the inventory script
has no options). `--reuse` in preprocessing reads existing processed layers;
it does not remove the need for the raw cohort and clinical inputs used to
define the sample set.

## Setup

```powershell
conda env create -f environment.yml
conda activate nsclc-mo
```

The environment targets Python 3.11. `requirements.txt` lists the runtime
packages and minimum versions; `environment.yml` installs the same
dependencies through conda and pip.

## Data and outputs

`data/sources.tsv` records source datasets and access terms.
`data/raw_manifest.tsv` records the URLs, file sizes, download timestamps,
and SHA-256 digests for the Xena files. Raw datasets are downloaded to
`data/raw/` and are excluded from Git; processed Parquet matrices are
recreated by step 3 and are also excluded. Small sample extracts and
tabular metadata remain available under `data/processed/`.

The pipeline writes its reports, tables, figures, model files, and cache to
`results/`. Knowledge-mining API caches are reproducibility aids, not source
data, and are excluded from Git. Dynamic external resources can change
between runs. `results/knowledge/knowledge_sources.tsv` records requests
made during a run, but cached responses do not add provenance entries; a
partial-refresh run can replace the file with an incomplete list, while a
cache-only run leaves an existing file unchanged. See
[`docs/provenance.md`](docs/provenance.md) for cache behavior and source-term
limitations.

Current saved outputs describe an 818-patient three-layer cohort
(451 LUAD, 367 LUSC), with 565 patients in the early-versus-late supervised
analysis (414 early, 151 late). The saved stage-model report gives a mean
patient-CV AUROC of 0.603 for MOFA factors and 0.609 for RNA alone. These
results are exploratory and should be interpreted with the limitations in
the reports, including the transductive MOFA fit and external-validation
sample size.

## Attribution and reuse

When reporting or redistributing outputs, acknowledge TCGA and the UCSC Xena
hub and identify the cohorts/datasets used; cite GEO series GSE37745, platform
GPL570, and Botling et al. (2013; DOI:
[10.1158/1078-0432.CCR-12-1139](https://doi.org/10.1158/1078-0432.CCR-12-1139)).
Also acknowledge Reactome, MSigDB Hallmark 2020 and Enrichr, Open Targets,
STRING, HGNC, and PubMed/NCBI E-utilities as applicable. Resource-specific
guidance and unresolved terms are documented in
[`docs/provenance.md`](docs/provenance.md).

Provider documentation identifies Reactome pathway data, Open Targets
Platform data, and HGNC data as CC0. Redistribution terms for the Xena-hosted
TCGA sources, GEO/GPL570, the Hallmark 2020 library served by Enrichr, and
STRING outputs still need confirmation.

Aggregate outputs are shared subject to the applicable terms of their source
data and resources; public access does not itself establish permission to
redistribute source data. No legal clearance is implied.

## Repository layout

```text
data/
  sources.tsv                 source and access manifest
  raw_manifest.tsv            downloaded-file checksums
  raw/                        downloaded inputs (not committed)
  processed/                  derived metadata and representative extracts
docs/
  objective.md                research question and analysis definitions
  provenance.md                data, environment, and reproduction notes
scripts/                       numbered pipeline steps
results/                       reports, result tables, figures, and models
environment.yml                conda environment
requirements.txt               Python package minimums
Makefile                       pipeline commands
```

## Scope and limitations

The main integration uses RNA, methylation, and copy-number views. RPPA is
excluded because its available cohort is smaller; mutation data are not
included in this Xena-based pipeline. Stage II and missing stage labels are
excluded from the binary stage-prediction task. The GEO transfer evaluates
only the RNA component of MOFA factors. See the analysis reports for
method-specific limitations.
