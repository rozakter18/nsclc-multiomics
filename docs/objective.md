# Research objective

## The question in one sentence

**Which molecular programmes, shared across several omics layers, separate
early-stage from late-stage NSCLC tumours — and are the genes driving those
programmes supported by known biology?**

## In plain terms

A tumour can be described at several levels at once:

| Layer | What it tells you |
|---|---|
| DNA copy number | Which genes are gained or lost |
| DNA methylation | Which genes are switched off or on epigenetically |
| RNA expression | Which genes are actually active |
| Protein (RPPA) | Which signalling proteins are present / activated; excluded from the main analysis because of cohort size |

Most studies look at one layer. A change observed coherently across multiple
layers can be more interpretable than a change in only one layer.

The pipeline does three things:

1. **Integrate** — use MOFA+ to find latent factors that capture variation
   across the three selected omics layers, with LUAD and LUSC as separate
   groups.
2. **Test** — check whether any factor tracks tumour stage (early vs late),
   using patient-level cross-validation and a simple baseline.
3. **Ground (knowledge mining)** — rank RNA-anchored genes by multi-omics
   support and evaluate pathway enrichment (Reactome and MSigDB Hallmark),
   interactions (STRING), Open Targets associations, and PubMed evidence.
   Evidence tags are *known*, *plausible*, *limited evidence*, or *ambiguous
   symbol*; limited evidence is not a claim of novelty.

## Why stage

- Clinically meaningful: stage drives treatment decisions and prognosis.
- Available for nearly all TCGA LUAD/LUSC tumours, so labels are open and
  well defined.
- Not trivially separable (unlike LUAD vs LUSC), so a result is informative.

## Expected outcome

Stage signal in bulk TCGA data is known to be modest. A weak predictive
result is compatible with the exploratory aim. The contribution is the
integrated factors and evidence-ranked candidate list, not a clinical
classifier.

## Operational label definition

| Label | AJCC pathologic stage |
|---|---|
| Early (0) | Stage I, IA, IB |
| Late (1) | Stage III, IIIA, IIIB, IV |
| Excluded | Stage II (intermediate), missing / discrepant |

The preprocessing script maps stage labels beginning with I to early, labels
III-IV as late, and retains stage II as intermediate for descriptive reports.
The binary prediction task excludes intermediate and missing labels. The
inventory script additionally maps a literal stage `0` to early, while the
preprocessing script does not. No stage-0 label appears in the current
processed metadata, but this parser difference should be resolved before
including any such samples in a future run.

## Decision rule after the sample inventory

| Inventory result | Decision |
|---|---|
| ≥ ~80 late-stage tumours with all layers | Proceed with early vs late |
| Late-stage count thin only because of RPPA | Drop RPPA, keep 3 layers |
| Still too thin | Unsupervised integration only; stage, histology and survival used to *annotate* factors, not as a prediction target |

## Leakage controls

- Splits by patient (first 12 characters of the TCGA barcode), never by sample.
- Primary tumours only (sample type `01`); normals excluded from the stage task.
- Feature selection, imputation, scaling, and regularization tuning are
  fitted inside each CV training fold.
- MOFA+ is fitted without stage labels on the full integration cohort before
  cross-validation; it is therefore transductive, though not label-informed.
- Tissue source site checked as a potential batch confounder.
- Histology (LUAD/LUSC) used as a stratification variable, not a feature.

## Inventory result and decision (17 Oct 2026)

Primary tumours, one per patient, Xena TCGA hub.

| Layer set | Patients | Early | Stage II | Late | Missing |
|---|---|---|---|---|---|
| RNA + methylation + CNV | 818 (LUAD 451, LUSC 367) | 414 | 245 | 151 | 8 |
| + RPPA | 315 | 164 | 81 | 65 | 5 |

- **Layers:** RNA (HiSeqV2), DNA methylation (450K), and CNV (GISTIC2
  thresholded).
- **RPPA excluded** from the main analysis: it removes ~60% of patients and over half of late-stage cases.
- **Mutation excluded:** no mutation dataset in the Xena TCGA hub cohort
  listing; CNV is the only genomics layer in this analysis.
- **Integration cohort:** 818 patients. **Supervised cohort:** 565 (early vs late, stage II and missing excluded).
- **Class imbalance** ~2.7 : 1 → stratified folds, class weights, AUROC and balanced accuracy.
- **LUSC late-stage n = 60** → used in grouped integration, not a separate LUSC classifier.
