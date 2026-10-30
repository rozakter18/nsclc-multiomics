# External validation report: GSE37745 (step 6)

- run date: 2026-10-18; platform GPL570; 196 samples
- fields used: stage = 'tumor stage', histology = 'histology', time = 'days to determined death status', event = 'dead', age = 'age', sex = 'gender'

## Pre-specified hypotheses (directions taken from TCGA)

Expected directions were specified from TCGA before validation testing. One-sided p-values test the expected direction and are nominal.

- **Factor2**: score lower in late (III-IV) than early (I) stage (TCGA AUC 0.423); survival HR per SD < 1 (TCGA HR 0.95, FDR 0.5)
- **Factor4**: score lower in late (III-IV) than early (I) stage (TCGA AUC 0.425); survival HR per SD < 1 (TCGA HR 0.84, FDR 0.011)
- **Factor5**: score higher in late (III-IV) than early (I) stage (TCGA AUC 0.565); survival HR per SD > 1 (TCGA HR 1.17, FDR 0.06)

## Transfer check in TCGA

- genes with RNA weights: 5000; measured in GSE37745: 4212 (84%)
- correlation in TCGA between the MOFA factor and the transferred score (all genes / only genes measured in the validation cohort):

| factor   |   r_all_genes |   r_overlap_genes |
|:---------|--------------:|------------------:|
| Factor2  |         0.935 |             0.937 |
| Factor4  |         0.884 |             0.886 |
| Factor5  |         0.946 |             0.946 |

Values near 1 mean the transferred score reproduces the factor; a low value means the validation result for that factor is not interpretable.

## Validation cohort

- LUAD/LUSC samples used: 172 ({'LUAD': 106, 'LUSC': 66})
- stage: {1: 110, 2: 34, 3: 24, 4: 4} -> early (I) 110, late (III-IV) 28
- survival: 172 with follow-up, 129 deaths

## Results

| factor   |   stage_auc |   n_late |   n_early |   stage_p_two_sided |   stage_p_one_sided | stage_concordant   |   stage_auc_LUAD |   stage_auc_LUSC |   os_hr_per_sd |   os_p_two_sided |   os_p_one_sided | os_concordant   |   os_n |   os_events |   os_hr_adj_stage |   os_p_adj_stage |
|:---------|------------:|---------:|----------:|--------------------:|--------------------:|:-------------------|-----------------:|-----------------:|---------------:|-----------------:|-----------------:|:----------------|-------:|------------:|------------------:|-----------------:|
| Factor2  |       0.427 |       28 |       110 |              0.233  |               0.116 | True               |            0.48  |            0.334 |          0.883 |           0.156  |           0.0778 | True            |    172 |         129 |             0.914 |           0.324  |
| Factor4  |       0.651 |       28 |       110 |              0.0141 |               0.993 | False              |            0.634 |            0.693 |          0.898 |           0.25   |           0.125  | True            |    172 |         129 |             0.861 |           0.114  |
| Factor5  |       0.501 |       28 |       110 |              0.989  |               0.495 | True               |            0.436 |            0.605 |          1.18  |           0.0628 |           0.0314 | True            |    172 |         129 |             1.17  |           0.0859 |

Stage AUC compares late (III-IV) with early (I); stage II is excluded, as in TCGA. Survival HRs are per SD of the factor score from Cox models stratified by histology. `os_hr_adj_stage` and `os_p_adj_stage` are from a model additionally adjusted for stage; this is covariate adjustment, not multiple-testing correction. One-sided p-values are derived from the two-sided Cox p-value using the expected TCGA direction.

- No stage association met the nominal directional threshold (one-sided p < 0.05).
- Exploratory, nominal directional survival result: Factor5 had HR 1.179, two-sided Cox p = 0.0628, and one-sided p = 0.0314 in the expected direction. This is not confirmed replication: the two-sided p-value is above 0.05, and the TCGA association had FDR 0.06.
- Validation p-values across the three factors and outcomes (stage, unadjusted survival, and stage-adjusted survival) were not multiplicity-adjusted. No multiplicity-adjusted validation p-value is reported.

Figures: results\validation\GSE37745/km_<factor>.png, results\validation\GSE37745/stage_<factor>.png

## Limitations

- Only the RNA part of each factor can be transferred; methylation and CNV contributions are not tested here.
- Microarray vs RNA-seq: within-cohort z-scoring removes scale differences, but platform-specific probe behaviour remains.
- One probe per gene (highest mean expression); genes absent from the platform are dropped (see transfer check).