# MOFA+ integration report (step 4)

- patients: 818 (LUAD 451, LUSC 367)
- views: rna 5000, methylation 5000, cnv 495
- settings: up to 15 factors, scale_views, groups centred, factors < 1% R2 dropped, seeds [1, 2, 3]

## Training and robustness

| seed | factors | final ELBO |
|---|---|---|
| 1 | 12 | -9,980,836 |
| 2 (best) | 13 | -9,936,063 |
| 3 | 13 | -9,944,566 |

Best |correlation| of each factor with any factor of the other runs (near 1 = reproducible). Factors below 0.8 in any run are marked unstable and excluded from conclusions:

|          |   max_abs_r_seed1 |   max_abs_r_seed3 | stable   |
|:---------|------------------:|------------------:|:---------|
| Factor1  |             1     |             0.996 | True     |
| Factor2  |             0.99  |             0.935 | True     |
| Factor3  |             0.993 |             0.966 | True     |
| Factor4  |             0.99  |             0.919 | True     |
| Factor5  |             0.994 |             0.987 | True     |
| Factor6  |             0.996 |             0.94  | True     |
| Factor7  |             0.976 |             0.856 | True     |
| Factor8  |             0.994 |             0.994 | True     |
| Factor9  |             0.978 |             0.892 | True     |
| Factor10 |             0.964 |             0.9   | True     |
| Factor11 |             0.935 |             0.652 | False    |
| Factor12 |             0.978 |             0.731 | False    |
| Factor13 |             0.333 |             0.825 | False    |

## Variance explained (%)

Total per view and histology:

| group   |   rna |   methylation |   cnv |
|:--------|------:|--------------:|------:|
| LUAD    |  37.3 |          43.6 |   9.2 |
| LUSC    |  37.8 |          36.8 |   7.9 |

Per factor (columns = view_histology):

| factor   |   r2_rna_LUAD |   r2_methylation_LUAD |   r2_cnv_LUAD |   r2_rna_LUSC |   r2_methylation_LUSC |   r2_cnv_LUSC |
|:---------|--------------:|----------------------:|--------------:|--------------:|----------------------:|--------------:|
| Factor1  |          1.29 |                 20.33 |          0.64 |          1.02 |                 13.88 |          0.54 |
| Factor2  |          7.41 |                  5    |          1.65 |          7.65 |                  4.28 |          1.13 |
| Factor3  |          0.51 |                  5.44 |          0.47 |          2.1  |                  6.05 |          0.85 |
| Factor4  |          2.02 |                  0.16 |          0.38 |          4.97 |                  3.98 |          2.48 |
| Factor5  |          4.24 |                  1.32 |          0.85 |          4.63 |                  0.99 |          0.34 |
| Factor6  |          4.04 |                  1.09 |          0.52 |          5.13 |                  0.64 |          0.14 |
| Factor7  |          0.39 |                  5.06 |          0.16 |          0.71 |                  4.88 |          0.23 |
| Factor8  |          4.81 |                  2.48 |          2.67 |          0.61 |                  0.15 |          0.27 |
| Factor9  |          2.06 |                  1.38 |          0.99 |          1.6  |                  0.88 |          0.71 |
| Factor10 |          2.64 |                  0.16 |          0.23 |          3.2  |                  0.17 |          0.43 |
| Factor11 |          3.08 |                  0.17 |          0.17 |          2.5  |                  0.17 |          0.25 |
| Factor12 |          3.69 |                  0.09 |          0.05 |          2.14 |                  0.03 |          0.04 |
| Factor13 |          1.12 |                  0.95 |          0.39 |          1.6  |                  0.69 |          0.45 |

## Factor annotation

Stage AUC: late vs early (0.5 = no association; <0.5 means lower in late stage). The SIGN of a MOFA factor is arbitrary, so AUC 0.30 is as strong as 0.70 - read the distance from 0.5. Site: excess eta2 over a within-histology permutation null. OS: Cox hazard ratio per SD, stratified by histology. FDR = Benjamini-Hochberg across factors. All associations are descriptive (no labels were used to train MOFA+).

| factor   | stable   |   stage_auc |   stage_auc_LUAD |   stage_auc_LUSC |   stage_fdr |   site_excess |   site_fdr |   os_hr_per_sd |   os_fdr |   age_rho |   sex_auc |
|:---------|:---------|------------:|-----------------:|-----------------:|------------:|--------------:|-----------:|---------------:|---------:|----------:|----------:|
| Factor1  | True     |       0.513 |            0.525 |            0.495 |       0.843 |         0.017 |      0.183 |          1.088 |    0.228 |     0.013 |     0.529 |
| Factor2  | True     |       0.423 |            0.402 |            0.452 |       0.03  |         0.046 |      0.022 |          0.953 |    0.499 |     0.174 |     0.404 |
| Factor3  | True     |       0.499 |            0.513 |            0.478 |       0.973 |         0.004 |      0.474 |          1.036 |    0.566 |    -0.118 |     0.447 |
| Factor4  | True     |       0.425 |            0.39  |            0.443 |       0.03  |         0.01  |      0.474 |          0.84  |    0.011 |     0.054 |     0.504 |
| Factor5  | True     |       0.565 |            0.587 |            0.532 |       0.048 |        -0.007 |      0.667 |          1.166 |    0.06  |     0.006 |     0.522 |
| Factor6  | True     |       0.433 |            0.405 |            0.473 |       0.047 |         0.068 |      0.022 |          0.957 |    0.499 |     0.044 |     0.453 |
| Factor7  | True     |       0.497 |            0.482 |            0.515 |       0.973 |         0.004 |      0.504 |          1.113 |    0.159 |     0.048 |     0.547 |
| Factor8  | True     |       0.519 |            0.525 |            0.491 |       0.788 |         0.006 |      0.474 |          1.089 |    0.228 |    -0.006 |     0.609 |
| Factor9  | True     |       0.485 |            0.479 |            0.491 |       0.837 |        -0.002 |      0.62  |          0.891 |    0.159 |    -0.157 |     0.479 |
| Factor10 | True     |       0.534 |            0.556 |            0.5   |       0.412 |         0.043 |      0.026 |          1.109 |    0.159 |     0.04  |     0.511 |
| Factor11 | False    |       0.426 |            0.419 |            0.436 |       0.03  |         0.177 |      0.022 |          0.924 |    0.228 |    -0.028 |     0.469 |
| Factor12 | False    |       0.509 |            0.509 |            0.512 |       0.869 |         0.038 |      0.026 |          0.992 |    0.892 |    -0.083 |     0.506 |
| Factor13 | False    |       0.556 |            0.557 |            0.564 |       0.087 |        -0.004 |      0.62  |          1.099 |    0.21  |    -0.137 |     0.556 |

- unstable across seeds (not interpreted): Factor11, Factor12, Factor13
- factors with a site effect (FDR < 0.05, excess > 0.05): Factor6, Factor11
- stage-associated, stable, no strong site effect (FDR < 0.05): **Factor2, Factor4, Factor5**
- stage-associated but stable with a site effect (interpret with caution): Factor6
- associated with overall survival (FDR < 0.05, stable): Factor4

Figures: results/mofa/variance_explained.png, results/mofa/factors_vs_stage.png