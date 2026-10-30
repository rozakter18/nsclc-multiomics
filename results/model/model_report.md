# Stage prediction report (step 5a)

- patients: 565 (151 late, 414 early); LUAD 335, LUSC 230
- tissue source sites: 64
- CV: 5-fold x 10 repeats; identical splits for every model
- stable MOFA factors used: Factor1, Factor2, Factor3, Factor4, Factor5, Factor6, Factor7, Factor8, Factor9, Factor10

## Cross-validated performance

AUROC 0.5 = chance. Balanced accuracy 0.5 = chance (corrects for the 2.7:1 class imbalance). Mean over repeats; sd = spread across repeats.

| model                   | cv      |   n_features |   auroc_mean |   auroc_sd |   bal_acc_mean |   auroc_LUAD |   auroc_LUSC |
|:------------------------|:--------|-------------:|-------------:|-----------:|---------------:|-------------:|-------------:|
| clinical (age+sex)      | patient |            2 |        0.546 |      0.012 |          0.534 |        0.534 |        0.567 |
| clinical (age+sex)      | site    |            2 |        0.538 |      0.021 |          0.53  |        0.53  |        0.552 |
| MOFA factors            | patient |           10 |        0.603 |      0.011 |          0.575 |        0.633 |        0.556 |
| MOFA factors            | site    |           10 |        0.603 |      0.01  |          0.578 |        0.636 |        0.553 |
| MOFA factors + clinical | patient |           12 |        0.596 |      0.013 |          0.567 |        0.622 |        0.557 |
| MOFA factors + clinical | site    |           12 |        0.594 |      0.007 |          0.575 |        0.623 |        0.552 |
| RNA only                | patient |         5000 |        0.609 |      0.015 |          0.574 |        0.625 |        0.583 |
| RNA only                | site    |         5000 |        0.605 |      0.008 |          0.57  |        0.623 |        0.573 |
| methylation only        | patient |         5000 |        0.558 |      0.013 |          0.553 |        0.582 |        0.527 |
| methylation only        | site    |         5000 |        0.555 |      0.01  |          0.545 |        0.567 |        0.537 |
| CNV only                | patient |          495 |        0.528 |      0.019 |          0.514 |        0.522 |        0.53  |
| CNV only                | site    |          495 |        0.518 |      0.016 |          0.51  |        0.518 |        0.519 |
| all omics concatenated  | patient |        10495 |        0.59  |      0.022 |          0.561 |        0.619 |        0.549 |
| all omics concatenated  | site    |        10495 |        0.585 |      0.011 |          0.558 |        0.603 |        0.557 |

## Paired comparison with the clinical baseline (patient CV)

Same splits for every model, so per-repeat differences are paired. `wins` = repeats in which the model beat the baseline.

| model                   |   delta_auroc_mean |   delta_sd | wins   |
|:------------------------|-------------------:|-----------:|:-------|
| MOFA factors            |              0.057 |      0.014 | 10/10  |
| MOFA factors + clinical |              0.051 |      0.013 | 10/10  |
| RNA only                |              0.064 |      0.017 | 10/10  |
| methylation only        |              0.012 |      0.017 | 6/10   |
| CNV only                |             -0.018 |      0.022 | 1/10   |
| all omics concatenated  |              0.045 |      0.025 | 9/10   |

## Label-permutation test (patient CV, repeat 1)

Labels shuffled 200 times; p = share of shuffles that reach the observed AUROC. Low-dimensional models only (fast).

| model                   |   observed_auroc |   null_mean |   null_95th |      p |
|:------------------------|-----------------:|------------:|------------:|-------:|
| clinical (age+sex)      |            0.544 |       0.494 |       0.549 | 0.0697 |
| MOFA factors            |            0.6   |       0.496 |       0.559 | 0.01   |
| MOFA factors + clinical |            0.593 |       0.501 |       0.551 | 0.01   |

## Stage-linked factors adjusted for age, sex and histology

Logistic regression of late stage on each factor (per SD). These factors were SELECTED for stage association in step 4 on the same patients, so p-values are optimistic; the question here is only whether each association survives adjustment.

| factor   | model                          |   n |   OR_per_SD |   CI_low |   CI_high |      p |
|:---------|:-------------------------------|----:|------------:|---------:|----------:|-------:|
| Factor2  | unadjusted                     | 565 |        0.77 |     0.63 |      0.92 | 0.0054 |
| Factor2  | adjusted for age+sex+histology | 545 |        0.79 |     0.65 |      0.96 | 0.019  |
| Factor4  | unadjusted                     | 565 |        0.76 |     0.63 |      0.92 | 0.0055 |
| Factor4  | adjusted for age+sex+histology | 545 |        0.76 |     0.62 |      0.92 | 0.0054 |
| Factor5  | unadjusted                     | 565 |        1.19 |     0.99 |      1.44 | 0.067  |
| Factor5  | adjusted for age+sex+histology | 545 |        1.16 |     0.96 |      1.41 | 0.13   |
| Factor6  | unadjusted                     | 565 |        0.82 |     0.68 |      0.99 | 0.039  |
| Factor6  | adjusted for age+sex+histology | 545 |        0.85 |     0.7  |      1.02 | 0.088  |

Figure: results/model/roc.png