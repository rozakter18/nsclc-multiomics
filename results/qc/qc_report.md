# QC report - preprocessing (step 3)

## 1. Cohort

- rna: 1017 primary-tumour patients
- methylation: 830 primary-tumour patients
- cnv: 1017 primary-tumour patients
- clinical: 1090 patients
- **all three layers + clinical: 818 patients**

Overall survival available for 800 patients (311 deaths)

| histology   |   early |   intermediate |   late |   missing |
|:------------|--------:|---------------:|-------:|----------:|
| LUAD        |     244 |            111 |     91 |         5 |
| LUSC        |     170 |            134 |     60 |         3 |

Supervised cohort (early vs late): 565 patients (414 early, 151 late)

## rna: reused data/processed/rna.parquet (5000 features; --reuse)

## methylation: reused data/processed/methylation.parquet (5000 features; --reuse)

## cnv: reused data/processed/cnv.parquet (495 features; --reuse)

## 5. QC

PCA on each processed layer (standardised features, NA -> feature mean for PCA only). eta2 = share of a PC's variance explained by a grouping.

Batch check: `eta2_site_resid` is the share of a PC's variance explained by tissue source site AFTER removing histology means. `site_null_mean` is the same quantity with site labels shuffled within histology (200 permutations); many small sites explain some variance by chance. `site_excess` = observed - null is the batch effect beyond chance; `site_p` is the permutation p-value.

### rna

| PC   | var_explained   |   eta2_histology |   eta2_stage |   eta2_site_resid |   site_null_mean |   site_excess |   site_p |
|:-----|:----------------|-----------------:|-------------:|------------------:|-----------------:|--------------:|---------:|
| PC1  | 15.4%           |            0.599 |        0.023 |             0.136 |            0.081 |         0.055 |    0.005 |
| PC2  | 7.2%            |            0.173 |        0.01  |             0.214 |            0.08  |         0.134 |    0.005 |
| PC3  | 4.8%            |            0.001 |        0.024 |             0.089 |            0.077 |         0.012 |    0.204 |
| PC4  | 3.7%            |            0.013 |        0.002 |             0.072 |            0.08  |        -0.008 |    0.706 |
| PC5  | 2.8%            |            0.006 |        0.005 |             0.125 |            0.081 |         0.045 |    0.005 |

### methylation

| PC   | var_explained   |   eta2_histology |   eta2_stage |   eta2_site_resid |   site_null_mean |   site_excess |   site_p |
|:-----|:----------------|-----------------:|-------------:|------------------:|-----------------:|--------------:|---------:|
| PC1  | 19.2%           |            0.211 |        0.001 |             0.089 |            0.077 |         0.011 |    0.189 |
| PC2  | 7.8%            |            0.124 |        0.012 |             0.111 |            0.081 |         0.031 |    0.035 |
| PC3  | 7.5%            |            0.252 |        0.005 |             0.082 |            0.08  |         0.002 |    0.453 |
| PC4  | 4.0%            |            0.117 |        0.007 |             0.102 |            0.084 |         0.017 |    0.129 |
| PC5  | 2.0%            |            0.069 |        0.004 |             0.09  |            0.077 |         0.012 |    0.154 |

### cnv

| PC   | var_explained   |   eta2_histology |   eta2_stage |   eta2_site_resid |   site_null_mean |   site_excess |   site_p |
|:-----|:----------------|-----------------:|-------------:|------------------:|-----------------:|--------------:|---------:|
| PC1  | 7.3%            |            0.473 |        0.012 |             0.082 |            0.08  |         0.002 |    0.423 |
| PC2  | 5.1%            |            0.058 |        0.004 |             0.066 |            0.079 |        -0.013 |    0.846 |
| PC3  | 4.1%            |            0.023 |        0.003 |             0.077 |            0.078 |        -0     |    0.512 |
| PC4  | 3.9%            |            0     |        0.004 |             0.077 |            0.077 |         0     |    0.512 |
| PC5  | 3.6%            |            0.001 |        0.005 |             0.097 |            0.077 |         0.021 |    0.08  |

### Per-patient checks

- patients with > 20% missing methylation or extreme RNA median (|z| > 4): 0
- flagged patients are reported, not removed; decide after review

## 6. Outputs

- rna: 818 patients x 5000 features -> data/processed/rna.parquet
- methylation: 818 patients x 5000 features -> data/processed/methylation.parquet
- cnv: 818 patients x 495 features -> data/processed/cnv.parquet
- metadata -> data/processed/metadata.tsv

## Parameters

- N_TOP = {'rna': 5000, 'methylation': 5000, 'cnv': 3000}
- RNA_MIN_EXPR = 1.0
- RNA_MIN_FRAC = 0.2
- METH_MAX_NA = 0.1
- CNV_MIN_ALT_FRAC = 0.05
- CNV_MAX_CORR = 0.9
- CNV_CANDIDATES = 8000
- N_PERM = 200