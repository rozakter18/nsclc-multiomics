# Knowledge-mining report (step 5b)

- factors mined: Factor2, Factor4, Factor5
- top 200 genes per factor for enrichment; top 100 for literature and curated evidence
- candidates are RNA-anchored: a gene must carry the factor in its own expression (scaled RNA weight >= 0.2); methylation and CNV add support (CNV discounted by region size)
- background (universe) for enrichment: genes measured in the RNA layer (the layer candidates are drawn from), not the whole genome
- literature queries include HGNC aliases and previous symbols

- Reactome: 2868 gene sets
- Hallmark: 50 gene sets
- Open Targets: 13288 targets for MONDO_0005233 non-small cell lung carcinoma
- network anchors: top 50 Open Targets NSCLC genes (e.g. EGFR, KRAS, ALK, ERBB2, ROS1, MET, DICER1, RET)

## Factor2

- stage AUC 0.423 (LUAD 0.402, LUSC 0.452); survival HR/SD 0.95 (FDR 0.5)
- variance explained: rna_LUAD 7.4%, methylation_LUAD 5.0%, cnv_LUAD 1.6%, rna_LUSC 7.7%, methylation_LUSC 4.3%, cnv_LUSC 1.1%
- RNA-supported genes: 1280 of 5000 measured

**Enriched gene sets** (FDR < 0.05: 5, of which 0 driven mainly by one genomic locus; top 10 shown)

| library   | term                                           |   overlap |   set_size |   expected |     fdr |   max_locus_share | positional_cluster   |
|:----------|:-----------------------------------------------|----------:|-----------:|-----------:|--------:|------------------:|:---------------------|
| Reactome  | Surfactant metabolism                          |        11 |         17 |       0.68 | 1.5e-09 |              0.11 | False                |
| Reactome  | Diseases associated with surfactant metabolism |         8 |         10 |       0.4  | 5.7e-08 |              0.17 | False                |
| Reactome  | Diseases of metabolism                         |        16 |         90 |       3.6  | 5.8e-05 |              0.07 | False                |
| Reactome  | Toll Like Receptor 4 (TLR4) Cascade            |         6 |         26 |       1.04 | 0.044   |              0.25 | False                |
| Reactome  | Metabolism of proteins                         |        27 |        343 |      13.72 | 0.044   |              0.04 | False                |
| Reactome  | Toll Like Receptor 2 (TLR2) Cascade            |         5 |         20 |       0.8  | 0.062   |              0.33 | False                |
| Reactome  | Toll Like Receptor TLR1:TLR2 Cascade           |         5 |         20 |       0.8  | 0.062   |              0.33 | False                |
| Reactome  | Disease                                        |        27 |        365 |      14.6  | 0.071   |              0.08 | False                |
| Reactome  | Toll-like Receptor Cascades                    |         6 |         32 |       1.28 | 0.076   |              0.25 | False                |
| Reactome  | Gastrulation                                   |         6 |         39 |       1.56 | 0.19    |              0.17 | False                |

**Evidence tags** (top 100 genes): limited evidence 55, known 23, plausible 20, ambiguous symbol 2

**Top 10 genes by multi-omics support**

| gene    | layers              |   data_support | rna_higher_in   |   pubmed_lung |   pubmed_progression |   ot_score | string_known_partners   | evidence         |   gap_score | aliases_searched                                      |
|:--------|:--------------------|---------------:|:----------------|--------------:|---------------------:|-----------:|:------------------------|:-----------------|------------:|:------------------------------------------------------|
| KNDC1   | rna+methylation     |          1.276 | early           |             0 |                    0 |      0     |                         | limited evidence |       0.934 | KIAA1768,bB439H18.3,FLJ25027,Very-KIND,v-KIND,RASGEF2 |
| CFTR    | rna+methylation+cnv |          1.244 | early           |            73 |                   30 |      0.101 |                         | known            |       0.163 | MRP7,ABC35,TNR-CFTR,dJ760C5.1,ABCC7                   |
| SFTA3   | rna+methylation+cnv |          1.221 | early           |            12 |                    8 |      0.078 |                         | plausible        |       0.388 | NANCI,SP-H,PAHRF,SFTPH                                |
| AFF3    | rna+methylation     |          1.139 | early           |             5 |                    5 |      0.436 |                         | plausible        |       0.556 | MLLT2-like,LAF4                                       |
| SCGB3A1 | rna+methylation+cnv |          1.122 | early           |            22 |                   15 |      0.09  |                         | limited evidence |       0.24  | UGRP2,HIN-1,HIN1,LU105,PnSP-2                         |
| PTPRN2  | rna+methylation+cnv |          1.121 | early           |             4 |                    1 |      0.014 |                         | limited evidence |       0.791 | KIAA0387,phogrin,ICAAR,IA-2beta                       |
| HOXB13  | rna+methylation     |          1.118 | late            |            11 |                    7 |      0.092 |                         | limited evidence |       0.408 |                                                       |
| IRX1    | rna+methylation+cnv |          1.099 | early           |             7 |                    7 |      0.06  |                         | limited evidence |       0.398 |                                                       |
| TCF21   | rna+methylation     |          1.095 | early           |            69 |                   34 |      0.08  |                         | known            |       0.061 | POD1,bHLHa23                                          |
| SFTPC   | rna+cnv             |          1.066 | early           |           100 |                   58 |      0.12  |                         | known            |      -0.051 | SP-C,PSP-C,SMDP2,BRICD6,SFTP2                         |

**Top 10 under-studied candidates** (highest gap score: strong data support, little progression literature)

| gene     | layers              |   data_support | rna_higher_in   |   pubmed_lung |   pubmed_progression |   ot_score | string_known_partners   | evidence         |   gap_score | aliases_searched                                      |
|:---------|:--------------------|---------------:|:----------------|--------------:|---------------------:|-----------:|:------------------------|:-----------------|------------:|:------------------------------------------------------|
| KNDC1    | rna+methylation     |          1.276 | early           |             0 |                    0 |      0     |                         | limited evidence |       0.934 | KIAA1768,bB439H18.3,FLJ25027,Very-KIND,v-KIND,RASGEF2 |
| PTPRN2   | rna+methylation+cnv |          1.121 | early           |             4 |                    1 |      0.014 |                         | limited evidence |       0.791 | KIAA0387,phogrin,ICAAR,IA-2beta                       |
| NECAB2   | rna+methylation+cnv |          1.007 | late            |             0 |                    0 |      0     |                         | limited evidence |       0.791 | EFCBP2                                                |
| SORCS2   | rna+methylation     |          1.005 | early           |             1 |                    1 |      0.002 |                         | limited evidence |       0.689 | KIAA1329                                              |
| TMEM145  | rna+methylation     |          0.915 | late            |             0 |                    0 |      0     |                         | limited evidence |       0.689 | FLJ90805                                              |
| NME5     | rna+methylation     |          0.908 | early           |             1 |                    0 |      0.007 |                         | limited evidence |       0.668 | nm23-H5,RSPH23                                        |
| MMP28    | rna+methylation+cnv |          1.019 | early           |             4 |                    2 |      0.007 |                         | limited evidence |       0.633 | MM28,EPILYSIN,MMP-28                                  |
| AFF3     | rna+methylation     |          1.139 | early           |             5 |                    5 |      0.436 |                         | plausible        |       0.556 | MLLT2-like,LAF4                                       |
| PPP1R16B | rna+methylation+cnv |          0.963 | early           |             2 |                    2 |      0.01  |                         | limited evidence |       0.551 | KIAA0823,TIMAP,ANKRD4                                 |
| PARP15   | rna+methylation+cnv |          1.015 | early           |             3 |                    3 |      0.005 |                         | limited evidence |       0.531 | FLJ40597,pART7,ARTD7,BAL3                             |

## Factor4

- stage AUC 0.425 (LUAD 0.390, LUSC 0.443); survival HR/SD 0.84 (FDR 0.011)
- variance explained: rna_LUAD 2.0%, methylation_LUAD 0.2%, cnv_LUAD 0.4%, rna_LUSC 5.0%, methylation_LUSC 4.0%, cnv_LUSC 2.5%
- RNA-supported genes: 777 of 5000 measured

**Enriched gene sets** (FDR < 0.05: 14, of which 8 driven mainly by one genomic locus; top 10 shown)

| library   | term                                         |   overlap |   set_size |   expected |     fdr |   max_locus_share | positional_cluster   |
|:----------|:---------------------------------------------|----------:|-----------:|-----------:|--------:|------------------:|:---------------------|
| Reactome  | Biological oxidations                        |        22 |        103 |       4.12 | 2.4e-08 |              0.27 | False                |
| Reactome  | Formation of the cornified envelope          |        16 |         74 |       2.96 | 4.8e-06 |              0.5  | True                 |
| Reactome  | Keratinization                               |        16 |         78 |       3.12 | 7.2e-06 |              0.5  | True                 |
| Reactome  | Phase I - Functionalization of compounds     |        13 |         54 |       2.16 | 1.4e-05 |              0.15 | False                |
| Reactome  | Metabolism                                   |        41 |        476 |      19.04 | 0.00011 |              0.15 | False                |
| Reactome  | Glucuronidation                              |         6 |         13 |       0.52 | 0.0004  |              1    | True                 |
| Reactome  | Aspirin ADME                                 |         8 |         27 |       1.08 | 0.0004  |              0.62 | True                 |
| Reactome  | Drug ADME                                    |        11 |         54 |       2.16 | 0.0004  |              0.55 | True                 |
| Reactome  | Phase II - Conjugation of compounds          |         9 |         44 |       1.76 | 0.0024  |              0.67 | True                 |
| Reactome  | Cytochrome P450 - arranged by substrate type |         8 |         36 |       1.44 | 0.0031  |              0.25 | False                |

**Evidence tags** (top 100 genes): plausible 42, limited evidence 40, known 18

**Top 10 genes by multi-omics support**

| gene     | layers              |   data_support | rna_higher_in   |   pubmed_lung |   pubmed_progression |   ot_score | string_known_partners   | evidence         |   gap_score | aliases_searched   |
|:---------|:--------------------|---------------:|:----------------|--------------:|---------------------:|-----------:|:------------------------|:-----------------|------------:|:-------------------|
| CALML3   | rna+methylation+cnv |          1.333 | early           |             7 |                    5 |      0.062 | KRAS                    | plausible        |        0.45 |                    |
| LRRC4    | rna+methylation+cnv |          1.301 | early           |             1 |                    1 |      0.025 |                         | limited evidence |        0.79 | NAG14,NGL-2        |
| FAM83C   | rna+methylation     |          1.155 | early           |             1 |                    1 |      0     |                         | limited evidence |        0.78 |                    |
| RAPGEFL1 | rna+methylation     |          1.028 | early           |             1 |                    0 |      0.007 |                         | limited evidence |        0.9  | Link-GEFII         |
| ADH7     | rna                 |          1     | early           |             5 |                    3 |      0.043 |                         | plausible        |        0.56 |                    |
| NTRK2    | rna+cnv             |          0.959 | early           |            97 |                   70 |      0.692 | EGFR,ERBB2,NTRK3,PIK3CA | known            |        0.02 | TRKB               |
| KRT7     | rna+methylation     |          0.957 | late            |           528 |                  390 |      0.061 |                         | known            |       -0.05 | CK7,CK-7,K2C7      |
| GLI2     | rna+methylation     |          0.953 | early           |           300 |                  176 |      0.097 |                         | known            |       -0.03 | THP2,HPE9,THP1     |
| NLGN1    | rna+methylation+cnv |          0.939 | early           |             2 |                    2 |      0.004 |                         | limited evidence |        0.61 | KIAA1070,NLG1      |
| RGMA     | rna+methylation+cnv |          0.923 | early           |             0 |                    0 |      0.009 |                         | plausible        |        0.84 |                    |

**Top 10 under-studied candidates** (highest gap score: strong data support, little progression literature)

| gene     | layers              |   data_support | rna_higher_in   |   pubmed_lung |   pubmed_progression |   ot_score | string_known_partners   | evidence         |   gap_score | aliases_searched            |
|:---------|:--------------------|---------------:|:----------------|--------------:|---------------------:|-----------:|:------------------------|:-----------------|------------:|:----------------------------|
| RAPGEFL1 | rna+methylation     |          1.028 | early           |             1 |                    0 |      0.007 |                         | limited evidence |        0.9  | Link-GEFII                  |
| RGMA     | rna+methylation+cnv |          0.923 | early           |             0 |                    0 |      0.009 |                         | plausible        |        0.84 |                             |
| LRRC4    | rna+methylation+cnv |          1.301 | early           |             1 |                    1 |      0.025 |                         | limited evidence |        0.79 | NAG14,NGL-2                 |
| FAM83C   | rna+methylation     |          1.155 | early           |             1 |                    1 |      0     |                         | limited evidence |        0.78 |                             |
| SOST     | rna+methylation     |          0.834 | early           |             0 |                    0 |      0.012 |                         | limited evidence |        0.76 | VBCH,DAND6                  |
| PCDH15   | rna+cnv             |          0.846 | early           |             1 |                    1 |      0.031 |                         | limited evidence |        0.66 | CDHR15,USH1F,DFNB23         |
| CYP26A1  | rna+methylation     |          0.836 | early           |             5 |                    1 |      0.022 |                         | plausible        |        0.64 | P450RAI,CP26,CYP26,P450RAI1 |
| NLGN1    | rna+methylation+cnv |          0.939 | early           |             2 |                    2 |      0.004 |                         | limited evidence |        0.61 | KIAA1070,NLG1               |
| ADH7     | rna                 |          1     | early           |             5 |                    3 |      0.043 |                         | plausible        |        0.56 |                             |
| MRAP2    | rna                 |          0.672 | early           |             0 |                    0 |      0.002 |                         | limited evidence |        0.54 | bA51G5.2,C6orf117           |

## Factor5

- stage AUC 0.565 (LUAD 0.587, LUSC 0.532); survival HR/SD 1.17 (FDR 0.06)
- variance explained: rna_LUAD 4.2%, methylation_LUAD 1.3%, cnv_LUAD 0.8%, rna_LUSC 4.6%, methylation_LUSC 1.0%, cnv_LUSC 0.3%
- RNA-supported genes: 1245 of 5000 measured

**Enriched gene sets** (FDR < 0.05: 12, of which 3 driven mainly by one genomic locus; top 10 shown)

| library   | term                                           |   overlap |   set_size |   expected |     fdr |   max_locus_share | positional_cluster   |
|:----------|:-----------------------------------------------|----------:|-----------:|-----------:|--------:|------------------:|:---------------------|
| Reactome  | Formation of the cornified envelope            |        16 |         74 |       2.96 | 9.7e-06 |              0.31 | False                |
| Reactome  | Keratinization                                 |        16 |         78 |       3.12 | 1.1e-05 |              0.31 | False                |
| Reactome  | Gap junction trafficking                       |         6 |         14 |       0.56 | 0.00083 |              0.5  | True                 |
| Reactome  | Gap junction trafficking and regulation        |         6 |         14 |       0.56 | 0.00083 |              0.5  | True                 |
| Reactome  | Gap junction assembly                          |         6 |         14 |       0.56 | 0.00083 |              0.5  | True                 |
| Reactome  | Developmental Biology                          |        33 |        412 |      16.48 | 0.0056  |              0.15 | False                |
| Reactome  | Diseases associated with surfactant metabolism |         4 |         10 |       0.4  | 0.029   |              0.25 | False                |
| Reactome  | Vesicle-mediated transport                     |        12 |        104 |       4.16 | 0.046   |              0.25 | False                |
| Reactome  | RET signaling                                  |         4 |         12 |       0.48 | 0.046   |              0.25 | False                |
| Reactome  | Membrane Trafficking                           |        10 |         78 |       3.12 | 0.046   |              0.3  | False                |

**Evidence tags** (top 100 genes): limited evidence 50, plausible 32, known 17, ambiguous symbol 1

**Top 10 genes by multi-omics support**

| gene    | layers              |   data_support | rna_higher_in   |   pubmed_lung |   pubmed_progression |   ot_score | string_known_partners   | evidence         |   gap_score | aliases_searched                           |
|:--------|:--------------------|---------------:|:----------------|--------------:|---------------------:|-----------:|:------------------------|:-----------------|------------:|:-------------------------------------------|
| ASCL1   | rna+methylation     |          1.255 | early           |           370 |                  222 |      0.051 | TP53                    | known            |       0     | ASH1,HASH1,bHLHa46                         |
| CTNND2  | rna+methylation     |          1.095 | early           |             6 |                    6 |      0.341 |                         | plausible        |       0.465 | NPRAP,GT24                                 |
| KHDRBS2 | rna+methylation+cnv |          1.07  | early           |             4 |                    3 |      0.026 |                         | limited evidence |       0.641 | SLM1,SLM-1,MGC26664,KHDRBS2-OT1,KHDRBS2-OT |
| SRCIN1  | rna+methylation+cnv |          1.064 | early           |            17 |                   11 |      0.108 |                         | plausible        |       0.232 | SNIP,p140Cap,KIAA1684                      |
| ATP8B3  | rna+methylation     |          1.045 | late            |             2 |                    2 |      0.004 |                         | limited evidence |       0.707 | ATPIK                                      |
| KRT16   | rna+cnv             |          1.031 | late            |            14 |                   11 |      0.107 |                         | plausible        |       0.212 | NEPPK                                      |
| KRT6A   | rna                 |          1     | late            |            44 |                   36 |      0.124 |                         | known            |       0.02  | CK6C,K6C,CK6D,K6D,KRT6D                    |
| LMO3    | rna+methylation+cnv |          0.979 | early           |            13 |                    8 |      0.065 |                         | limited evidence |       0.278 | Rhom-3,DAT1,RBTNL2                         |
| NOL4    | rna+cnv             |          0.963 | early           |             4 |                    2 |      0.007 |                         | limited evidence |       0.667 | NOLP,HRIHFB2255,CT125                      |
| KRT6B   | rna                 |          0.958 | late            |            10 |                    7 |      0.075 |                         | plausible        |       0.323 | KRTL1                                      |

**Top 10 under-studied candidates** (highest gap score: strong data support, little progression literature)

| gene    | layers              |   data_support | rna_higher_in   |   pubmed_lung |   pubmed_progression |   ot_score | string_known_partners         | evidence         |   gap_score | aliases_searched                           |
|:--------|:--------------------|---------------:|:----------------|--------------:|---------------------:|-----------:|:------------------------------|:-----------------|------------:|:-------------------------------------------|
| CALN1   | rna+cnv             |          0.884 | early           |             0 |                    0 |      0     |                               | limited evidence |       0.763 | CABP8                                      |
| CELF4   | rna+methylation+cnv |          0.883 | early           |             0 |                    0 |      0     |                               | limited evidence |       0.753 | BRUNOL4                                    |
| ATP8B3  | rna+methylation     |          1.045 | late            |             2 |                    2 |      0.004 |                               | limited evidence |       0.707 | ATPIK                                      |
| NOL4    | rna+cnv             |          0.963 | early           |             4 |                    2 |      0.007 |                               | limited evidence |       0.667 | NOLP,HRIHFB2255,CT125                      |
| DOK6    | rna+methylation+cnv |          0.834 | early           |             0 |                    0 |      0.027 | RET                           | plausible        |       0.662 | MGC20785,HsT3226,DOK5L                     |
| KHDRBS2 | rna+methylation+cnv |          1.07  | early           |             4 |                    3 |      0.026 |                               | limited evidence |       0.641 | SLM1,SLM-1,MGC26664,KHDRBS2-OT1,KHDRBS2-OT |
| TUBB2B  | rna+methylation+cnv |          0.93  | early           |             2 |                    2 |      0.613 | EML4,KIF5B,TUBA1B,TUBA1C,TUBB | known            |       0.626 | MGC8685,DKFZp566F223,bA506K6.1             |
| KRT6C   | rna                 |          0.908 | late            |             4 |                    2 |      0.058 |                               | plausible        |       0.596 | KRT6E                                      |
| HCG4    | rna+methylation     |          0.779 | late            |             0 |                    0 |      0     |                               | limited evidence |       0.53  | HCGIV-10,HCGIV.9,HCG4P10                   |
| NRXN3   | rna+methylation+cnv |          0.813 | early           |             1 |                    1 |      0.002 |                               | limited evidence |       0.525 | KIAA0743,C14orf60                          |

## Method notes and limitations

- Evidence tags: known = >= 20 lung-cancer progression papers or Open Targets score >= 0.5; plausible = in an enriched gene set (FDR < 0.05), a STRING partner (score >= 700) of a top-50 Open Targets NSCLC gene, or Open Targets score >= 0.1; limited evidence = none of these (not the same as novel - such genes may still have some literature and should be checked manually before any novelty claim).
- Aliases that equal another gene's official symbol (e.g. 'IRX-5' listed for IRX1) are excluded from literature queries.
- PubMed counts use title/abstract mentions of the gene symbol and its HGNC aliases/previous symbols (aliases shorter than 4 characters or ambiguous are dropped); symbols that are common words or very short are flagged 'ambiguous symbol' and excluded from novelty claims. Counts measure attention, not truth.
- Candidates are RNA-anchored; CNV features stand for several neighbouring genes, so CNV support is region-level and discounted by the square root of region size.
- Enriched gene sets whose overlapping genes mostly fall in one 2-Mb genomic bin are flagged `positional_cluster` (co-located gene families such as CYP4F or SERPINB can look like pathway enrichment) and are not used for 'plausible' tags.
- Methylation probes are assigned to their annotated genes; promoter vs gene-body effects are not distinguished.
- Candidates are hypotheses for follow-up, not validated biomarkers.

Outputs: results/knowledge/candidates.tsv, enrichment.tsv, knowledge_graph_edges.tsv, knowledge_graph.graphml, knowledge_sources.tsv, support_vs_literature.png