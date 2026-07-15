# Tumor Progression and Pan-Cancer Transcriptomic Analysis

## An Integrative Multi-Cancer Bioinformatics Framework for Molecular Network, Clinical, Single-Cell, and Therapeutic Validation

This repository contains the reproducible analysis workflow, processed results, figures, and supplementary materials associated with a multi-cancer transcriptomic study investigating molecular signatures associated with tumor progression and related cancer phenotypes across six cancer datasets.

The study integrates differential gene expression analysis, functional enrichment, protein–protein interaction network analysis, hub-gene prioritization, independent transcriptomic validation, pan-cancer integration, TCGA-based clinical validation, single-cell transcriptomic localization, and drug–gene interaction analysis.

---

## Study Overview

Six cancer datasets were analyzed using standardized transcriptomic workflows:

| Cancer Type | Discovery Dataset | Data Type |
|---|---|---|
| Breast cancer | GSE173839 | Microarray |
| Lung cancer | GSE33072 | Microarray |
| Ovarian cancer | GSE156699 | RNA-seq |
| Prostate cancer | GSE279730 | RNA-seq |
| Melanoma | GSE99898 | Microarray |
| Small-cell lung cancer (SCLC) | GSE188705 | RNA-seq |

Microarray datasets were analyzed using the `limma` framework, whereas RNA-seq datasets were analyzed using a DESeq2-based workflow.

Differentially expressed genes were standardized using:

- Benjamini–Hochberg false-discovery-rate correction
- Adjusted P value < 0.05
- Absolute log2 fold change ≥ 1.0

The resulting gene sets were subsequently evaluated using functional enrichment, interaction-network analysis, hub-gene prioritization, independent validation, clinical validation, single-cell localization, and therapeutic-relevance analysis.

---

## Analysis Workflow

The overall computational workflow consisted of the following major stages:

1. Transcriptomic data acquisition and preprocessing
2. Differential gene expression analysis
3. Standardization of DEG outputs across datasets
4. Quality-control and expression-pattern visualization
5. Gene Ontology and KEGG pathway enrichment analysis
6. STRING protein–protein interaction network analysis
7. Hub-gene identification using cytoHubba
8. Independent transcriptomic validation
9. Pan-cancer integrative analysis
10. Finalized hub-gene evidence integration
11. TCGA-based clinical validation
12. Candidate prioritization for single-cell analysis
13. TISCH-based single-cell localization
14. Drug–gene interaction analysis
15. Therapeutic target prioritization

The complete study workflow is summarized in Figure 1 of the manuscript.

---

## Repository Structure

```text
Tumor-Progression-PanCancer/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
│
├── Data/
│   ├── DEG_tables/
│   ├── Validation/
│   └── Metadata/
│
├── Scripts/
│   ├── 01_DESeq2.R
│   ├── 02_limma.R
│   ├── 03_GO_KEGG.R
│   ├── 04_STRING.R
│   ├── 05_cytoHubba.R
│   ├── 06_TCGA_validation.R
│   ├── 07_scRNA_validation.R
│   └── 08_DrugGene.R
│
├── Figures/
│
└── Supplementary/
    ├── Supplementary_Tables.xlsx
    └── Supplementary_Figures.pdf
