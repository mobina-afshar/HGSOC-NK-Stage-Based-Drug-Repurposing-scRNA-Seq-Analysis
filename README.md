# Stage-Based Drug Repurposing in High-Grade Serous Ovarian Cancer Using Single-Cell and Pseudobulk Transcriptomic Analysis of NK Cells

**Authors:** Mohadeseh Sadat Meftagh, Mobina Afshari Kaveh*

**Corresponding Author:** Mobina Afshari Kaveh

**Affiliations:**
1. Department of Genetics, Faculty of Sciences and Advanced Technologies in Biology, University of Science and Culture, Tehran, Iran.
2. Department of Genetics, Reproductive Biomedicine Research Center, Royan Institute for Reproductive Biomedicine, ACECR, Tehran, Iran.

---

## Overview

High-grade serous ovarian cancer (HGSOC) is an aggressive malignancy characterized by substantial molecular heterogeneity and progressive remodeling of the tumor microenvironment. Natural killer (NK) cells are important components of the tumor immune microenvironment and may undergo transcriptional changes during disease progression.

This project integrates single-cell RNA sequencing (scRNA-seq), patient-level pseudobulk differential expression analysis, pathway enrichment, and computational drug-repurposing approaches to investigate stage-associated transcriptional programs in NK-associated populations in HGSOC.

Using the publicly available dataset **GSE184880**, this study establishes a computational workflow extending from immune-cell characterization and cell-population selection to stage-associated pathway analysis and candidate compound prioritization.

The primary focus is the identification of transcriptional changes in NK-like cells between early- and late-stage HGSOC, followed by the computational prioritization and pharmacological annotation of candidate compounds.

## Research Objectives

- Characterize the immune-cell landscape of HGSOC tumors using single-cell transcriptomic data.
- Identify and characterize NK-associated cell populations within the tumor microenvironment.
- Investigate differences between early- and late-stage disease using patient-level pseudobulk differential expression analysis.
- Identify stage-associated biological pathways using Hallmark gene-set enrichment analysis.
- Prioritize candidate compounds with the potential to reverse stage-associated transcriptional signatures using the Connectivity Map (CMap).
- Investigate drug–gene interactions and pharmacological and clinical annotations for prioritized compounds using DGIdb and DrugBank.

## Study Design and Cell Populations

The study included seven HGSOC tumor samples from the GSE184880 dataset, classified into early- and late-stage disease groups.

Following immune-cell characterization and T/NK/ILC subclustering, two populations were selected for downstream investigation:

- **NK-like population (cluster 0):** The primary population for stage-associated differential expression, pathway analysis, and CMap-based drug prioritization.
- **ILC3 population (cluster 11):** Included as a complementary population for stage-associated differential expression and pathway-level investigation.

The NK-like designation represents an operational cluster annotation and should not be interpreted as evidence of a completely homogeneous or exclusively NK-cell population.

Patient-level pseudobulk profiles were generated using raw gene counts. The resulting analyses included:

| Population | Early-stage profiles | Late-stage profiles | Total |
|---|---:|---:|---:|
| NK-like | 3 | 4 | 7 |
| ILC3 | 3 | 3 | 6 |

Patients, rather than individual cells, were treated as the biological replicates in downstream differential expression analysis.

## Dataset

| Characteristic | Description |
|---|---|
| Dataset accession | GSE184880 |
| Data type | Single-cell RNA sequencing |
| Disease | High-grade serous ovarian cancer |
| Number of tumor samples | 7 |
| Disease-stage comparison | Early-stage vs. late-stage |
| Analysis platform | R |
| Main single-cell framework | Seurat |

The dataset is publicly available through the NCBI Gene Expression Omnibus:

**[GSE184880 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880)**

## Study Workflow

The computational workflow consisted of the following major steps:

1. **Data acquisition:** Retrieval of the HGSOC scRNA-seq dataset GSE184880.
2. **Quality control and preprocessing:** Assessment of cell quality, filtering of low-quality cells, and normalization using Seurat.
3. **Dimensionality reduction and clustering:** PCA, neighborhood graph construction, clustering, and UMAP visualization.
4. **Cell annotation:** Characterization of cell populations using canonical marker genes, supported by SingleR and CellTypist.
5. **T/NK/ILC subclustering:** Identification and annotation of immune-cell subpopulations.
6. **Population selection:** Selection of cluster 0 (NK-like) and cluster 11 (ILC3) for downstream investigation.
7. **Pseudobulk aggregation:** Aggregation of raw gene counts at the patient level for each selected population.
8. **Differential expression analysis:** Comparison of early- and late-stage samples using DESeq2.
9. **Pathway enrichment:** Preranked Hallmark GSEA using DESeq2 Wald statistics.
10. **Drug prioritization:** Querying the NK-like stage-associated transcriptional signature against CMap.
11. **Drug annotation:** Investigation of drug–gene interactions and pharmacological information using DGIdb and DrugBank.

## Methods

### 1. Single-Cell Processing and Quality Control

Single-cell RNA-seq data were processed using the Seurat framework. Quality control, normalization, variable-feature identification, dimensionality reduction, clustering, and UMAP visualization were performed to characterize the cellular landscape.

A total of **37,818 cells** were initially available, of which **29,390 cells** were retained following quality control, corresponding to approximately **77.7%** of the initial cell population.

### 2. Cell-State Annotation and Population Selection

Cell populations were characterized using canonical marker genes. SingleR and CellTypist were used as supporting annotation tools to assist cell-type identification.

Following global immune-cell characterization, T/NK/ILC subclustering was performed to investigate immune-cell heterogeneity in greater detail.

Two populations were selected for downstream stage-based analyses:

- **Cluster 0 (NK-like):** Selected as the primary population for pseudobulk differential expression, Hallmark GSEA, and CMap-based drug prioritization.
- **Cluster 11 (ILC3):** Selected as a complementary population for pseudobulk differential expression and pathway-level investigation.

The NK-like label is an operational annotation based on the analyzed cluster and does not imply that the population is necessarily composed exclusively of conventional NK cells.

### 3. Patient-Level Pseudobulk Analysis

To account for the biological replication structure of the dataset, raw gene counts were aggregated at the patient level.

Pseudobulk profiles were generated separately for the selected NK-like and ILC3 populations. Early-stage samples were compared with late-stage samples using DESeq2.

The final pseudobulk analysis included:

- NK-like: 7 patient-level profiles (3 early-stage and 4 late-stage).
- ILC3: 6 patient-level profiles (3 early-stage and 3 late-stage).

This approach treats patients, rather than individual cells, as the biological replicates and avoids treating individual cells from the same patient as independent biological observations.

### 4. Differential Expression Analysis

Differential expression analysis was performed using DESeq2 on patient-level pseudobulk raw counts.

The comparison was defined as **late-stage versus early-stage disease**, with early-stage samples serving as the reference group.

The resulting gene-level statistics, including log2 fold changes, Wald statistics, p-values, and adjusted p-values, were used for downstream pathway analysis and drug-signature construction.

### 5. Hallmark Gene Set Enrichment Analysis

Preranked Hallmark gene-set enrichment analysis was performed using the Wald statistics obtained from DESeq2.

All genes with valid ranking statistics were included in the ranking procedure rather than restricting the analysis to statistically significant differentially expressed genes.

Positive normalized enrichment scores (NES) indicate enrichment toward the late-stage group, whereas negative NES values indicate enrichment toward the early-stage group.

The analysis was conducted separately for NK-like and ILC3 populations.

For NK-like cells, the interferon-alpha response pathway was the only Hallmark pathway meeting the FDR-adjusted significance threshold. Additional pathways identified using nominal p-values were considered exploratory findings.

No FDR-significant Hallmark pathways were identified in the ILC3 analysis.

### 6. Connectivity Map-Based Drug Prioritization

The NK-like DESeq2 results were used to construct the CMap query signature based on Wald-statistic ranking.

The signature was derived independently of the GSEA leading-edge gene sets. Thus, CMap prioritization was based on the differential expression results rather than a subset of genes selected through pathway enrichment.

Compounds with connectivity scores **≤ −90** were prioritized for further investigation.

Negative connectivity scores indicate an inverse relationship between compound-associated perturbational expression profiles and the submitted query signature. Such compounds may represent candidates for further investigation based on their predicted transcriptional effects.

These scores represent computational prioritization and do not establish therapeutic efficacy, clinical suitability, or activity against HGSOC.

### 7. Drug–Gene and Pharmacological Annotation

Prioritized compounds were further investigated using two complementary resources:

- **DGIdb:** Used to retrieve drug–gene interaction records and their associated interaction annotations.
- **DrugBank:** Used to investigate pharmacological properties, mechanisms of action, regulatory information, and documented clinical applications.

Drug–gene interaction records were considered supporting annotation evidence rather than definitive proof of clinical approval or therapeutic relevance.

Regulatory status and clinical applications should be interpreted in the context of the specific drug records and, where necessary, confirmed using official regulatory labeling.

## Key Results

### 1. Single-Cell Landscape

- A total of 37,818 cells were initially available.
- Following quality control, 29,390 cells were retained for downstream analysis.
- Global clustering and immune-cell characterization were performed to investigate the cellular composition of HGSOC tumors.
- T/NK/ILC subclustering enabled the selection of NK-like and ILC3 populations for subsequent stage-associated analyses.

### 2. Patient-Level Differential Expression and Pathway Analysis

Patient-level pseudobulk analysis was performed separately for NK-like and ILC3 populations to investigate transcriptional differences between early- and late-stage disease.

Hallmark GSEA identified stage-associated biological programs in both populations.

In the NK-like population, **HALLMARK_INTERFERON_ALPHA_RESPONSE** was the only pathway meeting the FDR-adjusted significance threshold. Additional nominally enriched pathways were considered exploratory and should not be interpreted as FDR-significant findings.

No FDR-significant Hallmark pathways were identified in the ILC3 analysis.

### 3. CMap-Based Drug Prioritization

CMap analysis was performed using the NK-like stage-associated transcriptional signature.

Applying a connectivity-score threshold of **≤ −90**, seven candidate compounds were prioritized for further investigation:

- NSC-663284
- Phenethyl isothiocyanate
- LDN-193189
- Importazole
- BIIB021
- Mitomycin-C
- SN-38

These compounds represent computationally prioritized candidates based on their predicted inverse connectivity to the query signature. Their identification does not demonstrate therapeutic efficacy or clinical benefit in HGSOC.

### 4. Mitomycin-C: Pharmacological and Clinical Annotation

Among the prioritized compounds, mitomycin-C was selected for further pharmacological and clinical annotation.

DrugBank analysis characterized mitomycin-C as a DNA-alkylating chemotherapeutic agent and documented its clinical use in low-grade upper tract urothelial cancer (LG-UTUC).

This provides a clinically documented example of its use in a specific disease context. However, this annotation does not establish efficacy in HGSOC or demonstrate that mitomycin-C is effective in a particular HGSOC stage.

Further experimental investigation would be required to assess its relevance to HGSOC.

## Figures and Presentation Materials

The repository contains selected visualizations and presentation materials generated during the project.

### Presentation Materials

Available in the `Presentation_materials/` directory:

- Conference poster.
- Accepted conference abstract.

### Analysis Results

Selected results and visualizations are organized in the `Results/` directory:

- **Clustering:** Global cluster-size summaries.
- **CMap:** CMap results and prioritized candidate-compound visualization.
- **GSEA:** Hallmark enrichment results.
- **Metadata:** Sample-level metadata.
- **Pseudobulk:** Differential expression tables.
- **QC:** Quality-control visualizations before and after filtering.
- **RDS:** Selected serialized R objects and analysis results.
- **Umap:** Global and immune-cell subclustering visualizations.

## Repository Structure

```text
HGSOC-NK-Stage-Based-Drug-Repurposing-scRNA-Seq-Analysis/
│
├── Presentation_materials/
│   ├── Cancer Genomics Poster.jpg
│   └── Cancer genomics abstract.pdf
│
├── Results/
│   ├── Clustering/
│   │   └── Global_cluster_sizes.csv
│   │
│   ├── CMap/
│   │   ├── CMap_results.xlsx
│   │   └── CMap_Candidates.png
│   │
│   ├── GSEA/
│   │   └── Hallmark_GSEA_exploratory_nominalP05_NK_like.csv
│   │
│   ├── Metadata/
│   │   └── sample_metadata.csv
│   │
│   ├── Pseudobulk/
│   │   ├── DESeq2_Late_vs_Early_NK_like.csv
│   │   └── NK_like_DESeq2_all_results.csv
│   │
│   ├── QC/
│   │   ├── QC_after_filtering.png
│   │   └── QC_before_filtering.png
│   │
│   ├── RDS/
│   │   ├── NK_ILC3_StageAnalysis_Object.rds
│   │   ├── pseudobulk_counts.rds
│   │   └── DESeq2_Results_All_Populations.rds
│   │
│   └── Umap/
│       ├── Final_TNK_ILC_Umap.png
│       └── Global_umap_clusters.png
│
├── Scripts/
│   ├── 01_Project_Setup.R
│   ├── 02_Data_Import_and_QC.R
│   ├── 03_Preprocessing_and_Clustering.R
│   ├── 04_SingleR_Annotation.R
│   ├── 05_CellTypist_Annotation.R
│   ├── 06_TNK_ILC_Subclustering_and_Annotation.R
│   ├── 07_NK_ILC3_Selection_and_Pseudobulk.R
│   ├── 08_DESeq2_Stage_Analysis.R
│   ├── 09_Hallmark_GSEA.R
│   └── 10_CMap_and_Drug_Annotation.R
│   └── config.yaml
│
└── README.md
```

## Reproducibility and Data Availability

The analysis was performed in R using publicly available transcriptomic data and bioinformatics resources.

The original dataset can be accessed through NCBI GEO using accession **GSE184880**.

The repository provides the analysis scripts, configuration file, selected intermediate objects, and representative result files to document the computational workflow.

For reproducibility:

1. Clone or download this repository.
2. Download the original GSE184880 dataset from NCBI GEO.
3. Install the R packages and other software dependencies required by the analysis scripts.
4. Review and configure the project paths and parameters in `config.yaml`.
5. Ensure that the input data are placed in the expected project directories.
6. Execute the numbered scripts in their intended order, beginning with project setup.
7. Consult the corresponding subdirectories in `Results/` for generated outputs.

The numbered scripts are intended to represent the sequential analysis workflow. Exact reproduction may depend on the software environment, package versions, reference resources, and configuration parameters used during the original analysis.

## Limitations

This study is an exploratory computational investigation and has several limitations:

- The analysis includes only seven HGSOC tumor samples, limiting the number of independent biological replicates.
- Stage-associated findings should be interpreted cautiously because of the limited sample size and potential inter-patient heterogeneity.
- The NK-like population represents an operational cluster annotation and may retain cellular heterogeneity.
- CMap connectivity scores provide computational predictions and do not demonstrate drug efficacy, safety, or clinical benefit.
- Drug–gene interaction and pharmacological annotations depend on the coverage and evidence available in the referenced databases.
- The prioritization of mitomycin-C does not establish its therapeutic relevance to HGSOC.
- Independent validation using appropriate experimental models and larger patient cohorts is required before drawing translational conclusions.

## Conclusion

This project presents an integrated single-cell and patient-level pseudobulk transcriptomic framework for investigating stage-associated NK-like cell programs in HGSOC and prioritizing candidate compounds for further investigation.

By combining cell-level characterization, patient-level differential expression, Hallmark pathway enrichment, CMap-based compound prioritization, and pharmacological annotation, the study provides a computational strategy for exploring stage-informed drug-repurposing hypotheses in the HGSOC tumor microenvironment.

The findings provide a basis for further investigation rather than direct therapeutic recommendations.

## Conference Presentation

**3rd International Congress of Cancer Genomics and the First Provincial Congress of Cancer Genomics**

**Date:** 4–6 November 2026  
**Location:** Isfahan, Iran

**Poster title:**

*Stage-Based Drug Repurposing in High-Grade Serous Ovarian Cancer Using Single-Cell and Pseudobulk Transcriptomic Analysis of NK Cells.*

## Authors and Contributions

**Mohadeseh Sadat Meftagh**  
Department of Genetics, Faculty of Sciences and Advanced Technologies in Biology, University of Science and Culture, Tehran, Iran.

**Mobina Afshari Kaveh**  
Department of Genetics, Faculty of Sciences and Advanced Technologies in Biology, University of Science and Culture, Tehran, Iran; Department of Genetics, Reproductive Biomedicine Research Center, Royan Institute for Reproductive Biomedicine, ACECR, Tehran, Iran.

*Corresponding author: Mobina Afshari Kaveh.*

## Citation

If you use this repository, please cite the project as:

Meftagh MS, Afshari Kaveh M. **Stage-Based Drug Repurposing in High-Grade Serous Ovarian Cancer Using Single-Cell and Pseudobulk Transcriptomic Analysis of NK Cells.** Poster presented at the 3rd International Congress of Cancer Genomics and the First Provincial Congress of Cancer Genomics; 4–6 November 2026; Isfahan, Iran.

## Acknowledgments

We acknowledge the authors of the GSE184880 dataset and the developers and maintainers of the bioinformatics tools and databases used in this study, including Seurat, SingleR, CellTypist, DESeq2, MSigDB, Connectivity Map, DGIdb, and DrugBank.

---

**Disclaimer:** This repository documents computational research and candidate-compound prioritization. It does not provide clinical recommendations or establish the safety or efficacy of any compound for HGSOC.
