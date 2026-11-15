# Stage-Based Drug Repurposing in High-Grade Serous Ovarian Cancer Using Single-Cell and Pseudobulk Transcriptomic Analysis of NK Cells

**Authors:** Mohadeseh Sadat Meftagh, Mobina Afshari Kaveh*

**Affiliations:**
1. Department of Genetics, Faculty of Sciences and Advanced Technologies in Biology, University of Science and Culture, Tehran, Iran.
2. Department of Genetics, Reproductive Biomedicine Research Center, Royan Institute for Reproductive Biomedicine, ACECR, Tehran, Iran.

---

## Overview

High-grade serous ovarian cancer (HGSOC) is an aggressive malignancy characterized by substantial molecular heterogeneity and progressive remodeling of the tumor microenvironment. Natural killer (NK) cells are important components of the tumor immune microenvironment and may undergo transcriptional changes during disease progression.

This project integrates single-cell RNA sequencing (scRNA-seq), patient-level pseudobulk differential expression analysis, pathway enrichment, and computational drug-repurposing approaches to investigate stage-associated transcriptional programs in NK-associated populations in HGSOC.

The analysis uses the publicly available dataset **GSE184880** and follows a computational workflow from immune-cell characterization to candidate compound prioritization.

## Research Objectives

- Characterize the immune-cell landscape of HGSOC tumors using single-cell transcriptomic data.
- Identify NK-associated cell populations and characterize their transcriptional profiles.
- Investigate differences between early- and late-stage disease using patient-level pseudobulk analysis.
- Identify stage-associated biological pathways using Hallmark gene-set enrichment analysis.
- Prioritize compounds with the potential to reverse stage-associated transcriptional signatures using the Connectivity Map (CMap).
- Investigate drug–gene interaction information and pharmacological and clinical annotations for prioritized compounds using DGIdb and DrugBank.

## Study Workflow

The analysis followed these major steps:

1. **Data acquisition:** Retrieval of the HGSOC single-cell RNA-seq dataset GSE184880.
2. **Quality control and preprocessing:** Filtering of low-quality cells and normalization using Seurat.
3. **Dimensionality reduction and clustering:** PCA, neighborhood graph construction, clustering, and UMAP visualization.
4. **Cell annotation:** Identification of immune-cell populations using canonical marker genes, supported by SingleR and CellTypist.
5. **NK-associated population selection:** Selection of NK-associated populations for downstream stage-based analyses.
6. **Pseudobulk aggregation:** Aggregation of raw gene counts at the patient level.
7. **Differential expression:** Comparison of early- and late-stage samples using DESeq2.
8. **Pathway analysis:** Hallmark gene-set enrichment analysis to identify stage-associated biological programs.
9. **Drug prioritization:** Querying stage-associated transcriptional signatures against CMap.
10. **Drug annotation:** Investigation of drug–gene interactions using DGIdb and pharmacological and clinical information using DrugBank.

## Dataset

| Characteristic | Description |
|---|---|
| Dataset | GSE184880 |
| Data type | Single-cell RNA sequencing |
| Disease | High-grade serous ovarian cancer |
| Number of tumor samples | 7 |
| Comparison | Early-stage vs. late-stage disease |
| Analysis platform | R / Seurat |

The dataset is publicly available through the NCBI Gene Expression Omnibus:

[ GSE184880 – NCBI GEO ](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880)

## Methods

### 1. Single-Cell Processing and Quality Control

Single-cell RNA-seq data were processed using the Seurat framework. Quality control, normalization, variable-feature identification, dimensionality reduction, clustering, and UMAP visualization were performed to characterize the cellular landscape.

A total of **37,818 cells** were initially available, of which **29,390 cells** were retained following quality control, corresponding to approximately **77.7%** of the initial cell population.

### 2. Cell-State Annotation

Cell populations were characterized using canonical marker genes. SingleR and CellTypist were used as supporting annotation tools to assist cell-type identification.

NK-associated populations were selected for downstream analyses. The resulting annotations were used to investigate transcriptional differences among immune-cell populations and between disease stages.

### 3. Patient-Level Pseudobulk Analysis

To account for the biological replication structure of the dataset, raw gene counts were aggregated at the patient level.

Pseudobulk profiles were generated for the selected NK-associated population(s), and early-stage samples were compared with late-stage samples using DESeq2.

This approach treats patients, rather than individual cells, as the biological replicates.

### 4. Differential Expression and Hallmark GSEA

Differential expression analysis was performed using DESeq2.

The resulting gene-level statistics were used for Hallmark gene-set enrichment analysis to identify biological pathways associated with early- and late-stage disease.

The enrichment results were used to investigate stage-associated transcriptional programs, including immune-related and cellular signaling pathways.

### 5. Connectivity Map-Based Drug Prioritization

Stage-associated transcriptional signatures were queried against the Connectivity Map (CMap) to identify compounds with expression profiles predicted to oppose the queried disease-associated signature.

Compounds with connectivity scores **≤ −90** were prioritized for downstream investigation.

Negative connectivity scores indicate an inverse relationship between a compound-associated perturbational signature and the query signature. These scores represent computational prioritization and do not establish therapeutic efficacy.

### 6. Drug–Gene and Pharmacological Annotation

Prioritized compounds were further investigated using two complementary resources:

- **DGIdb:** Used to retrieve drug–gene interaction records and their associated interaction annotations.
- **DrugBank:** Used to investigate pharmacological properties, mechanisms of action, regulatory information, and documented clinical applications.

DGIdb drug–gene interaction records should not, by themselves, be interpreted as definitive evidence of FDA approval. Regulatory status should be confirmed using the relevant DrugBank records and, where necessary, official regulatory labeling.

## Key Results

### Single-Cell Landscape

- A total of 29,390 cells passed quality control.
- NK-associated populations represented a major component of the annotated immune compartment.
- UMAP visualization was used to display the annotated immune-cell landscape and the distribution of NK-associated populations.

### Stage-Associated Pathway Analysis

Patient-level pseudobulk analysis identified stage-associated transcriptional programs between early- and late-stage NK-associated cells.

Hallmark GSEA highlighted differences in biological pathway enrichment between the two disease-stage groups.

### CMap Drug Prioritization

CMap analysis identified multiple compounds with strong negative connectivity scores.

Using a prioritization threshold of ≤ −90, seven candidate compounds were retained for further investigation:

- NSC-663284
- Phenethyl isothiocyanate
- LDN-193189
- Importazole
- BIIB021
- Mitomycin-C
- SN-38

### Mitomycin-C Prioritization

Among the prioritized compounds, mitomycin-C was selected for further pharmacological and clinical annotation.

DrugBank analysis characterized mitomycin-C as a DNA-alkylating chemotherapeutic agent and documented its clinical use in low-grade upper tract urothelial cancer (LG-UTUC).

This provides a clinically documented example of its use in a specific disease context. It does not establish efficacy in HGSOC or demonstrate that mitomycin-C is effective in a particular HGSOC stage.

## Figures and Poster

The main figures generated for this project are available in the `figures/` directory.

The conference poster is available in:

`poster/HGSOC_Poster.pdf`

The poster was presented at the **3rd International Congress of Cancer Genomics and the First Provincial Congress of Cancer Genomics**, held on 4–6 November 2026.

## Repository Structure

```text
HGSOC-NK-Stage-Based-Drug-Repurposing-scRNA-Seq-Analysis/
│
├── Presentation_materials/
│   ├── Cancer Genomics Poster.jpg
│   └── Cancer Genomics abstract.pdf
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

*The directory and script names above describe the intended organization. Actual filenames may differ depending on the analysis files included in the repository.*

## Reproducibility and Data Availability

The analysis was performed in R using publicly available transcriptomic data and bioinformatics resources.

The original dataset can be accessed through NCBI GEO using accession **GSE184880**.

Analysis scripts and selected output files are provided in this repository to document the computational workflow.

For reproducibility:

1. Download the original dataset from GEO.
2. Install the R packages required by the corresponding analysis scripts.
3. Update input and output directory paths as needed.
4. Run the scripts in the order indicated by their filenames.
5. Consult the corresponding result folders for generated outputs.

Package versions and software requirements should be recorded alongside the scripts to support exact reproduction of the analysis.

## Limitations

This study is an exploratory computational investigation and has several limitations:

- The analysis includes seven HGSOC tumor samples, limiting the number of independent biological replicates.
- Stage-associated findings should be interpreted cautiously because of the limited sample size and potential inter-patient heterogeneity.
- CMap connectivity scores provide computational predictions and do not demonstrate drug efficacy or clinical benefit.
- Drug–gene interaction and pharmacological annotations depend on the coverage and evidence available in the referenced databases.
- The proposed mitomycin-C candidate requires further validation in appropriate experimental models before any therapeutic conclusions can be drawn.

## Conclusion

This project presents an integrated single-cell and pseudobulk transcriptomic framework for investigating stage-associated NK-cell programs in HGSOC and prioritizing candidate compounds for further investigation.

By combining cell-level characterization, patient-level differential expression, pathway enrichment, CMap analysis, and pharmacological annotation, the study provides a computational strategy for exploring stage-informed drug-repurposing hypotheses in the HGSOC tumor microenvironment.

## Conference Presentation

**3rd International Congress of Cancer Genomics and the First Provincial Congress of Cancer Genomics**

4–6 November 2026

Isfahan, Iran

**Poster title:** *Stage-Based Drug Repurposing in High-Grade Serous Ovarian Cancer Using Single-Cell and Pseudobulk Transcriptomic Analysis of NK Cells.*

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

**Disclaimer:** This repository documents computational research and candidate prioritization. It does not provide clinical recommendations or establish the safety or efficacy of any compound for HGSOC.
