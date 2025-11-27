# Multi-Omics Integration in Multiple Sclerosis (MS)

This project performs an integrative analysis of **RNA-seq** and **DNA methylation**
data from paired MS brain **lesion** and **NAWM** samples (N = 18).

We use a set of single-omics and multi-omics methods to identify:

- Differentially expressed genes (DEGs)
- Differentially methylated CpG sites (DMPs)
- Links between methylation and gene expression
- Latent factors (MOFA)
- Supervised multi-omics signatures (DIABLO)
- Patient clusters based on fused networks (SNF)
- Gene/CpG network modules (WGCNA)

---

## Graphical Abstract

> Below is the overall workflow summarizing the full analysis pipeline.

![Graphical abstract](graphical_abstract.png)

---

## Table of Contents

1. [Project Overview](#project-overview)  
2. [Data](#data)  
3. [mRNA Analysis](#mrna-analysis)  
4. [DNA Methylation Analysis](#dna-methylation-analysis)  
5. [MOFA – Unsupervised Integration](#mofa--unsupervised-integration)  
6. [DIABLO – Supervised Integration](#diablo--supervised-integration)  
7. [SNF – Similarity Network Fusion](#snf--similarity-network-fusion)  
8. [WGCNA – Co-expression/Co-methylation Modules](#wgcna--co-expressionco-methylation-modules)  
9. [How to Run](#how-to-run)  
10. [Team](#team)

---

## Project Overview

This repository integrates **transcriptomic (mRNA)** and **DNA methylation** data
from MS brain **lesion** and **NAWM** samples.

The goal is to:

- understand molecular differences between lesion and NAWM  
- identify epigenetic–transcriptomic signatures  
- discover patient subgroups and biological modules  
- compare different multi-omics integration methods  

---

## Data

- **Samples:** 18 patient-matched lesion / NAWM samples  
- **mRNA:** filtered, normalised expression matrices  
- **DNA methylation:** beta or M-values from Illumina array  
- **Clinical metadata:** stored in `data/clinical/`  

Raw full-size matrices are excluded using `.gitignore`.  
Small subsets used in analysis are stored under the `data/` directory.


## mRNA Analysis

This section describes the differential expression analysis of the **RNA-seq** data
comparing MS brain **lesion** vs **NAWM** samples.

---

### 🔹 1) Quality Control of Expression Values

Before running the differential expression analysis, the RNA-seq counts were 
filtered, normalised, and transformed (log2).  
The histogram below shows the distribution of processed expression values.

![mRNA QC](mRNA/00_mRNA_processed_counts_hist.png)

The distribution looks typical for RNA-seq data and shows no major artefacts.

---

### 🔹 2) Differential Expression (Lesion vs NAWM)

A differential expression analysis (e.g. DESeq2) identified genes significantly
up- or down-regulated in lesions compared to NAWM.

**Global volcano plot:**

![Volcano plot](mRNA/01_mRNA_lesion_vs_NAWM_volcano.png)

- Each dot = gene  
- X-axis = log2 fold change  
- Y-axis = −log10(p-value)  
- Red/blue points show significantly regulated genes

---

### 🔹 3) Annotated Volcano (Key DEGs)

We highlight biologically relevant DEGs (immune-related, myelin-related, signalling genes).

![Annotated volcano](mRNA/3_mRNA_lesion_vs_NAWM_volcano_annotated.png)

This plot makes it easier to identify potential biomarkers.

---

### 🔹 4) Heatmap of Significant DEGs

A heatmap of the top differentially expressed genes shows **clear separation**
between lesion and NAWM samples.

![DEG heatmap](mRNA/02_mRNA_lesion_vs_NAWM_heatmap.png)

Lesion and NAWM samples form distinct clusters, reflecting strong transcriptional differences.

---

### ✔️ Summary (mRNA)

- Lesion samples show strong transcriptional activation relative to NAWM  
- Clear gene-level signatures emerge (DEGs)  
- These genes are later integrated with methylation in MOFA, DIABLO, and SNF

---

## DNA Methylation Analysis

This section summarises the DNA methylation differences between **MS lesion** and **NAWM** brain tissue.  
We identify differentially methylated CpG sites (DMPs), visualise their patterns, and explore immune-cell involvement.

---

### 🔹 1) Quality Control of Methylation Values

Beta- and M-values were inspected to ensure correct distribution before filtering and normalisation.

![Methylation QC](methylation/pm_value.png)

The bimodal distribution is typical for Illumina methylation arrays and indicates good data quality.


---

### 🔹 2) Immune Cell Infiltration (CIBERSORT)

Epigenetically inferred immune-cell proportions show higher immune activation in lesion samples.

![Immune cells](methylation/immunecells_group.jpeg)

Lesions exhibit elevated T cells, macrophage activity, and innate immune signals—consistent with MS pathology.

---

### 🔹 3) Example Gene: PALMD (CpG + Expression)

PALMD shows multiple CpGs significantly altered between lesion and NAWM.  
These methylation changes correlate with its gene expression level.

![PALMD example](methylation/PALMDboxplot.jpeg)

This highlights a direct epigenetic–transcriptomic link relevant to MS lesion biology.

---

### ✔️ Summary (Methylation)

- Strong methylation signatures separate lesion vs NAWM  
- Many CpGs map to genes relevant to MS (e.g., PALMD)  
- Immune infiltration strongly increases in lesions  
- Results integrate later with MOFA, DIABLO, and WGCNA


## MOFA – Unsupervised Multi-Omics Integration

MOFA (Multi-Omics Factor Analysis) was used to learn shared latent factors that
capture coordinated variation between **mRNA expression** and **DNA methylation**.
These factors help summarise the major biological drivers separating **lesion**
and **NAWM** samples.

---

### 🔹 1) Variance Explained per View

MOFA learns latent factors (F1, F2, …) and estimates how much variation each
factor explains in each omic (mRNA vs methylation).

![MOFA variance](MOFA/plotvarianceF.png)

- **Factor 1** explains the largest proportion of methylation variance  
- **Factor 2/3** capture shared variation across both omics  
- Higher factors explain smaller residual structure  

This indicates strong coordinated structure between the two omics layers.

---

### 🔹 2) Factor Space (Lesion vs NAWM)

Plotting samples along the learned factors shows clear biological separation.

![MOFA scatter](MOFA/correlation between f1,3.png)

- Lesion and NAWM samples cluster apart in the MOFA factor space  
- Indicates the latent factors capture meaningful disease-related variation  

---

### 🔹 3) Top Features Driving the Factors

MOFA provides feature loadings showing which genes/CpGs contribute most to each factor.

**Top mRNA features (e.g., Factor 1):**

![Top mRNA features](MOFA/topfeatureF1mrna.png)

**Top CpG features (e.g., Factor 1 or 6):**

![Top CpG features](MOFA/topfeatureF6cpg.png)

These genes and CpGs form the core multi-omics signatures differentiating lesion vs NAWM.

---

### 🔹 4) Associations with Clinical Metadata

MOFA factors were tested against metadata such as tissue type, age and gender.

![MOFA metadata heatmap](MOFA/hheatmaP.png)

- Factors associated with **lesion/NAWM** represent biologically meaningful signals  
- Other factors capture patient-level heterogeneity  

---

### ✔️ Summary (MOFA)

- Learned latent factors capture shared variation across RNA + methylation  
- Factor 1 strongly separates lesion from NAWM  
- Multi-omics signatures involve genes (e.g., immune-related) and CpGs with coordinated changes  
- Metadata associations support biological interpretation  

