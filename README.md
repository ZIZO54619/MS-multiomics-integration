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

![Graphical abstract](<img width="2160" height="1215" alt="2" src="https://github.com/user-attachments/assets/be8caf8c-1619-40ca-b0ee-b966926fe9e6" />)

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
