# Multi-Omics Integration in Multiple Sclerosis (MS)
 
This repository contains an integrative analysis of **RNA-seq** and **DNA methylation**
data from paired MS brain **lesion** and **NAWM** samples (N = 18).

We apply several analysis steps:

- mRNA differential expression (lesion vs NAWM)
- DNA methylation analysis (DMPs, immune infiltration)
- **MOFA** – unsupervised multi-omics factor analysis
- **DIABLO** – supervised multi-omics integration
- **SNF** – similarity network fusion for patient clustering
- **WGCNA** – co-expression / co-methylation network modules

---

## Graphical abstract

![Graphical abstract](graphical_abstract.png)

---

## Table of contents

1. [Project overview](#project-overview)
2. [Data](#data)
3. [mRNA analysis](#mrna-analysis)
4. [DNA methylation analysis](#dna-methylation-analysis)
5. [MOFA – unsupervised integration](#mofa--unsupervised-integration)
6. [DIABLO – supervised integration](#diablo--supervised-integration)
7. [SNF – similarity network fusion](#snf--similarity-network-fusion)
8. [WGCNA – network modules](#wgcna--network-modules)
9. [How to run the analysis](#how-to-run-the-analysis)
10. [Team](#team)

---

## Project overview

In this project, we integrate **transcriptomic (mRNA)** and **DNA methylation** data
from MS brain **lesion** and **NAWM** samples from the same patients.

Our main goals are to:

- identify **differentially expressed genes (DEGs)** and **differentially methylated positions (DMPs)**  
- link CpG methylation to gene expression (e.g. PALMD, SLAIN1)  
- discover **shared latent factors** capturing multi-omics variation (MOFA)  
- build **supervised multi-omics signatures** that discriminate lesion vs NAWM (DIABLO)  
- cluster patients using a **fused similarity network** (SNF)  
- identify **co-expression / co-methylation modules** associated with lesion/NAWM and clinical traits (WGCNA)

---

## Data

- **Samples:** 18 paired MS brain samples (lesion / NAWM)  
- **mRNA:** RNA-seq expression matrices (after filtering and normalisation)  
- **DNA methylation:** Illumina array beta/M-values for matched samples  
- **Clinical:** phenotype / metadata tables (e.g. patient ID, lesion/NAWM, age, gender) stored under:

  ```text
  data/clinical/
