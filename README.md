# SchachtmanLab NGS data analysis method

## Data analysis method of 16 S rRNA amplicon
### Introduction
Root exudation is pivotal for plants to cope with the changing environment through interacting with the soil organisms. 16S rRNA amplicon next generation sequencing (NGS) is a powerful tool to understand how root exudate compound affect the mrobial composition in the root, rhizosphere and soil. Here we show the data analysis method developed by Dr. Peng Wang (wangplant@gmail.com) and Dr. Daniel Schachtman (daniel.schachtman@unl.edu) using different packages in R, combined with Qiime v1.9.1 and USEARCH v10. 

## Overview

**16 S rRNA data analysis consists of various steps by using different software or packages which are as follows:**

*Using 2 x 300 bp V4 region from the MiSeq platform as the example*

1. Quality filteration in USEARCH (version: 10.0.240) with the UPARSE pipeline. Primers and low quality reads are removed.
2. Dereplication in USEARCH. To find the unique sequences and create the input sequences of the 97% OTU clustering.
3. OTU clustering in USEARCH. 97% identity was used as the threshold.
4. Generation of OTU table in USEARCH.
5. Taxonomy assignment. Using the Ribosomal Database Project classifier (RDP) in QIIME's embeded python commands.
6. Removal of plastid and mitochondria in the OTU table by using the QIIME's embeded python commands.
7. alpha-diversity analysis using the QIIME's embeded python commands.
8. beta-diversity analysis (PCoA analysis with Bray-Curtis dissimilarity matrix) using the QIIME's embeded python commands.
9. Canonical Analysis of Principal coordinates (CAP) analysis. Using the *vegan* package in R and including the PERMONOVA analysis.
10. Co-occurrence networks analysis. Using the *SparCC* python command lines and *igraph* package in R. 


## Steps

1. Quality filteration in USEARCH
**Note: USEARCH is available online with full instructions (https://www.drive5.com/usearch/).**




```
cd 
```
