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


## Steps and command lines


### 1. Quality filteration in USEARCH


**Note: USEARCH is available online with full instructions (https://www.drive5.com/usearch/).**

After demultiplexing, the paired-end reads were merged with error correction using USEARCH. The maximum number of mismatches in the alignment was set at 10 base-pairs and minimum percentage of identity in the alignment was set to 80%. 

```
usearch -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastq_minmergelen 230 -fastq_maxmergelen 320 -fastq_pctid 80 -fastqout merged.fq
```
**Note: The parameters are set through referencing the USEARCH instruction manual.**


Remove the primers from the sequencing data to avoid substitutions in the primer sequences, which may be caused by the PCR reaction. 

```
usearch -fastx_truncate merge.fq -stripleft 19 -stripright 20 -fastqout stripped.fq
```

Filter sequencing data to remove the low-quality reads and keep high quality operational taxonomic unit (OTU) sequences.

```
usearch -fastq_filter stripped.fq -fastq_maxee 1.0 -fastaout filtered.fa
```

### 2. Dereplication to prepare for the OTU table generation

Perform the dereplication to identify the set of unique OTU sequences

```
usearch -fastx_uniques filtered.fa -fastaout uniques.fa -sizeout -relabel Uniq
```

### 3. OTU clustering with 97% threshold

```
usearch -cluster_otus uniques.fa -minsize 2 -otus otus.fa -relabel Otu
```

**This step also incorporates the removal of singletons from the clustered OTUs and removal of chimeras from sequencing data**

### 4. Generation of OTU table

```
usearch -usearch_global stripped.fq -db otus.fa -strand plus -id 0.97 -otutabout otutable.txt
```

**This command generates a table with the number of reads (counts) of all OTUs for each sample. The OTU table is used for downstream steps including differential abundance analyses and microbial diversity analyses**

### 5. Taxonomy assignment

 using the ribosomal database project classifier (RDP) by python command emmbedded in Qiime 
 
 ```
 assign_taxonomy.py -i otus.fa -m rdp -c 0.80
 ```
**Note: Install the MacQiime in the laptop if using OSX operating system**

### 6. Filteration. Removal of plastid and mitochondria in the OTU table. Additionally, OTUs that were not assigned at a Kingdom level RDP classification score of 0.8 were discarded


First, merge the RDP assignment information from the last step into the OTU table and make the "biom" format

```
biom convert -i otutable_rdp.txt -o otu_table_rdp.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
```

Second, remove the OTUs that were not assigned at a Kingdom level RDP classification, named as "Unclassified"

```
filter_taxa_from_otu_table.py -i otu_table_rdp.biom -o otu_table_rdp_no_unknow.biom -n Unclassified
```

Third, filter out the OTUs assigned to chloroplast and mitochondrion in the OTU table

```
filter_taxa_from_otu_table.py -i otu_table_rdp_no_unknow.biom -o otu_table_rdp_no_unknow_m_c.biom -n f__mitochondria,c__Chloroplast
```

### 7. alpha-diversity analysis


–First, Conduct a rarefaction analysis in QIIME v1.9.1


```
cd 
```
