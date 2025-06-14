# NGS data analysis method

## Data analysis method of 16 S rRNA amplicon
### Introduction
Root exudation is pivotal for plants to cope with the changing environment through interacting with the soil organisms. 16S rRNA amplicon next generation sequencing (NGS) is a powerful tool to understand how root exudate compound affect the mrobial composition in the root, rhizosphere and soil. Here we show the data analysis method developed by Dr. Peng Wang (wangplant@gmail.com) and Dr. Daniel Schachtman (daniel.schachtman@unl.edu) using different packages in R, combined with QIIME v1.9.1, QIIME v2 and USEARCH v10. 

### Overview

**16 S rRNA data analysis consists of various steps by using different software or packages which are as follows:**

*Using 2 x 300 bp V4 region from the MiSeq platform as the example*

1. Quality filteration in USEARCH (version: 10.0.240) with the UPARSE pipeline. Primers and low quality reads are removed.
2. Dereplication in USEARCH. To find the unique sequences and create the input sequences of the 97% OTU clustering.
3. OTU clustering in USEARCH. 97% identity was used as the threshold.
4. Generation of OTU table in USEARCH.
5. Taxonomy assignment. ([QIIME v2](https://docs.qiime2.org/2020.8/tutorials/moving-pictures/) or [RDP](http://qiime.org/scripts/assign_taxonomy.html))
6. Removal of plastid and mitochondria in the OTU table by using the QIIME's embeded python commands.
7. alpha-diversity analysis using the QIIME's embeded python commands.
8. beta-diversity analysis (PCoA analysis with Bray-Curtis dissimilarity matrix) using the QIIME's embeded python commands.
9. Canonical Analysis of Principal coordinates (CAP) analysis. Using the *vegan* package in R and including the PERMONOVA analysis.
10. Differential abundance analysis between wild type and RNAi lines
11. Co-occurrence networks analysis. Using the *SparCC* python command lines and *igraph* package in R. 


### Steps and command lines


### 1. Quality filteration in USEARCH


**Note:** USEARCH is available online with full instructions (https://www.drive5.com/usearch/), which is developed by Dr. Robert Edgar

After demultiplexing, the paired-end reads were merged with error correction using USEARCH. The maximum number of mismatches in the alignment was set at 10 base-pairs and minimum percentage of identity in the alignment was set to 80%. 

```
usearch -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastq_minmergelen 230 -fastq_maxmergelen 320 -fastq_pctid 80 -fastqout merged.fq
```
**Note:** The parameters are set through referencing the USEARCH instruction manual.


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

**Note:** This step also incorporates the removal of singletons from the clustered OTUs and removal of chimeras from sequencing data

### 4. Generation of OTU table

```
usearch -otutab stripped.fq -otus otus.fa -otutabout otutab.txt -mapout map.txt
```

**Note:** This command generates a table with the number of reads (counts) of all OTUs for each sample. The OTU table is used for downstream steps including differential abundance analyses and microbial diversity analyses

### 5. Taxonomy assignment

- using the ribosomal database project classifier (RDP) by python command emmbedded in QIIME v1.9.1
 
 ```
 assign_taxonomy.py -i otus.fa -m rdp -c 0.80
 ```
**Note:** Install the MacQiime in the laptop if using MAC OSX operating system

- using a pre-trained Naive Bayes classifier and the q2-feature-classifier plugin in QIIME v2 (recommended) [here is the link](https://docs.qiime2.org/2020.8/tutorials/feature-classifier/). To install the QIIME v2, [here is the link](https://docs.qiime2.org/2020.8/install/native/#install-qiime-2-within-a-conda-environment). *`conda install`* is recommended for MAC OSX.

To run the command in the QIIME v2 environment

```
conda activate qiime2-2020.8
```

To train the classifer and do the taxonomical assignment

```
qiime tools import --type 'FeatureData[Sequence]' --input-path 97_otus.fasta --output-path 97_otus.qza

qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 97_otu_taxonomy.txt --output-path ref-taxonomy.qza

qiime feature-classifier extract-reads --i-sequences 97_otus.qza --p-trunc-len 250 --p-min-length 100 --p-max-length 400 --o-reads ref-seqs.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza

qiime tools import --type 'FeatureData[Sequence]' --input-path ./new_otu_fa/otus.fa --output-path otus.qza

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads otus.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
```
**Note:** Greengenes 13_8 97% OTUs (reference sequences clustered at 97% sequence similarity) and its corresponding 97% taxonomy information are required to use at the same time. Please check the [link](https://docs.qiime2.org/2020.8/tutorials/feature-classifier/) with the title "*Obtaining and importing reference data sets*" for the detailed information


### 6. Filteration. Removal of plastid and mitochondria in the OTU table. Additionally, OTUs that were not assigned at a Kingdom level RDP classification score of 0.8 were discarded


- First, merge the RDP assignment information from the last step into the OTU table and make the "biom" format

```
biom convert -i otutable_rdp.txt -o otu_table_rdp.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
```

- Second, remove the OTUs that were not assigned at a Kingdom level RDP classification, named as "Unclassified"

```
filter_taxa_from_otu_table.py -i otu_table_rdp.biom -o otu_table_rdp_no_unknow.biom -n Unclassified
```

- Third, filter out the OTUs assigned to chloroplast and mitochondrion in the OTU table

```
filter_taxa_from_otu_table.py -i otu_table_rdp_no_unknow.biom -o otu_table_rdp_no_unknow_m_c.biom -n f__mitochondria,c__Chloroplast
```

### 7. alpha-diversity analysis


- First, Conduct a rarefaction analysis in QIIME v1.9.1, to estimate the lowest the sequencing depth that can capture the most different OTUs or richness to reaching the saturation 

```
biom summarize-table -i otu_table_rdp.biom >  otu_table_rdp.txt
```

and

```
single_rarefaction.py -i otu_table.biom -o otu_table_even100.biom -d 100
```

- Second, all the samples were rarefied to specific reads per sample for alpha-diversity analyses

```
alpha_rarefaction.py -i otu_table.biom -m Mapping_file.txt -o output_lpha_type_number -p shannon_index.txt -e raref_number -t all_phy_rooted_fasttree_midpoint.tre
```

**Note:** This analysis calculates diversity indices such as Shannon, Simpson, and Chao1


### 8. beta-diversity analysis

```
beta_diversity_through_plots.py -i otu_table.biom -o bd_type_number -p beta_d.txt -m Map.txt -e raref_number -t all_phy_rooted_fasttree_midpoint.tre
```

### 9. Canonical Analysis of Principal coordinates (CAP) analysis

- Conduct canonical analysis of principal coordinates (CAP) analysis using the capscale function by using the *`vegan`* package

- Permutational multivariate analysis of variance (PERMANOVA). Using the bray-curtis dissimilarity matrix generated in the step 8 to do the PERMANOVA analysis.

- Visulization by using *`ggplot2`* R package

**Note:** *please refer to `Transgenic_CAP.R`* in this github

### 10. Differential abundance analysis

```
differential_abundance.py -i otu_table_rdp_no_m_c_unknow.biom -o otu_table_rdp_no_m_c_unknow.txt -a DESeq2_nbinom -m Map.txt -c Genotype -x WT -y RNAi -d
```
**Note:** use the DESeq2 embedded in the Qiime v1.9.1

### 11. Co-occurrence networks analysis

- First, calculate the core OTU matrix in QIIME that can be used for following step

```
compute_core_microbiome.py -i otu_table.biom -o otu_table_core
```
**Note:** core-OTU is defined with more than 70% threshold indicating one OTU need to present in at least 70% samples

- Second, set up the python2.6.9 environment, which is required to run the SparCC

```
conda create --name SparCCEnv python=2.6.9
source activate SparCCEnv
conda install numpy=1.9.2
conda install pandas=0.16.2
```
- Third, download the *`SparCC_source_code`* file, inside of which has the python codes that are developed by Dr. Jonathan Friedman, and please check the [website](https://web.mit.edu/almlab/sparcc.html) for detail information on how to use it. The dataset in the folder *SparCC_data* can be for test.

- Fourth, navigate into the python code folder and change the mode of the file to make them executable

```
cd SparCC_source_code
chmod a+x *.py
ls -althr
```

- Fifth, make the correlation matrix by using the SparCC, spearman, pearson. Here we use the spearman as an example and actually are also used in the paper

    - construct the correlation matrix
    
    ```
    python SparCC.py ./AUG_WT_all.txt -i 20 --cor_file=./cor_AUG_WT_all.out -a spearman
    ```
    
    - pseudo p-values were calculated via a bootstrap procedure with at least 100 shuffles(1000 is strongly suggested but depending on the computer processing capacity) to determine the significance of the correlationship
    
    ```
    python MakeBootstraps.py ./AUG_WT_all.txt -n 100 -t permutation_#.txt -p ./pvals_WT_all/
    ```
    
    ```
    cd pvals_WT_all
    ```
    
    here is a loop for process all the data generated from last step
    ```
    for fq in permutation_*.txt; do python ../SparCC.py ./$fq -i 5 --cor_file=./perm_cor_$fq; done
    ```
    
    calculate the two tailed p-value 
    ```
    python ../PseudoPvals.py ../cor_AUG_WT_all.out ./perm_cor_#.txt 100 -o ./pvals.two_sided.txt -t two_sided
    ```
    
- Sixth, network visualization and ANCOVA stat analysis. Network analysis and visulization are done in R package *`igraph`*

**Note:** *please refer to `Network.R` and `ANCOVA.R`* in this github. The dataset in the folder *`Network_stats_data`*  can be download for the test to make graph.






**Example**: 2D network graph as shown in R studio for the rhizosphere core OTUs of wild type of transgenic sorghum in August

![Image of WT](https://github.com/SchachtmanLab/Transgenic-sorghum-sorgoleone/blob/master/2D_network.png)




**Example**: 3D network graph using the "rglplot" function in *`igraph`* R package for the rhizosphere core OTUs of wild type in August. You can hold the graph and turn to see it with different angles using rglplot

![Image of WT](https://github.com/SchachtmanLab/Transgenic-sorghum-sorgoleone/blob/master/3D_network.png)

  

