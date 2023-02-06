- [Feed-microbiome-host interactions in Atlantic salmon over life stages](#feed-microbiome-host-interactions-in-atlantic-salmon-over-life-stages)
  - [Getting Started](#getting-started)
    - [Step 1. Package dependencies in R](#step-1-package-dependencies-in-r)
    - [Step 2. Shell script](#step-2-shell-script)
    - [Step 3. Downstream analysis in R](#step-3-downstream-analysis-in-r)
- [Bugs](#bugs)


# Feed-microbiome-host interactions in Atlantic salmon over life stages
To meet future food demands, more efficient and sustainable animal production systems are needed. Given the crucial importance of the gut microbiota to animal (host) health and nutrition, selective enhancement of beneficial microbes via prebiotics may be a powerful approach for promoting farmed fish welfare and robustness. In this study, we employed three versions of a beta-mannan prebiotic that were fed to Atlantic salmon and explored the combined responses of both gut microbiota and the host from freshwater to seawater life stages. We have used weighted gene co-expression network analysis (WGCNA) of host RNA-seq and microbial 16S rRNA amplicon sequencing data to identify biological interactions between the gut ecosystem and the host. We observed several microbial and host modules through WGCNA that were significantly correlated with life stage, but not with diet. Microbial diversity was highest in the early life of Atlantic salmon and decreased over time. In particular, Lactobacillus and Paraburkholderia were the dominating genera and showed the highest correlation with host modules. Our findings demonstrate that salmon-microbiota interactions are mainly influenced by life stage, while further research is required to determine whether supplementation of selected prebiotics to diet can be used to modulate the salmon gut microbiota for improving host health and production sustainability.



## Getting Started
### Step 1. Package dependencies in R
Installing R packages can be done through various sources such as GitHub, the Comprehensive R Archive Network (CRAN), or by following the official website of the package.

To install a package from GitHub, use the devtools package and the `install_github()` function. For example, to install the "ggplot2" package from GitHub, run the following code:
```
library(devtools)
install_github("ggplot2")
```
To install a package from CRAN, use the `install.packages()` function. For example, to install the "dplyr" package from CRAN, run the following code:
```
install.packages("dplyr")
```
Finally, to install a package from its official website, download the package source code, and use the `install.packages()` function with the local file path as the argument. For example, to install the "reshape2" package from its official website, first download the source code, then run the following code:

```
install.packages("path/to/reshape2_package.tar.gz", repos = NULL, type = "source")
```
It is recommended to regularly update the installed packages to ensure compatibility and to benefit from new features and bug fixes.

### Step 2. Shell script
Primers were removed from the raw paired-end FASTQ files generated via MiSeq using “cutadapt”. Further, reads were analyzed by QIIME2 (qiime2-2021.8) pipeline through dada2 to infer the ASVs present and their relative abundances across the samples. For bed dust samples, using read quality scores for the dataset, forward and reverse reads were truncated at 280 bp and 260 bp, followed by trimming the 5′ end till 25 bp for both forward and reverse reads, respectively; other quality parameters used dada2 default values for both 16S rRNA gene sequencing. For 16S rRNA gene sequencing, taxonomy was assigned using a pre-trained Naïve Bayes classifier (Silva database, release 138, 99% ASV) were used.

The code is a shell script for processing paired-end sequencing data in order to perform a microbial analysis. It uses a combination of bash commands and QIIME2 (Quantitative Insights Into Microbial Ecology) commands.

```
ls -d -1 data/*_R1* > forward.txt
ls -d -1 data/*_R2* > reverse.txt

mkdir Trimmed
mkdir Trim.log

parallel -j 1 --xapply "cutadapt -g CCTAYGGGRBGCASCAG -G GGACTACHVGGGTWTCTAAT --pair-filter=any --discard-untrimmed -o $PWD/Trimmed/{1/.}.gz -p $PWD/Trimmed/{2/.}.gz {1} {2} &> $PWD/Trim.log/{1/.}.log" :::: forward.txt :::: reverse.txt

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format CasavaOneEightSingleLanePerSampleDirFmt --input-path Trimmed/ --output-path demux-paired-end.qza

qiime demux summarize --i-data demux-paired-end.qza --o-visualization step1_output.qzv

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza --o-table table --o-representative-sequences rep-seqs  --p-trunc-len-f 270 --p-trunc-len-r 250 --p-n-threads 4 --o-denoising-stats denoising-stats.qza --verbose &> dada2.log

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization step2_output.qzv
```
Creating Rooted Phylogenetic Tree
```
echo "Align rep-seqs and created rooted tree"
qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --p-n-threads 4
qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --p-n-threads 4
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
echo "Rooted tree has been created"
```
Taxonomy classification
```
echo "Classify rep-seqs using Silva132_99 classifier"
qiime feature-classifier classify-sklearn --i-reads rep-seqs.qza --i-classifier silva_99_341F806R_classifier.qza --p-n-jobs 6 --o-classification taxonomy.qza
echo "Classification done"
```
Export to biom file and import into phyloseq
```
echo "Export data from qiime"
qiime tools export  --input-path table.qza --output-path exported-feature-table
qiime tools export  --input-path taxonomy.qza --output-path exported-feature-table
echo "data has been exported"
change header of taxonomy file
```
sed -i "1 s/^.*$/#OTUID\ttaxonomy\tconfidence/" exported-feature-table/taxonomy.tsv
biom add-metadata -i exported-feature-table/feature-table.biom -o exported-feature-table/feature-table_taxonomy.biom --observation-metadata-fp exported-feature-table/taxonomy.tsv --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy
echo "BIOM file with OTU table and taxonomy saved as: exported-feature-table/feature-table_taxonomy.biom"
qiime tools export --input-path  rooted-tree.qza --output-path exported-feature-table/
echo "Rooted phylogenetic tree saved as: exported-feature-table/tree.nwk"
qiime tools export  --input-path merged_rep-seqs.qza --output-path exported-feature-table/
echo "Representative sequences saved as: exported-feature-table/dna-sequences.fasta"
```
Training the database
```
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
unzip Silva_132_release.zip
qiime tools import --type 'FeatureData[Sequence]' --input-path ./SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna  --output-path silva_99_seqs.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ./SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_all_levels.txt --output-path silva_99_consensus_taxonomy.qza
qiime feature-classifier extract-reads --i-sequences silva_99_seqs.qza --p-f-primer CCTAYGGGRBGCASCAG --p-r-primer GGACTACHVGGGTWTCTAAT --o-reads silva_99_341F806R.seqs.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva_99_341F806R.seqs.qza --i-reference-taxonomy silva_99_consensus_taxonomy.qza --o-classifier silva_99_341F806R_classifier.qza --verbose &> silva_99_341F806R_classifier.log
```

### Step 3. Downstream analysis in R
All downstream analyses were performed on this normalized ASVs table unless mentioned. We used two alpha diversity indices, i.e., observed richness and Shannon diversity index. Furthermore, beta diversity was calculated using weighted and unweighted UniFrac metric and visualized by principal coordinates analysis (PCoA). Alpha and beta diversity was calculated using phyloseq v1.38.0 and visualized with ggplot2 v3.3.5 in R v4.1.1. Comparison of community richness and diversity was assessed by the Kruskal-Wallis test between all the groups, and comparison between the two groups was done by Wilcoxon test with Benjamini-Hochberg FDR multiple test correction. Significance testing between the groups for beta diversity was assessed using permutational multivariate analysis of variance (PERMANOVA) using the “vegan” package.

Please follow the [link](README.html) to see all the downstream analysis done in R.


# Bugs
To inform us of any bugs, please open a new issue or send us an email to shashank.gupta@nmbu.no









