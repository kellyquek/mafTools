# mafTools - A package to process and visualize maf files. 

With advances in Cancer Genomics, maf format is being widley accepted and used to store variants detected. 
[The Cancer Genome Atlas](http://cancergenome.nih.gov) Project has seqenced over 30 different cancers with sample size of each cancer type being over 200. The [resulting data](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files) consisting of genetic variants is stored in the form of [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification). This package attempts to summerize such files either from TCGA or in house studies by providing various functions for plotting and parsing. 

#### Dependencies: 
Make sure you have following packages installed - plyr, rjson, reshape, ggplot2, RColorBrewer, [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap), gridExtra.

#### Installation: 
`library(devtools)`

`install_github(repo = "PoisonAlien/mafTools")`


Four main functions:

1. `maf_stats` Takes maf file as input and returns a list of 3 tables (number of variants per barcode, variants classified according to their type and classfication). Also plots a barplot of stacked barplot of variants per tumor sample barcode and a boxplot of distribution of different classes of variants. 
From TCGA [LAML](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files#TCGAMAFFiles-LAML:AcuteMyeloidLeukemia) study:

![alt tag](https://github.com/PoisonAlien/mafTools/blob/master/DATA/maf_stats.png)

2. `TiTv_stats` Takes maf file as input and classifies Single Nucleotide Variants into [Transtions and Transversions](http://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html). Returns a list tables with variants summerized into different classes of conversions. A boxplot is also plotted showing distribution of different conversion events.
From TCGA [LAML](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files#TCGAMAFFiles-LAML:AcuteMyeloidLeukemia) study: 

![alt tag](https://raw.githubusercontent.com/PoisonAlien/mafTools/master/DATA/tcga_aml_TiTv.tiff)

3. `oncoPlot` Takes maf file as input and draws a matrix similar to [oncoprint](http://www.cbioportal.org/faq.jsp#what-are-oncoprints) on [cBioPortal](http://www.cbioportal.org/index.do). This function uses excellent [oncoprint](https://github.com/jokergoo/ComplexHeatmap/blob/908b32ee4c495c74adfa077c967024a77c56b375/vignettes/oncoprint.R) script from [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) package by [Zuguang Gu](https://github.com/jokergoo), while taking care of format conversions.
From TCGA [LAML](https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files#TCGAMAFFiles-LAML:AcuteMyeloidLeukemia) study: 

![alt tag](https://github.com/PoisonAlien/mafTools/blob/master/DATA/oncoprint.png)

4. `oncotate` Takes variants as input and annotates them using borads oncotator api (http://www.broadinstitute.org/oncotator/). Output is a dataframe of annotated variants in maf format. Input should be a five column file with chr, start, end, ref_allele, alt_allele. Youd can also provide extra columns (such as tumor_sample_barcode, validation_status, etc), but only first five will be used, rest will attached to the resulting maf file.  
Note: Time consuming if input is huge.