# UKB_anc
Analysis scripts for UKB ancestry Nature Genetic paper

**Table of contents**

1. [Overview](#item-overview)
2. [Spatial mean plots](#item-spmean)
3. [Population structure plots](#item-popstr)
4. [Estimate ancestry specific allele frequency using EM based algorithm] (#item-emaf)
5. [Mean-centered Ancestry PGS construction](#item-mcpgs)

<!-- headings -->
<a id="item-overview"></a>
## Overview

Following the publication policy on Nature genetics, we deposit our scripts and some of the public available data for reproducibility, though we can't provide any individual level data as input for some of the scripts at this repo. However, we believe following the guidence listed here, you can reproduce the results using the scripts at this repo.

<a id="item-spmean"></a>
## Spatial mean plots

This R script "spatial_mean_ancestry_allele_freq_plot.R" was used for generating the main Figure 1b, Extended Data Figure 1, 3 and 5

Simply source the script file "spatial_mean_ancestry_allele_freq_plot.R" in R session and run the command like this:

```r
source('spatial_mean_ancestry_allele_freq_plot.R')

plot.ma.gbirl(n=1000,q=50,dir="Spatial_figure_out")
```

The input files for generating plots/results are described as belows:

- data/GB_IRELAND.rds: *sp* object of GB+Ireland map 
- data/new_boundary.rds: *sp* object of GB county boundary map data
- mapping.rds: *sp* individual level mapping data (mapping your application id to the UKB HRC imputed samples (same order)
- data/New_GB_boundaries.rds: *sp* object of GB boundary information
- withdraw_list.rds: individual level data, the withdraw list you received fro the UKB
- biobank_v2_eth.background.rds: individual level data, the self-reported ethnicity background in the UKB
- v2_487409.rds: individual level data, ACs matrix
- self_BI_487409.rds: individual level data, LOGI vector indicating if born in UK/Ireland
- v2_self_BI_487409.rds: individual level data, "sp.data.frame" object by mapping the ACs of each individual to the geographic coordinates. If the data frame is individual ancestral entropy or genotype, then this script can also be used to create spatial entropy plot in Extended data Figure 3 and regional allele frequency plot for Extended data Figure 5

In the main fuction "plot.ma.gbirl", user needs to provide three parameters: "*n*" is used for determing the dimension of pixel, e.g. *n*=1000 means the resolution of the figure will be 1000 * 1000 = 1,000,000; "*q*" is the parameter controlling the adaptive window of the Gaussian Kernal, by default we use *q*=50; "*dir*" is the output directory in which the figures will be generated. 


<a id="item-popstr"></a>
## Plot population structure in the UK/Ireland regions

This R script "UK_population_structure_plot.R" was used for generating the structure barplot for Main Figure 2a, Extended Data Figure 3

Simply source the script file "UK_population_structure_plot.R" in R session and run the command like this:

```r
source('UK_population_structure_plot.R')

barplot.gb(dir2="GB_region_ancestry_fig")
```
The input files for generating plots/results are described as belows:

- data/rastered_23GB_regions.rds: *sp* "*raster*" object of 23 GB+Ireland groups
- mapping.rds: *sp* individual level mapping data (mapping your application id to the UKB HRC imputed samples (same order)
- withdraw_list.rds: individual level data, the withdraw list you received fro the UKB
- biobank_v2_eth.background.rds: individual level data, the self-reported ethnicity background in the UKB
- v2_487409.rds: individual level data, ACs matrix
- self_BI_487409.rds: individual level data (UKB application specific), the IDs of UKB participants (in the same participants order as UKB HRC imputed data)
- v2_self_BI_487409.rds: individual level data, "sp.data.frame" object by mapping the ACs of each individual to the geographic coordinates. If the data frame is individual ancestral entropy or genotype, then this script can also be used to create spatial entropy plot in Extended data Figure 3 and regional allele frequency plot for Extended data Figure 5
- 426879_ind_birth_place.rds: individual level data, matched birth places for each of the White British/Irish participants
- Irish.rds: individual level data, ACs matrix for participants born in Ireland
- North_Irish.rds: individual level data, ACs matrix for participants born in North Ireland
- data/23_regions_color1.rds: colour assignment for each of the 23 UK+Ireland regions
- data/neighbour_region.rds: list object, for each of the 23 UK+Ireland regions, a sublist of regions neighbouring to that region

Just call the function "barplot.gb" and give the path of the directory for the output figures.  

<a id="item-emaf"></a>
## Estimate ancestry specific allele frequency using EM based algorithm
We used an *EM* based algorithm to estimate the allele frequency for ancestry regions (the ancestry regions are pre-defined). User need to provide
the following files to run the software:
* Ancestry matrix file contains a *N*x*A* matrix where *N* is the number of sample and *A* is the number of ancestries
* Genotype matrix file contains a *N*x*M* matrix where *N* is the number of (same) sample (as the Ancestry matrix) and *M* is the number of SNPs
* SNP annotation file contains at least one column called "SNP" include SNP IDs as columns in Genotype Matrix
* Output directory, if not specified, the software will try to create a folder under the directory where this script is running
* To imporve the performance, user can specify the number of threads to parallise the jobs, however, please bear in mind that the RAM may not be allocated properly so you job may be failed if you specify a big number.

To run the software, user should make sure "R/Rscript" has been installed on the research environment. To run the software, using the following command line:

```bash
./Run_EM_allele_frequency_estimation.sh -a PATH_TO_ancestry_file -g PATH_TO_genotype_file -s PATH_TO_SnpAnnotationFile -O PATH_TO_OUTPUT_DIR -p 4
```

If the input files are correct, the software will generate 2 files at "PATH_TO_OUTPUT_DIR" folder, the allele frequency file is called "PATH_TO_SnpAnnotationFile.freq.tab.gz"

<a id="item-mcpgs"></a>
## Mean-centered Ancestry PGS construction

As illustrated in our paper, we brought up an mean-centered ancestry PGS construction methods to deconvolute ordinary PGS into two independent ancestry PGSs (e.g. European and African PGS for European-African admixed populations).

Conventional approaches to construct such Ancestry PGSs will often lead to correlation between the PGS, and therefore make the results not robust. Our method, however, can correct the potential collinearity between the ancestry PGSs by introducing so-called *Mean-centering* technique to remove the correlation between the ancestry PGSs. 

To construct the mean-centered PGS, please use the function built in the script "generate_mean_centered_PGS.r". 

We would expect the user to run "hapmix" (diploid mode) (detail can be seen at "hapmix" paper) before using this script, and our script expect output from "hapmix"

Suppose user generated "hapmix" output at "/home/userA/hapmix_test" directory. "hapmix" sample id file is "/home/userA/data/sample_ind", "HAPMIX_MODE" for "hapmix" is "deploid" and admixed population is "ADMIXPOP"="EU_AF", then user can first open an R session and source the script as follows:

```r
source('generate_mean_centered_PGS.r')
``` 

Then user can run the command as follow to generate some important intermediate fileat specified folder "/home/userA/hapmix_PGS". User can also specify the number of cores to be used to speed up the job

```r
ri<-read.imp4("/home/userA/hapmix_test/data/sample_ind","/home/userA/hapmix_test/out/","EU_AF","diploid","/home/userA/hapmix_PGS",10)
```

Afterwards, user can call function "estimate_f" to generate the effect allele frequency between two populations:

```r
fq<-estimate_f()
```

To "mean-centered" genotypes (so as to PGS), the user should run the following function:

```r
mcg<-mean.center.geno()
```

By providing a table including the effect size estimation from other GWAS analysis, in which there must be a column called "SNP" and "BETA", user can generate their own ancestry PGS as follows():

```r
my.pgs<-cal.pgs('height_beta.rds')
```

Please bear in mind that we expect the input table should be in ".rds" format. If user's original file is plain text table, user can read it into R and use "saveRDS" function to convert the format of the original table.
