# UKB_anc
Analysis scripts for UKB ancestry Nature Genetic paper (accepted). If you find any of the scripts listed here is useful and you've used in your publication, please cite our preprint and journal paper as below:

**Preprint**:

<a id="ref-biorxiv">[1]</a> 
Hu, S.,et al. (2023). 
Leveraging fine-scale population structure reveals conservation in genetic effect sizes between human populations across a range of human phenotypes. 
Preprint at bioRxiv [https://doi.org/10.1101/2023.08.08.552281](https://doi.org/10.1101/2023.08.08.552281) (2023)

**Nature genetic**:

<a id="ref-ng">[2]</a> 
Hu, S.,et al. (2024).
Fine-scale population structure and widespread conservation of genetic effect sizes between human groups across traits.
Nature Genetics (accepted)

**Table of contents**

1. [Overview](#item-overview)
2. [Spatial mean plots](#item-spmean)
	1. [Define geographic regions for UK ancestries on the map](#item-spmean-bd)
	2. [Spatial mean ancestry plots](#item-spmean-acs)
	3. [Spatial mean entropy plots](#item-spmean-entrop)
	4. [Spatial mean AF plots](#item-spmean-af)
3. [Population structure plots](#item-popstr)
	1. [Structure barplots for 23 UK+Ireland Regions](#item-ukirl)
	2. [Structure barplots for UK cities](#item-ukcity)
4. [GWAS scripts](#item-gwas)
	1. [AC vs. PC](#item-gwas-acpc)
	2. [begenie-based gwas](#item-gwas-bgen)
	3. [boltlmm-based gwas](#item-gwas-bolt)
	4. [LD score regression](#item-gwas-ldsc)
4. [Estimate ancestry specific allele frequency by EM](#item-emaf)
5. [Mean-centered Ancestry PGS construction](#item-mcpgs)

<!-- headings -->
<a id="item-overview"></a>
## Overview

Following the publication policy on Nature genetics, we deposit our scripts and some of the public available data for reproducibility, though we can't provide any individual level data as input for some of the scripts at this repo. However, we believe following the guidence listed here, you can reproduce the results using the scripts at this repo.

<a id="item-spmean"></a>
## Spatial mean plots

This R scripts introduced in this section were used for generating the main Figure 1b, Extended Data Figure 1, 3 and 5

<a id="item-spmean-bd"></a>
### Define geographic regions for UK ancestries on the map

The R scrips used here were for defining the geographic regions of pre-defined GB/Ireland ancestry gruops in the ancestry pipeline, using the UKB participants born in the UK/Ireland

Source the R script file "define_boundary.r" in R session and run the command like this:

```r
source('define_boundary.r')

mb<-make.boundary()
```

Input files:

- make_boundary.rds: individual level data, a list of "*sp.data.frame*" object, with each of them the AC matrix for participants in one of the counties.

Output:

- "*mb*": a list of counties and their corresponding assigned ancestries.

As Cornwall is one county but in our ancestry pipeline, we found two clsuteres (groups) from the reference samples, so we need to define the geographic regions for these two groups: "Cornwall" and "Cornwall tip" by splitting the Cornwall county into two. To achieve this goal, we used the script "split_cornwall.r" to rasterise the map and assign ancestry to the "pixel" to define the boundaries of "Cornwall" and "Cornwall-tip" here:

```r
source('split_cornwall.r')

cornwall.list<-get.cornwall.pix(n=250,q=50)
```

The input files for generating plots/results are described as belows:

- data/GB_IRELAND.rds: *sp* object of GB+Ireland map
- data/GBR_adm2.rds: "sp" object of UK counties (secondary level)
- long_lat_ukbiobank.tsv: longtitude and lattitude of birth places of UKB participants
- biobank_v2_eth.background.rds: individual level data, the self-reported ethnicity background in the UKB
- biobank_v2_results.rds: individual level data, ACs matrix

Output:

"cornwall.list": is a list object including the rastered pixel for "Cornwall" and "Cornwall tip"

Finall, for visualisation purpose, we can rasterise the "*sp*" object of geographic boundaries for each of UK/Ireland groups using the script "raster_boundary.r" as below:

```r
source('raster_boundary')
```
After sourcing the script, you will end up with a R object called "mr", which is the rastered object for the 23 GB+Ireland gourps defined in the pipeline.

The input files for generating plots/results are described as belows:

- data/GB_IRELAND.rds: *sp* object of GB+Ireland map 
- data/new_boundary.rds: *sp* object of GB county boundary map data
- mapping.rds: *sp* individual level mapping data (mapping your application id to the UKB HRC imputed samples (same order)
- data/New_GB_boundaries.rds: *sp* object of GB boundary information
- v2_self_BI_487409.rds: individual level data, "sp.data.frame" object by mapping the ACs of each individual to the geographic coordinates. If the data frame is individual ancestral entropy or genotype, then this script can also be used to create spatial entropy plot in Extended data Figure 3 and regional allele frequency plot for Extended data Figure 5
- self_BI_counties_487409.rds: individual level data: county assignment for each of the UKB participants.

<a id="item-spmean-acs"></a>
### Spatial mean AC plots

The script used here is for Main Figure 1b and Extended Data Figure 1.

Simply source the script file "spatial_mean_ancestry_allele_freq_plot.R" in R session and run the command like this:

```r
source('spatial_mean_ancestry_allele_freq_plot.R')

plot.ma.gbirl(n=1000,q=50,dir="Spatial_AC_figure_out")
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

<a id="item-spmean-entrop"></a>
### Spatial mean Entropy plots

The script used here is for Extended Data Figure 3.

Simply source the script file "Entropy_plot_spatial.R" in R session and run the command like this:

```r
source('Entropy_plot_spatial.R')

plot.ma.entropy(n=1000,q=50,entropy.spdf=ukb.entropy,dir="Spatial_entropy_figure_out")
```

The input files for generating plots/results are described as belows:

- data/GB_IRELAND.rds: *sp* object of GB+Ireland map 
- data/UK_main_cities_England_Wales.rds.rds: *sp* object of main cities boundary in England + Wales
- data/Scotland/\*.rds: "sp" object of main cities of boundary in Scotland
- data/GBR_adm2.rds: "sp" object of UK counties (secondary level)
- long_lat_ukbiobank.tsv: longtitude and lattitude of birth places of UKB participants
- biobank_v2_eth.background.rds: individual level data, the self-reported ethnicity background in the UKB
- biobank_v2_results.rds: individual level data, ACs matrix

In the main fuction "plot.ma.entropy", user needs to provide three parameters: "*n*" is used for determing the dimension of pixel, e.g. *n*=1000 means the resolution of the figure will be 1000 * 1000 = 1,000,000; "*q*" is the parameter controlling the adaptive window of the Gaussian Kernal, by default we use *q*=50; "*dir*" is the output directory in which the figures will be generated. 

<a id="item-spmean-af"></a>
### Spatial mean AF plots

The script used here is for Extended Data Figure 5.

Simply source the script file "Spatial_regional_AF_raster.R" in R session and run the command like this:

```r
source('Spatial_regional_AF_raster.R')

plot.maf.gbirl(n=1000,q=50,dir="Spatial_AF_figure_out")
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
- data/noi_irl_region_af_mean.rds: Mean allele frequencies of 3 SNPs (including rs5743618 in Extended Data Figure 3) in Ireland and North Ireland.
- data/gwas_snp_regional_plot_af.rds: "*sp.data.frame*" object with the Data table is the genotype of each white British indviduals
- data/regions_af.rds: alelle frequencies of 3 SNPs (including rs5743618 in the Extended Data Figure 3)
- data/New_GB_boundaries.rds: *sp* object of GB boundary information

In the main fuction "plot.maf.gbirl", user needs to provide three parameters: "*n*" is used for determing the dimension of pixel, e.g. *n*=1000 means the resolution of the figure will be 1000 * 1000 = 1,000,000; "*q*" is the parameter controlling the adaptive window of the Gaussian Kernal, by default we use *q*=50; "*dir*" is the output directory in which the figures will be generated. 

<a id="item-popstr"></a>
## Plot population structure in the UK/Ireland regions

This R scripts described here were used for generating the structure barplot for Main Figure 2a, Extended Data Figure 3

<a id="item-ukirl"></a>
### Structure barplots for 23 UK+Ireland Regions

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

<a id="item-ukcity"></a>
### Structure barplots for UK cities

Simply source the script file "Entropy_plot_assign_individuals_to_cities.R" in R session and run the command like this:

```r
source('Entropy_plot_assign_individuals_to_cities.R')

make.figs()
```
The input files for generating plots/results are described as belows:

- v2_self_BI_487409.rds: individual level data, "sp.data.frame" object by mapping the ACs of each individual to the geographic coordinates. If the data frame is individual ancestral entropy or genotype, then this script can also be used to create spatial entropy plot in Extended data Figure 3 and regional allele frequency plot for Extended data Figure 5
- data/UK_main_cities_England_Wales.rds.rds: *sp* object of main cities boundary in England + Wales
- data/Scotland/\*.rds: "sp" object of main cities of boundary in Scotland
- data/GBR_adm2.rds: "sp" object of UK counties (secondary level)
- data/23_regions_color1.rds: colour assignment for each of the 23 UK+Ireland regions
- UKB_BI_ids_in_20_cities.rds: individual level data, the participants born in one of the 20 cities plotted

After call the main function "make.figs", the structure barplots of each 20 cities will be created at "GB_city_fig/" at the current working directory.

<a id="item-gwas"></a>
## GWAS scripts

The scripts in this section were used to run GWAS analysis and generate the input for the Main figure 3, Extended Data Figure 4,6 and Supplementary Figures 1-3

<a id="item-gwas-acpc"></a>
### AC vs. PC

The scripts used here were mainly for generating the Main Figure 3a-b and Supplemenetary Figure 1-2

<a id="item-emaf"></a>
## Estimate ancestry specific allele frequency by EM

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
