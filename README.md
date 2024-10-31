# UKB_anc
Analysis scripts for UKB ancestry Nature Genetics paper (accepted). If you find any of the scripts listed here is useful and you've used in your publication, please cite our preprint and journal paper (Nature Genetics) as below:

**Preprint**:

<a id="ref-biorxiv">[1]</a> 
Hu, S., et al. (2023). 
Leveraging fine-scale population structure reveals conservation in genetic effect sizes between human populations across a range of human phenotypes. 
Preprint at bioRxiv [https://doi.org/10.1101/2023.08.08.552281](https://doi.org/10.1101/2023.08.08.552281) (2023)

**Nature Genetic**:

<a id="ref-ng">[2]</a> 
Hu, S., et al. (2024).
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
	2. [bgenie-based GWAS](#item-gwas-bgen)
	3. [boltlmm-based GWAS](#item-gwas-bolt)
	4. [LD score regression](#item-gwas-ldsc)
	5. [GWAS peak plots](#item-gwas-peak)
	6. [Geno PC correlation](#item-gwas-gpc)
5. [Portablity of PGS](#item-pgs)
	1. [PGS calculation](#item-pgs-cal)
	2. [Run Hapmix](#item-pgs-hapmix)
	3. [Mean-centered Ancestry PGS construction](#item-pgs-mcpgs)
	4. [Plot European and African PGS effect size in bins](#item-pgs-bin)
	5. [Forest plot for effect size estimation for individual trait](#item-pgs-fst)
	6. [Simulation scripts and plots](#item-pgs-simu)
	7. [Trio PGS](#item-pgs-trio)
6. [Estimate ancestry specific allele frequency by EM](#item-emaf)

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

Simply source the script "AC_PC_pred.r" as follow and run the following functions in R session:

```r
source('AC_PC_pred.r')

### generate predicted PCs (by 127 ACs)
app<-ac.pred.pc()

### generate predicted ACs (by 140 PCs)
ppa<-pc.pred.ac()

### plot predicted PCs vs. real PCs
plot.ac.pred()

### plot predicted ACs vs. real ACs
plot.pc.pred()
```

Input files:

- v2_140PCs_487314.rds: individual level data, 140 PCs matrix across 487414 UKB participants
- v2_127ACs_487314.rds: individual level data, 127 ACs matrix across 487414 UKB participants

As UKB only released 40 PCs, the following github link will be useful for you to generate 140 PCs (highly close to the original UKB PCs with more than 99% correlation). [“pcapred.largedata”](https://github.com/danjlawson/pcapred.largedata) and [“pcapred”](https://github.com/danjlawson/pcapred) were used in our anlysis to generate the 140 UKB PCs.

Output files:

- "AC_pred_17_140PCs2.png": scatter plot from the function "plot.ac.pred()"
- "PC_pred_127AC2.png": scatter plot from the function "plot.pc.pred()"

<a id="item-gwas-bgen"></a>
### bgenie-based GWAS

Our gwas analysis was mostly conducted by ["bgenie"](https://jmarchini.org/software/#bgenie). The gwas results from "bgenie" were used in making the Main figure 3c-e.

Running the GWAS scripts also depends on the scheduler on the HPC (e.g. qsub or slurm. Here we just put the core scripts to call "bgenie": "run_bgenie_assoc.sh".

Plesase specify the input bash variables as follows to run the script.

```bash
BGENIEPATH=#PATH to the bgenie binary file
bgenfile=#PATH to the (imputed) bgen file
threads=#Number of threads for bgenie
Snpfile=#Subset of SNPs on which you run GWAS
phenofile=#phenotype file
covarfile=#covariate file
outfile=#output prefix for gwas summary statistics
```

To compare the performance between AC-corrected GWAS and PC-corrected GWAS, you only need to choose different covariate files containing either ACs columns or PCs columns (but not both) 

<a id="item-gwas-bolt"></a>
### boltlmm-based GWAS

We also conducted linear mixed model based GWAS by using software ["boltlmm"](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html). The gwas results from "boltlmm" were used in making the Supplementary Figure 3.

Here we use the script to call "bolt-lmm": "run_boltlmm.sh".

Plesase specify the input paths for bash variables as follows to run the script.

```bash
BOLT= # path to the executable bolt software
chipdir= #directory at which the plink model builds files are
modelbuildprefix= #plink bed/fam/bim prefix
removelistfile= #sample ids to be removed
impdir= #directory where the (imputed) bgen files are,bgen files are named as "ukb_imp_chr2_v3.bgen" and sample files are named as "ukb27960_imp_chr2_v3_s487324.sample"
dataodir= # output directory for bolt-lmm
phenoPath= # path to the phenotype file
covarPath= # path to the covariate file
ldscorePath= # path to the LDSCORE files
genmapPath= # path to the genetic map files
```

<a id="item-gwas-ldsc"></a>
### LD score regression

When we got the GWAS summary statistics, we applied the [ld score regression](https://github.com/bulik/ldsc) "LDSC" to compare the performance of correcting the gwas by ACs/PCs. We used the R script "run_ldsc.R" to wrap the LDSC software. We used this script to generate the input data for Extended Data Figure 4.

To use this script, user need to provide the following inputs files:

- PATH to the gwas summary statistics files as output from "bgenie". To save the storage space, each sumstat file only contain 4 "columns": beta, se, tstat and log10p
- PATH to the temporary directory: store the intermediate files from "munge" function in "ldsc" package
- PATH to the annotation file: Annotation file including the following columns: "vid": (chr:pos:ref:alt format). As the file is very big and not suitable to be deposited at github, you can find the Dropbox link here for the one used in the analysis <https://www.dropbox.com/scl/fi/ntn9lunw5rwuzy4ycgmcs/new_gwas_snps_anno_withaf.rds?rlkey=ewpeyf1ppil9sepufoia6xg04&dl=0> or generate your own. 
- sample size: sample size of the GWAS study

To run the script, just open a R session and run command as follows:

```r
source('run_ldsc.R')

run.all.ldsc()
```

Output of ldsc result will be directed to the temporary directory.

<a id="item-gwas-peak"></a>
### GWAS peak plots

The R scripts "peak_calling.r" and "gwas_plot_peak.r" used here were for generating the Main Figures 3c-e and Supplementary Figure 3. The script will select the peak/lead variants from either AC-corrected or PC-corrected GWAS.

To call the peak/lead variants either in AC-coreccted GWAS or PC-corrected GWAS from the output summary statistic data, user need to provide the following paths and assign them to the variables in R script:

- res.dir: the PATH to the directory where the sumstat files by AC-corrected/PC-corrected approaches stored
- var.anno.file: the PATH to the variants annotation file including two columns: "variant" (chr:pos:REF:ALT format) and "rsid"

Then user just simply run the following command in R:

```r
source('peak_calling.r')
# file.ac and file.pc are path to the AC-corrected GWAS and PC-corrected GWAS files, which have the same order of variants
# file.ac should include the following 2 columns: "ac.beta","ac.pv"
# file.pc should include the following 2 columns: "pc.beta","pc.pv" 
op<-output.peak(file.ac='gwas_ac.rds', file.pc='gwas_pc.rds')
```
The output "*op*" contains information of the peak variants for either AC-corrected or PC-corrected GWAS.

Then user can run the script "gwas_plot_peak.r" to plot the figure given the input from the output of "peak_calling.r", in this R script, we also provided the results from our analysis "data/98_traits_for_ldsc_EA.rds" as input to create the plots:

```r
source('gwas_plot_peak.r')
#tp variable expects input from the output of "op" from "peak_calling.r"
plot.phenos(dir.out='new_traits_peak/')
```
<a id="item-gwas-gpc"></a>
### Geno PC correlation

As we spotted 5 indpended signals were masked by the PC-corrected GWAS due in trait "waist circumference", we used this script to calculate the correlations between genotypes of these 5 SNPs and PCs. User need to provide raw genotype files and PC files, which are both individual level data, and run the R script as follows:

```r
source('calculate_var_PC_correlations.r'
cor.pc_geno<-pc.pred.geno(chr=15)
```
To generate the Extended figures 6, just plot the PC loadings against the genomic position of the variants.

<a id="item-pgs"></a>
## Portablity of PGS

The scripts described in this section were used for generating the main results in the portablity of PGS between ancestry groups. 

The central package used in the analysis is "ANCHOR", which can be accessed at our lab github: <https://github.com/MyersGroup/ANCHOR>.

Please note that the Main results in Main Figures 4-5, Extended Data Figures 7-10, and Supplemenary Figures 8-16 were generated from "ANCHOR" package, and scripts in this section were mainly used for PGS calculation, simulation and plots.  

<a id="item-pgs-cal"></a>
### PGS calculation

The scripts used here was mainly for calculating the "*normal*" polygenic score, including the results in the Main Figures 4-5, Extended Data figures 7-10 and Supplementary Figures 4-6 and 8-16

To construct the PGS, user need to first run the LD clumping. The R script "LD_clump.r" is for this purpose. To use this script, user needs to run the script as follows:

```r
source('LD_clump.r')
result.list<-run.all(p=0.05)
## if you need p-value threshold 0.01 or something else, simply change p as 0.01)
``` 

The following input files will be needed:

- res.dir: Path of directory where the gwas sumstats files stored
- clump.dir: A temporary directory called "PS_hapmap3_r2_0,1/" will be created at the current working directory, out put of LD clumping result will be created here
- var.anno.file: Path to the annotation file of the variants listed in the GWAS sumstats file. As the file is very big, use the Dropbox link to access the file we used in our analysis: <https://www.dropbox.com/scl/fi/ntn9lunw5rwuzy4ycgmcs/new_gwas_snps_anno_withaf.rds?rlkey=ewpeyf1ppil9sepufoia6xg04&dl=0> 
- hapmap3.snps: Annotation of Hapmap3 SNPs. When constructing the PGS, we restricted the list to Hapmap3 subset. User can use the Dropbox link here to access the file used in our analysis: <https://www.dropbox.com/scl/fi/3esbagk908ts355qmncwo/hapmap3_subset_SNPs_PS_filter_mean0.rds?rlkey=bmzxsryurxxn04hathe27yfg8&dl=0>

Given the above input files, user just need to call the "run.all" function and the plink based LD-clumping will be called internally and generate the clumped SNP list for each of the summary statistic data.

Once the LD-clumping is done, the next step is to calculate the PGS based on the clumped SNPs. Here we provided the R script "PGS_cal.r" to conduct this work

To use this script, user just need to run the script in the R seession as follows:

```r
source('PGS_cal.r')

all.pgs<-run.all.pgs()
```

The script will automatically collect the variants effect size and multiply to the individual genotype, and sum up the score across the genome to get the PGS.

Input files:

- snp.dir: Path to the directory where LD-clumped Effect size files stored
- geno.dir: Path to the directory of individual genotype files which have been converted by plink to human readable file
- out.dir: Path to the directory where the output individual PGS data is stored

The function "run.all.pgs" will return a list with each item is PGS of one the selected phenotypes calculated for the given samples 

<a id="item-pgs-hapmix"></a>
### Run Hapmix

We used "hapmix"<https://www.stats.ox.ac.uk/~myers/HapmixReleasev2/> to call the local ancestry. To run the Hapmix, user need to edit the configuration file ("example.par" in the released package). 

As the input of the "hapmix" should be haplotype data for both the reference samples (here we used 1000G CEU and YRI) and indviduals of which local ancestry to be called, and the site list and order of sites should be perfectly matched between reference and input samples, here we used two script to get the matched haplotype data for both reference and input.

We created a bundle of scripts to conduct this task of preparing input for "hapmix". In our analysis, we used 1000 Genome Phase3 phased haplotype data as reference panel. The raw data can be download using the link here <https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz>. As we only need European and West African subset of samples in our analysis, we also provide access to the processed CEU/YRI dataset using the Dropbox link as follows: <https://www.dropbox.com/scl/fo/bd415as5k46zqiq0y6pn6/AO1F9d9MTX1MBjikroPrvEU?rlkey=fxin4fzjsyh0me643ch90tmpf&dl=0>

User need to align the sites of the variants of input genotype in the same order as the 1000G reference panel. To achieve this, user can use the script "find_intersect_ukb_1kg.r" to find the intersect list of variants between 1000G and UKB:

```r
source('find_intersect_ukb_1kg.r')
rc<-run.chrs()
``` 

The script will call the function to find the intersected variants listt for all the 22 chromosomes.

Input files:

- ukb_snp_list: path to the variants list of the UKB genetic data
- 1kg_snp_list: path to the variants list of the 1000 G data

We provided the compressed site list of UKB <https://www.dropbox.com/scl/fo/9e1xsi4qgyncg4y6yrtip/AFlCMqjFxLCb1WdfcHVraBM?rlkey=6ex94zmh76liacst9m3alhh9k&dl=0> and 1000G here <https://www.dropbox.com/scl/fo/fiyj1157jy84wj8dmc829/AGxLvRaS-6JAsvnvCz2BRK8?rlkey=cm5c698p6gft7vulbhf0wjgrn&dl=0> used in our analysis.

Output:

The intersected site list can be seen at "tmp_dir/" at the current working directory. Here is the copy of which used in our analysis:<https://www.dropbox.com/scl/fo/v8dbaingvkcqtdodny06h/ABuU_gAKg5JfUSpRb-iFDcM?rlkey=k7x45ohhgp936077hpyp87i8i&dl=0>

Once we got intersected variants, the next thing is to retrieve the haplotype from each dataset. User can simply run script "get_1kg_line_number.sh" and "get_ukb_line_number.sh" to populate the matched line number

```bash
for i in {1..22};do
	bash get_1kg_line_number.sh $i
	bash get_ukb_line_number.sh $i
done 
```

These two scripts will create the row number of intersected variants in the original data file at "tmp_1kg_hap" and "tmp_ukb_hap", both at the current working directory, using the input files at "tmp_dir".

Next, using two R scripts "get_ukb_hap_rds.r" and "get_1kg_hap_rds.r" to extract the haplotype information from UKB samples and 1000G samples at intersected sites.

```r
source('get_ukb_hap_rds.r')
run(1:22)
```

and

```r
source('get_1kg_hap_rds.r')
ml<-mclapply(1:22,get.1kg.chr,mc.cores=4)
```

The two scripts will retrive the haplotypes from UKB and 1000G saperately and store the output at two directories: "tmp_ukb_hap" and "tmp_1kg_hap". 

Input files: (phased) vcf files for UKB (converted from chip array data) and 1000G. We provided "bcf" format files for 1000 G data but for UKB, we are not allowed to  provide such individual level data. 1000 data can be accessed here: <https://www.dropbox.com/scl/fo/byn7d8nfk2lw6gkeilttt/ADMDRJAoR5URMTOCYrINsRc?rlkey=z84ndildqo34nx928f3mfxiim&dl=0>. 

Output:

The two scripts will generate per chromosome haplotype data with variants matching the intersection list variants between UKB and 1000 G.

Usually people have more than 1000 samples to call their local ancestries, and current version of "hapmix" doesn't support scalling up the jobs well. Here, we used script "get_imp_site_list.r" to generate input for batch jobs (e.g. 100 samples per job for 1000 samples).

To use the script, simply type the following command in R session:

```r
source('get_imp_site_list.r')
read_chrs<-run.read.chr()
ukb_geno<-run.ukb.geno.chr()
kg_geno<-run.1kg.geno()
kg_hap<-run.1kg.hap()
ukb_sp_chr<-get.ukb.sp.chr()
```

Input files:

- ukb.sites.dir: Path to the folder where the annotation of the variants and individual genotype files are stored.
- GP3.sites.dir: Path to the folder where the 1000G site files stored. You can access the files we used (folder is called "1kg_snp_list") via Dropbox link: <https://www.dropbox.com/scl/fo/fiyj1157jy84wj8dmc829/AGxLvRaS-6JAsvnvCz2BRK8?rlkey=cm5c698p6gft7vulbhf0wjgrn&dl=0>

Output files:

- out.dir: Path to the output site list folder for input indvidiaul data, here we called it "imp_site_list/" and you can access via DB link: <https://www.dropbox.com/scl/fo/0lax966nmkyp3l7pc1sg4/AOtNiQ3y5lghRU03M5NTT_c?rlkey=gpmbvqp72t44lmixg5x7gsygl&dl=0>
- out.geno.dir: Path to the working directory where the "hapmix" will run on

User can use the script "create_hapmix_files_imp.r" to generate the aligned sites between the input haplotype and samples in the reference panel:

```r
source('create_hapmix_files_imp.r')
ra<-run.all()
``` 
Input files:

- tmp_ukb_hap: Path to the directory of UKB temporary haplotype files
- tmp_1kg_hap: Path to the drectory of 1000G temporary haplotypes files

Access link to the folder "tmp_1kg_hap" can be seen above. 

Output:

- hapmix_all_files: Path to the directory of aligned sites files


Running hapmix also required genetic map file for input variants. User can use the script "recomb_rate_imp.r" to generate the recombination map as follows:

```r
source('recomb_rate_imp.r')
rug<-run.gm()
raf<-run.rf()
```
Input files:

- genemap file: genetic map files. User can access the files used in our analysis via Dropbox link: <https://www.dropbox.com/scl/fo/baz2143vj0t3tclilfqu2/AHMFSPGI4dMTLiTYpUS_04w?rlkey=92n2wgnhemmm11qj98ucxix7h&dl=0>
- aligned site list: Path to the folder (called "imp_1kg_hap/") of aligned site list files, which can be accessed using the Dropbox link introduced above.

Output:

The script will generate the recombination rate files (by interpolation) at "imp_1kg_hap/" folder, which will be used when running "hapmix"

**Manual correction**. Two variants "rs12186596" at chromosome 5 and "rs4001921" at chromosome 19 will need to be corrected by script "correct_chr5_chr19_geno.r" (they were duplicated variants which will make "hapmix" terminated unexpectedly)

```r
source('_chr5_chr19_geno.r')
correct.data()
```

If user wanted to run "hapmix" in batches, they can call the R script "generate_hapmix_file.r" to split the input dataset into smaller batches, by running the command below:

```r
source('generate_hapmix_file.r')
run.split()
```

Input files:

- hapmix_all_files: output from script "create_hapmix_files_imp.r"

Output:

- splited individual genotypes files.

Given all the files prepared by the scripts above, user can genrate the ".par" file as configuration for "hapmix" and input forby using the script "prepare_par_file_imp.r"

```r
source('prepare_par_file_imp.r')
run.par()
```
Input files: all the files generated above

Output: Input files and configuration file for "hapmix" at the working directory.

Finally, user can run the hapmix by using the script "run_hapmix_imp.r"

```r
source('run_hapmix_imp.r')
run(1,length(file.seq))
```
The R script will internaly call another bash script: "run_hapmix_imp.sh" to run the (batch) "hapmix" jobs.

Input files: working directory including all the input files and configuration ".par" files

Output files: the same working directory.


<a id="item-pgs-mcpgs"></a>
### Mean-centered Ancestry PGS construction

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

<a id="item-pgs-bin"></a>
### Plot European and African PGS effect size in bins

To visualise the conservativeness of the causal effect for complex traits across individuals with different degree of ancestry, we binned the 
8003 African-British individuals in terms of their African Ancestry (estimated by hapmix), within each bin, we esimate the effect size for EPGS and APGS and plot the weighted mean of the $\beta_E$ and $\beta_A$ across traits.

We used R script "PGS_ancestry_bin_plot.r" to create result input for such bin plots in Main Figure 4c, Extended Data Figure 9, Supplementary Figures 9,12

In the R session, simply source the script and call the function as follows:

```r
source('PGS_ancestry_bin_plot.r')

eff.list<-get.all.dt(pv=0.05)
```
Input files:

- African_ancestry_8003.rds: individual level-data, inferred African ancestry by hapmix
- eu_af_cov.rds: individual level data, covariates files including age, sex and ACs
- eu_af_53_traits_pheno.rds: individual level data, phenotypes file, each phenotype per column
- eu_af_53_traits_pgs_pv0.05.rds: individual level data, pgs file, each PGS per column
- bootstrap_8003.mat.rds: individual level data, fixed (seed) bootstrap order for the 8003 individuals
- ObsEu_covar.rds: individual level data, covariates for Observed EU samples
- ObsEu_53_traits_pheno.rds: individual level data, phenotypes for Observed EU samples
- ObsEu_53_traits_pgs_pv0.05.rds: individual level data, PGSs for Observed EU samples
- ObsEu_bootstrap.rds: individual level data, fixed (seed) boostrap order for the Observed EU samples.

The output from the "eff.list" include the effect size estimation for EPGS and APGS with 95% CI across different ancestry bins.  

Given the output from the script above, user can use the script "plot_betaE_betaA_by_bins.R" to create figure with same style in Main Figure 4c, Extended Data Figure 9 and Supplementary Figure 9, 12, by simply source the script:

```r
source('plot_betaE_betaA_by_bins.R)
plot.all.panels()
```

Input:

- data/mean_af_anc_bins.rds: African Ancestry intervals and the mean African ancestry in each bin
- data/bin.ratio.rds: output from calling function "get.all.dt" from script "PGS_ancestry_bin_plot.r"


Output:

The "plot.all.panels" function will create a figure "ratio.bin.pdf" at the current working directory.

<a id="item-pgs-fst"></a>
### Forest plot for effect size estimation for individual trait #item-pgs-bin

After running "ANCHOR" on individual phenotypes, you may want to visualise the genetic correlation between e.g. European and African across the traits on which you run "ANCHOR". The script introduced here can be used to create figure shown in Main Figure 5 and Supplementary Figures 8, 11, 13-16.

To run this script "forest_ratio_plot.r" which was mainly used to create Main Figure 5, just simply run the command in R session as follows:

```r
source('forest_ratio_plot.r')
plot.multiple.ratio.mainfig.final(pfile="figure5_ratio_plot.pdf")

```  
Input files:

- data/figure5_plot_ratio_dt.rds: ratios with 95%CI for each phenotypes
- data/figure5_plot_ratio_anno.rds: annotation file for each phenotypes

Output:

After calling the function "plot.multiple.ratio.mainfig.final", a figure called "figure5_ratio_plot.pdf" will be created at the current working directory (In this example, it will produce figure 5).

<a id="item-pgs-simu"></a>
### Simulation scripts and plots

**Generate the simulated phenotypes**

The script used for phenotype simulation in GWAS analysis was created by Lino Ferreira @linoferreira <https://github.com/linoferreira>. To generate the simulated phenotypes, please run the R script "07c-simul-make-phen.R" with individual genotype data as input.

**Simulating the phenotypes of African ancestry individuals condition on European**

Given the simulated effect size of variants in European population, we can also simulated the effect size of variants in African population, conditioning on the effect size in European. 

Here the R script "simu_AF_effect_size.r" was used to simulate the effect size for African individuals, given different correlation parameters $\rho$:

```r
source('simu_AF_effect_size.r')
af.phenos<-all.pheno()
```

Input:

- file.dir: Path to the directory in which each of the file containg the effect size (BETA) simulated for each SNP in the European population (1st Column: variant IDs (chr:pos:ref:alt); 2nd Column: effect/alt allele; 3rd Column: BETA value

Output:

User specify the output directory and assign it to variable "out.dir"

The function "all.pheno" will return a list with each item the simulated phenotype for African participants.

**Plots based on simulation results**

After running the "ANCHOR" on simulated phenotypes, user can use the output result to make the plot for Main Figure 4b, Extended Data Figures 8-9, and Supplementary Figures 8-9 and 13-14.

The R script "figure4b.R" was used to create the plot in Main Figure 4b. Just source the script file as follows:

```r
source('figure4b.R')
plot.XY(pfile='figure4b_update.pdf')
```

Input file:

- data/combinned_bin_simu_and_real_traits53_pv0.05_pv1e4.rds: ratio estimation by "ANCHOR" for simulated traits and real tratis

Output figure:

The script will generate a figure named "figure4b_update.pdf" at the current working directory.

User can also use script "EDF8.r" to create figure for Extended Data Figure 8 without any input by using the command below:

```r
source('EDF8.r')

plot.res(pfile='ed8.pdf')
```

The script will create the output figure by calling funciton "plot.res".

User can also use R script "EDF9.r" to create the figure for Extended Data Figure 9 as follows:

```r
source('EDF9.r')
plot.all.panels(pfile='EDF9.pdf')
```

By calling the function "plot.all.panels", the script will create the Extended Data Figure 9 named as "EDF9.pdf".

Input files:

- data/mean_af_anc_bins.rds: African ancestry bins information
- data/simu_24_traits_ratio_bins.rds: output by running "ANCHOR" on 24 simulated traits.

<a id="item-pgs-trio"></a>
### Trio PGS

We identified 2 trios from the 8003 individuals with mixed African ancestry. 

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

