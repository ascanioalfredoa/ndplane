# Code to replicate Niche Divergence Plane analysis, submited to Ecological Monographs ![DOI in Zenodo badge](https://zenodo.org/badge/.svg "Zenodo DOI")

This repository contains data, instructions, and code for the paper
> Ascanio, A., Bracken, J.T., Stevens, H., and Jezkova, T. (20xx). *New theoretical and analytical framework for quantifying and classifying ecological niche differentiation* Submitted to Ecological Monographs. <DOI>

## How to cite

Please, cite this compendium as:
> Ascanio, A., Bracken, J.T., Stevens, H., and Jezkova, T. (20xx). *New theoretical and analytical framework for quantifying and classifying ecological niche differentiation - code companion* Submitted to Ecological Monographs. <DOI>

## How to download

You can download the compendium as a zip from from this URL:
<https://github.com/luismmontilla/biomass/archive/master.zip>

## Folder organization and Methods

This repository contains three main subdirectories. Each of these directories represents a section of our analysis: 1) Theoretical simulations on the Niche Divergence Plane; 2) Virtual Species simulations on the Niche Divergence Plane; and 3) Niche Divergence Plane application using _Ambystoma maculatum_ and _A. opacum_ data

### Theoretical simulations on the Niche Divergence Plane (/NDP_Simulations)

This directory contains the scripts to perform theoretical simulations following comparisons of 4-parameter Beta Function curves. These scripts are contained in /NDP_Simulations and should be able to run without any external data.


#### Scripts and order of execution

1. **Simulate_NDP_Space_Parallel.R**: First, this scripts creates a parameter space for the 4-parameter Beta functions curve. Second, it contains a parallelized for-loop that compares all possible Beta function curves from the previously created parameter space. Niche Divergence Plane indices are derived from these comparisons (Niche Exclusivity and Niche Similarity). Each Beta curve creates a CSV file as an output, containing the Niche Divergence Plane indices and the parameters of each curve used for comparison. R functions for estimating the shape of the Beta function curves are contained in **BetaFunctions.R**
2. **sim_merge_plots.R**: This scripts merges all CSV outputs from the previous script and combines them into a single data frame. This data frame is used to plot **Figure 3**

These simulations are not heavy themselves, but their number increases rapidly with the number of parameter values added to the 4-parameter Beta curves used for comparison. For these reasons, a set of R and Shell scripts was used to run these simulations in the RedHawk Cluster (slurm) from Miami University. These scripts can be found in /NDP_Simulations/Batch_version

### Virtual Species simulations on the Niche Divergence Plane (/NDP_VirtualSpecies)

This directory contains the scripts to perform virtual species simulations and comparisons using the Niche Divergence Plane, and further comparisons to previous measurements of Niche Similarity (Hellinger's I and Schoener's D). 

#### Data and settings

**Before running**: repository filesize limtis may prevent this script from executing, or it may take long to download. Before running, make sure that:

1. Worldclim data at a resolution of 5 arcminutes is downloaded in the folder /data/wc. If this data does not exist in the working directory, one of the scripts will download it.
2. Make sure you set your working directory in each script.

These simulations can be more computationally intensive than the ones in the previous section. For these reasons, we used a smaller parameter space for our Beta function curves and prepare the scripts to be automatically ran in batch (for parallel processing in the RedHawk cluster at Miami University). 

#### Scripts and order of execution

 Located in /NDP_VirtualSpecies/rscripts

1. **set_parameter_space.R**: Run this script first. This will ensure you have the environmental data from WorldClim and will also create the parameter space used to define the responses of the virtual species, based on Beta-function curves.
2. **virtualsp_simulations_batch_parallel_loop.R**: This script will produce virtual species based on Beta-function curves, provided a given set of parameter combinations (from the previous script), and the environmental data in which the virtual species responses are defined (BIO5 and BIO12). Make sure you set the working directory at the beginning of the script, and at the beginning of the parallel for-loop. This script relies on functions sourced from **BetaFunctions.R** to calculate the shape of every pair of virtual species response curves, and from **virtualsp_functions.R** to calculate Hellinger's I and Schoener's D niche similarities. Similar to the simulations in the previous section, this script produces a CSV for each response curve being compared to all others. 
3. **ggplot_virtualspecies_output.R**: This script merges the previous CSV outputs into one, and produces different exploratory plots, including **Figure 4** and supplementary **Figure S1**. It also performs the correlation assessment and statistical test that resulted in **Table 3**.


### Niche Divergence Plane application on _Ambystoma maculatum_ and _A. opacum_ (/NDP_Salamanders)

#### Data

**Before running**: Due to github limits on filesize, environmental information was not uploaded to the repository. Follow these instructions first:

1. Download the 19 bioclimatic variables from CHELSA and add them to /NDP_Salamanders/Data/CHELSA_v2_1 [https://chelsa-climate.org/bioclim/]
2. Download soil data was originally downloaded from the Gridded Soil Survey Database for the Conterminous United States (CONUS) [https://www.nrcs.usda.gov/resources/data-and-reports/gridded-national-soil-survey-geographic-database-gnatsgo]. 
	1. The repository contains an .RDS file with the already extracted soil varaiables for every pixel of the US. This file weights ~100 MB.
	2. You need to download the gNATSGO CONUS raster (~30x30 m cell size, ~6 GB) and add it to the /NDP_Salamanders/Data/NATSGO. The script was made to read that raster as a TIF, but any other raster file format should work.

Ambystoma occurrences were extracted from different sources, such as GBIF, HerpMapper, NA-HERPS, and personal communications with other researchers. Occurrences are located in /NDP_Salamanders/Data/Ambystoma

#### Scripts and order of execution

 Located in /NDP_Salamanders/Rscripts

1. **salamander_thinning.R**: With this scripts we thin the species occurrences to 1 in a 10 km radius. This scripts uses the function *thin_records()* located in the script **thin_records.R**. This function is a modified version of the thinning algorithm used in the R package spThin (Aiello-Lammens et al. 2015). From this script, we obtain a thinned occurrence layer only for _Ambystoma opacum_ and _A. maculatum_.
2. **sal_extract_env.R**: This script takes the previously thinned occurrences for both salamander species and extracts the relevant environmental variables for each of them. This scripts relies on having the CHELSA and NATSGO data.
3. **sal_NDP.R**: This scripts performs the Niche Divergence Plane analysis between both salamander species and produces **Figure 5** from the paper.
4. **sal_responses.R** This scripts produces supplementary **Figure S2**, containing all response curves used to construct the previous niche divergence plane.

## References

Aiello-Lammens, M. E., Boria, R. A., Radosavljevic,
  A. , Vilela, B. and Anderson, R. P. (2015). spThin:
  an R package for spatial thinning of species
  occurrence records for use in ecological niche
  models. Ecography, 38: 541-545. URL
  https://onlinelibrary.wiley.com/doi/10.1111/ecog.01132.