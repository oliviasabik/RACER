# RACER//Regional Association ComparER
# Introduction
RACER is a package that allows you to generate high quality regional association plots from genomic data using the function singlePlotRACER() and can stack two such plots using the function mirrorPlotRACER() such that they are mirrored across a common x-axis, enabling the direct visual comparison of two associations at the same locus. For example, in order to link the expression of a specific gene with changes in a phenotype it is common to identify colocalizing pairs of expression quantitative trait loci (eQTL) and quantitative trait loci (QTL) from genome-wide association studies (GWAS). The coloc package from Wallace, Giambartolomei, and Plagnol can be used to calculate the posterior probability of two association signals sharing a common genetic impetus, however the Mirror Plot was designed to provide a visual representation of the relationship between the two associations. RACER also incoroporates two helper functions: formatRACER() which helps you to format your data for plotting and ldRACER() which helps you to format LD information included in your dataset, or can use input rsID numbers to calculate LD using the 1000 genomes phase III data. 

## 1. Installation
In R, the following commands will install and load RACER:
```
install.packages("devtools") 
library(devtools) 
install_github("oliviasabik/RACER") 
library(RACER)
```
The package can also be installed by using the install command in R to install the package from the directory it was downloaded into:
```
install.packages("devtools") 
library(devtools)
install("{PATH}/RACER-master/")
library(RACER)
```
## 2. RACER vignette
The best way to get a feel for RACER is to refer to the [vignette](https://oliviasabik.github.io/RACER/articles/IntroToRACER.html) included in the package. 
