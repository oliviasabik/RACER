# genCom//mirrorplot
# Introduction
gencom is a package that allows you to generate high quality regional association plots from genomic data using the function gencomSinglePlot() and can stack two such plots using the function gencomMirrorPlot() such that they are mirrored across a common x-axis, enabling the direct visual comparison of two associations at the same locus. For example, in order to link the expression of a specific gene with changes in a phenotype it is common to identify colocalizing pairs of expression quantitative trait loci (eQTL) and quantitative trait loci (QTL) from genome-wide association studies (GWAS). The coloc package from Wallace, Giambartolomei, and Plagnol can be used to calculate the posterior probability of two association signals sharing a common genetic impetus, however the Mirror Plot was designed to provide a visual representation of the relationship between the two associations. Additionally, the gencomSinglePlot() function can be used to generate a regional association plot for one association data set. gencom also incoroporates two helper functions: gencomFormat() which helps you to format your data for plotting and gencomLD() which helps you to format LD information included in your dataset, or can use input rsID numbers to calculate LD using the 1000 genomes phase III data. 

## 1. Installation
In R, the following commands will install and load mirrorplot:
```
install.packages("devtools") 
library(devtools) 
install_github("oliviasabik/mirrorplot") 
library(mirrorplot)
```
The package can also be installed by using the install command in R to install the 
package from the directory it was downloaded into:
```
install.packages("devtools") 
library(devtools)
install("{PATH}/mirrorplot-master/")
library(mirrorplot)
```
## 2. genCom//mirrorplot vignette
The best way to get a feel for genCom/mirrorplot is to refer to the vignette included in the package. 
