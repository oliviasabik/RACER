setwd("~/Work/package_development/RACER")
library(devtools)
# install.packages("roxygen2")
library(roxygen2)
library(tidyverse)
#install.packages("ggpubr")
#library(ggpubr)
#install.packages("settings")
library(settings)
#devtools::create("mirrorplot")
#install_github("oliviasabik/mirrorplot")

setwd("~/Work/package_development/RACER/")
document()
install("../RACER/")
library(RACER)

devtools::check()

data("mark3_eqtl")
data("mark3_bmd_gwas")
head(mark3_bmd_gwas)
df1_f = formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3, pos_col = 4, p_col = 11,rs_col = 2)
head(mark3_eqtl)
df2_f = formatRACER(assoc_data = mark3_eqtl, chr_col = 10, pos_col = 11, p_col = 7,rs_col = 15)

head(df1_f)
df1_f_ld = ldRACER(assoc_data = df1_f, rs_col = 2, pops = "EUR", lead_snp = "rs11623869")
head(df2_f)
df2_f_ld = ldRACER(assoc_data = df2_f, rs_col = 15, pops = "EUR", lead_snp = "rs11623869")

scatterPlotRACER(assoc_data1 = df1_f_ld, assoc_data2 = df2_f_ld, chr = 14, name1 = "Mark3_GWAS", name2 = "Mark3_eQTL", region_start = 103750000, region_end = 104250000, ld_df = 1)






#### package down
library(pkgdown)
pkgdown::build_site()



source("https://bioconductor.org/biocLite.R")
biocLite("BiocStyle")

usethis::use_package("BiocStyle", "Suggests")
usethis::use_package("knitr", "Suggests")
usethis::use_package("rmarkdown", "Suggests")
usethis::use_package("utils", "Suggests")
usethis::use_data(mark3_eqtl, compress = "xz", overwrite = TRUE)
usethis::use_data(mark3_bmd_gwas, compress = "xz", overwrite = TRUE)
usethis::use_data(biomart_hg19, compress = "xz", overwrite = TRUE)
usethis::use_data(biomart_hg38, compress = "xz", overwrite = TRUE)
usethis::use_pkgdown()

data("mark3_eqtl")
head(mark3_eqtl)

library(rmarkdown)
render(input = "./vignettes/IntroToRACER.Rmd",output_format = "pdf_document", output_file = "IntroToRACER.pdf", output_dir = "./vignettes")


