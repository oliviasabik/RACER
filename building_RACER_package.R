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

devtools::build("../RACER")

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

df2_f_ld$LOG10P = rnorm(n = 9789, mean = 3, sd = 2)
scatterPlotRACER(assoc_data1 = df1_f_ld, assoc_data2 = df2_f_ld, chr = 14, name1 = "Mark3_GWAS", name2 = "Mark3_eQTL", region_start = 103750000, region_end = 104250000, ld_df = 1)

mirrorPlotRACER(assoc_data1 = df1_f_ld, assoc_data2 = df2_f_ld, chr = 14, name1 = "Mark3_GWAS", name2 = "Mark3_eQTL", plotby = "coord", start_plot = 103750000, end_plot = 104250000)

len1 = nchar(trunc(max(df1_f$LOG10P)))
len2 = nchar(trunc(max(df2_f$LOG10P)))



x = rnorm(n = 50, mean = 4, sd = 1)
y = rnorm(n = 50, mean = 10, sd = 2)
z = rnorm(n = 50, mean = 100, sd = 2)
scaleFUN0 <- function(x) sprintf("%.0f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)
scaleFUN2 <- function(x) sprintf("%.2f", x)
scaleFUN3 <- function(x) sprintf("%.3f", x)
scaleFUN4 <- function(x) sprintf("%.4f", x)

x2 = data.frame(x,x)
a = ggplot(data = x2, aes(x = x, y = x.1)) + geom_point() + scale_y_continuous(labels = scaleFUN3)
a
y2 = data.frame(x,y)
b = ggplot(data = y2, aes(x = x, y = y)) + geom_point() + scale_y_continuous(labels = scaleFUN2)
b
z2 = data.frame(x,z)
c = ggplot(data = z2, aes(x = x, y = z)) + geom_point() + scale_y_continuous(labels = scaleFUN1)
ggpubr::ggarrange(a,b,c, nrow = 3)

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



