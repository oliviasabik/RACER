setwd("~/Work/RACER/")
library(devtools)
# install.packages("roxygen2")
library(roxygen2)
library(tidyverse)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("settings")
library(settings)
#devtools::create("mirrorplot")
#install_github("oliviasabik/mirrorplot")

setwd("~/Work/RACER/")
document()
install("../RACER/")
library(RACER)
devtools::check()

usethis::use_data(mark3_eqtl, compress = "xz", overwrite = TRUE)
usethis::use_data(mark3_bmd_gwas, compress = "xz", overwrite = TRUE)
usethis::use_data(biomart_hg19, compress = "xz", overwrite = TRUE)
usethis::use_data(biomart_hg38, compress = "xz", overwrite = TRUE)



library(pkgdown)
pkgdown::build_site()

#toy_pvalues = read_csv("./toy_pvalues.csv")
# mirrorplot::single_plot_function(toy_pvalues)
load("../osteoblast_project/figure_data/b4galnt3_eqtl")
load("../osteoblast_project/figure_data/b4galnt3_gwas")
colnames(b4galnt3_eqtl) = c("RS_ID", "GENE_ID", "CHR", "CHR_POS", "ALLELE1",
                            "ALLELE0", "BUILD", "TSS_DIST", "MA_SA", "MA_CO",
                            "MAF", "P", "SLOPE", "SLOPE_SE", "LD")

colnames(b4galnt3_gwas) = c("RS_ID", "SNP_ID", "CHR", "CHR_POS", "ALLELE1",
                            "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P",
                            "N", "LD")

b4galnt3_eqtl$LOG10P = -log10(b4galnt3_eqtl$P)
b4galnt3_gwas$LOG10P = -log10(b4galnt3_gwas$P)
b4galnt3_gwas = b4galnt3_gwas[,-c(14,15)]

load("../osteoblast_project/figure_data/cadm1_eqtl")
load("../osteoblast_project/figure_data/cadm1_gwas")
colnames(cadm1_eqtl) = c("RS_ID", "GENE_ID", "CHR", "CHR_POS", "ALLELE1",
                            "ALLELE0", "BUILD", "TSS_DIST", "MA_SA", "MA_CO",
                            "MAF", "P", "SLOPE", "SLOPE_SE", "LD")

colnames(cadm1_gwas) = c("RS_ID", "SNP_ID", "CHR", "CHR_POS", "ALLELE1",
                            "ALLELE0", "A1FREQ", "INFO", "BETA", "SE", "P",
                            "N", "LD")

cadm1_eqtl$LOG10P = -log10(cadm1_eqtl$P)
cadm1_gwas$LOG10P = -log10(cadm1_gwas$P)

b4galnt3_gwas_no_ld = b4galnt3_gwas[,-13]
cadm1_gwas_no_ld = cadm1_gwas[,-13]



##### Testing gencom_format #####
b4galnt3_gwas_form = formatRACER(assoc_data = b4galnt3_gwas, chr_col = 3, pos_col = 4, log10p_col = 14, ld_col = 13)
b4galnt3_gwas_form = gencom_format(assoc_data = b4galnt3_gwas, pos_col = 4, log10p_col = 14, ld_col = 13)


colnames(b4galnt3_gwas_no_ld)
b4galnt3_gwas_no_ld_log10p = b4galnt3_gwas_no_ld[,-13]
B4galnt3_no_ld_formatted_A = mirrorplot::gencom_format(assoc_data = b4galnt3_gwas_no_ld, chr_col = 3, pos_col = 4, log10p_col = 13)
B4galnt3_no_ld_formatted_B = mirrorplot::gencom_format(assoc_data = b4galnt3_gwas_no_ld_log10p, chr_col = 3, pos_col = 4, p_col = 11)

head(b4galnt3_eqtl)
colnames(b4galnt3_eqtl)
B4galnt3_eqtl_form = mirrorplot::gencom_format(assoc_data = b4galnt3_eqtl, chr_col = 3, pos_col = 4, log10p_col = 16, ld_col = 15)
head(B4galnt3_eqtl_form)
# Success

##### Testing gencom_ld #####
B4galnt3_eqtl_form = B4galnt3_eqtl_form[,-c(15,17)]
head(B4galnt3_eqtl_form)
B4galnt3_eqtl_form_ld = gencom_ld(assoc_data = B4galnt3_eqtl_form, pops = "EUR", lead_snp = "rs6489548")
B4galnt3_eqtl_form_ld = gencom_ld(assoc_data = B4galnt3_eqtl_form, rs_col = 1, lead_snp = "rs6489548")
B4galnt3_eqtl_form_ld = gencom_ld(assoc_data = B4galnt3_eqtl_form, rs_col = 1, pops = "EUR")

head(B4galnt3_eqtl_form_ld)
#Success
data("mark3_bmd_gwas")
data("mark3_eqtl")
mark3_bmd_gwas_f = formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3, pos_col = 4, p_col = 11)
mark3_bmd_gwas_f_ld = ldRACER(assoc_data = mark3_bmd_gwas_f, rs_col = 2, pops = c("EUR"), lead_snp = "rs11623869")
singlePlotRACER(assoc_data = mark3_bmd_gwas_f_ld, chr = 14, build = "hg19", plotby = "coord", start_plot = 103500000, end_plot = 104500000)


#### package down
source("https://bioconductor.org/biocLite.R")
biocLite("BiocStyle")

usethis::use_package("BiocStyle", "Suggests")
usethis::use_package("knitr", "Suggests")
usethis::use_package("rmarkdown", "Suggests")
usethis::use_package("utils", "Suggests")
usethis::use_pkgdown()

data("mark3_eqtl")
head(mark3_eqtl)

library(rmarkdown)
render(input = "./vignettes/IntroToRACER.Rmd",output_format = "pdf_document", output_file = "IntroToRACER.pdf", output_dir = "./vignettes")


