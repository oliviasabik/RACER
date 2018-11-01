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

###### Testing Single Plot Funcation #####
single_plot_function(assoc_data = B4galnt3_eqtl_form_ld, chr = 12, plotby = "coord", start_plot = 470000, end_plot = 700000)

mirrorplot::single_plot_function(assoc_data = B4galnt3_no_ld_formatted, chr = 12, build = "hg19", plotby = "gene", gene_plot = "B4GALNT3")

mirrorplot::single_plot_function(assoc_data = B4galnt3_no_ld_formatted_B, chr =  12, plotby = "gene", gene_plot = "B4GALNT3")
#Success

###### Testing the mirrorplot function #####
head(b4galnt3_gwas_form)
head(B4galnt3_eqtl_form_ld)

mirror_plot_function(assoc_data1 = B4galnt3_eqtl_form_ld, assoc_data2 = b4galnt3_gwas_form, chr = 12, name1 = "eqtl", name2 = "gwas", plotby = "gene", gene_plot = "B4GALNT3")


data("mark3_bmd_gwas")


#### package down
source("https://bioconductor.org/biocLite.R")
biocLite("BiocStyle")

usethis::use_package("BiocStyle", "Suggests")
usethis::use_package("knitr", "Suggests")
usethis::use_package("rmarkdown", "Suggests")


