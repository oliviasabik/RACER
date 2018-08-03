# mirrorplot
## 1. Installation
In R, the following commands will install and load mirrorplot:

```
install.packages("devtools") 
library(devtools) 
install_github("oliviasabik/mirrorplot") 
library(mirrorplot)
```

## 2. Input
mirrorplot takes one or two data frames containing, at minimum
three columns,one containing the chromosome a variant is on (CHR), one 
containing the physical location of the variants (CHR_POS)
and one containing the -log10(p-value)s associated with theose variants.
This file can also have the rs_id values and/or the R2 information reflecting
the linkage disequalibrium of the SNPs based on the lead SNP. 
Here is an example file:
``` 
RS_ID	CHR_POS	LD	LOG10P 
rs1	580000	0.2	1.3		
rs2	585000	0.3	1.5  
rs3 	587000 	0.6 	8.0
rs4 	589000	0.8	1.4
```
Additionally, you are required to provide the chromosome you wish to plot, and the method
by which you are plotting the association data, either +/- 50kb of a snp, +/- 50kb of a snp,
or by a set of coordinates you provide. 

## 3. Options
(1) Plotting 
(2) LD information
(3) 

## 4. Examples
(1) Plotting two associations on chromosome 12 by coordinates specified in the command
with default LD information.
```
mirrorplot::mirror_plot_function(eqtl_snps, gwas_snps, chr = 12, plotby = "coord",
                                 "eQTL Data", "GWAS Data", 500000, 650000,
                                 pops = NULL, rs_id_1 = NULL, rs_id_2 = NULL) 
```
(2) Plotting two associations on chromosome 12 by coordinates specific in the command
with LD information from input. Even though SNP IDs were provided, the input had LD
information, so the input LD was used.
```
mirrorplot::mirror_plot_function(b4galnt3_eqtl, b4galnt3_gwas, chr = 12, plotby = "coord",
                                 "B4galnt3 eQTL", "eBMD GWAS", 500000, 650000,
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs6489548", rs_id_2 = "rs6489548")
```
(3) Plotting two associations on chromosome 11 by the location of gene CADM1, using input 
LD information.Even though SNP IDs were provided, the input had LD
information, so the input LD was used.
```
mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas, chr = 11, plotby = "gene", gene = "CADM1",
                                 "Cadm1 eQTL", "eBMD GWAS",
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs2509353", rs_id_2 = "rs2509353")
```
(4)Plotting two associations on chromosome 11 by the location of gene CADM1, using input 
LD information for association one, but calculating LD for association two using 1000 
genomes phase III for association dataset two using rs2509353 as the lead SNP, as the 
second dataset does not have LD information in the input. 
```
mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas_no_ld, chr = 11, plotby = "gene", gene = "CADM1",
                                 "Cadm1 eQTL", "eBMD GWAS",
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs2509353", rs_id_2 = "rs2509353")
```