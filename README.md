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
## 3. Examples
```
mirrorplot::mirror_plot_function(b4galnt3_eqtl, b4galnt3_gwas, chr = 12, plotby = "coord",
                                 "B4galnt3 eQTL", "eBMD GWAS", 500000, 650000,
                                 pops = NULL, rs_id_1 = NULL, rs_id_2 = NULL) 

mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas, chr = 11, plotby = "coord",
                                 "Cadm1 eQTL", "eBMD GWAS", 115390000, 115600000,
                                 pops = NULL, rs_id_1 = NULL, rs_id_2 = NULL)
 
mirrorplot::mirror_plot_function(b4galnt3_eqtl, b4galnt3_gwas, chr = 12, plotby = "coord",
                                 "B4galnt3 eQTL", "eBMD GWAS", 500000, 650000,
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs6489548", rs_id_2 = "rs6489548")

mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas, chr = 11, plotby = "snp",
                                 "Cadm1 eQTL", "eBMD GWAS",
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs2509353", rs_id_2 = "rs2509353")

mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas, chr = 11, plotby = "gene", gene = "CADM1",
                                 "Cadm1 eQTL", "eBMD GWAS",
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs2509353", rs_id_2 = "rs2509353")

mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas_no_ld, chr = 11, plotby = "gene", gene = "CADM1",
                                 "Cadm1 eQTL", "eBMD GWAS",
                                 pops = c("CEU","TSI","FIN","GBR","IBS"),
                                 rs_id_1 = "rs2509353", rs_id_2 = "rs2509353")
