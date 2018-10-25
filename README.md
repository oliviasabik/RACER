# mirrorplot
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
## 2. Input
mirror_plot() takes two data frames containing, at minimum
three columns, one containing the chromosome a variant is on (CHR), one 
containing the physical location of the variants (CHR_POS)
and one containing the -log10(p-value)s associated with theose variants (LOG10P).
This file can also have the rs_id values (RS_ID) to use to calculates linkage disequilibrium
and/or the R2 information reflecting the linkage disequalibrium of the SNPs based on 
the lead SNP. 

You can also plot one association using the single_plot() function. These have the same
required columns as described for for the mirror_plot() function above. 

Here is an example file:
``` 
RS_ID	CHR	CHR_POS	LD	LOG10P 
rs1	1	580000	0.2	1.3		
rs2	1	585000	0.3	1.5  
rs3	1	587000	0.6	8.0
rs4	1	589000	0.8	1.4
```
Additionally, you are required to provide:

(1) the chromosome you wish to plot ex: chr = 3

(2) the genome build of the data you are plotting, either "hg19" or "hg38"

(3) the methodby which you are plotting the association data  
- (1) either +/- 50kb of a gene, ex: plotby = "gene", gene_plot = "GENE_NAME",  
- (2) +/- 50kb of a snp, ex: plotby = "snp", snp_plot = "RS_#",  
- (3) or by a set of coordinates you provide, ex: plotby = "coord", start_plot = 100,000, plot_end = 103,000  

(4) the method by which you want to want to gather LD information  
- (1) either not include LD data in the plot (default), ex: ldby = "none",  
- (2) from the input data, ex: ldby = "input" and input data frame has a column named LD,  
- (3) or from the 1000 Genomes Phase III database, ex: ldby = "1000genomes", snp_ld_1 = "rs#", pops = "EUR"  
	If you want to calculate LD for the same SNP in both plots, just specify snp_ld_1, 
	but if you want a different lead SNP to be used for each plot, specify snp_ld_1 and snp_ld_2. 
	Finally, if you want to use 1000 Genomes to source LD information, you also need to specify 
	the populations you want to use to. List whatever subset of populations you wish. Pops include:
``` 
	(AFR) African
		(YRI) Yoruba in Ibadan, Nigera
		(LWK) Luhya in Webuye, Kenya
		(GWD) Gambian in Western Gambia
		(MSL) Mende in Sierra Leone
		(ESN) Esan in Nigera
		(ASW) Americans of African Ancestry in SW USA
		(ACB) African Carribbeans in Barbados
	(AMR) Ad Mixed American
		(MXL) Mexican Ancestry from Los Angeles, USA
		(PUR) Puerto Ricans from Puerto Rico
		(CLM) Colombians from Medellin, Colombia
		(PEL) Peruvians from Lima, Peru
	(EAS) East Asian
		(CHB) Han Chinese in Bejing, China
		(JPT) Japanese in Tokyo, Japan
		(CHS) Southern Han Chinese
		(CDX) Chinese Dai in Xishuangbanna, China
		(KHV) Kinh in Ho Chi Minh City, Vietnam
	(EUR) European
		(CEU) Utah Residents from North and West Europe
		(TSI) Toscani in Italia
		(FIN) Finnish in Finland
		(GBR) British in England and Scotland
		(IBS) Iberian population in Spain
	(SAS) South Asian
		(GIH) Gujarati Indian from Houston, Texas
		(PJL) Punjabi from Lahore, Pakistan
		(BEB) Bengali from Bangladesh
		(STU) Sri Lankan Tamil from the UK
		(ITU) Indian Telugu from the UK
``` 

## 3. Examples
(1) Plotting two associations on chromosome 12, plotted by coordinates specified in the command
with default LD information.
```
mirrorplot::mirror_plot_function(b4galnt3_eqtl, b4galnt3_gwas, chr = 12, plotby = "coord",
                                 "B4galnt3 eQTL", "eBMD GWAS", 
                                 start_plot = 500000, end_plot = 650000, ldby = "input")

```
(2) Plotting two associations on chromosome 12 by coordinates specific in the command
with LD information from the 1000 Genomes Phase III Database, calculated from the European
population, for SNP rs6489548. 
```
mirrorplot::mirror_plot_function(b4galnt3_eqtl, b4galnt3_gwas, chr = 12, plotby = "coord",
                                 start_plot = 500000, end_plot = 650000,
                                 "B4galnt3 eQTL", "eBMD GWAS",
                                 ldby = "1000genomes", pops = "EUR", snp_ld_1 = "rs6489548")
```
(3) Plotting two associations on chromosome 11 by the location of gene CADM1, with LD 
information from the 1000 Genomes Phase III Database, calculated from the Ad Mixed
populations, for SNP rs2509353.
```
mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas, chr = 11, 
				plotby = "gene", gene_plot = "CADM1", "Cadm1 eQTL", "eBMD GWAS",
                                ldby = "1000genomes", 
                                pops = c("MXL","PUR","CLM","PEL"), 
                                snp_ld_1 = "rs2509353")

```
(4)Plotting two associations on chromosome 11 by the location of gene CADM1, with no 
LD information included. 
```
mirrorplot::mirror_plot_function(cadm1_eqtl, cadm1_gwas, chr = 11, 
				plotby = "gene", gene = "CADM1","Cadm1 eQTL", "eBMD GWAS",
                                ldby = "none")
```
(5)Plotting one association on chromosome 11 by the location of gene CADM1, with no 
LD information included. All of the above options are also available in single_plot_function.
```
mirrorplot::single_plot_function(cadm1_gwas, chr = 11, 
				plotby = "gene", gene = "CADM1","eBMD GWAS", ldby = "none")
```
