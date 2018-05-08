# mirrorplot
## 1. Installation
In R or RStudio, these commands will install and load mirrorplot

```
install.packages("devtools") 
library(devtools) 
install_github("oliviasabik/mirrorplot") 
library(mirrorplot)
```

## 2. Input
mirrorplot takes one or two data frames containing, at minimum
two columns, one containing the physical location of the variants (CHR_POS)
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
## 3. Example
Under Construction
