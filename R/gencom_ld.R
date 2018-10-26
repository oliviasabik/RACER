#' Genomic Comparison Linkage Disequilibrium Calculation Function
#'
#' This group of functions allows you to creat a plot of -log10(P-values) of an association study by their genomic position, for example, the results of a GWAS or eQTL study. This function takes the rsID of a reference SNP and calculates LD for all other SNPs in the dataset using the 1000 Genomes Phase III Data. The input of the function should already have been formatted using gencom_format
#' @param assoc_data required. A dataframe produced by by gencom_format
#' @param rs_col required. numeric. index of column containing rsID numbers for SNPs
#' @param pops required. Populations used to calculate LD.
#' @param lead_snp required. Required if ldby = "1000genomes". snp used to calculate LD
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @examples
#' gencom_ld(assoc_data = assoc_data, rs_col = 5, pops = c("POP1", "POP2", etc...), lead_snp = "RS_ID")

gencom_ld <- function(assoc_data, rs_col, pops, lead_snp){

  if(missing(rs_col)){
    message("Please specify which column contains rsIDs.")
  }else if(missing(pops)){
    message("Please specify which 1000 Genomes populations to use to calculate LD.")
  }else if(missing(lead_snp)){
    message("Please specify which lead SNP to use to calculate LD.")
  }else{
    message("All inputs are go!")
  }

  `%>%` <- magrittr::`%>%`

  # read in, format, and filter data sets
  message("Reading in association data...")
  assoc_data <- as.data.frame(assoc_data)
  colnames(assoc_data)[rs_col] = "RS_ID"
  assoc_data$RS_ID = as.numeric(as.character(assoc_data$RS_ID))

  # calculate LD
  message(paste0("Calculating LD using", lead_snp, "..."))
  assoc_data$LD_BIN = 1
  assoc_data$LD_BIN = NA
  if(length(pops) == 1){
    ld_command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", lead_snp,
                                "&pop=",pops,"&r2_d=r2'")
    z = system(ld_command, intern = TRUE)
    z = as.data.frame(z)
    z = tidyr::separate(z, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                                 "col_9", "col_10"), sep = "\t")
    colnames(z) = z[1,]
    z = z[-1,]
    z = dplyr::select(z, RS_Number, R2)
    colnames(z) = c("RS_ID", "LD")
    assoc_data$LD = NA
    assoc_data = dplyr::select(assoc_data, -LD)
    assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
    assoc_data$LD = as.numeric(as.character(assoc_data$LD))
    assoc_data$LD_BIN <- cut(assoc_data$LD,
                        breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                        labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
    assoc_data$LD_BIN = as.character(assoc_data$LD_BIN)
    assoc_data$LD_BIN[is.na(assoc_data$LD_BIN)] <- "NA"
    assoc_data$LD_BIN = as.factor(assoc_data$LD_BIN)
    assoc_data$LD_BIN = factor(assoc_data$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
    }else if(length(pops) > 1){
      ld_command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", lead_snp,
                                "&pop=")
      for (i in 1:length(pops)){
        if (i < length(pops)){
          a = pops[i]
          ld_command = paste0(ld_command, a, "%2B")
          }else if(i == length(pops)){
            a = pops[i]
            ld_command = paste0(ld_command, a)
          }
        }
      ld_command = paste0(ld_command, "&r2_d=r2'")
      z = system(ld_command, intern = TRUE)
      z = as.data.frame(z)
      z = tidyr::separate(z, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                                 "col_9", "col_10"), sep = "\t")
      colnames(z) = z[1,]
      z = z[-1,]
      z = dplyr::select(z, RS_Number, R2)
      colnames(z) = c("RS_ID", "LD")
      assoc_data$LD = NA
      assoc_data = dplyr::select(assoc_data, -LD)
      assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
      assoc_data$LD = as.numeric(as.character(assoc_data$LD))
      assoc_data$LD_BIN <- cut(assoc_data$LD,
                          breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                          labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
      assoc_data$LD_BIN = as.character(assoc_data$LD_BIN)
      assoc_data$LD_BIN[is.na(assoc_data$LD_BIN)] <- "NA"
      assoc_data$LD_BIN = as.factor(assoc_data$LD_BIN)
      assoc_data$LD_BIN = factor(assoc_data$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
    }
  return(assoc_data)
}


