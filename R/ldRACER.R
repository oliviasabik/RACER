#' Calculating Linkage Disequilibrium Information for Regional Association ComparER
#'
#' This group of functions allows you to creat a plot of -log10(P-values) of an association study by their genomic position, for example, the results of a GWAS or eQTL study. This function takes the rsID of a reference SNP and calculates LD for all other SNPs in the dataset using the 1000 Genomes Phase III Data. This function takes a lead SNP, or finds the most significantly associated SNP in the input data set and use it as the lead SNP (auto_snp = TRUE). The input of the function should already have been formatted using formatRACER().
#' @param assoc_data required. A dataframe produced by by formatRACER()
#' @param rs_col required. numeric or character. index of column or name of column containing rsID numbers for SNPs
#' @param pops required. Populations used to calculate LD. Options can be found at the LD Link website.
#' @param lead_snp required, unless auto_snp = TRUE. Required if ldby = "1000genomes". snp used to calculate LD
#' @param auto_snp optional. default = FALSE, can be set to TRUE to calculate LD using the highest LOG10P SNP as the reference
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @examples
#' \dontrun{
#' data(mark3_bmd_gwas)
#' mark3_bmd_gwas_f = formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3, pos_col = 4, p_col = 11)
#' head(ldRACER(assoc_data = mark3_bmd_gwas_f, rs_col = 5, pops = c("EUR"), lead_snp = "rs11623869"))}

ldRACER <- function(assoc_data, rs_col, pops, lead_snp = NULL, auto_snp = FALSE){

  if(missing(rs_col)){
    stop("Please specify which column contains rsIDs.")
  }else if(missing(pops)){
    stop("Please specify which 1000 Genomes populations to use to calculate LD.")
  }else if(is.null(lead_snp) == TRUE && auto_snp == FALSE){
    stop("Please specify which lead SNP to use to calculate LD, or use auto SNP.")
  }else{
    message("All inputs are go!")
  }

  # read in, format, and filter data sets
  message("Reading in association data...")
  assoc_data <- as.data.frame(assoc_data)
  if(class(rs_col) == "numeric"){
    colnames(assoc_data)[rs_col] = "RS_ID"
  }else if(class(rs_col) == "character"){
    if((rs_col %in% colnames(assoc_data) == TRUE)){
      colnames(assoc_data)[which(colnames(assoc_data) == rs_col)] = "RS_ID"
    }else{
      stop("The rsID column you specified is not in the association data frame.")
    }
  }

  if(auto_snp == TRUE){
    lead_snp = assoc_data[which.max(assoc_data$LOG10P),"RS_ID"]
  }

  pot_pops = c("AFR", "AMR", "EAS", "EUR", "SAS", "YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB", "MXL", "PUR", "CLM", "PEL", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "GIH", "PJL", "BEB", "STU", "ITU")
  if(sum(pops %in% pot_pops == TRUE) == length(pops)){
    message("Populations selected.")
  }else{
    stop("The one or more of the populations you specified are not valid, please consult documentation for list of acceptable population codes.")
  }

  # calculate LD
  message(paste0("Calculating LD using ", lead_snp, "..."))
  assoc_data$LD_BIN = 1
  assoc_data$LD_BIN = NA
  if(length(pops) == 1){
    ld_command = paste0("https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", lead_snp,"&pop=",pops,"&r2_d=r2&token=c0f613f149ab")
    print(ld_command)
    z = as.data.frame(data.table::fread(ld_command))
    z = dplyr::select(z, "RS_Number", "R2")
    colnames(z) = c("RS_ID", "LD")
    assoc_data$LD = NA
    assoc_data = dplyr::select(assoc_data, -"LD")
    message(paste0("Merging input association data with LD..."))
    assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
    assoc_data$LD = as.numeric(as.character(assoc_data$LD))
    assoc_data$LD_BIN <- cut(assoc_data$LD,
                        breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                        labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
    assoc_data$LD_BIN = as.character(assoc_data$LD_BIN)
    assoc_data$LD_BIN[is.na(assoc_data$LD_BIN)] <- "NA"
    assoc_data$LD_BIN = as.factor(assoc_data$LD_BIN)
    assoc_data$LD_BIN = factor(assoc_data$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
    assoc_data$LABEL = NA
    assoc_data[which(assoc_data$RS_ID == lead_snp), which(colnames(assoc_data) == "LABEL")] = "LEAD"
    }else if(length(pops) > 1){
      ld_command = paste0("https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", lead_snp, "&pop=")
      for (i in 1:length(pops)){
        if (i < length(pops)){
          a = pops[i]
          ld_command = paste0(ld_command, a, "%2B")
          }else if(i == length(pops)){
            a = pops[i]
            ld_command = paste0(ld_command, a)
          }
        }
      ld_command = paste0(ld_command, "&r2_d=r2&token=c0f613f149ab")
      print(ld_command)
      z = as.data.frame(data.table::fread(ld_command))
      z = dplyr::select(z, "RS_Number", "R2")
      colnames(z) = c("RS_ID", "LD")
      assoc_data$LD = NA
      assoc_data = dplyr::select(assoc_data, -"LD")
      message(paste0("Merging input association data with LD..."))
      assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
      assoc_data$LD = as.numeric(as.character(assoc_data$LD))
      assoc_data$LD_BIN <- cut(assoc_data$LD,
                          breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                          labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
      assoc_data$LD_BIN = as.character(assoc_data$LD_BIN)
      assoc_data$LD_BIN[is.na(assoc_data$LD_BIN)] <- "NA"
      assoc_data$LD_BIN = as.factor(assoc_data$LD_BIN)
      assoc_data$LD_BIN = factor(assoc_data$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
      assoc_data$LABEL = NA
      assoc_data[which(assoc_data$RS_ID == lead_snp), which(colnames(assoc_data) == "LABEL")] = "LEAD"
    }
  return(assoc_data)
}


