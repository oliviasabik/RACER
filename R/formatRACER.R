#' Formating Data for Regional Association ComparER
#'
#' This group of functions allows you to creat a plot of -log10(P-values) of an association study by their genomic position, for example, the results of a GWAS or eQTL study. This function carries out the formatting necessary to generate a plot
#' @param assoc_data required. A dataframe that has columns containing the chromosome, physical position, and p-values or -log10(p-values) of the association, and can optionally have columns containing R2 information for LD in the region, or rsID numbers for the associated SNPs
#' @param chr_col required. numeric. index of column in assoc_data containing chromosome information
#' @param pos_col required. numeric. index of column in assoc_data containing genomic position information
#' @param log10p_col required, if no p_col specified. numeric. index of column in assoc_data containing -log10(p-value)s
#' @param p_col required, if no log10p_col column specified. numeric. index of column in assoc_data containing p-values
#' @param ld_col optional. numeric. Required if you want to use the LD data in your data set in your plot, index of column in assoc_data containing LD information, e.g. R2 or D' values
#' @param rs_col optional. numeric. Required if you want to use the use
#' ldRACER to pull LD information from the 1000 genomes phase III project,
#' or if you want to make a scatter comparison plot
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @examples
#' data(mark3_bmd_gwas)
#' formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3, pos_col = 4, p_col = 11, rs_col = 2)

formatRACER <- function(assoc_data, chr_col, pos_col, p_col=NULL, log10p_col=NULL, ld_col=NULL, rs_col = NULL){

  if(missing(chr_col)){
    stop("Please specify which column contains chromosome information.")
  }else if(missing(pos_col)){
    stop("Please specify which column contains genomic position information.")
  }else if(is.null(log10p_col) && is.null(p_col)){
    stop("Please specify which column contains p-values or -log10(p-values).")
  }else{
    message("All inputs are go!")
  }

  message("Formating association data...")
  colnames(assoc_data)[chr_col] = "CHR"
  colnames(assoc_data)[pos_col] = "POS"
  message("Processing -log10(p-values)...")
  if(!is.null(log10p_col)){
    colnames(assoc_data)[log10p_col] = "LOG10P"
  }else if(!is.null(p_col)){
    colnames(assoc_data)[p_col] = "P"
    assoc_data$LOG10P = -log10(assoc_data$P)
  }
  if(!is.null(rs_col)){
    colnames(assoc_data)[rs_col] = "RS_ID"
  }
  if(!is.null(ld_col)){
    message("Processing input LD information...")
    colnames(assoc_data)[ld_col] = "LD"
    assoc_data$LD = as.numeric(as.character(assoc_data$LD))
    assoc_data$LD_BIN <- cut(assoc_data$LD,
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                      labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
    assoc_data$LD_BIN = as.character(assoc_data$LD_BIN)
    assoc_data$LD_BIN[is.na(assoc_data$LD_BIN)] <- "NA"
    assoc_data$LD_BIN = as.factor(assoc_data$LD_BIN)
    assoc_data$LD_BIN = factor(assoc_data$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
  }

  # read in, format, and filter data sets
  message("Preparing association data...")
  assoc_data <- as.data.frame(assoc_data)
  assoc_data$POS = as.numeric(as.character(assoc_data$POS))
  assoc_data$LOG10P = as.numeric(as.character(assoc_data$LOG10P))
  assoc_data$CHR = as.numeric(as.character(assoc_data$CHR))

  return(assoc_data)
}

