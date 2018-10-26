#' Genomic Comparsion Formating Function
#'
#' This group of functions allows you to creat a plot of -log10(P-values) of an association study by their genomic position, for example, the results of a GWAS or eQTL study. This function carries out the formatting necessary to generate a plot
#' @param assoc_data required. A dataframe that has columns containing the chromosome, physical position, and p-values or -log10(p-values) of the association, and can optionally have columns containing R2 information for LD in the region, or rsID numbers for the associated SNPs
#' @param chr_col required. numeric. index of column in assoc_data containing chromosome information
#' @param pos_col required. numeric. index of column in assoc_data containing genomic position information
#' @param log10p_col required, if no p_col specified. numeric. index of column in assoc_data containing -log10(p-value)s
#' @param p_col required, if no log10p_col column specified. numeric. index of column in assoc_data containing p-values
#' @param rs_col optional. numeric. Required if you want to calculated LD using 1000 genomes for analysis, index of column containing rsID numbers for SNPs
#' @param ld_col optional. numeric. Required if you want to use the LD data in your data set in your plot, index of column in assoc_data containing LD information, e.g. R2 or D' values
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @examples
#' gencom_format(assoc_data = assoc_data, chr_col = 1, pos_col = 2, log10p_col = 3, ld_col = 4)

gen_comp_format <- function(assoc_data, chr_col, pos_col, p_col=NULL, log10p_col=NULL, ld_col=NULL){

  if(missing(chr_col)){
    message("Please specify which column contains chromosome information.")
  }else if(missing(pos_col)){
    message("Please specify which column contains genomic position information.")
  }else if(log10p_col == NULL & p_col == NULL){
    message("Please specify which column contains p-values or -log10(p-values).")
  }else{
    message("All inputs are go!")
  }

  message("Formating association data...")
  colnames(assoc_data)[chr_col] = "CHR"
  colnames(assoc_data)[pos_col] = "CHR_POS"
  message("Processing -log10(p-values)...")
  if(log10p_col != NULL){
    colnames(assoc_data)[log10p_col] = "LOG10P"
  }else if(p_col != NULL){
    colnames(assoc_data)[p_col] = "P"
    assoc_data$LOG10P = -log10(assoc_data$P)
  }
  if(ld_col != NULL){
    colnames(assoc_data)[ld_col] = "LD"}

  # read in, format, and filter data sets
  message("Preparing association data...")
  assoc_data <- as.data.frame(assoc_data)
  assoc_data$CHR_POS = as.numeric(as.character(assoc_data$CHR_POS))
  assoc_data$LOG10P = as.numeric(as.character(assoc_data$LOG10P))
  assoc_data$CHR = as.numeric(as.character(assoc_data$CHR))

  return(assoc_data)
}

