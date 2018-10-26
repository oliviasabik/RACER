#' Genomic Comparsion Formating Function
#'
#' This group of functions allows you to creat a plot of -log10(P-values) of an association study by their genomic position, for example, the results of a GWAS or eQTL study. This function carries out the formatting necessary to generate a plot
#' @param assoc_data required. A dataframe that has columns containing the chromosome, physical position, and p-values or -log10(p-values) of the association, and can optionally have columns containing R2 information for LD in the region, or rsID numbers for the associated SNPs.
#' @param chr_col required. numeric. index of column in assoc_data containing chromosome information
#' @param pos_col required. numeric. index of column in assoc_data containing genomic position information
#' @param log10p_col required. numeric. index of column in assoc_data containing -log10(p-value)s
#' @param rs_col optional. numeric. Required if you want to calculated LD using 1000 genomes for analysis, index of column containing rsID numbers for SNPs
#' @param ld_col optional. numeric. Required if you want to use the LD data in your data set in your plot, index of column in assoc_data containing LD information
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @import ggplot2
#' @examples
#' single_plot_function(assoc_data = assoc_data, chr = 1, plotby = "gene/snp/coord", x_plot = "GENE_NAME/RS_ID/START/END",
#'  ldby= "none/input/1000genomes", pops = c("POP1", "POP2", etc...), snp_ld = "RS_ID")

gen_comp_format <- function(assoc_data, chr_col, pos_col, p_col=NULL, log10p_col=NULL, rs_col=NULL, ld_col=NULL){

  if(missing(chr_col)){
    message("Please specify which column contains chromosome information.")
  }else if(missing(pos_col)){
    message("Please specify which column contains genomic position information.")
  }else if(log10p_col == NULL & p_col == NULL){
    message("Please specify which column contains p-values or -log10(p-values).")
  }else{
    message("All inputs are go.")
  }

  colnames(assoc_data)[chr_col] = "CHR"
  colnames(assoc_data)[pos_col] = "CHR_POS"
  colnames(assoc_data)[log10p_col] = "LOG10P"
  if(rs_col != NULL){
    colnames(assoc_data)[rs_col] = "RS_ID"}
  if(ld_col != NULL){
    colnames(assoc_data)[ld_col] = "LD"}

  # read in, format, and filter data sets
  message("Reading in association data")
  assoc_data <- as.data.frame(assoc_data)
  assoc_data$CHR_POS = as.numeric(as.character(assoc_data$CHR_POS))
  assoc_data$LOG10P = as.numeric(as.character(assoc_data$LOG10P))
  assoc_data$CHR = as.numeric(as.character(assoc_data$CHR))
  return(assoc_data)
}

