#' Single Regional Association Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values) of an association study
#' by their genomic position, for example, the results of a GWAS or eQTL study. Be sure
#' your input association data has been formatted using formatRACER. If you want to
#' include linkage disequilibirum data, you can use the ldRACER function to use the 1000
#' genomes project to calculate LD for your data.
#' @param assoc_data required. A dataframe that has been produced by formatRACER and has columns named CHR, POS
#' @param chr required. numeric. chromosome to plot
#' @param build optional. default = "hg19", can also optionally be set to "hg38", depending on the build of your input data
#' @param set optional. default = "protein_coding", however can be set to "all" to plot all RNAs in the genome
#' @param plotby required. "coord", "gene", or "snp". Which parameter to use to determine the reigon to be plotted.
#' @param gene_plot optional. Required if "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param snp_plot optional. Required if "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start_plot optional. Required if "coord" selected for plotby, then this will be lower bound of x axis
#' @param end_plot optional. Required if "coord" selected for plotby, then this will be upper bound of x axis
#' @param label_lead optional. default = FALSE, set = TRUE if you wish to add a label to your graph of the SNP used to calculate LD. If the SNP used to calculate LD is not in your data set, the SNP with the greatest -LOG10(P) will be labeled.
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' data(mark3_bmd_gwas)
#' mark3_bmd_gwas_f = formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3,
#' pos_col = 4, p_col = 11)
#' mark3_bmd_gwas_f_ld = ldRACER(assoc_data = mark3_bmd_gwas_f,
#' rs_col = 2, pops = c("EUR"), lead_snp = "rs11623869")
#' singlePlotRACER(assoc_data = mark3_bmd_gwas_f_ld, chr = 14,
#' build = "hg19", plotby = "coord", start_plot = 103500000, end_plot = 104500000)}

singlePlotRACER <- function(assoc_data, chr, build="hg19", set = "protein_coding", plotby, gene_plot = NULL, snp_plot = NULL, start_plot=NULL, end_plot = NULL, label_lead = FALSE){

  if(missing(assoc_data)){
    stop("Please provide a data set to plot.")
  }else if(missing(chr)){
    stop("Please specify which chromosome you wish to plot.")
  }else if(missing(plotby)){
    stop("Please specify the method by which you wish to plot.")
  }else if(plotby == "gene"){
    if(is.null(gene_plot)){
      stop("Please specify a gene to plot by.")
    }
  }else if(plotby == "snp"){
    if(is.null(snp_plot)){
      stop("Please specify a snp to plot by.")
    }
  }else if(plotby == "coord"){
    if(is.null(start_plot) | is.null(end_plot)){
      stop("Please specify start coordinate for plot.")
    }
  }else{
    message("All inputs are go.")
  }

  reqs = c("CHR", "POS", "LOG10P")
  cols = colnames(assoc_data)
  if(sum(reqs %in% cols) == 3){
  }else{stop("Association Data Set is missing a required column, please format your data set using formatRACER.R.")}

  reqs_2 = c("LD", "LD_BIN")
  if(sum(reqs_2 %in% cols) == 2){
  }else{message("Association Data Set is missing LD data, the resulting plot won't have LD information, but you can add it using the ldRACER.R function.")}

  `%>%` <- magrittr::`%>%`

  if(build == "hg38"){
    utils::data(hg38)
    chr_in = chr
    colnames(hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg38[hg38$CHR == chr_in,]
  }else if(build == "hg19"){
    utils::data(hg19)
    chr_in = chr
    colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "LENGTH", "GENE_NAME", "TYPE")
    gene_sub = hg19[hg19$CHR == chr_in,]
  }

  if(set == "protein_coding"){
    gene_sub = gene_sub[gene_sub$TYPE == "protein_coding",]
  }else{
    gene_sub = gene_sub
  }

  if(sum(is.null(plotby)) == 1){
    stop("Please specify a method by which to plot.")
  }

  if(sum(is.null(plotby)) == 0){
    message("Plotting by...")
    if((plotby == "coord") == TRUE){
      message("coord")
      start = start_plot
      end = end_plot
    }else if((plotby == "gene") == TRUE){
      message(paste("gene:",gene_plot))
      if(sum(is.null(gene_plot)) == 0){
        p = subset(gene_sub, gene_sub$GENE_NAME == gene_plot)
        start = min(p$TRX_START) - 500000
        end = max(p$TRX_END) + 500000
      }else{message("No gene specified.")}
    }else if((plotby == "snp") == TRUE){
      message(paste("snp",snp_plot))
      q = assoc_data[assoc_data$RS_ID == snp_plot,]
      w = q$POS
      w = as.numeric(as.character(w))
      start = w - 500000
      end = w + 500000}
  }

  # reading in gene data
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start-50000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end+50000))
  gene_sub = gene_sub[,c(3,4,6)]
  gene_sub = reshape2::melt(gene_sub,id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

  # read in, format, and filter data sets
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_in)
  in.dt = dplyr::filter(in.dt, POS > start) %>%
    dplyr::filter(POS < end)

  if(label_lead == TRUE){
    lsnp_row = which(in.dt$LABEL == "LEAD")
    label_data = in.dt[lsnp_row,]
    if(dim(label_data)[1] == 0){
      lsnp_row = in.dt[in.dt$LOG10P == max(in.dt$LOG10P),]
      label_data = lsnp_row[1,]
    }
  }
  # Generate plots
  message("Generating Plot")
  if("LD" %in% colnames(in.dt) && "LD_BIN" %in% colnames(in.dt)){
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
      ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                         hjust = -0.1,vjust = 0.3, size = 2.5) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::coord_cartesian(xlim = c(start,end), ylim = c(0,(max(gene_sub$y_value)+1)))

    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", y = "LOG10P", color = "LD_BIN")) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_manual(
        values = c("1.0-0.8" = "red", "0.8-0.6" = "darkorange1", "0.6-0.4" = "green1",
                   "0.4-0.2" = "skyblue1", "0.2-0.0" = "navyblue", "NA" = "grey"), drop = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Chromosome Position") + ggplot2::ylab("-log10(p-value)") +
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(min(in.dt$LOG10P),max(in.dt$LOG10P)))

  }else{
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
      ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                         hjust = -0.1,vjust = 0.3, size = 2.5) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::coord_cartesian(xlim = c(start,end), ylim = c(0,(max(gene_sub$y_value)+1)))

    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", y = "LOG10P")) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("Chromosome Position") +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(min(in.dt$LOG10P),max(in.dt$LOG10P)))
  }
  if(label_lead == TRUE){
    b = b + geom_point(data = label_data, aes_string(x = "POS", y = "LOG10P"), color = "purple")
    b = b + geom_text(data = label_data, aes_string(label = "RS_ID"), color = "black", size = 3, hjust = 1.25)
  }
  ggpubr::ggarrange(b, c, heights = c(3,1), nrow = 2, ncol = 1,
                    common.legend = TRUE, legend = "right")
}
