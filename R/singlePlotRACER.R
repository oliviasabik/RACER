#' Single Regional Association Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values) of an association study
#' by their genomic position, for example, the results of a GWAS or eQTL study. Sources
#' 1000K genomes phase III data for linkage disequilibrium calculations.
#' @param assoc_data required. A dataframe that has been produced by formatRACER and has columns named CHR, POS
#' @param chr required. numeric. chromosome to plot
#' @param build optional. default = "hg19", can also optionally be set to "hg38"
#' @param plotby required. "coord", "gene", or "snp". Which parameter to use to determine the reigon to be plotted.
#' @param gene_plot optional. Required if "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param snp_plot optional. Required if "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start_plot optional. Required if "coord" selected for plotby, then this will be lower bound of x axis
#' @param end_plot optional. Required if "coord" selected for plotby, then this will be upper bound of x axis
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

singlePlotRACER <- function(assoc_data, chr, build="hg19", plotby, gene_plot = NULL, snp_plot = NULL, start_plot=NULL, end_plot = NULL){

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
    utils::data(biomart_hg38)
    chr_in = chr
    colnames(biomart_hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
    gene_sub = subset(biomart_hg38, biomart_hg38$CHR == chr_in)
  }else if(build == "hg19"){
    utils::data(biomart_hg19)
    chr_in = chr
    colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
    gene_sub = subset(biomart_hg19, biomart_hg19$CHR == chr_in)
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
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start-5000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end+5000))
  myCol = paste0("desc(", "LENGTH)")
  gene_sub %>%
    dplyr::arrange_(.dots = c(myCol))
  gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID),]
  gene_sub = gene_sub[,c(3,4,5)]
  gene_sub = reshape2::melt(gene_sub,id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

  # read in, format, and filter data sets
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter_(in.dt, ~CHR == chr_in)
  in.dt = dplyr::filter_(in.dt, ~POS > start) %>%
    dplyr::filter_(~POS < end)

  # Generate plots
  message("Generating Plot")
  if("LD" %in% colnames(in.dt) && "LD_BIN" %in% colnames(in.dt)){
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
      ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                         hjust = -0.1,vjust = 0.3, size = 2.5) + ggplot2::xlim(start,end) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylim(0,(max(gene_sub$y_value)+1))

    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", y = "LOG10P", color = "LD_BIN")) +
      ggplot2::geom_point() + ggplot2::scale_colour_manual(
        values = c("1.0-0.8" = "red", "0.8-0.6" = "darkorange1", "0.6-0.4" = "green1",
                   "0.4-0.2" = "skyblue1", "0.2-0.0" = "navyblue", "NA" = "grey"), drop = FALSE) +
      ggplot2::theme_bw() + ggplot2::xlab("Chromosome Position") + ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt$LOG10P),max(in.dt$LOG10P))

    ggpubr::ggarrange(b, c, heights = c(4,1), nrow = 2, ncol = 1,
                      common.legend = TRUE, legend = "right")
  }else{
    c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
      ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                         hjust = -0.1,vjust = 0.3, size = 2.5) + ggplot2::xlim(start,end) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylim(0,(max(gene_sub$y_value)+1))

    b = ggplot2::ggplot(in.dt, ggplot2::aes_string(x = "POS", y = "LOG10P")) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("Chromosome Position") +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt$LOG10P),max(in.dt$LOG10P))

    ggpubr::ggarrange(b, c, heights = c(4,1), nrow = 2, ncol = 1,
                      common.legend = TRUE, legend = "right")
  }

}
