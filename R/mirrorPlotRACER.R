#' Mirror Plot -- Regional Association ComparER Plot
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data1 required. A dataframe that has columns named
#' POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values),
#' CHR, which contains the chromosome number, RS_ID, which
#' LD information. If no column named LD is in the input, LD will
#' be calculated from 1000 genomes phase III in relation to snp_ld_1
#'
#' @param assoc_data2 required. identifcal in format to assoc_data1
#' @param chr required. chromosome you wish to plot
#' @param name1 optional. name of association set 1
#' @param name2 optional. name of association set 2
#' @param plotby required. "coord", "gene", or "snp". Which parameter to use to
#' determine the reigon to be plotted.
#' @param gene_plot optional. If "gene" selected for plotby, then plot will be +/- 50kb of gene, should be a human gene symbol
#' @param snp_plot optional. If "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start_plot optional. If "coord" selected for plotby, then this will be lower bound of x axis
#' @param end_plot optional. If "coord" selected for plotby, then this will be upper bound of x axis


#' @keywords association plot
#' @export
#' @import ggplot2
#' @examples
#' data(mark3_bmd_gwas)
#' data(mark3_eqtl)
#'
#' mark3_bmd_gwas_f = RACER::formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3, pos_col = 4, p_col = 11)
#' mark3_eqtl_f = RACER::formatRACER(assoc_data = mark3_eqtl, chr_col = 10, pos_col = 11, p_col = 7)
#'
#' mark3_bmd_gwas_f_ld = RACER::ldRACER(assoc_data = mark3_bmd_gwas_f, rs_col = 2, pops = "EUR", lead_snp = "rs11623869")
#' mark3_eqtl_f_ld = RACER::ldRACER(assoc_data = mark3_eqtl_f, rs_col = 15, pops = "EUR", lead_snp = "rs11623869")
#' mirrorPlotRACER(assoc_data1 = mark3_bmd_gwas_f_ld, assoc_data2 = mark3_eqtl_f_ld, chr = 14, plotby = "gene", gene_plot = "MARK3")

mirrorPlotRACER <- function(assoc_data1, assoc_data2, chr, build = "hg19", name1="Association Dataset #1", name2="Association Dataset #2", plotby, gene_plot=NULL, snp_plot=NULL, start_plot=NULL, end_plot=NULL, ldby = "none", pops=NULL, snp_ld_1=NULL, snp_ld_2=NULL){
  reqs = c("CHR", "POS", "LOG10P")
  cols_1 = colnames(assoc_data1)
  cols_2 = colnames(assoc_data2)
  if(sum(reqs %in% cols_1) == 3){
  }else{stop("Association Data Set #1 is missing a required column.")}
  if(sum(reqs %in% cols_2) == 3){
  }else{stop("Association Data Set #2 is missing a required column.")}

  if(build == "hg38"){
    data(biomart_hg38)
    chr_in = chr
    colnames(biomart_hg38) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
    gene_sub = subset(biomart_hg38, biomart_hg38$CHR == chr_in)
  }else if(build == "hg19"){
    data(biomart_hg19)
    chr_in = chr
    colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
    gene_sub = subset(biomart_hg19, biomart_hg19$CHR == chr_in)
  }

  `%>%` <- magrittr::`%>%`

  if((sum(is.null(plotby)) == 0) == TRUE){
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
        }else{stop("No gene specified.")}
    }else if((plotby == "snp") == TRUE){
      message(paste("snp",snp_plot))
      q = subset(assoc_data1, RS_ID == snp_plot)
      w = q$POS
      w = as.numeric(as.character(w))
      start = w - 500000
      end = w + 500000}
  }else{
    stop("Please specify a parameter to plotby.")
  }

  # reading in gene data
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start-5000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end+5000))
  gene_sub = dplyr::arrange(gene_sub, desc(LENGTH))
  gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID),]
  gene_sub = gene_sub[,c(3,4,5)]
  gene_sub = reshape2::melt(gene_sub, id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

  # read in, format, and filter data sets
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_in)
  in.dt = dplyr::filter(in.dt, POS > start)%>%
    dplyr::filter(POS < end)

  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$POS = as.numeric(as.character(in.dt.2$POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  in.dt.2 = dplyr::filter(in.dt.2, CHR == chr_in)
  in.dt.2= dplyr::filter(in.dt.2, POS > start)%>%
    dplyr::filter(POS < end)

  # generate mirror plot
  message("Generating plot.")
  if("LD" %in% cols_1 && "LD_BIN" %in% cols_1){
    a = ggplot2::ggplot(data = in.dt, ggplot2::aes(x = POS, y = LOG10P, color = LD_BIN)) +
      ggplot2::geom_point() + ggplot2::scale_colour_manual(
        values = c("1.0-0.8" = "red", "0.8-0.6" = "darkorange1", "0.6-0.4" = "green1",
                   "0.4-0.2" = "skyblue1", "0.2-0.0" = "navyblue", "NA" = "grey"), drop = FALSE) +
      ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + ggplot2::ylab("-log10(p-value)") +
      ggplot2::scale_y_reverse() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                                axis.text.x=ggplot2::element_blank(),
                                axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlim(start,end) + ggplot2::ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))
  }else{
    message("No LD information for dataset #1.")
    a = ggplot2::ggplot(in.dt, ggplot2::aes(x = POS, y = LOG10P)) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::scale_y_reverse() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                                                  axis.text.x=ggplot2::element_blank(),
                                                  axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlim(start,end) + ggplot2::ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))
  }

  if("LD" %in% cols_2 && "LD_BIN" %in% cols_2){
    b = ggplot2::ggplot(data = in.dt.2, ggplot2::aes(x = POS, y = LOG10P, color = LD_BIN)) +
      ggplot2::geom_point() + ggplot2::scale_colour_manual(
        values = c("1.0-0.8" = "red", "0.8-0.6" = "darkorange1", "0.6-0.4" = "green1",
                   "0.4-0.2" = "skyblue1", "0.2-0.0" = "navyblue", "NA" = "grey"), drop = FALSE) +
      ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) +
      ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") +
      ggplot2::xlim(start,end) + ggplot2::ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P)) +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank())
    }else{
      b = ggplot2::ggplot(in.dt.2, ggplot2::aes(x = POS, y = LOG10P)) +
        ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) +
        ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") +
        ggplot2::xlim(start,end) + ggplot2::ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P)) +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank())
  }

    c = ggplot2::ggplot(gene_sub, ggplot2::aes(x = value, y = y_value)) +
      ggplot2::geom_line(ggplot2::aes(group = GENE_NAME), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes(x = value, y = y_value, label = GENE_NAME),
                hjust = -0.1,vjust = 0.3, size = 2.5) + ggplot2::xlim(start,end) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylim(0,(max(gene_sub$y_value)+1))

    ggpubr::ggarrange(a, b, c, heights = c(2,2,1), nrow = 3, ncol = 1,
                      common.legend = TRUE, legend = "right")

}


