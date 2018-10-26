#' Single Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values) of an association study
#' by their genomic position, for example, the results of a GWAS or eQTL study. Sources
#' 1000K genomes phase III data for linkage disequilibrium calculations.
#' @param assoc_data required. A dataframe that has columns named
#' CHR_POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values),
#' CHR, which contains the chromosome number, RS_ID, which
#' contains the rsID numbers of the SNPs, and LD, which contains
#' LD information. If you want to compute LD information for the dataset,
#' be sure to include rs_id numbers for all the SNPs in your input data.
#' @param chr_col required. numeric. index of column in assoc_data containing chromosome information
#' @param pos_col required. numeric. index of column in assoc_data containing genomic position information
#' @param log10p_col required. numeric. index of column in assoc_data containing -log10(p-value)s
#' @param chr required. numeric. chromosome to plot
#' @param plotby required. "coord", "gene", or "snp". Which parameter to use to
#' determine the reigon to be plotted.
#' @param gene_plot optional. Required if "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param snp_plot optional. Required if "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start_plot optional. Required if "coord" selected for plotby, then this will be lower bound of x axis
#' @param end_plot optional. Required if "coord" selected for plotby, then this will be upper bound of x axis
#' @param ldby required. default = "none", but can be "input", or "1000genomes"
#' @param rs_col optional. numeric. Required if ldby = "1000genomes", index of column containing rsID numbers for SNPs
#' @param ld_col optional. numeric. Required if ldby = "input", index of column in assoc_data containing LD information
#' @param pops optional. Required if ldby = "1000genomes". Populations used to calculate LD.
#' @param snp_ld optional. Required if ldby = "1000genomes". snp used to calculate LD
#'
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @import ggplot2
#' @examples
#' single_plot_function(assoc_data = assoc_data, chr = 1, plotby = "gene/snp/coord", x_plot = "GENE_NAME/RS_ID/START/END",
#'  ldby= "none/input/1000genomes", pops = c("POP1", "POP2", etc...), snp_ld = "RS_ID")

single_plot_function <- function(assoc_data, chr_col, pos_col, log10p_col, chr, build="hg19", plotby, gene_plot = NULL, snp_plot = NULL, start_plot=NULL, end_plot=NULL, ldby = "none", pops=NULL, snp_ld=NULL, rs_col=NULL, ld_col=NULL){

  if(missing(chr_col)){
    message("Please specify which column contains chromosome information.")
  }else if(missing(pos_col)){
    message("Please specify which column contains genomic position information.")
  }else if(missing(log10p_col)){
    message("Please specify which column contains genomic position information.")
  }else if(missing(chr)){
    message("Please specify which chromosome you wish to plot.")
  }else if(missing(plotby)){
    message("Please specify the method by which you wish to plot.")
  }else{
    message("All inputs are go.")
  }

  colnames(assoc_data)[chr_col] = "CHR"
  colnames(assoc_data)[pos_col] = "CHR_POS"
  colnames(assoc_data)[log10p_col] = "LOG10P"
  if(!missing(rs_col)){
    colnames(assoc_data)[rs_col] = "RS_ID"}
  if(!missing(ld_col)){
    colnames(assoc_data)[ld_col] = "LD"}

  `%>%` <- magrittr::`%>%`

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
      q = subset(assoc_data, RS_ID == snp_plot)
      w = q$CHR_POS
      w = as.numeric(as.character(w))
      start = w - 500000
      end = w + 500000}
  }

  if(sum(is.null(plotby)) == 1){
    message("Please specify a region to plot.")
  }

  # reading in gene data
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start-5000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end+5000))
  gene_sub = dplyr::arrange(gene_sub, desc(LENGTH))
  gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID),]
  gene_sub = gene_sub[,c(3,4,5)]
  gene_sub = reshape2::melt(gene_sub,id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

  # read in, format, and filter data sets
  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_in)
  in.dt = dplyr::filter(in.dt, CHR_POS > start) %>%
    dplyr::filter(CHR_POS < end)

  # calculate LD
  message("Calculating LD")
  if((sum(is.null(ldby)) == 0) == TRUE){
    if(ldby == "none"){
      message("Not including LD information in the plot.")
    }else if(ldby == "input"){
      if("LD" %in% colnames(in.dt)){
        message("Using LD info from the input data set.")
        in.dt$LD = as.numeric(as.character(in.dt$LD))
        in.dt$LD_BIN <- cut(in.dt$LD,
                            breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                            labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
        in.dt$LD_BIN = as.character(in.dt$LD_BIN)
        in.dt$LD_BIN[is.na(in.dt$LD_BIN)] <- "NA"
        in.dt$LD_BIN = as.factor(in.dt$LD_BIN)
        in.dt$LD_BIN = factor(in.dt$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
      }else{
        "No LD column in input. Please name your column of LD data in your input 'LD'."
      }
    }else if(ldby == "1000genomes"){
      message(paste0("Pulling LD data from 1000 Genomes Phase III data using ", snp_ld))
      in.dt$LD_BIN = 1
      in.dt$LD_BIN = NA
      if(!is.null(pops)){
        if(length(pops) == 1){
          ld_command_1 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", snp_ld,
                                "&pop=",pops,"&r2_d=r2'")
          z_1 = system(ld_command_1, intern = TRUE)
          z_1 = as.data.frame(z_1)
          z_1 = tidyr::separate(z_1, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                                 "col_9", "col_10"), sep = "\t")
          colnames(z_1) = z_1[1,]
          z_1 = z_1[-1,]
          z_1 = dplyr::select(z_1, RS_Number, R2)
          colnames(z_1) = c("RS_ID", "LD")
          in.dt$LD = NA
          in.dt = dplyr::select(in.dt, -LD)
          in.dt = merge(in.dt, z_1, by = "RS_ID", all.x = TRUE)
          in.dt$LD = as.numeric(as.character(in.dt$LD))
          in.dt$LD_BIN <- cut(in.dt$LD,
                              breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                              labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
          in.dt$LD_BIN = as.character(in.dt$LD_BIN)
          in.dt$LD_BIN[is.na(in.dt$LD_BIN)] <- "NA"
          in.dt$LD_BIN = as.factor(in.dt$LD_BIN)
          in.dt$LD_BIN = factor(in.dt$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
        }else if(length(pops) > 1){
          ld_command_1 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", snp_ld,
                                "&pop=")
          for (i in 1:length(pops)){
            if (i < length(pops)){
              a = pops[i]
              ld_command_1 = paste0(ld_command_1, a, "%2B")
            } else if(i == length(pops)){
              a = pops[i]
              ld_command_1 = paste0(ld_command_1, a)
            }
          }
          ld_command_1 = paste0(ld_command_1, "&r2_d=r2'")
          z_1 = system(ld_command_1, intern = TRUE)
          z_1 = as.data.frame(z_1)
          z_1 = tidyr::separate(z_1, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                                 "col_9", "col_10"), sep = "\t")
          colnames(z_1) = z_1[1,]
          z_1 = z_1[-1,]
          z_1 = dplyr::select(z_1, RS_Number, R2)
          colnames(z_1) = c("RS_ID", "LD")
          in.dt$LD = NA
          in.dt = dplyr::select(in.dt, -LD)
          in.dt = merge(in.dt, z_1, by = "RS_ID", all.x = TRUE)
          in.dt$LD = as.numeric(as.character(in.dt$LD))
          in.dt$LD_BIN <- cut(in.dt$LD,
                              breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                              labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
          in.dt$LD_BIN = as.character(in.dt$LD_BIN)
          in.dt$LD_BIN[is.na(in.dt$LD_BIN)] <- "NA"
          in.dt$LD_BIN = as.factor(in.dt$LD_BIN)
          in.dt$LD_BIN = factor(in.dt$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
        }else{
          "Please specify pops for LD calculation."
        }
      }
    }else{
      "Invalid LD input method. Please specify a valid LD input method."
    }
  }

  # Generate plots
  message("Generating Plot")
  if(ldby !="none"){
    c = ggplot2::ggplot(gene_sub, ggplot2::aes(x = value, y = y_value)) +
      ggplot2::geom_line(ggplot2::aes(group = GENE_NAME), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes(x = value, y = y_value, label = GENE_NAME),
                hjust = -0.1,vjust = 0.3, size = 2.5) + ggplot2::xlim(start,end) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylim(0,(max(gene_sub$y_value)+1))

    b = ggplot2::ggplot(in.dt, ggplot2::aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
      ggplot2::geom_point() + ggplot2::scale_colour_manual(
        values = c("1.0-0.8" = "red", "0.8-0.6" = "darkorange1", "0.6-0.4" = "green1",
                   "0.4-0.2" = "skyblue1", "0.2-0.0" = "navyblue", "NA" = "grey"), drop = FALSE) +
      ggplot2::theme_bw() + ggplot2::xlab("Chromosome Position") + ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt$LOG10P),max(in.dt$LOG10P))

    ggpubr::ggarrange(b, c, heights = c(4,1), nrow = 2, ncol = 1,
                      common.legend = TRUE, legend = "right")
  }else if(ldby == "none"){
    c = ggplot2::ggplot(gene_sub, ggplot2::aes(x = value, y = y_value)) +
      ggplot2::geom_line(ggplot2::aes(group = GENE_NAME), size = 2) + ggplot2::theme_bw() +
      ggplot2::geom_text(data = plot_lab, ggplot2::aes(x = value, y = y_value, label = GENE_NAME),
                hjust = -0.1,vjust = 0.3, size = 2.5) + ggplot2::xlim(start,end) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylim(0,(max(gene_sub$y_value)+1))

    b = ggplot2::ggplot(in.dt, ggplot2::aes(x = CHR_POS, y = LOG10P)) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("Chromosome Position") +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::xlim(start, end) + ggplot2::ylim(min(in.dt$LOG10P),max(in.dt$LOG10P))

    ggpubr::ggarrange(b, c, heights = c(4,1), nrow = 2, ncol = 1,
                      common.legend = TRUE, legend = "right")
  }

}

