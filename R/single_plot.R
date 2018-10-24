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
#' @param chr required. chromosome to plot
#' @param plotby required. "coord", "gene", or "snp". Which parameter to use to
#' determine the reigon to be plotted.
#' @param gene_plot optional. If "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param snp_plot optional. If "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start_plot optional. If "coord" selected for plotby, then this will be lower bound of x axis
#' @param end_plot optional. If "coord" selected for plotby, then this will be upper bound of x axis
#' @param ldby required. "none", "input", or "1000genomes"
#' @param pops optional. required if ldby = "1000genomes". Populations used
#' to calculate LD.
#' @param snp_ld optional. required if ldby = "1000genomes". snp used to calculate LD
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @import ggplot2
#' @examples
#' single_plot_function(assoc_data = assoc_data, chr = 1, plotby = "gene/snp/coord", x_plot = "GENE_NAME/RS_ID/START/END",
#'  ldby= "none/input/1000genomes", pops = c("POP1", "POP2", etc...), snp_ld = "RS_ID")

single_plot_function <- function(assoc_data, chr, plotby, gene_plot = NULL, snp_plot = NULL,
                                 start_plot=NULL, end_plot=NULL, ldby, pops=NULL, snp_ld=NULL){
  reqs = c("CHR", "CHR_POS", "LOG10P")
  cols = colnames(assoc_data)
  if(sum(reqs %in% cols) == 3){
    data = "good"
  }else{print("Association data set is missing a required column.")}

  `%>%` <- magrittr::`%>%`

  data(biomart_hg19)
  chr_in = chr
  colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
  gene_sub = subset(biomart_hg19, biomart_hg19$CHR == chr_in)

  if(sum(is.null(plotby)) == 0){
    print("Plotting by...")
    if((plotby == "coord") == TRUE){
      print("coord")
      start = start_plot
      end = end_plot
    }else if((plotby == "gene") == TRUE){
      print(paste("gene:",gene_plot))
            if(sum(is.null(gene_plot)) == 0){
              p = subset(gene_sub, gene_sub$GENE_NAME == gene_plot)
              start = min(p$TRX_START) - 500000
              end = max(p$TRX_END) + 500000
            }else{print("No gene specified.")}
    }else if((plotby == "snp") == TRUE){
      print(paste("snp",snp_plot))
      q = subset(assoc_data1, RS_ID == snp_plot)
      w = q$CHR_POS
      w = as.numeric(as.character(w))
      start = w - 500000
      end = w + 500000}
  }

  if(sum(is.null(plotby)) == 1){
    print("Please specify a region to plot.")
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
  print("Reading in association data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_in)
  in.dt = dplyr::filter(in.dt, CHR_POS > start) %>%
    dplyr::filter(CHR_POS < end)

  # calculate LD
  print("Calculating LD")
  if((sum(is.null(ldby)) == 0) == TRUE){
    if(ldby == "none"){
      print("Not including LD information in the plot.")
    }else if(ldby == "input"){
      if("LD" %in% colnames(in.dt)){
        print("Using LD info from the input data set.")
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
      print(paste0("Pulling LD data from 1000 Genomes Phase III data using ", snp_ld))
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
  print("Generating Plot")
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

