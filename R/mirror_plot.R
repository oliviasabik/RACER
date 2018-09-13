#' Mirror Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data1 required. A dataframe that has columns named
#' CHR_POS representing the position of the snp on
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
#' @param gene_plot optional. If "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param snp_plot optional. If "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start_plot optional. If "coord" selected for plotby, then this will be lower bound of x axis
#' @param end_plot optional. If "coord" selected for plotby, then this will be upper bound of x axis
#' @param ldby required. "none", "input", or "1000genomes"
#' @param pops optional. required if ldby = "1000genomes". Populations used
#' to calculate LD.
#' @param snp_ld_1 optional. required if ldby = "1000genomes". snp used to calculate LD
#' @param snp_ld_2 optional. required if ldby = "1000genomes", and you want to use a different snp
#'to calculate LD for the second dataset

#' @keywords association plot
#' @export
#' @examples
#' mirror_plot_function(assoc_data1, assoc_data2, chr, "name1", "name2", plotby = "coord", start_plot, end_plot,
#' ldby = "1000genomes", pops = c("CEU","TSI","FIN","GBR","IBS"), snp_ld_1 = "rs123456", snp_ld_2 = "rs123456")

mirror_plot_function <- function(assoc_data1, assoc_data2, chr, name1=NULL, name2=NULL,
                                 plotby, gene_plot=NULL, snp_plot=NULL, start_plot=NULL, end_plot=NULL,
                                 ldby, pops=NULL, snp_ld_1=NULL, snp_ld_2=NULL){
  reqs = c("CHR", "CHR_POS", "LOG10P")
  cols_1 = colnames(assoc_data1)
  cols_2 = colnames(assoc_data2)
  if(sum(reqs %in% cols_1) == 3){
  }else{print("Association Data Set #1 is missing a required column.")}
  if(sum(reqs %in% cols_2) == 3){
  }else{print("Association Data Set #2 is missing a required column.")}

  data(biomart_hg19)
  chr_in = chr
  colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
  gene_sub = subset(biomart_hg19, biomart_hg19$CHR == chr_in)

  `%>%` <- magrittr::`%>%`

  if((sum(is.null(plotby)) == 0) == TRUE){
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
  }else{
    print("Please specify a parameter to plotby.")
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
  print("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_in)
  in.dt = dplyr::filter(in.dt, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$CHR_POS = as.numeric(as.character(in.dt.2$CHR_POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  in.dt.2 = dplyr::filter(in.dt.2, CHR == chr_in)
  in.dt.2= dplyr::filter(in.dt.2, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  # determining input for LD info
  # NUMBER 1
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
      print(paste0("Pulling LD data from 1000 Genomes Phase III data using ", snp_ld_1))
      in.dt$LD_BIN = 1
      in.dt$LD_BIN = NA
      if(!is.null(pops)){
        if(length(pops) == 1){
          ld_command_1 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", snp_ld_1,
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
          ld_command_1 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", snp_ld_1,
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

  # NUMBER 2
  if((sum(is.null(snp_ld_2)) == 0) == FALSE){
    snp_ld_2 = snp_ld_1
  }else{
    print(snp_ld_1)
  }
  if((sum(is.null(ldby)) == 0) == TRUE){
    if(ldby == "none"){
      print("Not including LD information in the plot.")
    }else if(ldby == "input"){
      if("LD" %in% colnames(in.dt.2)){
        print("Using LD info from the input data set.")
        in.dt.2$LD = as.numeric(as.character(in.dt.2$LD))
        in.dt.2$LD_BIN <- cut(in.dt.2$LD,
                            breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                            labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
        in.dt.2$LD_BIN = as.character(in.dt.2$LD_BIN)
        in.dt.2$LD_BIN[is.na(in.dt.2$LD_BIN)] <- "NA"
        in.dt.2$LD_BIN = as.factor(in.dt.2$LD_BIN)
        in.dt.2$LD_BIN = factor(in.dt.2$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
      }else{
        "No LD column in input. Please name your column of LD data in your input 'LD'."
      }
    }else if(ldby == "1000genomes"){
      print(paste0("Pulling LD data from 1000 Genomes Phase III data using ", snp_ld_2))
      in.dt.2$LD_BIN = 1
      in.dt.2$LD_BIN = NA
      if(!is.null(pops)){
        if(length(pops) == 1){
          ld_command_2 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", snp_ld_2,
                                "&pop=",pops,"&r2_d=r2'")
          z_2 = system(ld_command_2, intern = TRUE)
          z_2 = as.data.frame(z_2)
          z_2 = tidyr::separate(z_2, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                                 "col_9", "col_10"), sep = "\t")
          colnames(z_2) = z_2[1,]
          z_2 = z_2[-1,]
          z_2 = dplyr::select(z_2, RS_Number, R2)
          colnames(z_2) = c("RS_ID", "LD")
          in.dt.2$LD = NA
          in.dt.2 = dplyr::select(in.dt.2, -LD)
          in.dt.2 = merge(in.dt.2, z_2, by = "RS_ID", all.x = TRUE)
          in.dt.2$LD = as.numeric(as.character(in.dt.2$LD))
          in.dt.2$LD_BIN <- cut(in.dt.2$LD,
                              breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                              labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
          in.dt.2$LD_BIN = as.character(in.dt.2$LD_BIN)
          in.dt.2$LD_BIN[is.na(in.dt.2$LD_BIN)] <- "NA"
          in.dt.2$LD_BIN = as.factor(in.dt.2$LD_BIN)
          in.dt.2$LD_BIN = factor(in.dt.2$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
        }else if(length(pops) > 1){
          ld_command_2 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", snp_ld_2,
                                "&pop=")
          for (i in 1:length(pops)){
            if (i < length(pops)){
              a = pops[i]
              ld_command_2 = paste0(ld_command_2, a, "%2B")
            } else if(i == length(pops)){
              a = pops[i]
              ld_command_2 = paste0(ld_command_2, a)
            }
          }
          ld_command_2 = paste0(ld_command_2, "&r2_d=r2'")
          z_2 = system(ld_command_2, intern = TRUE)
          z_2 = as.data.frame(z_2)
          z_2 = tidyr::separate(z_2, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                                 "col_9", "col_10"), sep = "\t")
          colnames(z_2) = z_2[1,]
          z_2 = z_2[-1,]
          z_2 = dplyr::select(z_2, RS_Number, R2)
          colnames(z_2) = c("RS_ID", "LD")
          in.dt.2$LD = NA
          in.dt.2 = dplyr::select(in.dt.2, -LD)
          in.dt.2 = merge(in.dt.2, z_2, by = "RS_ID", all.x = TRUE)
          in.dt.2$LD = as.numeric(as.character(in.dt.2$LD))
          in.dt.2$LD_BIN <- cut(in.dt.2$LD,
                              breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                              labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
          in.dt.2$LD_BIN = as.character(in.dt.2$LD_BIN)
          in.dt.2$LD_BIN[is.na(in.dt.2$LD_BIN)] <- "NA"
          in.dt.2$LD_BIN = as.factor(in.dt.2$LD_BIN)
          in.dt.2$LD_BIN = factor(in.dt.2$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
        }else{
          "Please specify pops for LD calculation."
        }
      }
    }else{
      "Invalid LD input method. Please specify a valid LD input method."
    }
  }

  # generate mirror plot
  print("Generating plot")
  if(ldby != "none"){
    a = ggplot2::ggplot(data = in.dt, ggplot2::aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
      ggplot2::geom_point() + ggplot2::scale_colour_manual(
        values = c("red", "darkorange1", "green1", "skyblue1", "navyblue", "grey")) +
      ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) + ggplot2::ylab("-log10(p-value)") +
      ggplot2::scale_y_reverse() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                                axis.text.x=ggplot2::element_blank(),
                                axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlim(start,end) + ggplot2::ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))


    b = ggplot2::ggplot(data = in.dt.2, ggplot2::aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
      ggplot2::geom_point() + ggplot2::scale_colour_manual(
        values = c("red", "darkorange1", "green1", "skyblue1", "navyblue", "grey")) +
      ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) +
      ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") +
      ggplot2::xlim(start,end) + ggplot2::ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P)) +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank())

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
  }else if(ldby == "none"){
    a = ggplot2::ggplot(in.dt, ggplot2::aes(x = CHR_POS, y = LOG10P)) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position")) +
      ggplot2::ylab("-log10(p-value)") +
      ggplot2::scale_y_reverse() + theme(axis.title.x=ggplot2::element_blank(),
                                axis.text.x=ggplot2::element_blank(),
                                axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlim(start,end) + ggplot2::ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))


    b = ggplot2::ggplot(in.dt.2, ggplot2::aes(x = CHR_POS, y = LOG10P)) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab(paste0("Chromosome ", chr_in, " Position (Mbp)")) +
      ggplot2::ylab("-log10(p-value)") + ggplot2::theme(legend.position = "bottom") +
      ggplot2::xlim(start,end) + ggplot2::ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P)) +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank())

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
}


