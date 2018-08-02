#' Mirror Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data1 required. A dataframe that has columns named
#' CHR_POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values),
#' CHR, which contains the chromosome number, RS_ID, which
#' contains the rsID numbers of the SNPs, and LD, which contains
#' LD information. If no column named LD is in the input, LD will
#' be calculated from 1000 genomes phase III in relation to rs_id_1.
#'
#' @param assoc_data2 required. identifcal in format to assoc_data1
#' @param chr required. chromosome you wish to plot
#' @param name1 optional. name of association set 1
#' @param name2 optional. name of association set 2
#' @param pops optional. required is rs_ids are specified. populations used
#' to calculate LD
#' @param rs_id_1 optional. lead snp used to calculate LD for first association dataset and used for coordinates
#' if plotby = "snp"
#' @param rs_id_2 optional. lead snp used to calculate LD for second association dataset
#' @param plotby optional. "coord", "gene", or "snp". Which parameter to use to
#' determine the reigon to be plotted.
#' @param gene optional. If "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param snp optional. If "snp" selected for plotby, then plot will be +/- 50kb of snp
#' @param start optional. If "coord" selected for plotby, then this will be lower bound of x axis
#' @param end optional. If "coord" selected for plotby, then this will be upper bound of x axis
#' @keywords association plot
#' @export
#' @examples
#' mirror_plot_function(assoc_data1, assoc_data2, "plot1", "plot2", 1000000, 110000000, pops = c("CEU","TSI","FIN","GBR","IBS"), rs_id_1 = "rs123456", rs_id_2 = "rs123456",  plotby = "coord")

mirror_plot_function <- function(assoc_data1, assoc_data2, chr, name1=NULL, name2=NULL, start=NULL, end=NULL, pops=NULL, rs_id_1=NULL, rs_id_2=NULL, plotby = NULL, gene = NULL){
  reqs = c("CHR", "CHR_POS", "LOG10P")
  cols_1 = colnames(assoc_data1)
  cols_2 = colnames(assoc_data2)
  if(sum(reqs %in% cols_1) == 3){
  }else{print("Association Data Set #1 is missing a required column.")}
  if(sum(reqs %in% cols_2) == 3){
  }else{print("Association Data Set #2 is missing a required column.")}

  data(biomart_hg19)
  chr_2 = chr
  chr_1 = chr
  colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
  gene_sub = subset(biomart_hg19, biomart_hg19$CHR == chr_2)

  if(sum(is.null(plotby)) == 0){
    print("Plotting by...")
    if((plotby == "coord") == TRUE){
      print("coord")
      start = start
      end = end
    }else if((plotby == "gene") == TRUE){
      print(paste("gene:",gene)
      if(sum(is.null(gene)) == 0){
        p = subset(gene_sub, gene_sub$GENE_NAME == gene)
        start = min(p$TRX_START) - 500000
        end = max(p$TRX_END) + 500000
      }else{print("No gene specified.")}
    }else if((plotby == "snp") == TRUE){
      print(paste("snp",rs_id_1))
      q = subset(assoc_data1, RS_ID == rs_id_1)
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
  gene_sub = arrange(gene_sub, desc(LENGTH))
  gene_sub = gene_sub[!duplicated(gene_sub$GENE_ID),]
  gene_sub = gene_sub[,c(3,4,5)]
  gene_sub = reshape2::melt(gene_sub)
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

  # read in, format, and filter data sets
  print("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_1)
  in.dt = dplyr::filter(in.dt, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$CHR_POS = as.numeric(as.character(in.dt.2$CHR_POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  in.dt.2 = dplyr::filter(in.dt.2, CHR == chr_2)
  in.dt.2= dplyr::filter(in.dt.2, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  # determining input for LD info
  if("LD" %in% colnames(in.dt) | "LD" %in% colnames(in.dt.2)){
    print("Using LD info from the input data set.")
  }else{
    print(paste0("Pulling LD data from 1000 Genomes Phase III data using", rs_id_1))
    if(!is.null(pops)){
      if(length(pops) == 1){
        ld_command_1 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id_1,
                              "&pop=",pops,"&r2_d=r2'")
        z_1 = system(ld_command_1, intern = TRUE)
        z_1 = as.data.frame(z_1)
        z_1 = tidyr::separate(z_1, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                               "col_9", "col_10"), sep = "\t")
        colnames(z_1) = z_1[1,]
        z_1 = z_1[-1,]
        z_1 = dplyr::select(z_1, RS_Number, LD)
        colnames(z_1) = c("RS_ID", "LD")
        assoc_data1$LD = NA
        assoc_data1 = dplyr::select(assoc_data1, -LD)
        assoc_data1 = merge(assoc_data1, z_1, by = "RS_ID", all.x = TRUE)
      }else if(length(pops) > 1){
        ld_command_1 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id_1,
                              "&pop=")
        for (i in 1:length(pops)){
          if (i < length(pops)){
            print(i)
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
        z_1 = dplyr::select(z_1, RS_Number, LD)
        colnames(z_1) = c("RS_ID", "LD")
        assoc_data1$LD = NA
        assoc_data1 = dplyr::select(assoc_data1, -LD)
        assoc_data1 = merge(assoc_data1, z_1, by = "RS_ID", all.x = TRUE)
      }
    }

    if(!is.null(pops)){
      if(length(pops) == 1){
        ld_command_2 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id_1,
                              "&pop=",pops,"&r2_d=r2'")
        z_2 = system(ld_command_2, intern = TRUE)
        z_2 = as.data.frame(z_2)
        z_2 = tidyr::separate(z_2, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                               "col_9", "col_10"), sep = "\t")
        colnames(z_2) = z_2[1,]
        z_2 = z_2[-1,]
        z_2 = dplyr::select(z_2, RS_Number, LD)
        colnames(z_2) = c("RS_ID", "LD")
        assoc_data2$LD = NA
        assoc_data2 = dplyr::select(assoc_data2, -LD)
        assoc_data2 = merge(assoc_data2, z_2, by = "RS_ID", all.x = TRUE)
      }else if(length(pops) > 1){
        ld_command_2 = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id_1,
                              "&pop=")
        for (i in 1:length(pops)){
          if (i < length(pops)){
            print(i)
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
        z_2 = dplyr::select(z_2, RS_Number, LD)
        colnames(z_2) = c("RS_ID", "LD")
        assoc_data2$LD = NA
        assoc_data2 = dplyr::select(assoc_data2, -LD)
        assoc_data2 = merge(assoc_data2, z_2, by = "RS_ID", all.x = TRUE)
      }
    }
  }

  #processing LD information
  in.dt$LD = as.numeric(as.character(in.dt$LD))
  in.dt$LD_BIN <- cut(in.dt$LD,
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                      labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
  in.dt$LD_BIN = as.character(in.dt$LD_BIN)
  in.dt$LD_BIN[is.na(in.dt$LD_BIN)] <- "NA"
  in.dt$LD_BIN = as.factor(in.dt$LD_BIN)
  in.dt$LD_BIN = factor(in.dt$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
  print(table(in.dt$LD_BIN))


  in.dt.2$LD = as.numeric(as.character(in.dt.2$LD))
  in.dt.2$LD_BIN <- cut(in.dt.2$LD,
                        breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                        labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
  in.dt.2$LD_BIN = as.character(in.dt.2$LD_BIN)
  in.dt.2$LD_BIN[is.na(in.dt.2$LD_BIN)] <- "NA"
  in.dt.2$LD_BIN = as.factor(in.dt.2$LD_BIN)
  in.dt.2$LD_BIN = factor(in.dt.2$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
  print(table(in.dt.2$LD_BIN))

  # generate mirror plot
  print("Generating plots")
  a = ggplot2::ggplot(data = in.dt, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("red", "darkorange1", "green1", "skyblue1", "navyblue", "grey")) +
    theme_bw() + xlab(paste0("Chromosome ", chr_1, " Position")) + ylab("-log10(p-value)") +
    scale_y_reverse() + theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank()) +
    theme(legend.position = "none") +
    xlim(start,end) + ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))


  b = ggplot2::ggplot(data = in.dt.2, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("red", "darkorange1", "green1", "skyblue1", "navyblue", "grey")) +
    theme_bw() + xlab(paste0("Chromosome ", chr_2, " Position (Mbp)")) + ylab("-log10(p-value)") +
    theme(legend.position = "bottom") +
    xlim(start,end) + ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P)) +
    theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank())

  c = ggplot(gene_sub, aes(x = value, y = y_value)) +
    geom_line(aes(group = GENE_NAME), size = 2) + theme_bw() +
    geom_text(data = plot_lab, aes(x = value, y = y_value, label = GENE_NAME),
              hjust = -0.1,vjust = 0.3, size = 2.5) + xlim(start,end) +
    theme(axis.title.y = element_text(color = "white", size = 28),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + xlab(paste0("Chromosome ", chr_1, " Position")) +
    ylim(0,(max(gene_sub$y_value)+1))

  ggpubr::ggarrange(a, b, c, heights = c(2,2,1), nrow = 3, ncol = 1,
                    common.legend = TRUE, legend = "right")
}


