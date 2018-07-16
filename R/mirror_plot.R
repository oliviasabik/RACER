#' Mirror Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data1 required. A dataframe that has columns named
#' CHR_POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values),
#' CHR, which contains the chromosome number, RS_ID, which
#' contains the rsID numbers of the SNPs, and R2, which contains
#' LD information. If you want to compute LD information for the dataset,
#' add a column to the input called LD, and fill it with NA, and include the
#' rs_id numbers for the lead SNPs you want to use to compute the LD.
#'
#' @param assoc_data2 required. identifcal in format to assoc_data1
#' @param name1 optional. name of association set 1
#' @param name2 optional. name of association set 2
#' @param start optional. lower bound on x axis
#' @param end optional. upper bound on x axis
#' @param pops optional. required is rs_ids are specified. populations used
#' to calculate LD
#' @param rs_id_1 optional. lead snp used to calculate LD for first association dataset
#' @param rs_id_2 optional. lead snp used to calculate LD for second association dataset
#'
#' @keywords association plot
#' @export
#' @examples
#' mirror_plot_function(assoc_data1, assoc_data2, "plot1", "plot2", 1000000, 110000000, pops = c("CEU","TSI","FIN","GBR","IBS"), rs_id1 = "rs123456", rs_id2 = "rs123456")

mirror_plot_function <- function(assoc_data1, assoc_data2, name1=NULL, name2=NULL, start=NULL, end=NULL, pops=NULL, rs_id_1=NULL, rs_id_2=NULL){
  print("Calculating LD")
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
      z_1 = dplyr::select(z_1, RS_Number, R2)
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
      z_1 = dplyr::select(z_1, RS_Number, R2)
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
      z_2 = dplyr::select(z_2, RS_Number, R2)
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
      z_2 = dplyr::select(z_2, RS_Number, R2)
      colnames(z_2) = c("RS_ID", "LD")
      assoc_data2$LD = NA
      assoc_data2 = dplyr::select(assoc_data2, -LD)
      assoc_data2 = merge(assoc_data2, z_2, by = "RS_ID", all.x = TRUE)
    }
  }

  # read in, format, and filter data sets
  print("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  chr_1 = min(in.dt$CHR)
  in.dt = dplyr::filter(in.dt, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$CHR_POS = as.numeric(as.character(in.dt.2$CHR_POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  chr_2 = min(in.dt.2$CHR)
  in.dt.2= dplyr::filter(in.dt.2, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  # process ld info
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

  # reading in gene data
  data(biomart_hg19)
  chr_2
  colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
  gene = subset(biomart_hg19, biomart_hg19$CHR == chr_2)
  gene = subset(gene, gene$TRX_START > (start-5000))
  gene = subset(gene, gene$TRX_END < (end+5000))
  gene = arrange(gene, desc(LENGTH))
  gene = gene[!duplicated(gene$GENE_ID),]
  gene = gene[,c(3,4,5)]
  gene = reshape2::melt(gene)
  gene$y_value = as.numeric(as.factor(gene$GENE_NAME))
  plot_lab = subset(gene, gene$variable == "TRX_END")
  c = ggplot(gene, aes(x = value, y = y_value)) +
    geom_line(aes(group = GENE_NAME), size = 2) + theme_bw() +
    geom_text(data = plot_lab, aes(x = value, y = y_value, label = GENE_NAME),
              hjust = -0.1,vjust = 0.3, size = 2.5) + xlim(start,end) +
    theme(axis.title.y = element_text(color = "white", size = 28),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + xlab(paste0("Chromosome ", chr_1, " Position")) +
    ylim(0,(max(gene$y_value)+1))
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
    #annotate("text", x=((min(in.dt$CHR_POS)+max(in.dt$CHR_POS))/5), y=(max(in.dt$LOG10P)-0.1), label = name1, size = 4) +
    xlim(start,end) + ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))


  b = ggplot2::ggplot(data = in.dt.2, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("red", "darkorange1", "green1", "skyblue1", "navyblue", "grey")) +
    theme_bw() + xlab(paste0("Chromosome ", chr_2, " Position (Mbp)")) + ylab("-log10(p-value)") +
    theme(legend.position = "bottom") +
    #annotate("text", x=((min(in.dt.2$CHR_POS)+max(in.dt.2$CHR_POS))/5), y=(max(in.dt.2$LOG10P)-0.1), label = name2, size = 4)+
    xlim(start,end) + ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P)) +
    theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank())

  #c = ggplot2::ggplot(genes, aes())
  ggpubr::ggarrange(a, b, c, heights = c(2,2,1), nrow = 3, ncol = 1,
                    common.legend = TRUE, legend = "right")
}


