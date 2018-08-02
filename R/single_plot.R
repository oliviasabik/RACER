#' Single Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values) of an association study
#' by their genomic position, for example, the results of a GWAS or eQTL study. Sources
#' 1000K genomes phase III data for linkage disequilibrium calculations.
#' @param assoc_data required. A dataframe that has columns named
#' CHR_POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values),
#' CHR, which contains the chromosome number, RS_ID, which
#' contains the rsID numbers of the SNPs, and R2, which contains
#' LD information. If you want to compute LD information for the dataset,
#' add a column to the input called LD, and fill it with NA, and include the
#' rs_id numbers for the lead SNPs you want to use to compute the LD.
#' @param chr required. chromosome to plot
#' @param rs_id optional. a character string. the rs_id number of the lead snp
#' used for LD calculations, set NULL if input has LD information
#' @param pops optional. required if an rs_id is specified. populations from the
#' 1000K Genomes dataset to use for LD calculations, set NULL if input has LD information
#' @param start required. starting chromosomal position of the association plot
#' @param end required. ending chromosomal position of the association plot
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @examples
#' single_plot_function(assoc_data = assoc_data, rs_id = "rs00000001", pops = c("CEU","TSI","FIN","GBR","IBS"), start = 1000000, end = 1230000)

single_plot_function <- function(assoc_data, chr, start=NULL, end=NULL, rs_id=NULL, pops=NULL){
  print("Calculating LD")
  if(!is.null(pops)){
    if(length(pops) == 1){
      ld_command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id,
                          "&pop=",pops,"&r2_d=r2'")
      z = system(ld_command, intern = TRUE)
      z = as.data.frame(z)
      z = tidyr::separate(z, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                         "col_9", "col_10"), sep = "\t")
      colnames(z) = z[1,]
      z = z[-1,]
      z = dplyr::select(z, RS_Number, R2)
      colnames(z) = c("RS_ID", "LD")
      assoc_data$LD = NA
      assoc_data = dplyr::select(assoc_data, -LD)
      assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
    }else if(length(pops) > 1){
      ld_command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id,
                        "&pop=")
      for (i in 1:length(pops)){
        if (i < length(pops)){
          print(i)
          a = pops[i]
          ld_command = paste0(ld_command, a, "%2B")
        } else if(i == length(pops)){
          a = pops[i]
          ld_command = paste0(ld_command, a)
          }
        }
      ld_command = paste0(ld_command, "&r2_d=r2'")
      z = system(ld_command, intern = TRUE)
      z = as.data.frame(z)
      z = tidyr::separate(z, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                         "col_9", "col_10"), sep = "\t")
      colnames(z) = z[1,]
      z = z[-1,]
      z = dplyr::select(z, RS_Number, R2)
      colnames(z) = c("RS_ID", "LD")
      assoc_data$LD = NA
      assoc_data = dplyr::select(assoc_data, -LD)
      assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
    }
  }
  print("Formatting Data")
  in.dt <- as.data.frame(assoc_data)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$LD = as.numeric(as.character(in.dt$LD))
  in.dt$LD_BIN <- cut(in.dt$LD,
                   breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                   labels=c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))
  in.dt$LD_BIN = as.character(in.dt$LD_BIN)
  in.dt$LD_BIN[is.na(in.dt$LD_BIN)] <- "NA"
  in.dt$LD_BIN = as.factor(in.dt$LD_BIN)
  in.dt = dplyr::filter(in.dt, CHR == chr)
  in.dt = dplyr::filter(in.dt, CHR_POS > start)
  in.dt = dplyr::filter(in.dt, CHR_POS < end)
  print("Generating Plot")

  # reading in gene data
  data(biomart_hg19)
  colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
  dat_chr = in.dt$CHR[1]
  gene = subset(biomart_hg19, biomart_hg19$CHR == dat_chr)
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
          axis.ticks.y = element_blank()) + xlab(paste0("Chromosome ", dat_chr, " Position")) +
    ylim(0,(max(gene$y_value)+1))

  b = ggplot2::ggplot(data = in.dt, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("navyblue", "skyblue1", "green1", "darkorange1", "red", "grey")) +
    theme_bw() + xlab("Chromosome Position") + ylab("-log10(p-value)") +
    xlim(start, end) + ylim(min(in.dt$LOG10P),max(in.dt$LOG10P))

  ggpubr::ggarrange(b, c, heights = c(4,1), nrow = 2, ncol = 1,
                    common.legend = TRUE, legend = "right")
}

