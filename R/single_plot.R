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
#' add a column to the input called LD, and fill it with NA, and include the
#' rs_id numbers for the lead SNPs you want to use to compute the LD.
#' @param chr required. chromosome to plot
#' @param rs_id optional. a character string. the rs_id number of the lead snp
#' used for LD calculations, and can be used to set the plotting coordinates
#' @param pops optional. required if an rs_id is specified. populations from the
#' 1000K Genomes dataset to use for LD calculations, set NULL if input has LD information
#' @param plotby required. "coord", "gene", or "snp". Which parameter to use to
#' determine the reigon to be plotted.
#' @param gene optional. If "gene" selected for plotby, then plot will be +/- 50kb of gene
#' @param start optional. If "coord" selected for plotby, then this will be lower bound of x axis
#' @param end optional. If "coord" selected for plotby, then this will be upper bound of x axis
#' @keywords association plot, gwas, linkage disequilibrium.
#' @export
#' @examples
#' single_plot_function(assoc_data = assoc_data, chr = 1, rs_id = "rs00000001", pops = c("CEU","TSI","FIN","GBR","IBS"), start = 1000000, end = 1230000)

single_plot_function <- function(assoc_data, chr, plotby, start=NULL, end=NULL, rs_id=NULL, pops=NULL, gene=NULL){
  reqs = c("CHR", "CHR_POS", "LOG10P")
  cols = colnames(assoc_data)
  if(sum(reqs %in% cols) == 3){
  }else{print("Association Data Set is missing a required column.")}

  data(biomart_hg19)
  chr_in = chr
  colnames(biomart_hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "GENE_NAME", "LENGTH")
  gene_sub = subset(biomart_hg19, biomart_hg19$CHR == chr_in)

  if(sum(is.null(plotby)) == 0){
    print("Plotting by...")
    if((plotby == "coord") == TRUE){
      print("coord")
      start = start
      end = end
    }else if((plotby == "gene") == TRUE){
      print(paste("gene:",gene))
            if(sum(is.null(gene)) == 0){
              p = subset(gene_sub, gene_sub$GENE_NAME == gene)
              start = min(p$TRX_START) - 500000
              end = max(p$TRX_END) + 500000
            }else{print("No gene specified.")}
    }else if((plotby == "snp") == TRUE){
      print(paste("snp",rs_id))
      q = subset(assoc_data1, RS_ID == rs_id)
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
  in.dt <- as.data.frame(assoc_data)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter(in.dt, CHR == chr_in)
  in.dt = dplyr::filter(in.dt, CHR_POS > start)%>%
    dplyr::filter(CHR_POS < end)

  # calculate LD
  print("Calculating LD")
  if("LD" %in% colnames(in.dt) | "LD" %in% colnames(in.dt)){
    print("Using LD info from the input data set.")
  }else{if(!is.null(pops)){
    if(length(pops) == 1){
      ld_command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=", rs_id,
                          "&pop=",pops,"&r2_d=r2'")
      z = system(ld_command, intern = TRUE)
      z = as.data.frame(z)
      z = tidyr::separate(z, 1, into = c("col_1", "col_2", "col_3", "col_4", "col_5", "col_6", "col_7", "col_8",
                                         "col_9", "col_10"), sep = "\t")
      colnames(z) = z[1,]
      z = z[-1,]
      z = dplyr::select(z, RS_Number, LD)
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
      colnames(z) = c("RS_ID", "LD")
      assoc_data$LD = NA
      assoc_data = dplyr::select(assoc_data, -LD)
      assoc_data = merge(assoc_data, z, by = "RS_ID", all.x = TRUE)
    }
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

  # Generate plots
  print("Generating Plots.")
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

