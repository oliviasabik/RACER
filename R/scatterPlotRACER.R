#' Scatter Plot -- Regional Association ComparER Plot
#'
#' This function allows you to create a scatter plot of -log10(P-values).
#' @param assoc_data1 required. A dataframe that has columns named
#' POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values),
#' CHR, which contains the chromosome number, RS_ID,
#' which contains rsIDs for the SNPs in the data frame. If you
#' have processed your data using formatRACER() these columns will
#' be properly labeled
#'
#' @param assoc_data2 required. identifcal in format to assoc_data1
#' @param chr required. chromosome you wish to plot
#' @param name1 optional. name of association set 1
#' @param name2 optional. name of association set 2
#' @param region_start start coordinates on chr to be compared
#' @param region_end end coordinates on the chr to be compared
#' @param ld_df data frame containing the LD data to use to color the plot
#' @param label optional. If TRUE, will add a label to a the maximum combined LOG10P of the plot
#'
#' @keywords association plot
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' data("mark3_eqtl")
#' data("mark3_bmd_gwas")
#' df1_f = formatRACER(assoc_data = mark3_bmd_gwas, chr_col = 3, pos_col = 4, p_col = 11,rs_col = 2)
#' df2_f = formatRACER(assoc_data = mark3_eqtl, chr_col = 10, pos_col = 11, p_col = 7,rs_col = 15)
#' df1_f_ld = ldRACER(assoc_data = df1_f, rs_col = 2, pops = "EUR", lead_snp = "rs11623869")
#' df2_f_ld = ldRACER(assoc_data = df2_f, rs_col = 15, pops = "EUR", lead_snp = "rs11623869")
#' scatterPlotRACER(assoc_data1 = df1_f_ld, assoc_data2 = df2_f_ld, chr = 14,
#' name1 = "Mark3_GWAS", name2 = "Mark3_eQTL",
#' region_start = 103750000, region_end = 104250000, ld_df = 1)}

scatterPlotRACER <- function(assoc_data1, assoc_data2, chr, name1="Association Dataset #1", name2="Association Dataset #2", region_start, region_end, ld_df = NULL, label = FALSE){
  reqs = c("CHR", "POS", "LOG10P", "RS_ID")
  cols_1 = colnames(assoc_data1)
  cols_2 = colnames(assoc_data2)
  if(sum(reqs %in% cols_1) == 4){
  }else{stop("Association Data Set #1 is missing a required column.")}
  if(sum(reqs %in% cols_2) == 4){
  }else{stop("Association Data Set #2 is missing a required column.")}

  `%>%` <- magrittr::`%>%`

  message("Reading in association data")
  in.dt <- as.data.frame(assoc_data1)
  in.dt$POS = as.numeric(as.character(in.dt$POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$CHR = as.numeric(as.character(in.dt$CHR))
  in.dt = dplyr::filter_(in.dt, ~CHR == chr)
  in.dt = dplyr::filter_(in.dt, ~POS > region_start)%>%
    dplyr::filter_(~POS < region_end)

  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt.2$POS = as.numeric(as.character(in.dt.2$POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$CHR = as.numeric(as.character(in.dt.2$CHR))
  in.dt.2 = dplyr::filter_(in.dt.2, ~CHR == chr)
  in.dt.2= dplyr::filter_(in.dt.2, ~POS > region_start)%>%
    dplyr::filter_(~POS < region_end)

  if(ld_df > 0){
    if(ld_df == 1){
      in.dt.final = dplyr::select(in.dt, "RS_ID", "LOG10P", "LD", "LD_BIN")
      colnames(in.dt.final) = c("RS_ID", "LOG10P1", "LD", "LD_BIN")
      in.dt.2.final = dplyr::select(in.dt.2, "RS_ID", "LOG10P")
      colnames(in.dt.2.final) = c("RS_ID", "LOG10P2")
    }else if(ld_df == 2){
      in.dt.2.final = dplyr::select(in.dt.2, "RS_ID", "LOG10P", "LD", "LD_BIN")
      colnames(in.dt.2.final) = c("RS_ID", "LOG10P2", "LD", "LD_BIN")
      in.dt.final = dplyr::select(in.dt, "RS_ID", "LOG10P")
      colnames(in.dt.final) = c("RS_ID", "LOG10P1")
    }
  }else{
    in.dt.final = dplyr::select(in.dt, "RS_ID", "LOG10P")
    colnames(in.dt.final) = c("RS_ID", "LOG10P1")
    in.dt.2.final = dplyr::select(in.dt.2, "RS_ID", "LOG10P")
    colnames(in.dt.2.final) = c("RS_ID", "LOG10P2")
  }
  df_plot = merge(in.dt.final, in.dt.2.final, by = 'RS_ID')

  lab.in = df_plot[which.max(df_plot$LOG10P1 + df_plot$LOG10P2),]
  print(lab.in)

  message("Generating plot.")
  if(ld_df > 0){
  ggplot2::ggplot(df_plot, aes_string(x = "LOG10P1", y = "LOG10P2", color = "LD_BIN")) +
    ggplot2::geom_point() + ggplot2::xlab(paste0("LOG10P for ", name1)) +
    ggplot2::ylab(paste0("LOG10P for ", name2)) + ggplot2::theme_bw() + ggplot2::scale_colour_manual(
      values = c("1.0-0.8" = "red", "0.8-0.6" = "darkorange1", "0.6-0.4" = "green1",
                 "0.4-0.2" = "skyblue1", "0.2-0.0" = "navyblue", "NA" = "grey"), drop = FALSE) +
      ggplot2::geom_point(data = lab.in, color = "purple") +
      geom_text(data = lab.in, aes(label = RS_ID), color = "black", size = 3, hjust = 1.25)
  }else{
    ggplot2::ggplot(df_plot, aes_string(x = "LOG10P1", y = "LOG10P2")) +
      ggplot2::geom_point() + ggplot2::xlab(paste0("LOG10P for ", name1)) +
      ggplot2::ylab(paste0("LOG10P for ", name2)) + ggplot2::theme_bw() +
      ggplot2::geom_point(data = lab.in, color = "purple") +
      geom_text(data = lab.in, aes(label = RS_ID), color = "black", size = 3, hjust = 1.25)
  }

}

