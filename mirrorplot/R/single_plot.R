#' Single Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data a data frmae that has columns named CHR_POS representing
#' the position of the snp on the chromosome, LOG10P which contains the
#' -log10(P-values), and optional columns CHR, which contains the chromosome
#' number, RS_ID, which contains the rsID numbers of the SNPs, and R2, which
#' contains LD information.
#' @keywords association plot
#' @export
#' @examples
#' single_plot_function()

single_plot_function <- function(assoc_data){
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
  ggplot2::ggplot(data = in.dt, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("navyblue", "skyblue1", "green1", "darkorange1", "red", "grey")) +
    theme_bw() + xlab("Chromosome Position") + ylab("-log10(p-value)")
}

