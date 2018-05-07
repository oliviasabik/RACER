#' Mirror Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data1 object in your current R environment that has
#' columns named CHR_POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values), and optional
#' columns CHR, which contains the chromosome number, RS_ID, which
#' contains the rsID numbers of the SNPs, and R2, which contains LD information.
#' @param assoc_data2 identifcal in format to 
#' @keywords association plot
#' @export
#' @examples
#' base_plot_function()

mirror_plot_function <- function(assoc_data1, assoc_data2){
  # read and process association data set 1
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
  
  # read and process association data set 2
  in.dt.2 <- as.data.frame(assoc_data2)
  in.dt$CHR_POS = as.numeric(as.character(in.dt$CHR_POS))
  in.dt$LOG10P = as.numeric(as.character(in.dt$LOG10P))
  in.dt$LD = as.numeric(as.character(in.dt$LD))
  in.dt$LD_BIN <- cut(in.dt$LD,
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                      labels=c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))
  in.dt$LD_BIN = as.character(in.dt$LD_BIN)
  in.dt$LD_BIN[is.na(in.dt$LD_BIN)] <- "NA"
  in.dt$LD_BIN = as.factor(in.dt$LD_BIN)
  
 # generate mirror plot
  ggplot2::ggplot(data = in.dt, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("navyblue", "skyblue1", "green1", "darkorange1", "red", "grey")) +
    theme_bw() + xlab("Chromosome Position") + ylab("-log10(p-value)")
}

