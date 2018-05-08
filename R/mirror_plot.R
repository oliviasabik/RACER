#' Mirror Plot Function
#'
#' This function allows you to creat a plot of -log10(P-values).
#' @param assoc_data1 object in your current R environment that has
#' columns named CHR_POS representing the position of the snp on
#' the chromosome, LOG10P which contains the -log10(P-values), and optional
#' columns CHR, which contains the chromosome number, RS_ID, which
#' contains the rsID numbers of the SNPs, and R2, which contains LD information.
#' @param assoc_data2 identifcal in format to assoc_data1
#' @param name1 name of association set 1
#' @param name2 name of association set 2
#' @param x1 lower bound on x axis
#' @param x2 upper bound on x axis
#' @keywords association plot
#' @export
#' @examples
#' base_plot_function()

mirror_plot_function <- function(assoc_data1, assoc_data2, name1, name2, x1, x2){
  # read and process association data set 1
  in.dt <- as.data.frame(assoc_data1)
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
  in.dt.2$CHR_POS = as.numeric(as.character(in.dt.2$CHR_POS))
  in.dt.2$LOG10P = as.numeric(as.character(in.dt.2$LOG10P))
  in.dt.2$LD = as.numeric(as.character(in.dt.2$LD))
  in.dt.2$LD_BIN <- cut(in.dt.2$LD,
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                      labels=c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))
  in.dt.2$LD_BIN = as.character(in.dt.2$LD_BIN)
  in.dt.2$LD_BIN[is.na(in.dt.2$LD_BIN)] <- "NA"
  in.dt.2$LD_BIN = as.factor(in.dt.2$LD_BIN)

  # filter out region requested
  in.dt = dplyr::filter(in.dt, CHR_POS > x1)%>%
    dplyr::filter(CHR_POS < x2)
  in.dt.2= dplyr::filter(in.dt.2, CHR_POS > x1)%>%
    dplyr::filter(CHR_POS < x2)

 # generate mirror plot
  a = ggplot2::ggplot(data = in.dt, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("navyblue", "skyblue1", "green1", "darkorange1", "red", "grey")) +
    theme_bw() + xlab("Chromosome Position") + ylab("-log10(p-value)") +
    scale_y_reverse() + theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank()) +
    theme(legend.position = "none") +
    #annotate("text", x=((min(in.dt$CHR_POS)+max(in.dt$CHR_POS))/5), y=(max(in.dt$LOG10P)-0.1), label = name1, size = 4) +
    xlim(min(in.dt$CHR_POS),max(in.dt$CHR_POS)) + ggtitle(paste0("Mirror Plot of ", name1, " and ", name2 ))


  b = ggplot2::ggplot(data = in.dt.2, aes(x = CHR_POS, y = LOG10P, color = LD_BIN)) +
    geom_point() + scale_colour_manual(
      values = c("navyblue", "skyblue1", "green1", "darkorange1", "red", "grey")) +
    theme_bw() + xlab("Chromosome Position (Mbp)") + ylab("-log10(p-value)") +
    theme(legend.position = "bottom") +
    #annotate("text", x=((min(in.dt.2$CHR_POS)+max(in.dt.2$CHR_POS))/5), y=(max(in.dt.2$LOG10P)-0.1), label = name2, size = 4)+
    xlim(min(in.dt.2$CHR_POS),max(in.dt.2$CHR_POS)) + ylim(min(in.dt.2$LOG10P),max(in.dt.2$LOG10P))


  ggpubr::ggarrange(a, b, nrow = 2, ncol = 1,
                    common.legend = TRUE, legend = "right", heights = 1)
}


