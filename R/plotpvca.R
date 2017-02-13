#' A function for plotting a result from the PVCA function
#'
#' This function takes as input a proportion vector calculated from the 'pvca' function and generate a bar plot.
#'
#' @param pvca.res The output from the 'pvca' function
#' @param title  The title of the bar plot
#'
#' @return p A ggplot2 object
#'
#' @export
#'

PlotPVCA <- function(pvca.res, title){
  plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
  p <- ggplot2::ggplot(plot.dat, aes(x=eff, y=prop))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::geom_bar(stat="identity", fill="steelblue", colour="steelblue")
  p <- p + ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4)
  p <- p + ggplot2::scale_x_discrete(limits=names(pvca.res))
  p <- p + ggplot2::scale_y_continuous(limits = c(0,1))
  p <- p + ggplot2::labs(x= "Effects", y= "Weighted average proportion variance")
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(plot.background = element_blank() ,panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank() ,panel.border = element_blank(), panel.background = element_blank())
  p <- p + ggplot2::theme(axis.line = element_line(color = 'black'))
  p <- p + ggplot2::theme(axis.title.x = element_text(size = 15, vjust=-0.5))
  p <- p + ggplot2::theme(axis.title.y = element_text(size = 15, vjust= 1.0))
  p <- p + ggplot2::theme(axis.text = element_text(size = 12))
  p <- p + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust= 0.5, hjust=1))
  p
}

