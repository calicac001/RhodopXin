#' Visualize Helix Alignments
#'
#' A function that takes the output of createHelixAlignments() and create visual
#' representations of the alignments.
#'
#' @param all_pwa a list of AAStringSet objects containing the combined pairwise
#' alignments of each template's helix. Each index in the list corresponds to the
#' helix where that alignment was made e.g. the AAStringSet at index 1 corresponds
#' to the combined pairwise alignments of all the query sequences against Helix 1
#'
#' @return a combined visual of all the helix alignments
#'
#' @references ref here
#'
#' @export
#'
#' @importFrom aplot as.patchwork
#' @importFrom ggplot2 theme
#' @importFrom ggmsa ggmsa geom_seqlogo geom_msa
#' @importFrom patchwork wrap_plots plot_layout
visualizeHelixAlignments <- function(all_pwa){
  plots <- list()
  for(i in 1:length(all_pwa)){
    plot <- ggmsa::ggmsa(all_pwa[[i]], seq_name = TRUE, color = "Hydrophobicity")  +
      ggmsa::geom_seqlogo(color = "Hydrophobicity") +
      ggmsa::geom_msaBar()

    plots[[i]] <- aplot::as.patchwork(plot) +
      ggplot2::theme(plot.margin = unit(c(0,50,50,0), "pt"))
  }
  num_rows <- ceiling(length(all_pwa)/2)
  combined <- patchwork::wrap_plots(plots) +
    patchwork::plot_layout(ncol = 2, nrow = num_rows, heights = rep(1, num_rows))

  return(combined)
}
