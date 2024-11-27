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

#' Map Helix Alignment to Structure
#'
#' @param template_ranges the dataframe returned by createHelixAlignments()
#'
#' @param rcsb_id the template rhodopsin to map to
#'
#' @param query_num to highlight
#'
#' @return r3dmol object
#'
#' @export
#'
#' @importFrom dplyr %>% filter
#' @importFrom r3dmol r3dmol m_add_model m_set_style m_sel m_style_cartoon m_zoom_to m_render
visualizeHelixMapping <- function(template_ranges, template, rcsb_id, query_num){
  query <- paste("Query", query_num)

  mapping_df <- templateMapping(template_ranges = template_ranges,
                                template = template,
                                rcsb_id = rcsb_id) %>%
    dplyr::filter(Query == query)

  struct <- getPDBstruct(rcsb_id = rcsb_id)


  viewer <- r3dmol::r3dmol() %>%
    r3dmol::m_add_model(r3dmol::m_bio3d(struct)) %>% # Load the PDB file
    r3dmol::m_set_style(style = r3dmol::m_style_cartoon()) # Cartoon representation


  for (i in 1:nrow(mapping_df)){
    viewer <- viewer %>%
      r3dmol::m_set_style(
        sel = r3dmol::m_sel(resi = mapping_df$mapped_start[i]:mapping_df$mapped_end[i]),
        style = r3dmol::m_style_cartoon(color = "blue")
      )
  }

  viewer <- viewer %>%
    r3dmol::m_zoom_to() %>%
    r3dmol::m_render()

  return(viewer)
}
# [END]
