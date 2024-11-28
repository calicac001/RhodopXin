#' Visualize Helix Alignments
#'
#' A function that takes the output of createHelixAlignments() and create visual
#' representations of the alignments. This combines all of the helices alignments
#' in one plot. For each helix plot, we have the helix sequence of the query aligned with
#' the query sequences. The left hand side of this plot shows which residue number
#' this alignment comes from as well as which query. Additionally, a sequence logo
#' and a consensus plot of the alignment is present on the top and bottom of the plot.
#'
#' @param all_pwa the all_pwa output of createHelixAlignments which is a list of
#' AAStringSet objects containing the combined pairwise alignments of each template's
#' helix.
#'
#' @return a combined visual of all the helix alignments in one plot.
#'
#' @examples
#' # Create the pairwise alignments first
#' template <- template_rhodopsins[1]
#' sequences <- sample_rhodopsins
#' rcsb_id <- "1QHJ"
#' pwa_results <- createHelixAlignments(template = template, sequences = sequences, rcsb_id = rcsb_id)
#'
#' # Visualize the alignments
#' all_pwa <- pwa_results$all_pwa
#' result <- visualizeHelixAlignments(all_pwa = all_pwa)
#'
#' # Note that output may look distorted/squished when run in a markdown file,
#' # best to view it in the shiny app of this package. If not, adjust the
#' fig.height or the R chunk.
#'
#' @references
#' H. Wickham. _ggplot2: Elegant Graphics for Data Analysis_. Springer-Verlag New York, 2016.
#'
#' L Zhou, T Feng, S Xu, F Gao, TT Lam, Q Wang, T Wu, H Huang, L Zhan, L Li, Y Guan, Z Dai, G
#' Yu. _ggmsa: a visual exploration tool for multiple sequence alignment and associated data_.
#' Bioinformatics. 2022, 23(4):bbac222. 10.1093/bib/bbac222
#'
#' Pedersen T (2024). _patchwork: The Composer of Plots_. R package version 1.3.0,
#' <https://CRAN.R-project.org/package=patchwork>.
#'
#' Yu G (2024). _aplot: Decorate a 'ggplot' with Associated Information_. R package version
#' 0.2.3, <https://CRAN.R-project.org/package=aplot>.
#'
#' @export
#'
#' @importFrom aplot as.patchwork
#' @importFrom ggplot2 theme
#' @importFrom ggmsa ggmsa geom_seqlogo geom_msa
#' @importFrom patchwork wrap_plots plot_layout
visualizeHelixAlignments <- function(all_pwa){
  # Create list to store plots in
  plots <- list()

  # something buggy with ggmsa, warning keeps appearing but has no effect on output
  suppressWarnings({
  # Iterate through each alignment i.e. each AAStringSet object
  for(i in 1:length(all_pwa)){
    # Create the alignment plot, sequence logo and consensus plot
    plot <- ggmsa::ggmsa(all_pwa[[i]], seq_name = TRUE, color = "Hydrophobicity")  +
      ggmsa::geom_seqlogo(color = "Hydrophobicity") +
      ggmsa::geom_msaBar()

    # Append to the list, convert to patchwork so they can be combined,
    # adjust margins for spacing
    plots[[i]] <- aplot::as.patchwork(plot) +
      ggplot2::theme(plot.margin = unit(c(0,50,50,0), "pt"))
  }
  # Display the plots in cols of 2 so count number of rows
  num_rows <- ceiling(length(all_pwa)/2)

  # Combine plots using patchwork in columns of 2
  combined <- patchwork::wrap_plots(plots) +
    patchwork::plot_layout(ncol = 2, nrow = num_rows, heights = rep(1, num_rows))
  })
  return(combined)
}

#' Map Helix Alignment to Structure
#'
#' Given the alignments produced from createHelixAlignments, map the conserved
#' residues in to the 3D structure of the given template rhodopsin.
#'
#' @param template_ranges the dataframe template_ranges returned by createHelixAlignments()
#'
#' @param rcsb_id tthe RCSB accession code of the template rhodopsin
#'
#' @param query_num the number of query which to highligh in the structure. Follows
#' the order in the FASTA file
#'
#' @return r3dmol object is the viewer with the loaded 3d structure of the template
#' and the aligned residues of the query highlighted.
#'
#' @examples
#' # Create the pairwise alignments first
#' template <- template_rhodopsins[1]
#' sequences <- sample_rhodopsins
#' rcsb_id <- "1QHJ"
#' pwa_results <- createHelixAlignments(template = template, sequences = sequences, rcsb_id = rcsb_id)
#'
#' # Visualize the alignments
#' template_ranges <- pwa_results$template_ranges
#' result <- visualizeHelixMapping(template_ranges = template_ranges, template = template, rcsb_id = rcsb_id, query_num = 1)
#'
#' @references
#' Su W, Johnston B (2024). _r3dmol: Create Interactive 3D Visualizations of Molecular Data_.
#' R package version 0.2.0, commit 0eefc1c2c430fff0122ae71ae8c9cc4e8912190c,
#' <https://github.com/swsoyee/r3dmol>.
#'
#' Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data_
#' _Manipulation_. R package version 1.1.4, <https://CRAN.R-project.org/package=dplyr>.
#'
#' @export
#'
#' @importFrom dplyr %>% filter
#' @importFrom r3dmol r3dmol m_add_model m_set_style m_sel m_style_cartoon m_zoom_to m_render m_add_res_labels
visualizeHelixMapping <- function(template_ranges, template, rcsb_id, query_num){
  # Determine which query to highlight
  query <- paste("Query", query_num)

  # Get the mapping dataframe and filter by query
  mapping_df <- templateMapping(template_ranges = template_ranges,
                                template = template,
                                rcsb_id = rcsb_id) %>%
    dplyr::filter(Query == query)

  # Get 3d structure of template
  struct <- getPDBstruct(rcsb_id = rcsb_id)

  # Initialize viewer and load template structure
  viewer <- r3dmol::r3dmol() %>%
    r3dmol::m_add_model(r3dmol::m_bio3d(struct)) %>% # Load the PDB file
    r3dmol::m_set_style(style = r3dmol::m_style_cartoon()) # Cartoon representation

  # For each mapped alignment, color this position in the structure
  for (i in 1:nrow(mapping_df)){
    viewer <- viewer %>%
      r3dmol::m_set_style(
        sel = r3dmol::m_sel(resi = mapping_df$mapped_start[i]:mapping_df$mapped_end[i]),
        style = r3dmol::m_style_cartoon(color = "#e95420")
      )
  }

  # zoom in the structure and render
  viewer <- viewer %>%
    r3dmol::m_zoom_to() %>%
    r3dmol::m_render()

  return(viewer)
}

# [END]
