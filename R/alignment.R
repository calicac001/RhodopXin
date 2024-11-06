#' PW Align Helices with Query Sequences
#'
#' @param template AAStringSet of template rhodopsin to get helices from
#'
#' @param sequences AAStringSet of query sequences to align to helices
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' Aboyoun P, Gentleman R (2024). _pwalign: Perform pairwise sequence alignments_.doi:10.18129/B9.bioc.pwalign <https://doi.org/10.18129/B9.bioc.pwalign>, R package version 1.2.0, <https://bioconductor.org/packages/pwalign>.
#'
#' @export
#' @import Biostrings
#' @import dplyr
#' @import pwalign

createPWAlignments <- function(template, sequences){
  helices_seq <- helixSegments(template = template)

  for (i in seq_along(helices_seq)){
    for (j in seq_along(sequences)){
      print(paste("Alignment with Helix", LETTERS[i], "and Query sequence", j))
      print(pwalign::pairwiseAlignment(helices_seq[[i]], sequences[[j]], type = "local"))
    }
  }
}

#' Get Helix Segments of Sequence
#'
#' @param template
#'
#' @return dataframe
#'
#' #' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @export
#' @import Biostrings
helixSegments <- function(template){
  template_id <- substr(names(template)[[1]], 1, 4)
  helices_df <- findHelices(template_id)

  segments <- list()
  for (i in 1:nrow(helices_df)) {
    start_pos <- helices_df$start[i]
    end_pos <- helices_df$end[i]

    segment <- subseq(template[[1]], start=start_pos, end=end_pos)

    segments[[i]] <- Biostrings::AAStringSet(as.character(segment))
  }

  return(segments)
}
