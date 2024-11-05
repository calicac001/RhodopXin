#' Load Amino Acid Sequence Input
#'
#' A function that takes in a string of an amino acid corresponding to a rhodopsin
#' or a FASTA file path containing the sequence information and creates an AAString
#' object out of it
#'
#' @param inputStr A string input corresponding to a a FASTA file path
#' containing the amino acid sequences of multiple rhodopsins.
#'
#' @return AAStringSet
#'
#' @examples
#' # Example 1:
#'
#' # Define the path to the fasta file (this one is included in the package)
#' fastaPath <- system.file("extdata", "template.fasta", package = "RhodopXin")
#' seq <- loadSequence(fastaPath)
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @export
#' @import Biostrings, dplyr, pwalign
#'
createPWAlignments <- function(template, sequences){
  helices_seq <- helixSegments(template = template)

  for (i in seq_along(helices_seq)){
    for (j in seq_along(sequences)){
      print(paste0("Alignment with Helix ", i, "and Query sequence ", j))
      print(pwalign::pairwiseAlignment(helices_seq[[i]], sequences[[j]], type = "local"))
    }
  }
}

#' Find Helices Positions Given RCSB ID
#'
#' Adapted directly from rPDBapi's documentation
#'
#' @param rcsb_id
#'
#' @return dataframe
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
