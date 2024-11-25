#' Sequence Alignment of Template Rhodopsin to the Query Sequences
#'
#' A function that performs pairwise alignment between each of the template
#' rhodopsin's helices and each of the query sequences, allowing users to see
#' the conservation of residues across a given helix.
#'
#' This function mimics the template-based alignment performed in the del Val
#' et al. (2014) study using T-Coffee which is an alignment algorithm that allows
#' PDB structures to be used as a 'template' when aligning. Since this is not
#' available in R and there are no similar R packages that can implement it, this
#' function performs a similar approach by identifying the helices in a given
#' structure and performing pairwise alignments to each one with the sequence queries.
#'
#' @param template AAStringSet of template rhodopsin to get helices from
#'
#' @param sequences AAStringSet of query sequences to align to the template
#' rhodopsin's helices
#'
#' @return a list of AAStringSet objects containing the combined pairwise alignments
#' of each template's helix. Each index in the list corresponds to the helix where
#' that alignment was made e.g. the AAStringSet at index 1 corresponds to the combined
#' pairwise alignments of all the query sequences against Helix 1
#'
#' @references
#' Aboyoun P, Gentleman R (2024). _pwalign: Perform pairwise sequence alignments_.
#' doi:10.18129/B9.bioc.pwalign <https://doi.org/10.18129/B9.bioc.pwalign>, R
#' package version 1.2.0, <https://bioconductor.org/packages/pwalign>.
#'
#' del Val, C., Royuela-Flor, J., Milenkovic, S., & Bondar, A.-N. (2014).
#' Channelrhodopsins: A bioinformatics perspective. Biochimica et Biophysica Acta
#' (BBA) - Bioenergetics, 1837(5), 643–655. https://doi.org/10.1016/j.bbabio.2013.11.005
#'
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @export
#'
#' @importFrom DECIPHER AlignSeqs
#' @importFrom pwalign pairwiseAlignment alignedSubject start end subject


createHelixAlignments <- function(template, sequences){
  # Get the AA sequence for each helix present in the template
  helices_seq <- helixSegments(template = template)

  # Generate a vector for query names
  query_names <- c(paste("Query", seq_along(1:length(sequences))))

  # Extend the dataframe returned by findHelices to include the ranges where each
  # of the query sequences aligned with the helix of the template
  ranges <- findHelices(substr(names(template)[[1]], 1, 4))

  # Add a column for each query sequence that will hold the ranges where they
  # are aligned with the template helix
  for(q in query_names){
    ranges[[q]] <- NA
  }

  # Create a list that will hold of the pairwise alignments which will be
  # returned by the function
  all_pwa <- list()

  # Iterate thru each helix and each query sequence to perform pairwise alignment
  for (i in seq_along(helices_seq)){
    # store the alignments for this helix in a vector
    alignments <- helices_seq[[i]]

    for (j in seq_along(sequences)){
      # perform pairwise alignment with the library pwalign
      alignment <- pwalign::pairwiseAlignment(helices_seq[[i]], sequences[[j]], type = "local")

      # Extract the part of the query sequence where there is alignment
      query_aligned <- pwalign::alignedSubject(alignment)

      # Add the aligned segment of the query to the alignments vector
      alignments <- c(alignments, query_aligned)

      # Get the position where this alignment starts in the query sequence
      start <- pwalign::start(pwalign::subject(alignment))

      # Get the position where this alignment ends in the query sequence
      end <- pwalign::end(pwalign::subject(alignment))

      # Add the start and end positions of the alignment in the ranges dataframe
      # Add 2 in the column to skip the start and positions of the helix
      ranges[i, j+2] <- paste0(start, ":", end)
    }

    # Combine the helix sequence with all the aligned query sequences into a
    # singular AAStringSet object
    combined_alignment <- DECIPHER::AlignSeqs(alignments)

    # Add names to the AAStringSet object as this got lost during alignment
    names(combined_alignment) <- c(paste0("Helix", i, "-", ranges[i, "start"], ":", ranges[i, "end"]),
                                   paste(query_names, "-", ranges[i, 3:(length(sequences)+2)]))

    # Append this pairwise alignment to the all_pwa list
    all_pwa[[i]] <- combined_alignment
  }

  return(all_pwa)
}

#' Get the Sequence for Each Helix of the Template
#'
#' A function that returns a list of AAString objects containing the sequence
#' info of the given template rhodopsin's helices. The order in the list
#' corresponds to the order the helices were found in the original sequence.
#'
#' @param template AAStringSet of template rhodopsin to get helices from
#'
#' @return a list of AAStringSet objects containing the sequences of the template
#' rhodopsin's helices.
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @importClassesFrom Biostrings AAStringSet
#' @importFrom Biostrings subseq

helixSegments <- function(template){
  template_id <- substr(names(template)[[1]], 1, 4)
  helices_df <- findHelices(template_id)

  segments <- list()
  for (i in 1:nrow(helices_df)) {
    start_pos <- helices_df$start[i]
    end_pos <- helices_df$end[i]

    segment <- Biostrings::subseq(template[[1]], start=start_pos, end=end_pos)

    segments[[i]] <- Biostrings::AAStringSet(as.character(segment))
  }

  return(segments)
}
