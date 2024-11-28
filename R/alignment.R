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
#' @param rcsb_id the RCSB accession code for the template rhodopsin to be used
#'
#' @return returns a list containing two objects:
#' \describe{
#'  \item{all_pwa}{a list of AAStringSet objects containing the combined pairwise alignments
#'  of each template's helix. Each index in the list corresponds to the helix where
#'  that alignment was made e.g. the AAStringSet at index 1 corresponds to the combined
#'  pairwise alignments of all the query sequences against Helix 1}
#'  \item{template_ranges}{a dataframe containing information on which positions
#'  of the template sequence the queries where aligned to. It has the following
#'  format:
#'    \describe{
#'        \item{rownames}{Name of the helix e.g. Helix 1}
#'        \item{start}{start position of helixin template sequence}
#'        \item{end}{end position of helix in template sequence}
#'        \item{Query #}{Have a column for each query containing the ranges for
#'        which they have alignment to the template. Ranges are in the string format
#'        'start:end'.}
#'    }
#'  }
#' }
#'
#' @examples
#' template <- template_rhodopsins[1]
#' sequences <- sample_rhodopsins
#' rcsb_id <- "1QHJ"
#' results <- createHelixAlignments(template = template, sequences = sequences, rcsb_id = rcsb_id)
#'
#' @references
#' Aboyoun P, Gentleman R (2024). _pwalign: Perform pairwise sequence alignments_.
#' doi:10.18129/B9.bioc.pwalign <https://doi.org/10.18129/B9.bioc.pwalign>, R
#' package version 1.2.0, <https://bioconductor.org/packages/pwalign>.
#'
#' del Val, C., Royuela-Flor, J., Milenkovic, S., & Bondar, A.-N. (2014).
#' _Channelrhodopsins: A bioinformatics perspective_. Biochimica et Biophysica Acta
#' (BBA) - Bioenergetics, 1837(5), 643–655. https://doi.org/10.1016/j.bbabio.2013.11.005
#'
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient _
#' _manipulationof biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' Wright ES (2016). _Using DECIPHER v2.0 to Analyze Big Biological _
#' _Data in R._ The R Journal, *8*(1), 352-359.
#'
#' @export
#'
#' @importFrom DECIPHER AlignSeqs
#' @importFrom pwalign pairwiseAlignment alignedSubject start end subject
createHelixAlignments <- function(template, sequences, rcsb_id){
  # Validate the given rcsb_id. See structure.R for this function
  validateRcsbId(rcsb_id = rcsb_id)

  # Get the AA sequence for each helix present in the template
  helices_seq <- helixSequences(template = template, rcsb_id = rcsb_id)

  # Generate a vector for query names
  query_names <- c(paste("Query", seq_along(1:length(sequences))))

  # Extend the dataframe returned by findHelices to include the ranges where each
  # of the query sequences aligned with the helix of the template
  subject_ranges <- findHelices(rcsb_id = rcsb_id)

  # Create a dataframe to store which positions in the template each query is
  # aligned to. will be used to map the alignments in the structure
  template_ranges <- subject_ranges

  # Add a column for each query sequence that will hold the subject_ranges where they
  # are aligned with the template helix
  for(q in query_names){
    subject_ranges[[q]] <- NA
    template_ranges[[q]] <- NA
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
      query_aligned <- pwalign::subject(alignment)

      # Add the aligned segment of the query to the alignments vector
      alignments <- c(alignments, query_aligned)

      # Get the position where this alignment starts in the query sequence
      query_start <- pwalign::start(pwalign::subject(alignment))

      # Get the position where this alignment ends in the query sequence
      query_end <- pwalign::end(pwalign::subject(alignment))

      # Add the start and end positions of the alignment in the subject_ranges dataframe
      # Add 2 in the column to skip the start and positions of the helix
      subject_ranges[i, j+2] <- paste0(query_start, ":", query_end)

      # Get the position where this alignment starts in the template sequence
      template_start <- pwalign::start(pwalign::pattern(alignment))

      # Get the position where this alignment ends in the template sequence
      template_end <- pwalign::end(pwalign::pattern(alignment))

      # Add the start and end positions of the alignment in the subject_ranges dataframe
      template_ranges[i, j+2] <- paste0(template_start, ":", template_end)
    }

    # Combine the helix sequence with all the aligned query sequences into a
    # singular AAStringSet object
    combined_alignment <- DECIPHER::AlignSeqs(alignments)

    # Add names to the AAStringSet object as this got lost during alignment
    names(combined_alignment) <- c(paste0("Helix ", i, " - ", subject_ranges[i, "start"], ":", subject_ranges[i, "end"]),
                                   paste(query_names, "-", subject_ranges[i, 3:(length(sequences)+2)]))

    # Append this pairwise alignment to the all_pwa list
    all_pwa[[i]] <- combined_alignment
  }

  return(list(all_pwa = all_pwa, template_ranges = template_ranges))
}

#' Get the Sequence for Each Helix of the Template
#'
#' A function that returns a list of AAString objects containing the sequence
#' info of the given template rhodopsin's helices. The order in the list
#' corresponds to the order the helices were found in the original sequence.
#'
#' @param template AAStringSet of template rhodopsin to get helices from
#'
#' @param rcsb_id the RCSB accession code for the template rhodopsin to be used
#'
#' @return a list of AAStringSet objects containing the sequences of the template
#' rhodopsin's helices.
#'
#' @examples
#' template <- template_rhodopsins[1]
#' rcsb_id <- "1QHJ"
#' results <- helixSequences(template = template, rcsb_id = rcsb_id)
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @importClassesFrom Biostrings AAStringSet
#' @importFrom Biostrings subseq
helixSequences <- function(template, rcsb_id){
  # Get the dataframe containing the start and end positions of every helix in
  # the template
  helices_df <- findHelices(rcsb_id)

  # Initialize a list that will contain all the sequence of the helices
  helices_seq <- list()

  # Loop thru each helix int the template
  for (i in 1:nrow(helices_df)) {
    start_pos <- helices_df$start[i] # Get the start position of helix
    end_pos <- helices_df$end[i]     # Get the end position of the helix

    # Given the positions of the sequence, extract it from the template sequence
    # using the subseq function from BioString
    seq <- Biostrings::subseq(template[[1]], start=start_pos, end=end_pos)

    # Put the extracted sequence as an AAStringSet object and store in the list
    helices_seq[[i]] <- Biostrings::AAStringSet(as.character(seq))
  }

  # Return the list of helices sequences
  return(helices_seq)
}

#' Process Template Ranges Dataframe for Mapping
#'
#' The rannges returned by template_ranges of createHelixAlignments need to be processed
#' because they do not correspond to the absolute ranges during alignment i.e. taking
#' into account the whole length of the sequence and not just the helix it was aligned to.
#' Additionally, the resolved structure of the template may be missing some residues
#' so the absolute positions to the full sequence may again nedd to be adjusted
#' to take into account the missing residues in the structure.
#'
#' This function processes the dataframe so it reflects the mapping positions of the
#' alignments to the template's resolved structure.
#'
#' @param template_ranges the dataframe template_ranges returned by createHelixAlignments()
#'
#' @param template AAStringSet of template rhodopsin to get helices from
#'
#' @param rcsb_id the RCSB accession code for the template rhodopsin used
#'
#' @return a processed template_ranges dataframe with the following columns:
#' \describe{
#'  \item{Helix}{name of the helix where this alignment occurs}
#'  \item{Query}{name of query this alignment belongs to}
#'  \item{mapped_start}{the residue start positions mapped to the resolved structure}
#'  \item{mapped_end}{the residue end positions mapped to the resolved structure}
#' }
#'
#' @examples
#' # Create the alignments first
#' template <- template_rhodopsins[1]
#' sequences <- sample_rhodopsins
#' rcsb_id <- "1QHJ"
#' results <- createHelixAlignments(template = template, sequences = sequences, rcsb_id = rcsb_id)
#'
#' # Run the mapping
#' template_ranges <- results$template_ranges
#' results <- templateMapping(template_ranges = template_ranges, template = template, rcsb_id = rcsb_id)
#'
#' @references
#' Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package version 3.2.1,
#' <https://CRAN.R-project.org/package=tibble>.
#'
#' Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data
#' Manipulation_. R package version 1.1.4, <https://CRAN.R-project.org/package=dplyr>.
#'
#' Wickham H, Henry L (2023). _purrr: Functional Programming Tools_. R package version 1.0.2,
#' <https://CRAN.R-project.org/package=purrr>.
#'
#' Wickham H, Vaughan D, Girlich M (2024). _tidyr: Tidy Messy Data_. R package version 1.3.1,
#' <https://CRAN.R-project.org/package=tidyr>.
#'
#' @importFrom purrr map2_lgl
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr %>% mutate select starts_with
templateMapping <- function(template_ranges, template, rcsb_id){
  # Adjust for numbering in alignment and correct residue position
  mapping_df <- template_ranges %>%
    tibble::rownames_to_column(var = "HelixName") %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("Query"),
                        names_to = "Query",
                        values_to = "Range") %>%
    dplyr::mutate(relative_start = as.numeric(sapply(strsplit(Range, ":"), "[", 1)),
                  relative_end = as.numeric(sapply(strsplit(Range, ":"), "[", 2)),
                  absolute_start = start + relative_start - 1,
                  absolute_end = start + relative_end - 1
                  ) %>%
    dplyr::select(HelixName, Query, absolute_start, absolute_end)

  # Adjustment for resolved structure, get the mapping dataframe
  resolved_df <- resolvedMapDf(template = template, rcsb_id = rcsb_id)

  mapping_df <- mapping_df %>%
    mutate(map = purrr::map2_lgl(absolute_start, absolute_end, function(start, end) {
      all(resolved_df$resolved[start:end] != "-")
    }))

  mapping_df$map[is.na(mapping_df$map)] <- FALSE

  mapping_df <- mapping_df %>%
    mutate(mapped_start = ifelse(map, resolved_df$resolved[absolute_start], NA),
           mapped_end = ifelse(map, resolved_df$resolved[absolute_end], NA)) %>%
    filter(map == TRUE) %>%
    select(HelixName, Query, mapped_start, mapped_end)

  return(mapping_df)
}

#' Generate Mapping Dataframe of Resolved Structure
#'
#' This function creates a dataframe containing the the residue numbers in the full
#' sequence of the template and the corresponding residue numbers in the resolved
#' structure when missing residues are taken into account.
#'
#' @param template AAStringSet of template rhodopsin to get helices from
#'
#' @param rcsb_id the RCSB accession code for the template rhodopsin used
#'
#' @return a dataframe with two columns:
#' \describe{
#'    \item{full}{the numbering of residues in the full sequence of the template}
#'    \item{resolved}{the corresponding numbering of positions in the resolved
#'    structure taking into account the missing residues}
#' }
#'
#' @examples
#' template <- template_rhodopsins[1]
#' rcsb_id <- "1QHJ"
#'
#' result <- resolvedMapDf(template = template, rcsb_id = rcsb_id)
#'
#' @references
#' Aboyoun P, Gentleman R (2024). _pwalign: Perform pairwise sequence alignments_.
#' doi:10.18129/B9.bioc.pwalign <https://doi.org/10.18129/B9.bioc.pwalign>, R package version
#' 1.2.0, <https://bioconductor.org/packages/pwalign>.
#'
#' Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D. (2006).
#' _Bio3D: An R package for the comparative analysis of protein structures_.
#' <http://thegrantlab.org/bio3d/>
#'
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation of_
#' _biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @importClassesFrom Biostrings AAStringSet
#' @importFrom bio3d pdbseq
#' @importFrom pwalign pairwiseAlignment
resolvedMapDf <- function(template, rcsb_id){
  # Full sequence as given by FASTA file
  full <- as.character(template[[1]])

  # Get residues in resolved structure from pdb file
  pdb_struct <- getPDBstruct("3UG9")
  resolved <- paste0(bio3d::pdbseq(pdb_struct), collapse = "")
  resolved_set <- Biostrings::AAStringSet(resolved)

  # Align the full sequence and resolved sequence to identify whats missing
  alignment <- pwalign::pairwiseAlignment(template[1], resolved_set)

  full_seq <- strsplit(as.character(alignment@pattern), "")[[1]]
  resolved_seq <- as.character(alignment@subject)

  # Find the positions of non-'-' characters
  non_dash_indices <- which(strsplit(resolved_seq, "")[[1]] != "-")

  # Create a vector with the adjusted indices
  adjusted <- sapply(seq_along(strsplit(resolved_seq, "")[[1]]), function(i) {
    if (strsplit(resolved_seq, "")[[1]][i] != "-") {
      return(which(non_dash_indices == i))
    } else {
      return("-")
    }
  })

  # create the mapping dataframe
  mapping <- data.frame(full = 1:length(full_seq), resolved = adjusted)

  return(mapping)
}

# [END]
