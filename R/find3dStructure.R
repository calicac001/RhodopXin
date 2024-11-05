#' Find 3D Structure Given Amino Acid Sequence
#'
#' @param aaSeq <AAString> The amino acid sequence that we want to find the 3D structure of
#'
#' @examples
#' fastaPath <- system.file("extdata", "rcsb_pdb_8GI8.fasta", package = "RhodopXin")
#' seq <- loadSequence(fastaPath)
#' struct_df <- find3dStructure(seq)
#'
#' @export
#'
#' @import rPDBapi

find3dStructure <- function(aaSeq){
  aaSeqChar <- as.character(aaSeq)
  searchOp <- rPDBapi::SequenceOperator(aaSeqChar, sequence_type = "PROTEIN",
                                        identity_cutoff = 1)
  ids <- rPDBapi::perform_search(search_operator = searchOp)
  structDataframe <- data.frame(rcsb_id = ids)

  for (id in ids){
    info <- rPDBapi::get_info(id)
    structDataframe$title[structDataframe$rcsb_id == id] <- info$struct$title
    structDataframe$keywords[structDataframe$rcsb_id == id] <- info$struct_keywords$text
    structDataframe$paper[structDataframe$rcsb_id == id] <- info$rcsb_primary_citation$title
  }

  return(structDataframe)
}

#' Display 3D structure Given RCSB ID
#'
#' Adapted directly from rPDBapi's documentation
#'
#' @param rcsb_id
#'
#' @return Viewer
#'
#' @export
#'
#' @import rPDBapi, r3dmol

loadStructure <- function(rcsb_id){
  # Retrieve and parse a PDB structure
  pdb_path <- get_pdb_file(rcsb_id, filetype = "pdb", save = TRUE)

  # Visualize the tertiary structure using r3dmol
  viewer <- r3dmol() %>%
    m_add_model(pdb_path$path, format = "pdb") %>% # Load the PDB file
    m_set_style(style = m_style_cartoon()) %>% # Cartoon representation
    m_zoom_to()

  return(viewer)
}
