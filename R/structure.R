#' Find 3D Structure Given Amino Acid Sequence
#'
#' @param aaSeq An AAString object of the rhodopsin sequence that we want to find the 3D structure of
#'
#' @return dataframe with information on structures
#'
#' @examples
#' #fastaPath <- system.file("extdata", "template.fasta", package = "RhodopXin")
#' #seq <- loadSequence(fastaPath)
#' #struct_df <- find3dStructure(seq)
#'
#' @references
#' Korkmaz S, Yamasan B (2024). _rPDBapi: A Comprehensive Interface for Accessing the Protein
#' Data Bank_. R package version 2.1.1, <https://CRAN.R-project.org/package=rPDBapi>.
#'
#' @export
#'
#' @importFrom rPDBapi SequenceOperator perform_search get_info

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
#' @param rcsb_id desc
#'
#' @return Viewer
#'
#' @references
#' Korkmaz S, Yamasan B (2024). _rPDBapi: A Comprehensive Interface for Accessing the Protein Data Bank_. R package version 2.1.1, <https://CRAN.R-project.org/package=rPDBapi>.
#'
#' Su W, Johnston B (2021). _r3dmol: Create Interactive 3D Visualizations of Molecular Data_. Rpackage version 0.1.2, <https://CRAN.R-project.org/package=r3dmol>.
#'
#' Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023). _dplyr: A Grammar of Data Manipulation_. <https://dplyr.tidyverse.org>.
#'
#' @export
#'
#' @importFrom rPDBapi get_pdb_file
#' @importFrom r3dmol r3dmol m_add_model m_set_style m_zoom_to
#' @importFrom dplyr %>%

loadStructure <- function(rcsb_id){
  # Retrieve and parse a PDB structure
  pdb_path <- rPDBapi::get_pdb_file(rcsb_id, filetype = "pdb", save = TRUE)

  # Visualize the tertiary structure using r3dmol
  viewer <- r3dmol::r3dmol() %>%
    r3dmol::m_add_model(pdb_path$path, format = "pdb") %>% # Load the PDB file
    r3dmol::m_set_style(style = r3dmol::m_style_cartoon()) %>% # Cartoon representation
    r3dmol::m_zoom_to()

  return(viewer)
}


#' Find Helices Positions Given RCSB ID
#'
#' @param rcsb_id desc
#'
#' @return dataframe
#'
#' @references
#' Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D. (2006). _Bio3D: An R package for the comparative analysis of protein structures_. <http://thegrantlab.org/bio3d/>
#'
#' @export
#'
#' @importFrom bio3d read.pdb

findHelices <- function(rcsb_id){
  pdb_struct <- bio3d::read.pdb(rcsb_id)
  helices_df <- data.frame(start = pdb_struct$helix$start,
                           end = pdb_struct$helix$end)

  rownames(helices_df) <-  paste("Helix", 1:nrow(helices_df))

  return(helices_df)
}
