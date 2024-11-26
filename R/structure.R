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
#' A function that finds all the helices in a pdb structure given by an rcsb id
#' and then return a dataframe for the start and end positions of each helix.
#'
#' @param rcsb_id the RCSB accession code to find the helices from
#'
#' @return dataframe with row names corresponding to the helices found in the
#' structure i.e. 'Helix 1' ... 'Helix n' and two columns: 'start' for the start
#' position of the helix in the sequence and 'end' for the end position
#'
#' @references
#' Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D. (2006).
#' _Bio3D: An R package for the comparative analysis of protein structures_.
#' <http://thegrantlab.org/bio3d/>
#'
#' @export
#'
#' @importFrom bio3d read.pdb

findHelices <- function(rcsb_id){
  # Check if the given rcsb_id is valid
  validateRcsbId(rcsb_id = rcsb_id)

  # Retrieve the pdb structure
  pdb_struct <- getPDBstruct(rcsb_id)

  # Create a dataframe with the start and end positions for each helix
  helices_df <- data.frame(start = pdb_struct$helix$start,
                           end = pdb_struct$helix$end)

  # Add row names for the dataframe
  rownames(helices_df) <-  paste("Helix", 1:nrow(helices_df))

  return(helices_df)
}

#' Validate Given RCSB ID
#'
#' A function that validates the rcsb id input given by the user so that the PDB
#' structure can be fetched in the RCSB database
#'
#' @param rcsb_id the RCSB accession code given by the user
#'
#' @return invisible null if the rcsb id is valid, otherwise throw an appropriate error
#'
#' @references
#' Korkmaz S, Yamasan B (2024). _rPDBapi: A Comprehensive Interface for Accessing_
#' _the Protein Data Bank_. R package version 2.1.1,
#' <https://CRAN.R-project.org/package=rPDBapi>.
#'
#' @importFrom rPDBapi get_info

validateRcsbId <- function(rcsb_id){
  # Check if rcsb_id is a character
  if (!is.character(rcsb_id)) {
    stop("Invalid parameter type: rcsb_id must be of type character")
  }

  # Check that only one input is given
  if (length(rcsb_id) != 1){
    stop("More than 1 input: provide only one string for the rcsb_id")
  }

  # Check that the rcsb_id has length 4
  if (nchar(rcsb_id) != 4){
    stop("Incorrect length: rcsb_id must only have a length of 4")
  }

  # Catch errors from  rPDBapi, they deal with non-existent ids so no need to check it
  tryCatch(
    {
      rPDBapi::get_info(rcsb_id)
      return(invisible(NULL))
    },
    error = function(e) {
        stop(e$message)
    }
  )
}

#' Retrieve PDB file from RCSB
#'
#' A function to retrieve the pdb file given an rcsb id. This function also sets
#' up a cache to handle where the pdb files are stored.
#'
#' @param rcsb_id the RCSB accession code given by the user
#'
#' @return structure
#'
#' @references
#' Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D. (2006).
#' _Bio3D: An R package for the comparative analysis of protein structures_.
#' <http://thegrantlab.org/bio3d/>
#'
#' OpenAI. (2024). ChatGPT (November 2024 version). Retrieved from
#' https://chat.openai.com
#'
#' @importFrom bio3d read.pdb get.pdb
getPDBstruct <- function(rcsb_id){
  # Define the temporary directory and file path
  temp_path <- tempdir()
  pdb_file <- file.path(temp_path, paste0(rcsb_id, ".pdb"))

  if (file.exists(pdb_file)) {
    message("Using existing PDB file from tempdir: ", pdb_file)

    # Load the PDB file from the temporary directory
    pdb <- read.pdb(pdb_file)
  } else {
    message("Downloading PDB file and saving to tempdir: ", pdb_file)

    # Download and save the PDB file in the temporary directory
    bio3d::get.pdb(rcsb_id, path = temp_path)
    pdb <- bio3d::read.pdb(rcsb_id)
  }

  return(pdb)
}
