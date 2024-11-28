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
#' @examples
#' result <- findHelices(rcsb_id = "3UG9")
#'
#' @references
#' Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D. (2006).
#' _Bio3D: An R package for the comparative analysis of protein structures_.
#' <http://thegrantlab.org/bio3d/>
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
#' @examples
#' result <- validateRcsbId("3UG9")
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
#' A function to retrieve the pdb file given an rcsb id. Stores the pdb file in
#' a temporary directory. The package deals with all this so the user does not have
#' to manually manage and download the pdb files themselves.
#'
#' @param rcsb_id the RCSB accession code given by the user
#'
#' @return pdb class (see ?rPDBapi::read.pdb for specific info on components)
#'
#' @examples
#' result <- getPDBstruct("3UG9")
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

  # Kept getting a warning on file being redownloaded so asked ChatGPT for help
  # on how to fix warning and work with temporary files
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
