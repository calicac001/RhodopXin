#' Load Amino Acid Sequence Input from File
#'
#' A function that takes in a FASTA file path containing the sequence information
#' and creates an AAStringSet object out of it.
#'
#' @param filePath A string input corresponding to a a FASTA file path
#' containing the amino acid sequences of multiple rhodopsins.
#'
#' @return AAStringSet of the sequences in the FASTA file
#'
#' @examples
#' # Define the path to the fasta file (this one is included in the package)
#' fastaPath <- system.file("extdata", "rhodopsins.fasta", package = "RhodopXin")
#' seq <- loadSequence(fastaPath)
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @export
#' @importFrom Biostrings readAAStringSet
loadSequence <- function(filePath){
  # Check if filePath is a character
  if (!is.character(filePath)) {
    stop("Invalid parameter type: filePath must be of type character")
  }

  # Check that only one input is given
  if (length(filePath) != 1){
    stop("More than 1 input: provide only one string for the file path")
  }

  # Check if the file format is valid
  if (!grepl("\\.(fa|fasta|txt)$", filePath)){
    stop("Incorrect file format: must be one of '.fa', '.fasta', '.txt'")
  }

  # Check if the file exists
  if (!file.exists(filePath)){
    stop("Invalid file: specified path to the fasta file does not exist")
  }

  # Create AAStringSet object given the input file
  seq_set <- Biostrings::readAAStringSet(filePath)

  if (length(seq_set) != 0){
    return(seq_set)
  } else {
    stop("Empty file: no sequences detected")
  }
}

#' Load Sequence Given an RCSB PDB ID
#'
#' A function that takes in a RCSB pdb id and returns an AAStringSet object of
#' the structure's sequence
#'
#' @param rcsb_id the RCSB accession code for the template rhodopsin to get the
#' sequence of
#'
#' @return AAStringSet of the sequence with the given rcsb_id
#'
#' @examples
#' seq <- loadFromRCSB("3UG9")
#'
#' @references
#' Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D. (2006).
#' _Bio3D: An R package for the comparative analysis of protein structures_.
#' <http://thegrantlab.org/bio3d/>
#'
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
#' @export
#'
#' @importFrom Biostrings AAStringSet
#' @importFrom bio3d aa321
loadFromRCSB <- function(rcsb_id){
  # Get the structure given the rcsb_id. This function also validates input so
  # no need to do here. See structure.R for this function.
  pdb_struct <- getPDBstruct(rcsb_id = rcsb_id)

  # Extract the sequence
  print(pdb_struct)

  # Convert sequence from 3-letter to 1-letter code then combine to one string
  # since the aa321 function returns a vector of characters
  sequence <- paste0(bio3d::aa321(pdb_struct$seq), collapse = "")

  # Create AAStringSet object given the input file
  seq_set <- Biostrings::AAStringSet(sequence)
  names(seq_set) <- rcsb_id

  return(seq_set)
}

# [END]
