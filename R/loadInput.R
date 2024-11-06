#' Load Amino Acid Sequence Input
#'
#' A function that takes in a FASTA file path containing the sequence information
#' and creates an AAStringSet object out of it
#'
#' @param filePath A string input corresponding to a a FASTA file path
#' containing the amino acid sequences of multiple rhodopsins.
#'
#' @return AAStringSet
#'
#' @examples
#' # Example 1:
#' # Using a fasta file as input
#'
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
#' @import Biostrings
#'
loadSequence <- function(filePath){
  # Check if filePath is a character
  if (!is.character(filePath)) {
    stop("Invalid parameter type: filePath must be of type character")
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

  return(seq_set)
}


