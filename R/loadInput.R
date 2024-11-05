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
loadSequence <- function(inputStr){
  # Check if inputStr is a character
  if (!is.character(inputStr)) {
    stop("Invalid parameter type: inputStr must be of type character")
  }

  if (!grepl("\\.(fa|fasta|txt)$", inputStr)){
    stop("Incorrect file format: must be one of '.fa', '.fasta', '.txt'")
  }

  if (!file.exists(inputStr)){
    stop("Invalid file: specified path to the fasta file does not exist")
  }

  seq_set <- Biostrings::readAAStringSet(inputStr)

  return(seq_set)
}


