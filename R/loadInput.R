#' Load Amino Acid Sequence Input
#'
#' A function that takes in a string of an amino acid corresponding to a rhodopsin
#' or a FASTA file path containing the sequence information and creates an AAString
#' object out of it
#'
#' @param inputStr <Character> A string input corresponding to the amino acid
#' sequence or a path to a FASTA file containing the amino acid sequence.
#' @param inputType <Character> The type of input provided. Accepts either "fasta"
#' for a FASTA file format ('.fa', '.fasta', '.txt') or "string" for an amino
#' acid sequence.
#'
#' @return AAString
#'
#' @examples
#' # Example 1:
#' # Using an amino acid string as input
#' sampleAA <- "MPFYDSRPPEGWPKGSINDMDYPLLGSICAVCCVFVAGSGIWMLYRLDLGMGYSCKPYKSGRA"
#' seq <- loadSequence(sampleAA)
#'
#' # Example 2:
#' # Using a fasta file as input
#'
#' # Define the path to the fasta file (this one is included in the package)
#' fastaPath <- system.file("extdata", "rcsb_pdb_8GI8.fasta", package = "RhodopXin")
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
loadSequence <- function(inputStr, inputType = "fasta"){
  # Check if inputStr is a character
  if (!is.character(inputStr)) {
    stop("Invalid parameter type: inputStr must be of type character")
  }

  # Check if inputType is one of two given options
  if (!(inputType %in% c("fasta", "string"))) {
    stop("Invalid parameter: inputType must be one of 'fasta' or 'string'")
  } else if (inputType == "fasta"){
    validateFile(inputStr)
    seq_set <- Biostrings::readAAStringSet(inputStr)
    seq <- seq_set[[1]]
  } else if (inputType == "string") {
    capitalized <- toupper(inputStr)
    validateString(capitalized)
    seq <- Biostrings::AAString(capitalized)
  }

  return(seq)
}

#' Validate Input File
#'
#' Check whether the input file given above is of the correct format and that
#' such file exists
#'
#'  @param inputStr <Character> A string to the path of an input file
#'
#'  @export

validateFile <- function(inputStr){
  if (!grepl("\\.(fa|fasta|txt)$", inputStr)){
    stop("Incorrect file format: must be one of '.fa', '.fasta', '.txt'")
  }

  if (!file.exists(inputStr)){
    stop("Invalid file: specified path to the fasta file does not exist")
  }
}

#' Validate Input String
#'
#' Check whether the input sring given above only contains valid amino acids
#'
#'  @param inputStr <Character> An amino acid string
#'
#'  @export

validateString <- function(inputStr){
  if (!all(strsplit(inputStr, NULL)[[1]] %in% Biostrings::AA_STANDARD)) {
    stop("Invalid amino acid string: must contain only valid amino acid letters.")
  }
}
