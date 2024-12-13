% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadInput.R
\name{loadSequence}
\alias{loadSequence}
\title{Load Amino Acid Sequence Input from File}
\usage{
loadSequence(filePath)
}
\arguments{
\item{filePath}{A string input corresponding to a a FASTA file path
containing the amino acid sequences of multiple rhodopsins.}
}
\value{
AAStringSet of the sequences in the FASTA file
}
\description{
A function that takes in a FASTA file path containing the sequence information
and creates an AAStringSet object out of it.
}
\examples{
# Define the path to the fasta file (this one is included in the package)
fastaPath <- system.file("extdata", "rhodopsins.fasta", package = "RhodopXin")
seq <- loadSequence(fastaPath)

}
\references{
Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). \emph{Biostrings: Efficient manipulation
of biological strings}. doi:10.18129/B9.bioc.Biostrings
\url{https://doi.org/10.18129/B9.bioc.Biostrings}, R package version 2.73.1,
\url{https://bioconductor.org/packages/Biostrings}.
}
