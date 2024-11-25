#' Selected Examples of Rhodopsin Variants
#'
#' @description This is the sample AAStringSet output produced by running
#' loadSequence() with the rhodopsins.fasta file located in inst/extdata. It
#' contains three rhodopsins selected from the del Val et al. study (2014) which
#' originally looked at 11 channelrhodopsins in a bioinformatic perspective. They
#' included the UniProt entries for the sequences which was then used to fetch the
#' sequencing data to construct the rhodopsins.fasta file.
#'
#' It contains the amino acid sequence of the rhodopsins with the following UniProt
#' entries:
#' \enumerate{
#'    \item{B4Y103 - Volvox carter Channelrhodopsin-1}
#'    \item{F8UVI5 - Mesostigma viride Channelrhodopsin-1}
#'    \item{Q93WP2 - Chlamydomonas reinhardtii Archeal-type opsin 1}
#' }
#'
#' @format An AAStringSet object of length 3 with the sequences of the rhodopsins.
#' \describe{
#'    \item{width}{Length of the rhodopsin's sequence}
#'    \item{seq}{Amino acid sequence of the rhodopsin}
#'    \item{names}{the FASTA header containing the UniProt accession as well as
#'    the entry name of the rhodopsin}
#' }
#'
#' @source del Val, C., Royuela-Flor, J., Milenkovic, S., & Bondar, A.-N. (2014).
#' Channelrhodopsins: A bioinformatics perspective. Biochimica et Biophysica Acta
#' (BBA) - Bioenergetics, 1837(5), 643–655. https://doi.org/10.1016/j.bbabio.2013.11.005
#'
#' UniProt Consortium. (2023). Acop1—Archaeal-type opsin 1—Chlamydomonas reinhardtii
#' (Chlamydomonas smithii) | UniProtKB. Retrieved November 24, 2024, from
#' https://www.uniprot.org/uniprotkb/Q93WP2/entry
#'
#' UniProt Consortium. (2023). Channelopsin 1—Mesostigma viride (Green alga) | UniProtKB.
#' Retrieved November 24, 2024, from https://www.uniprot.org/uniprotkb/F8UVI5/entry
#'
#' UniProt Consortium. (2023). Channelrhodopsin-1—Volvox carteri f. Nagariensis |
#' UniProtKB. Retrieved November 24, 2024, from https://www.uniprot.org/uniprotkb
#' /B4Y103/entry
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
"sample_rhodopsins"


#' Selected Microbial Rhodopsins to be Used as Templates
#'
#' @description This is the sample AAStringSet output produced by running
#' loadSequence() with the templates.fasta file located in inst/extdata. It
#' contains two template rhodopsins selected from the del Val et al. study (2014):
#' a bacteriorhodopsin and a chimera for channelrhodopsin 1 and 2.They
#' included the RCSB PDB entries for these sequences which was then used to fetch
#' their sequencing data to construct the templates.fasta file. Additionally, three
#' more types of microbial rhodopsins (halorhodopsin, proteorhodopsin, and xanthorhodopsin)
#' were included to cover a variety of differnet types.
#'
#' It contains the amino acid sequence of the rhodopsins with the following RCSB
#' PDB entries in the same order as in the AAStringSet:
#' \enumerate{
#'    \item{1QHJ - Bacteriorhodopsin, Halobacterium salinarum}
#'    \item{3UG9 - Channelrhodopsin-1/2 Chimera, Chlamydomonas reinhardtii}
#'    \item{3A7K - Halorhodopsin, Natronomonas pharaonis}
#'    \item{4JQ6 - Proteorhodopsin, MED12BPR Uncultured Medium}
#'    \item{3DDL - Xanthorhodopsin, Salinibacter ruber}
#' }
#'
#' @format An AAStringSet object of length 3 with the sequences of the rhodopsins.
#' \describe{
#'    \item{width}{Length of the rhodopsin's sequence}
#'    \item{seq}{Amino acid sequence of the rhodopsin}
#'    \item{names}{the FASTA header containing the UniProt accession as well as
#'    the entry name of the rhodopsin}
#' }
#'
#' @source
#' Belrhali, H., Nollert, P., Royant, A., Menzel, C., Rosenbusch, J. P., Landau,
#' E. M., & Pebay-Peyroula, E. (1999). Protein, lipid and water organization in
#' bacteriorhodopsin crystals: A molecular view of the purple membrane at 1.9 Å
#' resolution. Structure, 7(8), 909–917. https://www.rcsb.org/structure/1QHJ
#'
#' del Val, C., Royuela-Flor, J., Milenkovic, S., & Bondar, A.-N. (2014).
#' Channelrhodopsins: A bioinformatics perspective. Biochimica et Biophysica Acta
#' (BBA) - Bioenergetics, 1837(5), 643–655. https://doi.org/10.1016/j.bbabio.2013.11.005
#'
#' Kato, H. E., Zhang, F., Yizhar, O., Ramakrishnan, C., Nishizawa, T., Hirata,
#' K., Ito, J., Aita, Y., Tsukazaki, T., Hayashi, S., Hegemann, P., Maturana, A.
#' D., Ishitani, R., Deisseroth, K., & Nureki, O. (2012). Crystal structure of
#' the channelrhodopsin light-gated cation channel. Nature, 482(7385), 369–374.
#' https://www.rcsb.org/structure/3UG9
#'
#' Kouyama, T., Kanada, S., Takeguchi, Y., Narusawa, A., Murakami, M., & Ihara,
#' K. (2010). Crystal Structure of the Light-Driven Chloride Pump Halorhodopsin
#' from Natronomonas pharaonis. Journal of Molecular Biology, 396(3), 564–579.
#' https://www.rcsb.org/structure/3A7K
#'
#' Luecke, H., Schobert, B., Stagno, J., Imasheva, E. S., Wang, J. M., Balashov,
#' S. P., & Lanyi, J. K. (2008). Crystallographic structure of xanthorhodopsin,
#' the light-driven proton pump with a dual chromophore. Proceedings of the National
#' Academy of Sciences, 105(43), 16561–16565. https://www.rcsb.org/structure/3DDL
#'
#' Ran, T., Ozorowski, G., Gao, Y., Sineshchekov, O. A., Wang, W., Spudich, J.
#' L., & Luecke, H. (2013). Cross-protomer interaction with the photoactive site
#' in oligomeric proteorhodopsin complexes. Acta Crystallographica Section D:
#' Biological Crystallography, 69(10), 1965–1980.https://www.rcsb.org/structure/4JQ6
#'
#'
#' @references
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). _Biostrings: Efficient manipulation
#' of biological strings_. doi:10.18129/B9.bioc.Biostrings
#' <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
#' <https://bioconductor.org/packages/Biostrings>.
#'
"template_rhodopsins"

# [END]
