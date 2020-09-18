#' Get.AlignmentLength: get the total length of an alignment in the input fasta file
#'
#' This function returns the number of bp found in the input fasta file
#' @param string.PathToFastaFile String defining the path to the input fasta file
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return numeric.Total_AlignmentLength Number of base pairs in the input fasta file
#' @export
#' @examples
#'
#' #####################################
#' # Read example chromosome alignment #
#' #####################################
#' String.Path_ExampleChromosomeAlignment <-  system.file("extdata", "Example_Chr10.fasta", package="PhyloWGA")
#'
#' ########################
#' # Get alignment length #
#' ########################
#' numeric.AlignmentLength <- Get.AlignmentLength(string.PathToFastaFile = String.Path_ExampleChromosomeAlignment)


#######################
# Get.AlignmentLength #
#######################
Get.AlignmentLength <- function(string.PathToFastaFile){

  #################################################
  # Set string for getting locus length using awk #
  #################################################
  string.SystemCommand_GetLocusLength <- "awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next} { seqlen = seqlen +length($0)}END{print seqlen} {n++} n>seqlen {exit}' XXX"
  string.SystemCommand_GetLocusLength <- gsub(pattern = "XXX", replacement = string.PathToFastaFile, x = string.SystemCommand_GetLocusLength)
  handle.Print_CommandResults <- system(command = string.SystemCommand_GetLocusLength, intern = T)[[2]]
  numeric.Total_AlignmentLength <- as.numeric(handle.Print_CommandResults)

  string.PrintSummary <- "Found an alignment of XXXbp in fasta file YYY"
  string.PrintSummary <- gsub(pattern = "XXX", replacement = numeric.Total_AlignmentLength, x = string.PrintSummary)
  string.PrintSummary <- gsub(pattern = "YYY", replacement = string.PathToFastaFile, x = string.PrintSummary)
  print(string.PrintSummary)

  return(numeric.Total_AlignmentLength)
}
