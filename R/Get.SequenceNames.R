#' Get.SequenceNames: get a vector of the sequence names provided in a fasta file
#'
#' This function returns a vector of sequence names
#' @param string.PathToFastaFile String defining the path to the input fasta file
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return vector.SequenceNames Vector containing the set of sequence names obtained from fasta file
#' @export
#' @examples
#'
#'
#'
#' #####################################
#' # Read example chromosome alignment #
#' #####################################
#' String.Path_ExampleChromosomeAlignment <-  system.file("extdata", "Example_Chr10.fasta", package="PhyloWGA")
#'
#' ########################
#' # Get alignment length #
#' ########################
#' vector.SequenceNames <- Get.SequenceNames(string.PathToFastaFile = String.Path_ExampleChromosomeAlignment)



##############################
# Function_Get_SequenceNames #
##############################
Get.SequenceNames <- function(string.PathToFastaFile){

  #######################################################
  # Set string for getting sequence names back with awk #
  #######################################################
  string.SystemCommand_GetSequenceNames <- 'grep ">" XXX'
  string.SystemCommand_GetSequenceNames <- gsub(pattern = "XXX", replacement = string.PathToFastaFile, string.SystemCommand_GetSequenceNames)
  handle.Results_SequenceNames <- system(command = string.SystemCommand_GetSequenceNames, intern = T)
  numeric.NumberOfSequences <- length(handle.Results_SequenceNames)
  handle.Splited_Results_SequenceNames <- strsplit(x = handle.Results_SequenceNames, split = ">", fixed = T)
  vector.SequenceNames <- rep(NA,numeric.NumberOfSequences)

  for (i in 1:numeric.NumberOfSequences){
    vector.SequenceNames[i] <- handle.Splited_Results_SequenceNames[[i]][2]
  }


  string.Print_Summary <- "Found XXX sequences in fasta file YYY"
  string.Print_Summary <- gsub(pattern = "XXX", replacement = numeric.NumberOfSequences, x = string.Print_Summary)
  string.Print_Summary <- gsub(pattern = "YYY", replacement = string.PathToFastaFile, x = string.Print_Summary)
  print(string.Print_Summary)

  return(vector.SequenceNames)


}
