#' Extract.WindowAlignment: function to extract an alignment matrix given an input vector of sequences, coordinates, and a path to a fasta file
#'
#' This function returns matrix consisting of an extracted alignment
#' @param string.PathToFastaFile String defining the path to the input fasta file
#' @param vector.SequenceNames Vector of sequence names to be extracted
#' @param vector.Coordinates_ExperimentalLocus Vector containing two coordinates for the window (start and end)
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return matrix.Window_Alignment Matrix containing the alignment for the input window coordinates from the fasta file
#' @export
#' @examples
#'
#' ################
#' # Load depends #
#' ################
#' library(PhyloWGA)
#' library(ape)
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
#'
#' ######################
#' # Get sequence names #
#' ######################
#' vector.SequenceNames <- Get.SequenceNames(string.PathToFastaFile = String.Path_ExampleChromosomeAlignment)
#'
#' ##################################
#' # Define experimental parameters #
#' ##################################
#' matrix.ExperimentalParams <- Define.PhyloWGA_Experiment(numeric.WindowSize = 1000, numeric.StepSize = 1000, numeric.TotalLength = numeric.AlignmentLength)
#'
#' ############################
#' # Extract window alignment #
#' ############################
#' matrix.WindowAlignment <- Extract.WindowAlignment(string.PathToFastaFile = String.Path_ExampleChromosomeAlignment,
#'                                                   vector.SequenceNames = vector.SequenceNames,
#'                                                   vector.Coordinates_ExperimentalLocus = matrix.ExperimentalParams[5,])
#'
#'
#' write.dna(x = matrix.WindowAlignment, file = '~/Desktop/ExampleWindow.fasta', colsep = "", format = "fasta")
#'

###########################
# Extract.WindowAlignment #
###########################
Extract.WindowAlignment <- function(string.PathToFastaFile, vector.SequenceNames, vector.Coordinates_ExperimentalLocus){


  #############################
  # Summarize input arguments #
  #############################
  options(scipen=9999999)
  numeric.NumberOfSequences <- length(vector.SequenceNames)
  numeric.WindowStartPosition <- vector.Coordinates_ExperimentalLocus[1]
  numeric.WindowEndPosition <- vector.Coordinates_ExperimentalLocus[2]
  numeric.WindowSize <- length(numeric.WindowStartPosition:numeric.WindowEndPosition)

  ######################################
  # Define matrix for window alignment #
  ######################################
  matrix.Window_Alignment <- matrix(nrow = numeric.NumberOfSequences, ncol = numeric.WindowSize)
  rownames(matrix.Window_Alignment) <- vector.SequenceNames

  ##############################
  # Loop through each sequence #
  ##############################
  for (i in 1:numeric.NumberOfSequences){

    #####################
    # Get sequence name #
    #####################
    string.SequenceName_i <- vector.SequenceNames[i]

    ##################################
    # Define samtools command string #
    ##################################
    string.Command_SamTools <- "samtools faidx ZZZ AAA:BBB-CCC"
    string.Command_SamTools <- gsub(pattern = "AAA", replacement = string.SequenceName_i, x = string.Command_SamTools)
    string.Command_SamTools <- gsub(pattern = "BBB", replacement = numeric.WindowStartPosition, x = string.Command_SamTools)
    string.Command_SamTools <- gsub(pattern = "CCC", replacement = numeric.WindowEndPosition, x = string.Command_SamTools)
    string.Command_SamTools <- gsub(pattern = "ZZZ", replacement = string.PathToFastaFile, x = string.Command_SamTools)

    ####################
    # Execute samtools #
    ####################
    handle.Results_ExecuteSamtools <- system(command = string.Command_SamTools, intern = T)
    numeric.Length_Results_ExecuteSamtools <- length(handle.Results_ExecuteSamtools)
    vector.Sequence_AlignmentData <- vector()

    for (j in 2:numeric.Length_Results_ExecuteSamtools){
      vector.SequenceData_Chunk_j <- strsplit(x = handle.Results_ExecuteSamtools[j], split = "")[[1]]
      vector.Sequence_AlignmentData[(length(vector.Sequence_AlignmentData)+1):((length(vector.Sequence_AlignmentData))+length(vector.SequenceData_Chunk_j))] <- vector.SequenceData_Chunk_j
    }

    #########################
    # Process sequence data #
    #########################
    vector.Sequence_AlignmentData <- gsub(pattern = "*", replacement = "N", x = vector.Sequence_AlignmentData, fixed = T)
    matrix.Window_Alignment[string.SequenceName_i,] <- vector.Sequence_AlignmentData



  }
  return(matrix.Window_Alignment)
}
