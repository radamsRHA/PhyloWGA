#' Prune.MissingDataTaxa: function to pruned matrix (high-missing data sequences removed)
#'
#' This function returns a matrix that has been pruned of sequences with missing data above a given threshold
#' @param matrix.WindowAlignment Matrix of input chromosomal window alignment
#' @param numeric.MissingDataThreshold Numeric (proportion) of missing data permited of any one sequence
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return matrix.Pruned_WindowAlignment Matrix of window alignment with sequences with % missing data above the threshold removed
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
#' ####################################################
#' # Prune taxa that are only missin data from matrix #
#' ####################################################
#' matrix.Pruned_WindowAlignment <- Prune.MissingDataTaxa(matrix.WindowAlignment = matrix.WindowAlignment, numeric.MissingDataThreshold = 0.999)
#'


#########################
# Prune.MissingDataTaxa #
#########################
Prune.MissingDataTaxa <- function(matrix.WindowAlignment, numeric.MissingDataThreshold){

  #############################
  # Summarize input arguments #
  #############################
  numeric.WindowSize <- length(matrix.WindowAlignment[1,])
  vector.SequenceNames <- rownames(matrix.WindowAlignment)
  numeric.NumberOfSequences <- length(vector.SequenceNames)
  vector.PassingSequenceNames <- vector()

  ##############################################
  # Loop through sequences in window alignment #
  ##############################################
  for (i in 1:numeric.NumberOfSequences){

    ################################
    # Extract sequence information #
    ################################
    string.SequenceName_i <- vector.SequenceNames[i]
    vector.Sequence_NucleotideString <- matrix.WindowAlignment[string.SequenceName_i,]

    numeric.Count_N_Symbols <- length(vector.Sequence_NucleotideString[vector.Sequence_NucleotideString=="N"])
    numeric.Count_Gap_Symbols <- length(vector.Sequence_NucleotideString[vector.Sequence_NucleotideString=="-"])

    if (((numeric.Count_N_Symbols+numeric.Count_Gap_Symbols)/numeric.WindowSize) <=  numeric.MissingDataThreshold){
      vector.PassingSequenceNames[(length(vector.PassingSequenceNames)+1)] <- string.SequenceName_i

    }
  }

  #############################
  # Define post-pruned matrix #
  #############################
  matrix.Pruned_WindowAlignment <- matrix(nrow = length(vector.PassingSequenceNames), ncol = numeric.WindowSize)
  rownames(matrix.Pruned_WindowAlignment) <- vector.PassingSequenceNames

  ###################################
  # Loop throughed pruned sequences #
  ###################################
  if (length(vector.PassingSequenceNames) >= 3){
    for (j in 1:length(vector.PassingSequenceNames)){
      matrix.Pruned_WindowAlignment[vector.PassingSequenceNames[j],] <- matrix.WindowAlignment[vector.PassingSequenceNames[j],]
    }
    return(matrix.Pruned_WindowAlignment)
  }
}
