#' Execute.PhylogeneticInference_iqtree: function to estimate a phylogenetic tree model for a window with iqtree
#'
#' This function returns a list containing parameter estimates and a phylogenetic tree estimate
#' @param matrix.WindowAlignment Matrix of the alignment for a window
#' @param string.PathParentDir Path used for estimating trees with iqtree
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return List List containing (1) list of parameter estimates and (2) phylogenetic tree estimate
#' @export
#' @examples
#'
#'
#'
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
#' ##############################
#' # Estimate model with iqtree #
#' ##############################
#' Execute.PhylogeneticInference_iqtree(matrix.WindowAlignment = matrix.Pruned_WindowAlignment, string.PathParentDir = '~/Desktop/', string.Commands_iqtree = "-m GTR+G")
#'


########################################
# Execute.PhylogeneticInference_iqtree #
########################################

Execute.PhylogeneticInference_iqtree <- function(matrix.WindowAlignment, string.PathParentDir, string.Commands_iqtree){

  ######################################################################
  # Define parent dir used for phylogenetic tree inference with iqtree #
  ######################################################################
  string.Path_PhylogenyInfernce_ParentDir = paste0(string.PathParentDir, '/PhylogeneticInference_iqtree')
  unlink(string.Path_PhylogenyInfernce_ParentDir, recursive = T)
  dir.create(string.Path_PhylogenyInfernce_ParentDir, showWarnings = T, recursive = T)

  ###################################################
  # Write window alignment to fasta file for iqtree #
  ###################################################
  string.PathFastaFile <- paste0(string.Path_PhylogenyInfernce_ParentDir, '/WindowAlignment.fa')
  write.FASTA(x = as.DNAbin(matrix.WindowAlignment), file = string.PathFastaFile)

  #############################
  # Define command for iqtree #
  #############################
  setwd(string.Path_PhylogenyInfernce_ParentDir)
  string.Command_IqTree <- paste0("iqtree -safe -s WindowAlignment.fa ", string.Commands_iqtree)
  handle.Results_IqTree <- system(command = string.Command_IqTree, intern = T)

  #####################################
  # Read the results file from iqtree #
  #####################################
  setwd(string.PathParentDir)
  string.PathOutputFile <- paste0(string.Path_PhylogenyInfernce_ParentDir, '/WindowAlignment.fa.iqtree')
  list.ParameterEstimates_Window <- Read.OutputFile_iqtree(string.Path_OutputFile_iqtree = string.PathOutputFile)
  handle.PhylogeneticTree_Estimate <- read.tree(file = paste0(string.Path_PhylogenyInfernce_ParentDir, "/WindowAlignment.fa.treefile"))

  return(list(list.ParameterEstimates_Window = list.ParameterEstimates_Window,
              handle.PhylogeneticTree_Estimate = handle.PhylogeneticTree_Estimate))


}
