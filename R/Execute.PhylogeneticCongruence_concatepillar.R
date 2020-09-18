#' Execute.PhylogeneticCongruence_concatepillar: function to execute concatepillar on a pair of alignments
#'
#' This function returns a list of results from concatepillar
#' @param matrix.WindowAlignment_01 Matrix of alignment for the first window
#' @param matrix.WindowAlignment_02 Matrix of alignment for the second window
#' @param string.PathParentDir Path to parent dir for all analyses
#' @param numeric.NumberOfCores Number of cores
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return List List of results from concatepillar
#' @export
#' @examples
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
#' matrix.WindowAlignment_01 <- Extract.WindowAlignment(string.PathToFastaFile = String.Path_ExampleChromosomeAlignment,
#'                                                   vector.SequenceNames = vector.SequenceNames,
#'                                                   vector.Coordinates_ExperimentalLocus = matrix.ExperimentalParams[5,])
#'
#'
#' ####################################################
#' # Prune taxa that are only missin data from matrix #
#' ####################################################
#' matrix.Pruned_WindowAlignment_01 <- Prune.MissingDataTaxa(matrix.WindowAlignment = matrix.WindowAlignment_01, numeric.MissingDataThreshold = 0.999)
#'
#' ############################
#' # Extract window alignment #
#' ############################
#' matrix.WindowAlignment_02 <- Extract.WindowAlignment(string.PathToFastaFile = String.Path_ExampleChromosomeAlignment,
#'                                                   vector.SequenceNames = vector.SequenceNames,
#'                                                   vector.Coordinates_ExperimentalLocus = matrix.ExperimentalParams[6,])
#'
#'
#' ####################################################
#' # Prune taxa that are only missin data from matrix #
#' ####################################################
#' matrix.Pruned_WindowAlignment_02 <- Prune.MissingDataTaxa(matrix.WindowAlignment = matrix.WindowAlignment_02, numeric.MissingDataThreshold = 0.999)
#'
#' ##############################
#' # Estimate model with iqtree #
#' ##############################
#' Execute.PhylogeneticCongruence_concatepillar(matrix.Pruned_WindowAlignment_01, matrix.Pruned_WindowAlignment_02, string.PathParentDir, numeric.NumberOfCores)
#'



################################################
# Execute.PhylogeneticCongruence_concatepillar #
################################################
Execute.PhylogeneticCongruence_concatepillar <- function(matrix.WindowAlignment_01, matrix.WindowAlignment_02, string.PathParentDir, numeric.NumberOfCores){

  ###################################################################################
  # Define parent dir used for phylogenetic congruence inference with concatepillar #
  ###################################################################################
  string.Path_PhylogenyCongruence_ParentDir = paste0(string.PathParentDir, '/PhylogeneticCongruence_concatepillar')
  unlink(string.Path_PhylogenyCongruence_ParentDir, recursive = T)
  dir.create(string.Path_PhylogenyCongruence_ParentDir, showWarnings = T, recursive = T)

  ###################################################
  # Write window alignment to fasta file for iqtree #
  ###################################################
  string.PathFastaFile <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/WindowAlignment_01.seq')
  write.FASTA(x = as.DNAbin(matrix.WindowAlignment_01), file = string.PathFastaFile)
  string.PathFastaFile <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/WindowAlignment_02.seq')
  write.FASTA(x = as.DNAbin(matrix.WindowAlignment_02), file = string.PathFastaFile)

  #############################
  # Define command for iqtree #
  #############################
  setwd(string.Path_PhylogenyCongruence_ParentDir)
  string.Command_ConcatePillar <- paste0("concaterpillar.py -t -m GTR -c ", numeric.NumberOfCores)
  handle.Command_ConcatePillar <- system(command = string.Command_ConcatePillar, intern = T)

  #####################################
  # Read the results file from iqtree #
  #####################################
  setwd(string.PathParentDir)
  string.PathOutputFile <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/results.ccp')
  handle.OutputFile <- readLines(string.PathOutputFile)
  string.FinalLine <- handle.OutputFile[length(handle.OutputFile)]

  if (length(grep(pattern = "2 genes", x = string.FinalLine)) == 1){
    numeric.GeneCounter <- 2
  }
  if (length(grep(pattern = "1 genes", x = string.FinalLine)) == 1){
    numeric.GeneCounter <- 1
  }

  #####################
  # Summarize Results #
  #####################
  if (numeric.GeneCounter == 2){

    ###########################
    # Read concatenated fasta #
    ###########################
    string.ConcatenatedFasta <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/set000.seq')
    handle.ConcatenatedPhylip <- read.dna(file = string.ConcatenatedFasta, format = "sequential", as.matrix = T)
    handle.ConcatenatedPhylip <- as.character(as.matrix(handle.ConcatenatedPhylip))
    string.Path_SupergeneWindow_tree <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/set000.tre')
    handle.SupergeneWindow_tree <- read.tree(string.Path_SupergeneWindow_tree)

    return(list(boo.Concate = TRUE,
                matrix.ConcatenatedWindows = handle.ConcatenatedPhylip,
                handle.SupergeneWindow_tree = handle.SupergeneWindow_tree))

  }
  if (numeric.GeneCounter == 1){

    #########################
    # Read individual trees #
    #########################
    string.Path_Tree_Window_01 <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/set000.tre')
    handle.Window_tree_01 <- read.tree(file = string.Path_Tree_Window_01)
    string.Path_Tree_Window_02 <- paste0(string.Path_PhylogenyCongruence_ParentDir, '/set001.tre')
    handle.Window_tree_02 <- read.tree(file = string.Path_Tree_Window_02,)


    return(list(boo.Concate = FALSE,
                matrix.WindowAlignment_01 = matrix.WindowAlignment_01,
                matrix.WindowAlignment_02 = matrix.WindowAlignment_02,
                handle.Window_tree_01 = handle.Window_tree_01,
                handle.Window_tree_02 = handle.Window_tree_02))

  }
}
