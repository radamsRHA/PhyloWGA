#' Read.OutputFile_iqtree: function to extract a list of parameter estimates from an iqtree outputfile
#'
#' This function returns a list of the parameter estimates
#' @param string.Path_OutputFile_iqtree String defining path to the output iqtree file
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return List List of containing estimates of the GTR+G model
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
#' Execute.PhylogeneticInference_iqtree(matrix.WindowAlignment = matrix.Pruned_WindowAlignment, string.PathParentDir = '~/Desktop/')
#'
#' ####################
#' # Read output file #
#' ####################
#' Read.OutputFile_iqtree(string.Path_OutputFile_iqtree = '~/Desktop/PhylogeneticInference_iqtree/WindowAlignment.fa.iqtree')
#'

##########################
# Read.OutputFile_iqtree #
##########################
Read.OutputFile_iqtree <- function(string.Path_OutputFile_iqtree){

  ######################################
  # Read in the outputfile from iqtree #
  ######################################
  handle.ReadLines_OutputFile_iqtree <- readLines(string.Path_OutputFile_iqtree)

  #############################
  # Summarize param estimates #
  #############################
  #vector.ParameterEstimates <- vector()
  #vector.NamesParameterEstimates <- vector()
  
  vector.ParameterEstimates <- rep("-", 13)
  names(vector.ParameterEstimates) <- c("r_AC", "r_AG", "r_AT", "r_CG", "r_CT", "r_GT", "pi_A", "pi_C", "pi_G", "pi_T", "alpha", "Invar", "Model")
  
  #####################
  # Relative rate A-C #
  #####################
  numeric.LineWithAC <- grep(pattern = "A-C", x = handle.ReadLines_OutputFile_iqtree)
  if (length(numeric.LineWithAC) >= 1){
    string.LineWithAC <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithAC]
    numeric.AC <- as.numeric(strsplit(x = string.LineWithAC, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "r_AC"] <- numeric.AC
  }

  #####################
  # Relative rate A-G #
  #####################
  numeric.LineWithAG <- grep(pattern = "A-G", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithAG) >= 1){
    string.LineWithAG <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithAG]
    numeric.AG <- as.numeric(strsplit(x = string.LineWithAG, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "r_AG"] <- numeric.AG
  }


  #####################
  # Relative rate A-T #
  #####################
  numeric.LineWithAT <- grep(pattern = "A-T", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithAT) >= 1){
    string.LineWithAT <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithAT]
    numeric.AT <- as.numeric(strsplit(x = string.LineWithAT, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "r_AT"] <- numeric.AT
  }


  #####################
  # Relative rate C-G #
  #####################
  numeric.LineWithCG <- grep(pattern = "C-G", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithCG) >= 1){
    string.LineWithCG <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithCG]
    numeric.CG <- as.numeric(strsplit(x = string.LineWithCG, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "r_CG"] <- numeric.CG
  }

  #####################
  # Relative rate C-T #
  #####################
  numeric.LineWithCT <- grep(pattern = "C-T", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithCT) >= 1){
    string.LineWithCT <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithCT]
    numeric.CT <- as.numeric(strsplit(x = string.LineWithCT, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "r_CT"] <- numeric.CT
  }

  #####################
  # Relative rate G-T #
  #####################
  numeric.LineWithGT <- grep(pattern = "G-T", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithGT) >= 1){
    string.LineWithGT <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithGT]
    numeric.GT <- as.numeric(strsplit(x = string.LineWithGT, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "r_GT"] <- numeric.GT
  }


  ##############################
  # Equilibrium frequency of A #
  ##############################
  numeric.LineWithA <- grep(pattern = "pi(A)", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithA) >= 1){
    string.LineWithA <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithA]
    numeric.A <- as.numeric(strsplit(x = string.LineWithA, split = "=", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "pi_A"] <- numeric.A
  }

  ##############################
  # Equilibrium frequency of C #
  ##############################
  numeric.LineWithC <- grep(pattern = "pi(C)", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithC) >= 1){
    string.LineWithC <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithC]
    numeric.C <- as.numeric(strsplit(x = string.LineWithC, split = "=", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "pi_C"] <- numeric.C
  }


  ##############################
  # Equilibrium frequency of G #
  ##############################
  numeric.LineWithG <- grep(pattern = "pi(G)", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithG) >= 1){
    string.LineWithG <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithG]
    numeric.G <- as.numeric(strsplit(x = string.LineWithG, split = "=", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "pi_G"] <- numeric.G
  }
  
  ##############################
  # Equilibrium frequency of T #
  ##############################
  numeric.LineWithT <- grep(pattern = "pi(T)", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if (length(numeric.LineWithT) >= 1){
    string.LineWithT <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithT]
    numeric.T <- as.numeric(strsplit(x = string.LineWithT, split = "=", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "pi_T"] <- numeric.T
  }
  
  ###############
  # Gamma rates #
  ###############
  numeric.LineWithAlpha <- grep(pattern = "Gamma shape alpha:", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if(length(numeric.LineWithAlpha) >= 1){
    string.LineWithAlpha <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithAlpha]
    numeric.Alpha <- as.numeric(strsplit(x = string.LineWithAlpha, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "alpha"] <- numeric.Alpha
  }
  
  ###################
  # Invariant sites #
  ###################
  numeric.LineWithInvar <- grep(pattern = "Proportion of invariable sites", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if(length(numeric.LineWithInvar) >= 1){
    string.LineWithInvar <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithInvar]
    numeric.Invar <- as.numeric(strsplit(x = string.LineWithInvar, split = ":", fixed = T)[[1]][2])
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "Invar"] <- numeric.Invar
  }

  ##############
  # Best Model #
  ##############
  numeric.LineWithBestModel <- grep(pattern = "Model of substitution", x = handle.ReadLines_OutputFile_iqtree, fixed = T)
  if(length(numeric.LineWithBestModel) >= 1){
    string.LineWithBestModel <- handle.ReadLines_OutputFile_iqtree[numeric.LineWithBestModel]
    string.BestModel <- strsplit(x = string.LineWithBestModel, split = ":", fixed = T)[[1]][2]
    
    vector.ParameterEstimates[names(vector.ParameterEstimates) == "Model"] <- string.BestModel
  }

  ######################
  # Set results vector #
  ######################
  
  return(vector.ParameterEstimates)
#  return(list(numeric.A = numeric.A,
#              numeric.C = numeric.C,
#              numeric.G = numeric.G,
#              numeric.T = numeric.T,
#              numeric.AC = numeric.AC,
#              numeric.AG = numeric.AG,
#              numeric.AT = numeric.AT,
#              numeric.CG = numeric.CG,
#              numeric.CT = numeric.CT,
#              numeric.GT = numeric.GT,
#              numeric.Alpha =numeric.Alpha))
}
