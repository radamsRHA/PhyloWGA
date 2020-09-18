#' Organize.ParallelPhyloWGA: function to partition and organize subdirectories for parallel WGA analyses
#'
#' This function returns a XXX
#' @param numeric.NumberSubsets Number of subsets for the original alignment and analyses
#' @param string.PathToFastaFile String defining the path to the input fasta file
#' @param string.PathParentDir Path to parent dir for all analyses
#' @param numeric.WindowSize Size of each chromosomal window
#' @param numeric.StepSize Spacing of each chromosomal window
#' @param numeric.NumberOfCores Number of cores
#' @keywords phylogenetics, whole genome, species tree, coalescent theory, parallel
#' @return XXX XXX
#' @export
#' @examples
#'
#'


#############################
# Organize.ParallelPhyloWGA #
#############################
Organize.ParallelPhyloWGA <- function(numeric.NumberSubsets, string.PathToFastaFile, string.PathParentDir, string.Analysis, numeric.WindowSize, string.Commands_iqtree, numeric.StepSize, numeric.NumberOfCores){
  
  #############################################
  # Define parent dir used for the experiment #
  #############################################
  string.Path_ParenDir_Organize_PhyloWGA = paste(string.PathParentDir, '/Organize_PhyloWGA_', Sys.Date(), sep = "")
  unlink(string.Path_ParenDir_Organize_PhyloWGA, recursive = T)
  dir.create(string.Path_ParenDir_Organize_PhyloWGA, showWarnings = T, recursive = T)
  
  #########################
  # Summarize input fasta #
  #########################
  numeric.Fasta_TotalLength <- Get.AlignmentLength(string.PathToFastaFile = string.PathToFastaFile)
  vector.Fasta_SequenceNames <- Get.SequenceNames(string.PathToFastaFile = string.PathToFastaFile)
  
  #################
  # Subset params #
  #################
  numeric.SubsetSize <- numeric.Fasta_TotalLength/numeric.NumberSubsets
  numeric.SubsetSize <- floor(numeric.SubsetSize)
  matrix.Experimental_WindowCoordinates <- Define.PhyloWGA_Experiment(numeric.WindowSize = numeric.SubsetSize, numeric.StepSize = numeric.SubsetSize, numeric.TotalLength = numeric.Fasta_TotalLength)
  numeric.NumberOfWindows <- length(matrix.Experimental_WindowCoordinates[,1])
  
  ######################################
  # Loop through subsets in experiment #
  ######################################
  for (i in 1:numeric.NumberSubsets){
    
    ##################
    # Dir for subset #
    ##################
    string.Path_SubsetDir = paste0(string.Path_ParenDir_Organize_PhyloWGA, "/Subset_", i, '/')
    unlink(string.Path_SubsetDir, recursive = T)
    dir.create(string.Path_SubsetDir, showWarnings = T, recursive = T)
    
    ##################################
    # Extract coordinates for window #
    ##################################
    vector.Coordinates_Window_i <- matrix.Experimental_WindowCoordinates[i,]
    
    ################################
    # Extract alignment for window #
    ################################
    matrix.Alignment_Window_i <- Extract.WindowAlignment(string.PathToFastaFile = string.PathToFastaFile,
                                                         vector.SequenceNames = vector.Fasta_SequenceNames,
                                                         vector.Coordinates_ExperimentalLocus = vector.Coordinates_Window_i)
    string.PathFastaFile <- paste0(string.Path_SubsetDir, '/Subset', i, '.fa')
    write.FASTA(x = as.DNAbin(matrix.Alignment_Window_i), file = string.PathFastaFile)
    
    #########################
    # Set file for R script #
    #########################
    string.Script_R <- paste0(string.Path_SubsetDir, '/_Analysis_WGA_', i, '.R')
    
    if (string.Analysis == "Chromo.Phylome"){
      
      string.Commands_R <-  "
library(PhyloWGA)
library(phylotools)
library(ape)


handle.RESULTS_ZZZ <- Chromo.Phylome(string.PathParentDir = 'AAA',
                                           string.PathToFastaFile = 'BBB',
                                           numeric.WindowSize = CCC,
                                           numeric.StepSize = DDD,
                                           string.Commands_iqtree = 'EEE')
      "
      string.Commands_R <- gsub(pattern = "ZZZ", replacement = i, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "AAA", replacement = string.Path_SubsetDir, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "BBB", replacement = paste0(string.Path_SubsetDir,'/Subset', i, '.fa'), x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "CCC", replacement = numeric.WindowSize, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "DDD", replacement = numeric.StepSize, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "EEE", replacement = string.Commands_iqtree, x = string.Commands_R)
      #string.Commands_R <- gsub(pattern = "FFF", replacement = paste0(string.Path_SubsetDir,'/OUTPUT', i, '.txt'), x = string.Commands_R)
      write(x = string.Commands_R, file = string.Script_R)
      
    }
    
    if (string.Analysis == "Chromo.Crawl"){
      
      string.Commands_R <-  "
library(PhyloWGA)
library(phylotools)
library(ape)


Chromo.Crawl(string.PathParentDir = 'AAA',
                         numeric.WindowSize = CCC,
                         numeric.StepSize = DDD,
                         string.PathToFastaFile = 'BBB', 
                         numeric.NumberOfCores = 'OOO')
      "
      string.Commands_R <- gsub(pattern = "AAA", replacement = string.Path_SubsetDir, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "BBB", replacement = paste0(string.Path_SubsetDir,'/Subset', i, '.fa'), x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "CCC", replacement = numeric.WindowSize, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "DDD", replacement = numeric.StepSize, x = string.Commands_R)
      string.Commands_R <- gsub(pattern = "OOO", replacement = numeric.NumberOfCores, x = string.Commands_R)
      write(x = string.Commands_R, file = string.Script_R)
      
      
    }
  }
}