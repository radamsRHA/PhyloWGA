#' Chromo.Phylome: function to infer a chromosome-specific set of phylogenetic trees
#'
#' This function returns a list of the phylogenetic tree models, one for each window set by the user-defined window and step size input parameters
#' @param string.PathParentDir Path to parent dir for all analyses
#' @param string.PathToFastaFile Path to fasta file
#' @param numeric.WindowSize Size of each chromosomal window
#' @param numeric.StepSize Spacing of each chromosomal window
#' @param string.Commands_iqtree String used for nucleotide substitution model for IqTree. Leave as "" for model selection
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return Results Results writted to output files
#' @export
#' @examples
#'
#'
#'
#'
#' #####################################
#' # Read example chromosome alignment #
#' #####################################
#' String.Path_ExampleChromosomeAlignment <-  system.file("extdata", "Example_Chr10.fasta", package="PhyloWGA")
#'
#'
#' #########################
#' # Conduct.ChromoPhylome #
#' #########################
#' handle.RESULTS <- Chromo.Phylome(string.PathParentDir = '~/Desktop/',
#'                                           string.PathToFastaFile = String.Path_ExampleChromosomeAlignment,
#'                                           numeric.WindowSize = 1000,
#'                                           numeric.StepSize = 1000,
#'                                           string.Commands_iqtree = "-m GTR+G")


##################
# Chromo.Phylome #
##################
Chromo.Phylome <- function(string.PathParentDir, string.PathToFastaFile, numeric.WindowSize, numeric.StepSize, string.Commands_iqtree){

  #############################################
  # Define parent dir used for the experiment #
  #############################################
  string.Path_ParenDir_PPD_Experiment = paste(string.PathParentDir, '/Construct.ChromoPhylome_', Sys.Date(), sep = "")
  unlink(string.Path_ParenDir_PPD_Experiment, recursive = T)
  dir.create(string.Path_ParenDir_PPD_Experiment, showWarnings = T, recursive = T)

  ####################################
  # File to write checkpoint results #
  ####################################
  string.Path_File_CheckPoints_Results <- paste0(string.Path_ParenDir_PPD_Experiment, '/Results_CheckPoint_w', numeric.WindowSize, "_s", numeric.StepSize, ".txt")
  write(x = "Window_i\tStart_i\tEnd_i\tPercentN_i\tPercentGap_i\tNumberSequences_i\tPercent_GC\tr_AC\tr_AG\tr_AT\tr_CG\tr_CT\tr_GT\tpi_A\tpi_C\tpi_G\tpi_T\talpha\tInvar\tModel\tMeanBranchLengths\tSumBranchLengths\tSDBranchLengths\tTree", file = string.Path_File_CheckPoints_Results)

  #########################
  # Summarize input fasta #
  #########################
  numeric.Fasta_TotalLength <- Get.AlignmentLength(string.PathToFastaFile = string.PathToFastaFile)
  vector.Fasta_SequenceNames <- Get.SequenceNames(string.PathToFastaFile = string.PathToFastaFile)

  ##################################
  # Define experimental parameters #
  ##################################
  matrix.Experimental_WindowCoordinates <- Define.PhyloWGA_Experiment(numeric.WindowSize = numeric.WindowSize, numeric.StepSize = numeric.StepSize, numeric.TotalLength = numeric.Fasta_TotalLength)
  numeric.NumberOfWindows <- length(matrix.Experimental_WindowCoordinates[,1])

  ################################################
  # List to collect estimated models from iqtree #
  ################################################
  list.EstimatedModels_iqtree <- list()

  ######################################
  # Loop through windows in experiment #
  ######################################
  for (i in 1:numeric.NumberOfWindows){

    if (i == 1){
      numeric.StartTime <- Sys.time()
    }
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

    ######################################
    # Prune out missing data from window #
    ######################################
    matrix.Pruned_Alignment_Window_i <- Prune.MissingDataTaxa(matrix.WindowAlignment = matrix.Alignment_Window_i, numeric.MissingDataThreshold = 0.95)

    if (length(matrix.Pruned_Alignment_Window_i) >= 1){

      ################################################
      # Build subdirectory for window-based analyses #
      ################################################
      string.Path_Window_i <- paste0(string.Path_ParenDir_PPD_Experiment, '/Analyses_Window_', i)
      unlink(string.Path_Window_i, recursive = T)
      dir.create(string.Path_Window_i, showWarnings = T, recursive = T)

      ###########################################
      # Estimate phylogenetic trees with iqtree #
      ###########################################
      list.Estimated_PhylogeneticTreeModel <- Execute.PhylogeneticInference_iqtree(matrix.WindowAlignment = matrix.Pruned_Alignment_Window_i,
                                                                                   string.PathParentDir = string.Path_Window_i,
                                                                                   string.Commands_iqtree = string.Commands_iqtree)

      #####################
      # Summarize results #
      #####################
      numeric.NumberOfSequences <- length(matrix.Pruned_Alignment_Window_i[,1])
      numeric.PercentGap <- length(matrix.Pruned_Alignment_Window_i[matrix.Pruned_Alignment_Window_i=="-"])/(length(matrix.Pruned_Alignment_Window_i))
      numeric.PercentN <- length(matrix.Pruned_Alignment_Window_i[matrix.Pruned_Alignment_Window_i=="N"])/(length(matrix.Pruned_Alignment_Window_i))
      numeric.Total <- length(matrix.Pruned_Alignment_Window_i) - (numeric.PercentGap*length(matrix.Pruned_Alignment_Window_i)) - (numeric.PercentN*length(matrix.Pruned_Alignment_Window_i))
      numeric.PercentGC <- (length(matrix.Pruned_Alignment_Window_i[matrix.Pruned_Alignment_Window_i=="G"]) +  length(matrix.Pruned_Alignment_Window_i[matrix.Pruned_Alignment_Window_i=="C"]))/numeric.Total

      list.Estimated_PhylogeneticTreeModel$numeric.NumberOfSequences <- numeric.NumberOfSequences
      list.Estimated_PhylogeneticTreeModel$numeric.PercentGap <- numeric.PercentGap
      list.Estimated_PhylogeneticTreeModel$numeric.PercentN <- numeric.PercentN
      list.Estimated_PhylogeneticTreeModel$numeric.PercentGC <- numeric.PercentGC
      list.EstimatedModels_iqtree[[i]] <- list.Estimated_PhylogeneticTreeModel
      list.Estimated_PhylogeneticTreeModel_ModelParams <- list.Estimated_PhylogeneticTreeModel$list.ParameterEstimates_Window

      ###################################################
      # Append to vector for writing checkpoint results #
      ###################################################
      #vector.WindowModels <- c(i, vector.Coordinates_Window_i[1], vector.Coordinates_Window_i[2],
      #                         numeric.PercentN, numeric.PercentGap, numeric.NumberOfSequences, numeric.PercentGC,
      #                         list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.AC, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.AG, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.AT,
      #                         list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.CG, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.CT, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.GT,
      #                         list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.A, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.C, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.G, list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.T,
      #                         list.Estimated_PhylogeneticTreeModel_ModelParams$numeric.Alpha, mean(list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate$edge.length), sum(list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate$edge.length), sd(list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate$edge.length), write.tree(phy = list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate, file = ""))
      
     vector.WindowModels <- c(i, vector.Coordinates_Window_i[1], vector.Coordinates_Window_i[2],
                              numeric.PercentN, numeric.PercentGap, numeric.NumberOfSequences, numeric.PercentGC,
                              list.Estimated_PhylogeneticTreeModel_ModelParams, 
                              mean(list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate$edge.length), sum(list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate$edge.length), sd(list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate$edge.length), write.tree(phy = list.Estimated_PhylogeneticTreeModel$handle.PhylogeneticTree_Estimate, file = ""))
      

      string.Results_Window_i <- paste(vector.WindowModels, collapse = "\t")

      write(x = string.Results_Window_i, file = string.Path_File_CheckPoints_Results, append = T)

      #########################################
      # Change into dir and bzip2 the results #
      #########################################
      setwd(string.Path_Window_i)
      string.Command_Bzip2 <- "find . -type f -exec bzip2 -9 {} +"
      system(string.Command_Bzip2, intern = T)
      setwd(string.PathParentDir)
      
      ########################
      # Estimate time to run #
      ########################
      if (i == 1){
        
        numeric.EndTime <- Sys.time()
        numeric.TimeToRun <- numeric.EndTime-numeric.StartTime

        print(gsub(pattern = "XXX", replacement = round(digits = 6, numeric.TimeToRun/60), x = "Time to run the first window: XXX minutes"))
        print(gsub(pattern = "XXX", replacement = round(digits = 6, numeric.TimeToRun*numeric.NumberOfWindows/60), x = "total (ROUGH) estimated time to completiong: XXX minutes"))
      }
    }
    
    print(gsub(pattern = "XXX", replacement = i/numeric.NumberOfWindows*100, x = "Completed XXX% of windows"))
  }

  print("Returning results of Chromo.Phylome in a list (one item for each window containing a sublist of window estimates)")
  return(list.EstimatedModels_iqtree)
}
