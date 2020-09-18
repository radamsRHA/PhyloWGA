#' Chromo.Crawl: (original version) function to crawl over an alignment of chromosomes while attempting to concatenate adjacent genomic windows
#'
#' This function writes the results to an output file. Additionally, this function provides a directory where the resulting supergenes are placed.
#' @param string.PathParentDir Path to parent dir for all analyses
#' @param string.PathToFastaFile Path to fasta file
#' @param numeric.WindowSize Size of each chromosomal window
#' @param numeric.StepSize Spacing of each chromosomal window
#' @param numeric.NumberOfCores Number of cores
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return Results Results written to outputfiles
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
#' #########################
#' # Conduct.ChromoPhylome #
#' #########################
#' Chromo.Crawl(string.PathParentDir = '~/Desktop/',
#'                          numeric.WindowSize = 1000,
#'                          numeric.StepSize = 1000,
#'                          string.PathToFastaFile = String.Path_ExampleChromosomeAlignment)
#'

################
# Chromo.Crawl #
################
Chromo.Crawl <- function(string.PathParentDir, numeric.WindowSize, numeric.StepSize, string.PathToFastaFile, numeric.NumberOfCores){
  
  #############################################
  # Define parent dir used for the experiment #
  #############################################
  string.Path_ParenDir_ChromoCrawl_Experiment = paste(string.PathParentDir, '/Conduct.ChromoPhyloCrawl_', Sys.Date(), sep = "")
  unlink(string.Path_ParenDir_ChromoCrawl_Experiment, recursive = T)
  dir.create(string.Path_ParenDir_ChromoCrawl_Experiment, showWarnings = T, recursive = T)
  
  #####################################################
  # Define parent dir to collect supergene alignments #
  #####################################################
  string.Path_SupergeneAlignments = paste0(string.Path_ParenDir_ChromoCrawl_Experiment, '/SupergeneAlignments/')
  unlink(string.Path_SupergeneAlignments, recursive = T)
  dir.create(string.Path_SupergeneAlignments, showWarnings = T, recursive = T)
  
  ####################################
  # File to write checkpoint results #
  ####################################
  string.Path_File_CheckPoints_Results <- paste0(string.Path_ParenDir_ChromoCrawl_Experiment, '/Results_CheckPoint_w', numeric.WindowSize, "_s", numeric.StepSize, ".txt")
  write(x = "Supergene\tStart_i\tEnd_i\tNumberWindows\tTrees", file = string.Path_File_CheckPoints_Results)
  
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
  
  ##############################
  # Count windows crawled over #
  ##############################
  numeric.Count_Crawler_Windows <- 0
  numeric.Count_Supergenes <- 0
  
  ######################
  # Crawl over windows #
  ######################
  while (numeric.Count_Crawler_Windows < numeric.NumberOfWindows){
    
    ##########################
    # Extract current window #
    ##########################
    numeric.Count_Crawler_Windows <- numeric.Count_Crawler_Windows + 1
    vector.Coordinates_Window_i <- matrix.Experimental_WindowCoordinates[numeric.Count_Crawler_Windows,]
    matrix.Alignment_Window_i <- Extract.WindowAlignment(string.PathToFastaFile = string.PathToFastaFile,vector.SequenceNames = vector.Fasta_SequenceNames,vector.Coordinates_ExperimentalLocus = vector.Coordinates_Window_i)
    matrix.Pruned_Alignment_Window_i <- Prune.MissingDataTaxa(matrix.WindowAlignment = matrix.Alignment_Window_i, numeric.MissingDataThreshold = 0.95)
    vector.Pruned_Alignment_Window_i_Names <- rownames(matrix.Pruned_Alignment_Window_i)
    
    ########################################
    # Check if pruned window stille exists #
    ########################################
    if (length(matrix.Pruned_Alignment_Window_i) > 1 & length(matrix.Pruned_Alignment_Window_i[,1])> 3){
      
      #########################
      # Found a new supergene #
      #########################
      numeric.Count_Supergenes <- numeric.Count_Supergenes + 1
      numeric.Count_Concatenated_Windows <- 1
      string.Path_Supergene_i = paste0(string.Path_ParenDir_ChromoCrawl_Experiment, '/Supergene_', numeric.Count_Supergenes)
      unlink(string.Path_Supergene_i, recursive = T)
      dir.create(string.Path_Supergene_i, showWarnings = T, recursive = T)
      
      #######################################
      # Count adjacent windows crawler over #
      #######################################
      numeric.Count_AdjacentWindows <- numeric.Count_Crawler_Windows
      
      ####################################
      # Crawl over next adjacent windows #
      ####################################
      while (numeric.Count_AdjacentWindows < numeric.NumberOfWindows){
        
        ########################
        # Extract nexxt window #
        ########################
        matrix.Alignment_Window_j <- vector.Pruned_Alignment_Window_j_Names <- vector.IntersectWindowSeqNames <- NA
        numeric.Count_AdjacentWindows <- numeric.Count_AdjacentWindows + 1
        numeric.Count_Crawler_Windows <- numeric.Count_Crawler_Windows + 1
        vector.Coordinates_Window_j <- matrix.Experimental_WindowCoordinates[(numeric.Count_AdjacentWindows),]
        matrix.Alignment_Window_j <- Extract.WindowAlignment(string.PathToFastaFile = string.PathToFastaFile,vector.SequenceNames = vector.Fasta_SequenceNames,vector.Coordinates_ExperimentalLocus = vector.Coordinates_Window_j)
        matrix.Pruned_Alignment_Window_j <- Prune.MissingDataTaxa(matrix.WindowAlignment = matrix.Alignment_Window_j, numeric.MissingDataThreshold = 0.95)
        vector.Pruned_Alignment_Window_j_Names <- rownames(matrix.Pruned_Alignment_Window_j)
        vector.IntersectWindowSeqNames <- unique(c(intersect(vector.Pruned_Alignment_Window_j_Names, vector.Pruned_Alignment_Window_i_Names),intersect(vector.Pruned_Alignment_Window_i_Names, vector.Pruned_Alignment_Window_j_Names)))
        
        #######################################
        # Check if pruned window still exists #
        #######################################
        if (length(matrix.Pruned_Alignment_Window_j) > 1  & length(matrix.Pruned_Alignment_Window_j[,1])> 3 & length(vector.IntersectWindowSeqNames) >= 4){
          
          
          ####################################
          # Found an adjacent window to test #
          ####################################
          #string.Path_Supergene_i_Window_j = paste0(string.Path_Supergene_i, '/Supergene_', numeric.Count_Supergenes, "_Window_", numeric.Count_AdjacentWindows)
          
          if (numeric.Count_Concatenated_Windows == 1){
            string.Path_Supergene_i_Window_j <- NA
            string.Path_Supergene_i_Window_j = paste0(string.Path_Supergene_i, '/Supergene_', numeric.Count_Supergenes, "_Window_", (numeric.Count_Crawler_Windows-1), "_Window_", numeric.Count_AdjacentWindows)
            unlink(string.Path_Supergene_i_Window_j, recursive = T)
            dir.create(string.Path_Supergene_i_Window_j, showWarnings = T, recursive = T)
          }
          
          if (numeric.Count_Concatenated_Windows > 1){
            string.Path_Supergene_i_Window_j <- NA
            string.Path_Supergene_i_Window_j = paste0(string.Path_Supergene_i, '/Supergene_', numeric.Count_Supergenes, "_Window_", numeric.Count_AdjacentWindows)
            unlink(string.Path_Supergene_i_Window_j, recursive = T)
            dir.create(string.Path_Supergene_i_Window_j, showWarnings = T, recursive = T)
          }

          ####################################################
          # Conduct test of combinability with CONCATEPILLAR #
          ####################################################
          handle.Results_CONCATEPILLAR <- Execute.PhylogeneticCongruence_concatepillar(matrix.WindowAlignment_01 = matrix.Pruned_Alignment_Window_i,
                                                                                       matrix.WindowAlignment_02 = matrix.Pruned_Alignment_Window_j,
                                                                                       string.PathParentDir = string.Path_Supergene_i_Window_j,
                                                                                       numeric.NumberOfCores = numeric.NumberOfCores)
          
          ##################################
          # Evaluate CONCATEPILLAR results #
          ##################################
          if (handle.Results_CONCATEPILLAR$boo.Concate == T){
            
            ###################
            # Add to counters #
            ###################
            matrix.Pruned_Alignment_Window_i <- handle.Results_CONCATEPILLAR$matrix.ConcatenatedWindows
            numeric.Count_Concatenated_Windows <- numeric.Count_Concatenated_Windows + 1
            
            print(gsub(pattern = "XXX", replacement = numeric.Count_Crawler_Windows/numeric.NumberOfWindows*100, x = "Crawled across XXX% of windows"))
            
            
          }
          if (handle.Results_CONCATEPILLAR$boo.Concate == F){
            
            ############################################
            # Add to counter because supergene is done #
            ############################################
            numeric.Count_Crawler_Windows <- (numeric.Count_AdjacentWindows-1)
            handle.Results_CONCATEPILLAR$handle.SupergeneWindow_tree <- handle.Results_CONCATEPILLAR$handle.Window_tree_01
            handle.Results_CONCATEPILLAR$matrix.ConcatenatedWindows <- handle.Results_CONCATEPILLAR$matrix.WindowAlignment_01
            
            ##################################
            # Summarize start and end points #
            ##################################
            numeric.SupergeneStart <- vector.Coordinates_Window_i[1]
            numeric.SupergeneEnd <- matrix.Experimental_WindowCoordinates[(numeric.Count_AdjacentWindows-1),2]
            
            ############
            # ADD HERE #
            ############
            if (numeric.Count_AdjacentWindows < numeric.NumberOfWindows){
              vector.SupergeneResults <- c(numeric.Count_Supergenes, numeric.SupergeneStart, numeric.SupergeneEnd, numeric.Count_Concatenated_Windows, write.tree(phy = handle.Results_CONCATEPILLAR$handle.SupergeneWindow_tree, file = ""))
              write(x = paste(vector.SupergeneResults, collapse = "\t"), file = string.Path_File_CheckPoints_Results, append = T)
              write.FASTA(x = as.DNAbin(handle.Results_CONCATEPILLAR$matrix.ConcatenatedWindows), file = paste0(string.Path_SupergeneAlignments, '/Supergene_', numeric.Count_Supergenes, '.fasta'))
            }
            if (numeric.Count_AdjacentWindows == numeric.NumberOfWindows){
              
              vector.SupergeneResults <- c(numeric.Count_Supergenes, numeric.SupergeneStart, numeric.SupergeneEnd, numeric.Count_Concatenated_Windows, write.tree(phy = handle.Results_CONCATEPILLAR$handle.SupergeneWindow_tree, file = ""))
              write(x = paste(vector.SupergeneResults, collapse = "\t"), file = string.Path_File_CheckPoints_Results, append = T)
              write.FASTA(x = as.DNAbin(handle.Results_CONCATEPILLAR$matrix.ConcatenatedWindows), file = paste0(string.Path_SupergeneAlignments, '/Supergene_', numeric.Count_Supergenes, '.fasta'))
              
              
              numeric.LastSupergeneStart <- matrix.Experimental_WindowCoordinates[(numeric.Count_AdjacentWindows),1]
              numeric.LastSupergeneEnd <- matrix.Experimental_WindowCoordinates[(numeric.Count_AdjacentWindows),2]
              vector.SupergeneResults <- c(numeric.Count_Supergenes+1, numeric.LastSupergeneStart, numeric.LastSupergeneEnd, 1, write.tree(phy = handle.Results_CONCATEPILLAR$handle.Window_tree_02, file = ""))
              write(x = paste(vector.SupergeneResults, collapse = "\t"), file = string.Path_File_CheckPoints_Results, append = T)
              write.FASTA(x = as.DNAbin(handle.Results_CONCATEPILLAR$matrix.ConcatenatedWindows), file = paste0(string.Path_SupergeneAlignments, '/Supergene_', numeric.Count_Supergenes, '.fasta'))
            
            }
            
            
          
            break
          }
        }
      }
      
      
      ########################
      ## Summarize supergene #
      ########################
      ##vector.SupergeneResults <- c(numeric.Count_Supergenes, vector.Coordinates_Window_i[1], vector.Coordinates_Window_j[2], numeric.Count_Concatenated_Windows, write.tree(phy = handle.Results_CONCATEPILLAR$handle.SupergeneWindow_tree, file = ""))
      #vector.SupergeneResults <- c(numeric.Count_Supergenes, numeric.SupergeneStart, numeric.SupergeneEnd, numeric.Count_Concatenated_Windows, write.tree(phy = handle.Results_CONCATEPILLAR$handle.SupergeneWindow_tree, file = ""))
      #write(x = paste(vector.SupergeneResults, collapse = "\t"), file = string.Path_File_CheckPoints_Results, append = T)
      #write.FASTA(x = as.DNAbin(handle.Results_CONCATEPILLAR$matrix.ConcatenatedWindows), file = paste0(string.Path_SupergeneAlignments, '/Supergene_', numeric.Count_Supergenes, '.fasta'))
      
      #########################################
      # Change into dir and bzip2 the results #
      #########################################
      setwd(string.Path_Supergene_i)
      string.Command_Bzip2 <- "find . -type f -exec bzip2 -9 {} +"
      system(string.Command_Bzip2, intern = T, ignore.stdout = T, ignore.stderr = T)
      setwd(string.PathParentDir)
      
      #########################################
      # Change into dir and bzip2 the results #
      #########################################
      options(warn = -1)
      setwd(string.Path_SupergeneAlignments)
      string.Command_Bzip2 <- "find . -type f -exec bzip2 -9 {} +"
      system(string.Command_Bzip2, intern = T, ignore.stdout = T, ignore.stderr = T)
      options(warn = 0)
      setwd(string.PathParentDir)
      
    }
  }
}
