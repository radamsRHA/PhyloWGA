#' Define.PhyloWGA_Experiment: function to generate a matrix of window coordinates for a given alignment
#'
#' This function returns a matrix of the window coordinates for an alignment
#' @param numeric.WindowSize Number of bp for each window, window length
#' @param numeric.StepSize Number of bp seperating each window
#' @param numeric.TotalLength Total number of bp for the entire alignment
#' @keywords phylogenetics, whole genome, species tree, coalescent theory,
#' @return matrix.WindowCoordinates Matrix of k x 2 dimensions, each row indicating the start and end coordinates for the respective window
#' @export
#' @examples
#'
#'
#' ###################################
#' # Specify experimental parameters #
#' ###################################
#' matrix.ExperimentalParams <- Define.PhyloWGA_Experiment(numeric.WindowSize = 1000, numeric.StepSize = 500, numeric.TotalLength = 10000)
#'


#####################################
# Define.PhyloWGA_Experiment #
#####################################
Define.PhyloWGA_Experiment <- function(numeric.WindowSize, numeric.StepSize, numeric.TotalLength){

  ##############################
  # Summarize Window Positions #
  ##############################
  #vector.WindowStartPositions <- seq(from = 1, to = (numeric.TotalLength-numeric.WindowSize), by = numeric.StepSize)
  #vector.WindowStartPositions <- seq(from = 1, to = (numeric.TotalLength), by = numeric.StepSize)

  if (numeric.TotalLength %% numeric.WindowSize == 0){
    vector.WindowStartPositions <- seq(from = 1, to = (numeric.TotalLength), by = numeric.StepSize)
  }
  if (numeric.TotalLength %% numeric.WindowSize != 0){
    vector.WindowStartPositions <- seq(from = 1, to = (numeric.TotalLength-numeric.WindowSize), by = numeric.StepSize)
  }

  vector.WindowEndPositions <- vector.WindowStartPositions + numeric.WindowSize - 1
  numeric.NumberOfWindows <- length(vector.WindowStartPositions)

  ################################################
  # Make matrix that contains window coordinates #
  ################################################
  matrix.WindowCoordinates <- matrix(nrow = numeric.NumberOfWindows, ncol = 2)
  rownames(matrix.WindowCoordinates) <- paste0("Window_", 1:numeric.NumberOfWindows)
  colnames(matrix.WindowCoordinates) <- c("Coordinate_Start", "Coordinate_End")
  matrix.WindowCoordinates[,1] <- vector.WindowStartPositions
  matrix.WindowCoordinates[,2] <- vector.WindowEndPositions

  ############################
  # Print summary of results #
  ############################
  string.PrintSummary <- "Generated coordinates for k = AAA chromsome windows each of size BBB with start sites spaced apart by CCC sites"
  string.PrintSummary <- gsub(pattern = "AAA", replacement = numeric.NumberOfWindows, x = string.PrintSummary)
  string.PrintSummary <- gsub(pattern = "BBB", replacement = numeric.WindowSize, x = string.PrintSummary)
  string.PrintSummary <- gsub(pattern = "CCC", replacement = numeric.StepSize, x = string.PrintSummary)
  print(string.PrintSummary)

  return(matrix.WindowCoordinates)

}
