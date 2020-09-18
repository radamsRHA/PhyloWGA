---
# PhyloWGA: phylogenetic analyses and interrogation of whole genome alignments (WGAs)
**NOTE See the file https://github.com/radamsRHA/PhyloWGA/blob/master/PhyloWGA.pdf for more detailed instructions**

## Installing R package PhyloWGA from github
The R package PhyloWGA is freely available to download and distribute from github <https://github.com/radamsRHA/PhyloWGA/>. To install and load PhyloWGA, you must first install the R package `devtools`, Additionally, make sure the most updated version of R is installed 

```
install.packages("devtools")
```
Now using devtools we can install `PhyloWGA` from github:

```
library(devtools)
install_github("radamsRHA/PhyloWGA")
library(PhyloWGA) # Load package ThetaMater
```

Additionally, a number of dependencies are required by `PhyloWGA`, including the following:
```
library(PhyloWGA)
library(ape)
library("phylotools")
```

`PhyloWGA` also requires the following dependencies to be installed and available on your path `$PATH`:

* samtools: http://www.htslib.org (`samtools` must be in your path)
* concatepillar: http://leigh.geek.nz/software.shtml (`concaterpillar.py` must be in your path, can be found also at github.com/radamsRHA/PhyloWGA/inst/extdata/concaterpillar-1.8a.zip)
* iqtree: http://www.iqtree.org (`iqtree` must be in your path, can be found also at github.com/radamsRHA/PhyloWGA/inst/extdata/iqtree-1.6.12-MacOSX.zip)
* raxml:must be version 7.3! (used by concatepillar.py, `raxmlHPC` must be in your path, can be found also at github.com/radamsRHA/PhyloWGA/inst/extdata/raxmlHPC.zip)

To begin using `PhyloWGA` try using the examples associated with each function and the file `PhyloWGA.pdf` provided with the R package. 

## Example: Phylogenetic analyses along windows of a chromosome alignment 

In this first example, we will use the function `Chromo.Phylome` to infer a chromosome-specific set of gene trees (i.e., "chromo phylome") for the example segment. 

First, let's load `PhyloWGA` and its dependencies and find the path to the example chromosome alignment 


```
################
# Load depends #
################
library(PhyloWGA)
library(ape)
library("phylotools")

#####################################
# Read example chromosome alignment #
#####################################
String.Path_ExampleChromosomeAlignment <-  system.file("extdata", "Example_Chr10.fasta", package="PhyloWGA")
```

Now, let's use `Chromo.Phylome` to infer a set of gene trees from this alignment using a very fine scale (window size = step size = 100bp) analyses for this example.`Chromo.Phylome` will proceed along the alignment in windows of 100bp using a GTR+G model. 

```
####################################
# Conduct chromo phylome inference #
####################################
handle.Chr1_w100_s100_Phylome <- Chromo.Phylome(string.PathParentDir = '~/Desktop/', 
                                                          numeric.WindowSize = 100, 
                                                          numeric.StepSize = 100, 
                                                          string.PathToFastaFile = String.Path_ExampleChromosomeAlignment, 
                                                          string.Commands_iqtree= "-m GTR+G")
```

After running `Chromo.Phylome`, check the output directory on the Desktop to see the results, including the file `Results_CheckPoint_w100_s100.txt` that summarizes the run. Also, the object `handle.Chr1_w100_s100_Phylome` contains the results of the analysis. 



## Example: construct a set of statistically-justified and genome-informed supergenes 

The function `Chromo.Crawl` will "crawl" along the input alignment and apply the likelihood ratio test of CONCATEPILLAR to infer model-based supergenes (i.e., genomic windows that have been concated together due to evidence of a shared tree). Let's run `Chromo.Crawl` on our example dataset:

```
########################
# Conduct chromo crawl #
########################
handle.Chr1_w100_s100_Crawler <- Chromo.Crawl(string.PathParentDir = '~/Desktop/', 
                         numeric.WindowSize = 100, 
                         numeric.StepSize = 100, 
                         string.PathToFastaFile = String.Path_ExampleChromosomeAlignment, 
                         numeric.NumberOfCores = 2)
```

The results of `Chromo.Crawl` will be found on the parent directory placed on the Desktop, and include a directory of supergene alignments (one for each supegene), and a results file that summarizes the analysis `Results_CheckPoint_w10000_s10000.txt`. These can then be used to further dissect evidence of phylogenetic conflict, and the resulting supergene alignments (found in the subdirectory `SupergeneAlignments`) can also be used for downstream species tree inference, and so forth. 


## Example: setting up PhyloWGA for parallel analyses

The function `Organize.ParallelPhyloWGA` can be used to organize a directory that contains a set of subdirectories, each of which contain a slice of a WGA and an R script for conducting a particular PhyloWGA analyses. Batch scripts can then be run on each subdirectory. For example, one can provide sbatch scripts for each subdirectory for analyses on a cluser. 

```
###############################################
# Organize PhyloWGA for Chromo.Crawl analysis #
###############################################
Organize.ParallelPhyloWGA(numeric.NumberSubsets = 2, 
                          string.PathToFastaFile = String.Path_ExampleChromosomeAlignment, 
                          string.PathParentDir = '~/Desktop/', 
                          string.Analysis = "Chromo.Crawl", 
                          numeric.WindowSize = 1000,
                          numeric.StepSize = 1000, 
                          string.Commands_iqtree = "", 
                          numeric.NumberOfCores = 2)

```

## Example: conducting Chromo.Phylome on a customized set of loci

The function `Chromo.Phylome.Custom` allows users to infer locus-specific phylogenetic trees given a user-specified matrix of locus coordinates. You can provide a 2-columns matrix with the number of rows equal to the number of loci that will be analyzed. See the arguments `matrix.WindowCoordinates = matrix.LocusCoordinates` below. 

```
LocusCoordinates <- matrix(nrow = 3, ncol = 2)
LocusCoordinates[1,] <- c(1,500)
LocusCoordinates[2,] <- c(600,1200)
LocusCoordinates[3,] <- c(2000,5000)

Chromo.Phylome.Custom(string.PathParentDir = "~/Desktop/", 
                      string.PathToFastaFile = String.Path_ExampleChromosomeAlignment, 
                      matrix.WindowCoordinates = matrix.LocusCoordinates,
                      string.Commands_iqtree = "")
```



