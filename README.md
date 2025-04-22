# ScITree
ScITree: R package for scalable and robust mechanistic integration of epidemiological and genomic data for transmission tree inference 
  
This repository contains the R package (developed with Rcpp) to fit the ScITree method.

### Required packages

The following packages are required to be installed in order to install and run the example code for ScITree. However, they should be loaded automatically upon installing ScITree:

- ape
- coda
- grDevices
- Rcpp (>= 1.0.0)
- utils
- BH 

### To build and load the package via RStudio:

- Load the ScITree.Rproj file, then build the package in the "build" tab and install.

### To build the package from source:

- Download the "ScITree_1.1.0.tar.gz" file
- **In R** , from the directory in which you saved the .tar.gz package:

```R
install.packages("ScITree_1.1.0.tar.gz",type="source",repos=NULL)
```

Some example scripts, which will run a short analysis on your machine and summarize its results, are in the ./test/ directory. You may source() the "simulation_inference_analysis.R" file after altering the base directory at the head of the script.
