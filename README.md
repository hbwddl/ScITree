# ScTree
ScTree: R package for scalable and robust mechanistic integration of epidemiological and genomic data for transmission tree inference 
  
This repository contains the R package (developed with Rcpp) to fit the ScTree method.

### To build and load the package via RStudio:

- Load the ScITree.Rproj file, then build the package in the ``build'' tab and install.

### To build the package from source:

- Download the "ScITree_1.1.0.tar.gz" file
- **In R** , from the directory in which you saved the .tar.gz package:

```R
install.packages("ScITree_1.1.0.tar.gz",type="source",repos=NULL)
```


An example script, which will run the inference on your own machine, is in vignettes/
