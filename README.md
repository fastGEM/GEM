GEM
----------------------


## Installation

```R
## install dependent packages
if(!require("methods")){
    install.packages("methods")
}
if(!require("ggplot2")){
    install.packages("ggplot2")
}
if(!require("Rcpp")){
    install.packages("Rcpp")
}

## install gem package
install.packages("devtools")
devtools::install_github("fastGEM/GEM")
```

## Launch the GUI

```R
GEM_GUI()
```
