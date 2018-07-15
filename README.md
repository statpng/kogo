# kogo
2018 KOGO statistical genetics workshop
(6. Multi-omics Data Analysis by Hokeun Sun)

## Installation
You need to install:

1. `R` from [__CRAN__](https://cran.r-project.org/)
2. `devtools` package
3. development tools (`rtools` for _Windows users_ and `Xcode Command Line Tools` for _macOS users_)

### 
```
## Install R packages 'devtools' and 'kogo'
## The 'kogo' package contains 'glmnet', 'SGL', 'gglasso'
install.packages("devtools")
library(devtools)
install_github("statsun78/kogo")

## Install R packages 'pclgot' and 'VennDiagram'
## For Window users, 'pclogit' requires the program 'Rtools'.
install_github("statsun78/pclogit")
install.packages("VennDiagram")

## Install Bioconductor R packages
source("http://www.bioconductor.org/biocLite.R")
biocLite("golubEsets")
biocLite("bit")
biocLite("AnnotationDbi")
biocLite("hu6800.db")
biocLite("qvalue")

## Load the seven packages below. If there are no errors, you are all set.
library(kogo)
library(pclogit)
library(VennDiagram)
library(golubEsets)
library(hu6800.db)
library(AnnotationDbi)
library(qvalue)
```
