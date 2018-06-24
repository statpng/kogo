# kogo
2018 KOGO statistical genetics workshop 

## Installation
```
## Install R packages 'devtools' and 'kogo'
## The 'kogo' package contains 'glmnet', 'SGL', 'gglasso'
install.packages("devtools")
library(devtools)
install_github("statsun78/kogo")

## Install R packages 'pclgot' and 'VennDiagram'
## For window users, 'pclogit' requires the program 'Rtools'.
install_github("statsun78/pclogit")
install.packages("VennDiagram")

## Install Bioconductor R packages
source("http://www.bioconductor.org/biocLite.R")
biocLite("golubEsets")
biocLite("bit")
biocLite("AnnotationDbi")
biocLite("hu6800.db")
biocLite("qvalue")

## Load seven R packages. If there are no errors occurs, you are all set.
library(kogo)
library(pclogit)
library(VennDiagram)
library(golubEsets)
library(hu6800.db)
library(AnnotationDbi)
library(qvalue)
```
