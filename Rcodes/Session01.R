##### p.18
set.seed(12345)
x <- rnorm(400)
y <- 0.5*x + 0.5*x^2 + 0.3*x^3 + rnorm(400, sd=2)
data <- data.frame(x=x, y=y)

##### p.19
dg <- seq(1, 10)
MSE <- NULL
for (k in 1:length(dg)) {
    u <- sample(rep(seq(2), length=length(y)))
    tran <- which(u==1)
    g <- lm(y ~ poly(x, k), subset=tran)
    MSE[k] <- mean((y - predict(g, data))[-tran]^2)
}
plot(dg, MSE, type="b", col=2, xlab="Degree of Polynomial",
     ylab="Mean Squared Error", ylim=c(2,7), lwd=2, pch=19)

##### p.20
set.seed(4321)
K <- 10
dg <- seq(1, 10)
MSE1 <- matrix(0, length(dg), K)
for (i in 1:K) {
    u <- sample(rep(seq(2), length=length(y)))
    for (k in 1:length(dg)) {
        tran <- which(u==1)
        g <- lm(y ~ poly(x, k), subset=tran)
        MSE1[k, i] <- mean((y - predict(g, data))[-tran]^2)
    }
}
matplot(dg, MSE1, type="l", xlab="Degree of Polynomial", lty=1,
        ylab="Mean Squared Error", ylim=c(2,7))
aMSE1 <- apply(MSE1, 1, mean)
plot(dg, aMSE1, type="b", col=2, xlab="Degree of Polyenomial",
     ylab="Mean Squared Error", ylim=c(2,7), lwd=2, pch=19)

##### p.26
N <- 10 ## number of simulation replications
K <- 10 ## K-fold cross validation
dg <- seq(1, 10); n <- length(y)
CVE <- matrix(0, length(dg), N)
set.seed(1010)
for (j in 1:N) {
    u <- sample(rep(seq(K), length=n))
    for (i in 1:length(dg)) {
        MSE <- NULL
        for (k in 1:K) {
            tran <- which(u!=k)
            g <- lm(y ~ poly(x, i), subset=tran)
            MSE[k] <- mean((y - predict(g, data))[-tran]^2)
        }
    CVE[i, j] <- mean(MSE)
    }
}
matplot(dg, CVE, type="l", xlab="Degree of Polynomial", lty=1,
        ylab="Mean Squared Error", ylim=c(2,7), main="10-fold CV")

##### p.27
## Make a plot matrix (2 rows and 2 columns)
par(mfrow=c(2,2))
matplot(dg, MSE1, type="l", xlab="Degree of Polynomial", lty=1,
        ylab="MSE", ylim=c(3.5,6), main="Validation-set")
matplot(dg, CVE, type="l", xlab="Degree of Polynomial", lty=1,
        ylab="MSE", ylim=c(3.5,6), main="10-fold CV")
plot(dg, aMSE1, type="b", col=2, xlab="Degree of Polynomial",
     ylab="MSE", ylim=c(3.5,6), lwd=2, pch=19,
     main="Averaged Validation-set")
aCVE <- apply(CVE, 1, mean)
plot(dg, aCVE, type="b", col="darkblue", xlab="Degree of Polynomial",
     ylab="MSE", ylim=c(3.5,6), lwd=2, pch=19,
     main="Averaged 10-fold CV")

##### p.28
set.seed(1234)
n <- 100; x <- rnorm(n)
y <- 0.5*x + 0.5*x^2 + 0.3*x^3 + rnorm(n, sd=2)
data <- data.frame(x=x, y=y)
N <- 10; K <- 5; dg <- seq(1, 6)
KCV <- matrix(0, length(dg), N)
for (j in 1:N) {
    u <- sample(rep(seq(K), length=n))
    for (i in 1:length(dg)) {
        MSE <- NULL
        for (k in 1:K) {
            tran <- which(u!=k)
            g <- lm(y ~ poly(x, i), subset=tran)
            MSE[k] <- mean((y - predict(g, data))[-tran]^2)
        }
    KCV[i, j] <- mean(MSE)
    }
}
matplot(dg, KCV, type="l", xlab="Degree of Polynomial", lty=1,
        ylab="Mean Squared Error", main="5-fold CV for small samples")

##### p.29
## Install R packages `devtools' and `kogo'
## The `kogo' package contains `glmnet', `SGL', `gglasso'
install.packages("devtools")
library(devtools)
install_github("statsun78/kogo")

## Install R packages `pclgot' and `VennDiagram'
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

##### p.30
## Load the seven packages below.
## If there are no errors, you are all set.
library(kogo)
library(pclogit)
library(VennDiagram)
library(golubEsets)
library(hu6800.db)
library(AnnotationDbi)
library(qvalue)


