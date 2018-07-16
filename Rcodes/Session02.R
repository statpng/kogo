##### p.6
## Install a bioconductor package 'golubEsets'
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("golubEsets")
#library(golubEsets)
data(Golub_Merge)

## Golub data
?Golub_Merge
dim(Golub_Merge)
varLabels(Golub_Merge)

## Covariates
Golub_Merge$ALL.AML
Golub_Merge$Gender

## Expression data
dim(exprs(Golub_Merge))
head(exprs(Golub_Merge))
tail(exprs(Golub_Merge))
hist(exprs(Golub_Merge), nclass=50, col="orange")

##### p.7
## Install a package 'glmnet'
#install.packages('glmnet')
#library(glmnet)

## Generate x (genomic data) and y (outcome)
x <- scale(t(exprs(Golub_Merge)))
hist(x, nclass=50, col="orange")
dim(x)
y <- as.numeric(Golub_Merge$ALL.AML)
length(y)

## Fit a ridge regression between x and y
g <- glmnet(x, y, alpha=0, family="binomial")
plot(g, "lambda")
g$lambda
g$df
dim(g$beta)

##### p.10
## Fit a lasso
g2 <- glmnet(x, y, family="binomial")
g2 <- glmnet(x, y, alpha=1, family="binomial")
plot(g2, "lambda")
g2$lambda
log(g2$lambda)
dim(g2$beta)

## Total number of selected genes
g2$df
apply(g2$beta, 2, function(x) sum(x!=0))

## Coefficient estimates for selected genes
k <- 50
wh <- which(g2$beta[, k]!=0)
g2$beta[wh, k]
abline(v=log(g2$lambda[k]), lwd=2, lty=2, col="red")

##### p.13
## Random split: set the seed number
set.seed(1234)

## 5-fold CV for Lasso
gs <- cv.glmnet(x, y, alpha=1, family="binomial", nfolds=5)
plot(gs)
names(gs)

## Optimal lambda for minimum CVE
gs$lambda[which.min(gs$cvm)]
gs$lambda.min
log(gs$lambda.min)

## One-standard-error rule
wh <- min(which(gs$cvm < gs$cvup[which.min(gs$cvm)]))
gs$lambda[wh]
gs$lambda.1se
log(gs$lambda.1se)

##### p.14
## Re-fit the lasso with two selected lambda values
fit1 <- glmnet(x, y, family="binomial", lambda=gs$lambda.min)
fit2 <- glmnet(x, y, family="binomial", lambda=gs$lambda.1se)

## The number of selected genes.
fit1$df
fit2$df

## Location of selected genes
w1 <- which(as.matrix(fit1$beta) != 0)
w2 <- which(as.matrix(fit2$beta) != 0)

## Coefficient estimates and their probe id.
data.frame(probeID=rownames(fit1$beta)[w1], beta=fit1$beta[w1])
data.frame(probeID=rownames(fit2$beta)[w2], beta=fit2$beta[w2])

##### p.15
## Independent two-sample t-test
pval <- apply(x, 2, function(t) t.test(t[y==1], t[y==2])$p.val)
b.pval <- p.adjust(pval, method="bonferroni")
wp <- which(b.pval < 0.05)
length(wp)

## Wilcoxon rank sum test
pval2 <- apply(x, 2, function(t) wilcox.test(t ~ y)$p.val)
b.pval2 <- p.adjust(pval2, method="bonferroni")
wp2 <- which(b.pval2 < 0.05)
length(wp2)

## Overlapped genes selected by lambda.min and t-test
intersect(w1, wp)
colnames(x)[intersect(w1, wp)]

## Overlapped genes selected by lambda.min and Wilcoxon test
intersect(w1, wp2)
colnames(x)[intersect(w1, wp2)]

##### p.24
## Consider three covariates (non-genetic effects)
varLabels(Golub_Merge)
Golub_Merge$BM.PB
Golub_Merge$PS
Golub_Merge$Source
z <- cbind(Golub_Merge$BM.PB, Golub_Merge$PS, Golub_Merge$Source)

## Create new x-matrix (3 covariates + 7129 genes)
xz <- cbind(z, x)
dim(xz)
c(ncol(x), ncol(xz))
colnames(xz)[1:3] <- c("BM.PM", "PS", "Source")

## Create a penalty factor.
## Covariates (non-penalized): 0 and genes (penalized): 1
pp <- c(rep(0, ncol(z)), rep(1, ncol(x)))

##### p.25
set.seed(1234)

## Fix group ID for 5-fold cross-validation
## Balance group ID between ALL(y=1) and AML(y=2)
fold1 <- sample(1:5, size=sum(y==1), replace=TRUE)
fold2 <- sample(1:5, size=sum(y==2), replace=TRUE)
foldid <- rep(0, length(y))
foldid[y==1] <- fold1
foldid[y==2] <- fold2

## Find the optimal alpha first
alp <- seq(0, 1, 0.05)
cvm <- NULL
for (j in 1:length(alp)) {
    g <- cv.glmnet(xz, y, alpha=alp[j], family="binomial",
    foldid=foldid, penalty.factor=pp)
    cvm[j] <- min(g$cvm)
}
wh <- which.min(cvm)
alp[wh]

##### p.26
## Find the optimal lambda
gp <- cv.glmnet(xz, y, alpha=alp[wh], family="binomial",
foldid=foldid, penalty.factor=pp)

## Re-fit the elastic net with two selected lambda values
fit3 <- glmnet(xz, y, alpha=alp[wh], family="binomial",
               penalty.factor=pp, lambda=gp$lambda.min)
fit4 <- glmnet(xz, y, alpha=alp[wh], family="binomial",
               penalty.factor=pp, lambda=gp$lambda.1se)

## The number of selected genes.
c(fit3$df, fit4$df)

## Location of selected genes
w3 <- which(as.matrix(fit3$beta) != 0)
w4 <- which(as.matrix(fit4$beta) != 0)

## Coefficient estimates and their probe id.
data.frame(probeID=rownames(fit3$beta)[w3], beta=fit3$beta[w3])
data.frame(probeID=rownames(fit4$beta)[w4], beta=fit4$beta[w4])


