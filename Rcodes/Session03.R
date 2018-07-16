##### p.6
## Install a package 'gglasso'
#install.packages('gglasso')

## Load gglasso package
#library(gglasso)

## Load bardet data set
data(bardet)
str(bardet)
dim(bardet$x)
length(bardet$y)

## Define group index
grp <- rep(1:20, each=5)

## Fit group lasso for a quantitative outcome
g <- gglasso(x=bardet$x, y=bardet$y, group=grp, loss="ls")
names(g)
dim(g$beta); g$df; g$lambda

##### p.7
## Make a plot matrix (1 row, 3 columns)
par(mfrow=c(1,3))

## Plot the coefficients against the log-lambda sequence
plot(g)

## Plot group norm against the log-lambda sequence
plot(g, group=TRUE)

## Plot against the lambda sequence
plot(g, log.l=FALSE)

## Coefficients at particular lambda values
coef(g, s=c(g$lambda[1], g$lambda[5], g$lambda[10], g$lambda[20]))

## Compute group norms of 20 groups
gnorm <- function(x) sqrt(sum(x^2))
gn <- apply(g$beta, 2, function(t) tapply(t, grp, gnorm))
apply(gn, 2, function(t) sum(t!=0))

##### p.8
set.seed(12345)

## 5 fold cross-validation to pick up the optimal lambda
cv1 <- cv.gglasso(x=bardet$x, y=bardet$y, group=grp, loss="ls",
                  nfolds=5)
plot(cv1)

## Two tuning parameter values
c(cv1$lambda.min, cv1$lambda.1se)
log(c(cv1$lambda.min, cv1$lambda.1se))

## Refit with the selected lambdas
fit1 <- gglasso(x=bardet$x, y=bardet$y, group=grp, loss="ls",
                lambda=cv1$lambda.min)
fit2 <- gglasso(x=bardet$x, y=bardet$y, group=grp, loss="ls",
                lambda=cv1$lambda.1se)
cbind(fit1$beta, fit2$beta)
cbind(fit1$df, fit2$df)

##### p.9
## Load colon data set
data(colon)
str(colon)
dim(colon$x)
table(colon$y)

## Define group index
grp <- rep(1:20, each=5)

## Fit group lasso for a binary outcome
g2 <- gglasso(x=colon$x, y=colon$y, group=grp, loss="logit")
names(g2)
dim(g2$beta); g2$df; g2$lambda

## Plot the coefficients in three different ways
par(mfrow=c(1,3))
plot(g2)
plot(g2, group=TRUE)
plot(g2, log.l=FALSE)

##### p.10
set.seed(123)

## 5 fold cross-validation to pick up the optimal lambda
cv2 <- cv.gglasso(x=colon$x, y=colon$y, group=grp, loss="logit",
                  nfolds=5)
plot(cv2)

## Two tuning parameter values
c(cv2$lambda.min, cv2$lambda.1se)
log(c(cv2$lambda.min, cv2$lambda.1se))

## Refit group lasso with the selected lambda values
fit3 <- gglasso(x=colon$x, y=colon$y, group=grp, loss="logit",
                lambda=cv2$lambda.min)
fit4 <- gglasso(x=colon$x, y=colon$y, group=grp, loss="logit",
                lambda=cv2$lambda.1se)
cbind(fit3$beta, fit4$beta)
cbind(fit3$df, fit4$df) 

##### p.13
## Install a package 'SGL'
#install.packages('SGL')

## Load gglasso package
#library(SGL)

## Define group index
grp <- rep(1:20, each=5)

## Fit SGL for a quantitative outcome
data1 <- list(x = bardet$x, y = bardet$y)
gg <- SGL(data1, grp, alpha = 0.95, type = "linear")
names(gg)
dim(gg$beta)
gg$lambdas

## Pathwise solution to beta along with lambdas
matplot(log(gg$lambdas), t(gg$beta), type="l", lty=1,
        xlab="log.lambda", ylab="Coefficients")

##### p.14
set.seed(54321)

## Cross validation to pick up the optimal lambda
cv3 <- cvSGL(data1, grp, alpha=0.95, type = "linear", nfold=5)
plot(cv3)

## Create lambda.min
lambda.min <- cv3$lambda[which.min(cv3$lldiff)]
log(lambda.min)
abline(v=log(lambda.min), lty=2, lwd=2, col="red")

## Create lambda.1se
cvup <- cv3$lldiff + cv3$llSD/sqrt(5)
wh <- min(which(cv3$lldiff < cvup[which.min(cv3$lldiff)]))
lambda.1se <- cv3$lambda[wh]
log(lambda.1se)
abline(v=log(lambda.1se), lty=2, lwd=2, col="blue")        

##### p.15
## Refit SGL with the selected lambda values
gfit1 <- SGL(data1, grp, alpha = 0.95, type = "linear",
             lambda=lambda.min)
gfit2 <- SGL(data1, grp, alpha = 0.95, type = "linear",
             lambda=lambda.1se)

## Coefficients for two lambda values
coef <- cbind(gfit1$beta, gfit2$beta)
colnames(coef) <- c("lambda.min", "lambda.1se")
rownames(coef) <- paste("V", seq(1:100), sep="")
coef

## The number of selected genes
apply(coef, 2, function(t) sum(t!=0))

##### p.20
## Install the packages 'devtools' and 'pclogit'
#install.packages('devtools')
#library(devtools)
#install_github("statsun78/pclogit")

## Load pclogit package
#library(pclogit)

## Define group size
## Do not confused with group index grp
## grp <- rep(1:20, each=5)
gsize <- rep(5, 20)

## Change y coding to 0-1 values
y <- colon$y
y[y==-1] <- 0

## Fit pclogit for a binary outcome using ring and f.con networks
gp1 <- pclogit(colon$x, y, alpha=0.1, group=gsize, type="ring")
gp2 <- pclogit(colon$x, y, alpha=0.1, group=gsize, type="fcon")

##### p.21
## The number of selected genes
cbind(gp1$df, gp2$df)

## Pathwise solution to beta along with lambda
par(mfrow=c(1,2))
matplot(log(gp1$lambda), t(gp1$beta), type="l", lty=1,
        xlab="log.lambda", ylab="Coefficients", main="Ring network")
matplot(log(gp2$lambda), t(gp2$beta), type="l", lty=1,
        xlab="log.lambda", ylab="Coefficients", main="F.con network")
set.seed(1234)

## Selection probabilities
sp1 <- sel.pclogit(colon$x, y, alpha=0.1, group=gsize, type="ring")
sp2 <- sel.pclogit(colon$x, y, alpha=0.1, group=gsize, type="fcon")

## SP ranking of genes
SP <- cbind(sp1$maxsel, sp2$maxsel)
colnames(SP) <- c("ring.genes", "ring.sp", "fcon.genes", "fcon.sp")
SP

##### p.22
## Golub data
#library(golubEsets)
data(Golub_Merge)
annotation(Golub_Merge)

## Install annotation packages for golub data
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("bit")
#biocLite("AnnotationDbi")
#biocLite("hu6800.db")
#library(hu6800.db)
#library(AnnotationDbi)

## Convert probe ID to gene symbol
symb <- mget(featureNames(Golub_Merge), env=hu6800SYMBOL)
sum(!is.na(unlist(symb)))
dim(exprs(Golub_Merge))        

##### p.23
## The number of genes
ww <- which(is.na(unlist(symb)))
genes <- unlist(symb)[-ww]
c(length(genes), length(unique(genes)))

## List genes by alphabetical order
oo <- order(genes)
genes <- genes[oo]
genes[1:10]

## Arrange gene expressions corresponding to gene names
x <- scale(t(exprs(Golub_Merge)))
x <- x[,-ww]
x <- x[,oo]
dim(x)

## Generate a phenotype outcome
y <- as.numeric(Golub_Merge$ALL.AML)
y[y==2] <- -1        

##### p.24
## Generate group index
grp <- as.numeric(factor(genes))

## Perform cross-validation for group lasso
set.seed(123)
gls <- cv.gglasso(x, y, group=grp, loss="logit", nfolds=5)
gls1 <- gglasso(x, y, group=grp, loss="logit", lambda=gls$lambda.min)
gls2 <- gglasso(x, y, group=grp, loss="logit", lambda=gls$lambda.1se)

## The result of selected genes
cbind(gls1$df, gls2$df)
ww1 <- which(gls1$beta!=0)
ww2 <- which(gls2$beta!=0)
data.frame(genes=genes[ww1], beta=gls1$beta[ww1])
data.frame(genes=genes[ww2], beta=gls2$beta[ww2])

##### p.25
## Generate group sizes
gsize <- as.numeric(table(grp))

## Recode the phenotype outcome
ry <- y
ry[y==-1] <- 0

## Compute selection probabilities
set.seed(111)
sp1 <- sel.pclogit(x, ry, alpha=0.1, group=gsize, type="ring")
sp2 <- sel.pclogit(x, ry, alpha=0.1, group=gsize, type="fcon")

## Display top 20 genes ranked by selection probabilities
top <- 20
data.frame(sp1$maxsel[1:top,], genes=genes[sp1$maxsel[1:top,1]])
data.frame(sp2$maxsel[1:top,], genes=genes[sp2$maxsel[1:top,1]])

