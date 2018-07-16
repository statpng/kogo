##### p.7
## Open the ovarian cancer data
#library(kogo)
data(ovarian27k)
?ovarian27k
str(ovarian27k)
names(ovarian27k)

## `genename' contains Illumina ID and gene symbols
head(ovarian27k$genename)

## `group' includes group sizes of 12,770 genes
gsize <- ovarian27k$group
c(length(gsize), sum(gsize))
table(gsize)

## Find the gene name that has the most number of CpG sites
ww <- which.max(gsize)
cc <- cumsum(gsize)[(ww-1):ww]
ovarian27k$genename[(cc[1]+1):cc[2], c(2,3)]

##### p.8
## Phenotype outcome
y <- ovarian27k$y
y
length(y)
table(y)

## DNA methylation beta values
x <- ovarian27k$x
dim(x)
summary(as.numeric(x))
hist(x, nclass=50, col="orange")
set.seed(25)
par(mfrow=c(5,5))
s <- sample(1:ncol(x), 25)
for (i in 1:25) {
    si <- s[i]
    hist(x[,si], nclass=20, col="orange", xlab="",
    main=paste("CpG site = ", s[i]))
}

##### p.9
## The p-values of independent two-sample t-test
pval <- apply(x, 2, function(t) t.test(t[y==0], t[y==1])$p.val)

## Bonferroni adjustment
b.pval <- p.adjust(pval, method="bonferroni")
wp1 <- which(b.pval < 0.05)
oo1 <- order(b.pval[wp1])
length(wp1)
gn1 <- ovarian27k$genename[wp1[oo1], c(2,3)]
data.frame(gn1, p.values=b.pval[wp1[oo1]])

## Benjamini and Hochberg adjustment
fdr <- p.adjust(pval, method="BH")
wp2 <- which(fdr < 0.05)
oo2 <- order(fdr[wp2])
length(wp2)
gn2 <- ovarian27k$genename[wp2[oo2], c(2,3)]
data.frame(gn2, fdr=fdr[wp2[oo2]])

##### p.10
#library(qvalue)

## False discover rate control using qvalue
qval <- qvalue(pval)$qvalues
wp3 <- which(qval < 0.05)
oo3 <- order(qval[wp3])
length(wp3)
gn3 <- ovarian27k$genename[wp3[oo3], c(2,3)]
data.frame(gn3, q.values=qval[wp3[oo3]])

## Manhattan plot
cc <- rainbow(length(gsize)); ss <- rep(1:length(gsize), gsize)
plot(-log10(pval), type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab=expression(paste("-log"[10],"(P)")),
     main="Independent two sample T-test")
abline(h=-log10(max(pval[wp1])), col="red", lty=2, lwd=2 )
abline(h=-log10(max(pval[wp2])), col="darkblue", lty=3, lwd=2)
abline(h=-log10(max(pval[wp3])), col="darkgreen", lty=4, lwd=2)
legend("topleft", c("5% Bonferroni","5% BH","5% qvalue"), cex=1.5,
       lty=c(2,3,4), col=c("red", "darkblue", "darkgreen"), lwd=2)

##### p.11
set.seed(1111)

## Fix alpha=0.1 and perform 5-fold cross-validation
alpha=0.1
g1 <- cv.glmnet(x, y, alpha=alpha, family="binomial", nfolds=5)
plot(g1)

## Re-fit the elastic net with two selected lambda values
fit11 <- glmnet(x, y, alpha=alpha, family="binomial",
                lambda=g1$lambda.min)
fit12 <- glmnet(x, y, alpha=alpha, family="binomial",
                lambda=g1$lambda.1se)

## The number of selected CpG sites.
c(fit11$df, fit12$df)

## Location of selected CpG sites
w11 <- which(as.matrix(fit11$beta)!= 0)
w12 <- which(as.matrix(fit12$beta)!= 0)

##### p.12
## List genes from the largest absolute value of beta
o11 <- order(abs(fit11$beta[w11]), decreasing=TRUE)
o12 <- order(abs(fit12$beta[w12]), decreasing=TRUE)

## Selected CpG sites and their corresponding genes
gn11 <- ovarian27k$genename[w11[o11], c(2,3)]
gn12 <- ovarian27k$genename[w12[o12], c(2,3)]
data.frame(gn11, beta=fit11$beta[w11[o11]])
data.frame(gn12, beta=fit12$beta[w12[o12]])

## Coefficient plot over entire CpG sites
par(mar=c(5,5,4,3))
plot(abs(fit11$beta), type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab=expression(paste("|",hat(beta),"|")),
     main="Elastic Net using lambda.min")
plot(abs(fit12$beta), type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab=expression(paste("|",hat(beta),"|")),
     main="Elastic Net using lambda.1se")
     
##### p.13
## Generate group index
grp <- rep(1:length(gsize), gsize)

## Perform cross-validation for group lasso
set.seed(1111)
ry <- y
ry[y==0] <- -1
g2 <- cv.gglasso(x, ry, group=grp, loss="logit", nfolds=5)
plot(g2)

## Re-fit the elastic net with two selected lambda values
fit21 <- gglasso(x, ry, group=grp, loss="logit", lambda=g2$lambda.min)
fit22 <- gglasso(x, ry, group=grp, loss="logit", lambda=g2$lambda.1se)

## The number of selected CpG sites
cbind(fit21$df, fit22$df)

## Location of selected CpG sites
w21 <- which(fit21$beta!=0)
w22 <- which(fit22$beta!=0)     

##### p.14
## List genes from the largest absolute value of beta
o21 <- order(abs(fit21$beta[w21]), decreasing=TRUE)
o22 <- order(abs(fit22$beta[w22]), decreasing=TRUE)

## Selected CpG sites and their corresponding genes
gn21 <- ovarian27k$genename[w21[o21], c(2,3)]
gn22 <- ovarian27k$genename[w22[o22], c(2,3)]
data.frame(gn21, beta=fit21$beta[w21[o21]])
data.frame(gn22, beta=fit22$beta[w22[o22]])

## Coefficient plot over entire CpG sites
par(mar=c(5,5,4,3))
plot(abs(fit21$beta), type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab=expression(paste("|",hat(beta),"|")),
     main="Group Lasso using lambda.min")
plot(abs(fit22$beta), type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab=expression(paste("|",hat(beta),"|")),
     main="Group Lasso using lambda.1se")     

##### p.15
set.seed(1212)

## Compute selection probabilities
g3 <- sel.pclogit(x, y, alpha=0.1, group=gsize, type="ring")
g4 <- sel.pclogit(x, y, alpha=0.1, group=gsize, type="fcon")

## Display top 50 genes ranked by selection probabilities
top <- 50
gn31 <- ovarian27k$genename[g3$maxsel[1:top,1], c(2,3)]
gn41 <- ovarian27k$genename[g4$maxsel[1:top,1], c(2,3)]
data.frame(gn31, SP=g3$maxsel[1:top,2])
data.frame(gn41, SP=g4$maxsel[1:top,2])

## Selection probability plot over entire CpG sites
plot(g3$beta[,5], type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab="Selection probability",
     main="Network-based regularization using a ring network")
plot(g4$beta[,5], type="p", pch=20, col=cc[ss], xlab="CpG sites",
     ylab="Selection probability",
     main="Network-based regularization using a f.con network")
     
##### p.16
## 12 CpG sites from T-test using qvalue
wh1 <- sort(wp3)
length(wh1)

## 46 CpG sites from Elastic net (E.net) using lambda.1se
wh2 <- sort(w12)
length(wh2)

## 54 CpG sites from Group lasso (Glasso) using lambda.min
wh3 <- sort(w21)
length(wh3)

## Top 50 CpG sites from a ring network
wh4 <- sort(g3$maxsel[1:top, 1])
length(wh4)

## Top 50 CpG sites from a f.con network
wh5 <- sort(g4$maxsel[1:top, 1])
length(wh5)

##### p.17
#library(VennDiagram)

## Venn Diagram of 5 methods
vd1 <- venn.diagram(list(T.test=wh1, E.net=wh2, GLasso=wh3,
                    Ring=wh4,F.con=wh5), fill=2:6, alpha=0.25,
                    filename=NULL,cat.cex=1.5, cex=1.5, margin=0.1)
x11(100, 100)
grid.draw(vd1)

### Genes detected by all of 5 methods
all <- unique(c(wh1, wh2, wh3, wh4, wh5))
mat <- matrix(0, length(all), 5)
for (i in 1:5) {
    mat[,i] <- as.numeric(all %in% get(paste0("wh",i)))
}
colnames(mat) <- c("T.test", "E.net", "GLasso", "Ring", "F.con")
rownames(mat) <- ovarian27k$genename[all,3]
wt <- apply(mat, 1, sum)
names(which(wt==5))     

##### p.18
## Venn Diagram of 4 methods without T.test
vd2 <- venn.diagram(list(E.net=wh2, GLasso=wh3,
                    Ring=wh4,F.con=wh5), fill=2:5, alpha=0.25,
                    filename=NULL,cat.cex=1.5, cex=1.5, margin=0.1)
x11(100, 100)
grid.draw(vd2)

## Genes detected by all of 4 regularization methods
wt2 <- apply(mat[,-1], 1, sum)
names(which(wt2==4))

## Genes detected by at least 3 regularization methods
names(which(wt2>=3))
