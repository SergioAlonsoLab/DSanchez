# data analysis for David Sanchez
# the data has been previously extensively analyzed by VHIO genomics platform

setwd("/imppc/labs/mplab/salonso/DSanchez/")

library(umap)
library(magrittr)
library(gplots)

data0 <- read.delim("data/results/Normalized.Filtered.csv",sep=";",nrows = 138750,dec = ",")
rownames(data0) <- data0$X
data0$X <- NULL

m <- regexpr("_.+_",colnames(data0))
patient <- regmatches(colnames(data0),m) %>% gsub("_","",.)

clinical0 <- gdata::read.xls("data/Macrophages clinical data.xlsx",nrows=22)
rownames(clinical0) <- clinical0$ID

genes0 <- read.csv("data/results/Annotations.AllGenes.csv",sep="\t",nrows = 138750)
table(genes0$Probe %in% row.names(data0)) 

# select only probes with gene information

inGenes <- genes0$Probe[which(genes0$Symbol!="")] 
inGenes <- intersect(inGenes,rownames(data0))
rownames(genes0) <- genes0$Probe

data1 <- data0[inGenes,] %>% limma::normalizeBetweenArrays()
rownames(data1) <- genes0[inGenes,"Symbol"]

boxplot(data1,las=2,pch=19)

type1 <- factor(substr(colnames(data0),1,3))
type2 <- factor(substr(colnames(data0),5,7))

colSums(is.na(data1))

locPalette <- c("lightblue3","darkblue")
disPalette <- c("green3","orange3")

par(mfrow=c(2,2),mar=c(5,5,2,2))

pca1 <- prcomp(t(data1),scale. = T)
plot(pca1$x[,1:2],col=locPalette[type1],pch=19,main=c("PCA - tissue location"))
plot(pca1$x[,1:2],col=disPalette[type2],pch=19,main=c("PCA - disease status"))

summary(glm(type1 ~ .,as.data.frame(pca1$x[,c(1,2,4)]),family=binomial))
plot(pca1$x[,c(1,4)],col=locPalette[type1],pch=19,main=c("PCA - tissue location"))

summary(glm(type1 ~ pca1$x[,1] + pca1$x[,2],family="binomial"))
summary(lm(pca1$x[,1] ~ type1))
summary(lm(pca1$x[,2] ~ type1))

summary(lm(pca1$x[,1] ~ type2))
summary(lm(pca1$x[,2] ~ type2))

# There are no significant differences in global expression according to disease vs not diseased


umap1 <- umap(scale(t(data1)))
plot(umap1$layout,col=locPalette[type1],pch=19,xlab="UMAP1",ylab="UMAP2",main="UMAP - tissue location")
identify(umap1$layout[,1],umap1$layout[,2],rownames(umap1$layout))
points(umap1$layout["SUB.UNH_550174_1",],col="red")

plot(umap1$layout,col=disPalette[type2],pch=19,xlab="UMAP1",ylab="UMAP2",main="UMAP - disease status")

par(mfrow=c(1,2),mar=c(12,2,1,1))
boxplot(data1,col=locPalette[type1],las=2)
boxplot(data1,col=disPalette[type2],las=2)



# Lineal models with and without interactions

lm0 <- apply(data1,1,function(x) lm(x ~ type1 + type2))
lm1 <- apply(data1,1,function(x) lm(x ~ type1 * type2))

sapply(lm1,function(x) {
  summary(x) %>% coef %>% .[4,4]
}) -> interaction.pval

table(interaction.pval < 0.01)

candidates <- which(interaction.pval < 0.05)

sapply(names(candidates),function(x) {
  
  boxplot(data1[x,] ~ type2 + type1,main=x,ylab="Expression (log)",xlab=NULL)
  means <- tapply(data1[x,],list(type2,type1),mean)
  segments(c(1,3),means[c(1,3)],c(2,4),means[c(2,4)])
  points(1:4,means,pch=21,cex=1.5,bg=locPalette[c(1,1,2,2)])
})





as.data.frame((sort(interaction.pval)))




var1 <- apply(data1,1,var) 
relvar1 <- apply(data1,1,function(x) var(x)/mean(x))


selected <- order(relvar1,decreasing = T)[1:500]

pca2 <- prcomp(t(data1[selected,]),scale. = T)
plot(pca2$x[,1:2],pch=19,col=type1)
plot(pca2$x[,1:2],pch=19,col=type2)

umap2 <- umap(t(data1[selected,]))
plot(umap2$layout,pch=19,col=type1,xaxt="n",yaxt="n",xlab="UMAP1",ylab="UMAP2")
title("UMAP projection using 500 most variable genes")
identify(umap2$layout,labels = colnames(data0))

sapply(patient,function(x) {
  
  columns <- grep(x,colnames(data1))
  if(length(columns)==2) {
    dpatient <- data1[,columns]
    plot(dpatient)
  }
  
})

par(mar=c(10,5,3,1))

data1 %>% t %>% scale %>% dist %>% hclust %>% as.dendrogram %>% plot
data1 %>% t %>% dist %>% hclust %>% plot


# clustering

hclust(dist(t(scale(data1[selected,])))) %>% plot

heatmap.2(data1)



# find genes most associated with glucose levels

myf <- function(x) lm(x ~ type1 * type2)     

lmodels <- apply(data1,1,myf) 

names(lmodels) <- genes1$Symbol

interactions <- sapply(lmodels,function(x) coef(summary(x))[4,c(1,4)])

interactions <- t(interactions) %>% as.data.frame()

plot(interactions$Estimate,-log10(interactions$`Pr(>|t|)`))
abline(h=2)

identify(interactions$Estimate,-log10(interactions$`Pr(>|t|)`),rownames(data1))

plot(interactions$Estimate,-log10(p.adjust(interactions$`Pr(>|t|)`,method="fdr")))
p.a

glucose <- sapply(lmodels,function(x) coef(summary(x))[3,c(1,4)])


plot(t(interactions),log="y")
plot(t(glucose),log="y")
qqplot(runif(1000),interactions[2,],log="xy")
abline(0,1)
which.min(interactions[2,])
which.min(glucose[2,])

plot(unlist(data1[5768,]) ~ clinical0[patient,"Glu"],col=type1)


plot(clinical0$IMC,clinical0$Glu)
cor.test(clinical0$IMC,clinical0$Diabetes,method = "spear")

plot(clinical0$Glu ~ clinical0$Diabetes)
