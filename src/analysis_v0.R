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

colnames(data1) %>% gsub(".","_",.,fixed=T) %>% strsplit(.,"_") %>% 
  unlist %>% matrix(.,byrow=T,ncol=4) %>% as.data.frame -> sampleData

colnames(sampleData) <- c("location","status","patient","batch")

boxplot(data1,las=2,pch=19) # check the normalization

location <- factor(sampleData$location)
status <- factor(sampleData$status)

sum(is.na(data1)) # there are no NAs

# color palettes

locPalette <- c("lightblue3","darkblue")
disPalette <- c("green3","red")

# PCA analysis using all genes

par(mfrow=c(2,2),mar=c(5,5,2,2))

pca1 <- prcomp(t(data1),scale. = T)
plot(pca1$x[,1:2],col=locPalette[location],pch=19,main=c("PCA - tissue location"))
plot(pca1$x[,1:2],col=disPalette[status],pch=19,main=c("PCA - disease status"))
summary(glm(location ~ .,as.data.frame(pca1$x[,c(1,2)]),family=binomial))
summary(lm(pca1$x[,1] ~ location))
summary(lm(pca1$x[,2] ~ location))

summary(glm(status ~ .,as.data.frame(pca1$x[,c(1,2)]),family=binomial))
summary(lm(pca1$x[,1] ~ status))
summary(lm(pca1$x[,2] ~ status))


# There are no significant differences in global expression according to disease vs not diseased


umap1 <- umap(scale(t(data1)))
plot(umap1$layout,col=locPalette[location],pch=19,xlab="UMAP1",ylab="UMAP2",main="UMAP - tissue location")
identify(umap1$layout[,1],umap1$layout[,2],rownames(umap1$layout))

plot(umap1$layout,col=disPalette[status],pch=19,xlab="UMAP1",ylab="UMAP2",main="UMAP - disease status")

par(mfrow=c(1,2),mar=c(12,4,3,1))
boxplot(data1,col=locPalette[location],las=2)
boxplot(data1,col=disPalette[status],las=2)



# Lineal models with and without interactions

# Healthy vs Unhealthy
# In ALL
lm0 <- apply(data1,1,function(x) lm(x ~ status))

# In SUB
lm1 <- apply(data1[,location=="SUB"],1,function(x) lm(x ~ status[location=="SUB"]))

# In VIS
lm2 <- apply(data1[,location=="VIS"],1,function(x) lm(x ~ status[location=="VIS"]))

# 2-factor models
# Without interactions
lm3 <- apply(data1,1,function(x) lm(x ~ status + location))

# With interactions
lm4 <- apply(data1,1,function(x) lm(x ~ status * location))



# Get coefficients, t-statistic, or P-values from a list of linear models

getSummaryValues <- function(x,what=c("Estimate","Std.Error","t value","pValue")) {
  what <- match.arg(what)
  if(what=="pValue") what <- "Pr(>|t|)"
  sapply(x,function(x) summary(x) %>% coef %>% .[,what]) %>% t %>% as.data.frame
}


# generate a table with pvalues


pValues <- data.frame(Gene=names(lm0),
                       P_global=getSummaryValues(lm0,"p")[,2],
                       P_sub=getSummaryValues(lm1,"p")[,2],
                       P_vis=getSummaryValues(lm2,"p")[,2],
                       P_status=getSummaryValues(lm3,"p")[,2],
                       P_location=getSummaryValues(lm3,"p")[,3],
                       P_interaction=getSummaryValues(lm4,"p")[,4])

pValues[,2:7] %>% signif(.,6) -> pValues[,2:7]
write.csv(pValues,file = "output/pvalues.csv",quote = F,row.names = F)


# generate a table with t-statistic values

tValues <- data.frame(Gene=names(lm0),
                      t_global=getSummaryValues(lm0,"t")[,2],
                      t_sub=getSummaryValues(lm1,"t")[,2],
                      t_vis=getSummaryValues(lm2,"t")[,2],
                      t_status=getSummaryValues(lm3,"t")[,2],
                      t_location=getSummaryValues(lm3,"t")[,3],
                      t_interaction=getSummaryValues(lm4,"t")[,4])


tValues[,2:7] %>% signif(.,6) -> tValues[,2:7]
write.csv(tValues,file = "output/tvalues.csv",quote = F,row.names = F)

# boxplot with location and status ----

myBoxPlot <- function(gene) {
  yrange <- range(data1[gene,])
  prange <- yrange + c(0,diff(yrange)/2)
  boxplot(data1[gene,] ~ status + location,las=1,ylab="Expression (log10)",pch=19,cex=.5,col=disPalette,ylim=prange,xaxt="n",xlab=NULL)
  axis(side = 1,1:4,labels=paste(c("HEA","UNH","HEA","UNH"),c("SUB","SUB","VIS","VIS"),sep="\n"),padj = .5)
  means <- tapply(data1[gene,],list(status,location),mean,na.rm=T)
  segments(c(1,3),means[c(1,3)],c(2,4),means[c(2,4)],lwd=2)
  points(1:4,means,pch=23,bg="yellow",cex=2)
  segments(2.5,prange[1]/2,2.5,yrange[1]+diff(yrange)*1.1)
  segments(0,yrange[1]+diff(yrange)*1.1,5,yrange[1]+diff(yrange)*1.1)
  title(ifelse(is.numeric(gene),rownames(data1)[gene],gene))
  text(1.5,yrange[1] + diff(yrange)*1.05,sprintf("P=%1.2g",pValues[gene,"P_sub"]))
  text(3.5,yrange[1] + diff(yrange)*1.05,sprintf("P=%1.2g",pValues[gene,"P_vis"]))
  text(2.5,yrange[1] + diff(yrange)*1.3,sprintf("Pglobal=%1.2g\nPstatus=%1.2g\nPlocation=%1.2g\nPinter=%1.2g",
                                               pValues[gene,"P_all"],
                                               pValues[gene,"P_status"],
                                               pValues[gene,"P_location"],
                                               pValues[gene,"P_interaction"]))
  
}



# generate some graphs with different types of genes ----


pdf("output/status.pdf",8,12) 
par(mfrow=c(4,3),mar=c(4,4,2,1))
for(i in order(pValues$P_status)[1:24]) myBoxPlot(i)
dev.off()

pdf("output/subcutaneous.pdf",8,12) 
par(mfrow=c(4,3),mar=c(4,4,2,1))
for(i in order(pValues$P_sub)[1:24]) myBoxPlot(i)
dev.off()

pdf("output/visceral.pdf",8,12) 
par(mfrow=c(4,3),mar=c(4,4,2,1))
for(i in order(pValues$P_vis)[1:24]) myBoxPlot(i)
dev.off()

pdf("output/interactions.pdf",8,12) 
par(mfrow=c(4,3),mar=c(4,4,2,1))
for(i in order(pValues$P_interaction)[1:24]) myBoxPlot(i)
dev.off()

# scraps ----


lm4.pvals[order(lm4.pvals$`statusUNH:locationVIS`),][1:10,]


table(lm1.pvals$statusUNH < 0.05,lm2.pvals$statusUNH < 0.05)

plot(lm1.pvals[,2],lm2.pvals[,2],log="xy")

lm0 <- apply(data1,1,function(x) lm(x ~ type1 + type2))
lm1 <- apply(data1,1,function(x) lm(x ~ type1 * type2))

sapply(lm0,function(x) {
  summary(x) %>% coef %>% .[,4]
}) %>% t %>% as.data.frame -> no.interaction.pvals


sapply(lm1,function(x) {
  summary(x) %>% coef %>% .[,4]
}) %>% t %>% as.data.frame -> interaction.pvals

plot(no.interaction.pvals$type1VIS,interaction.pvals$type1VIS,log="xy",pch=".")

candidates <- which(interaction.pvals$`type1VIS:type2UNH` < 0.05)

table(p.adjust(interaction.pvals$`type1VIS:type2UNH`) < 0.05)


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

