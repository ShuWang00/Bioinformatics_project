# Version info: R 2.14.1, Biobase 2.15.3, GEOquery 2.23.2, limma 3.10.1
# R scripts generated  Wed Apr 2 03:03:43 EDT 2014

################################################################
#   Differential expression analysis with limma
library(Biobase)
source("http://bioconductor.org/biocLite.R")
library(biocLite)
library(GEOquery)
library(limma)

biocLite("rgl")

## move to a project folder and create a result folder in it 
basedir <- "/Users/wangshu/Downloads/GSE37648_RAW"
setwd(basedir)

resfolder <- "GSE37648.results"
if (! file.exists(resfolder) ){
  dir.create(resfolder, showWarnings = FALSE, recursive = FALSE, mode = "0777")
  Sys.chmod(resfolder, mode = "0777", use_umask=TRUE)
}

#load (once only) series and platform data from GEO
gset <- getGEO("GSE37648", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1

#Arrange the (normalized) data and apply log2 transformation
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples (Simvastatin excluded)
sml <- c("G0","G0","G0","G1","G1","G1","X","X","X","X","X","X");

# eliminate samples marked as "X" (Simvastatin not used here)
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#set up the DE contrasts and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)

# the differential expression contrast is set here to "LPS" versus "Control"
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))


#Save the resulting table to file for reuse
filename <- "GSE37648.data.tsv"
outfile <- paste(basedir, resfolder, filename, sep="/")
write.table(tT, file=outfile, 
            quote = FALSE, dec=",", sep="\t", col.names = NA, row.names = T)


#Principal Component Analysis (PCA) of the 6 samples
# copy the ex object created above and transpose it to get rows as samples and columns as probes
pca.data <- ex

# rename samples
grps <- c(rep("Cont",3), rep("LPS",3))
grpcol <- c(rep("blue",3), rep("red",3))
colnames(pca.data) <- paste(grps, colnames(pca.data), sep="-")

# remove NAs
pca.data <- na.omit(as.matrix(pca.data))

# transpose
pca.data <- t(pca.data)

# inspect pca.data
#dim(pca.data)
#str(pca.data)

# preview first 5 probe-sets
pca.data[,1:5]

# compute PCA
pca <- prcomp(pca.data, scale=TRUE)

# identify variance in components
summary(pca)

# the first 2 component group 64% of the total variance
# the first 3 component group 79% of the total variance
# the first 4 component group 90% of the total variance

# components #1 and #2
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca.data), cex=0.75)

# show other component pairs for the example
# components #1 and #3
plot(pca$x[,1], pca$x[,3], xlab="PCA1", ylab="PCA3", main="PCA for components 1&3", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,3], rownames(pca.data), cex=0.75)

# components #2 and #3
plot(pca$x[,2], pca$x[,3], xlab="PCA2", ylab="PCA3", main="PCA for components 2&3", type="p", pch=10, col=grpcol)
text(pca$x[,2], pca$x[,3], rownames(pca.data), cex=0.75)


#heatmap
source("http://bioconductor.org/biocLite.R")
biocLite("rafalib")
library(oligo)
storedFiles <- basename(as.character(pData(gset)$supplementary_file))

smpls <- as.character(pData(gset)$title)
cells <- factor(sub("2_(\\w+)_[ab]", "\\1", smpls))
replicate <- factor(sub("2_\\w+_([ab])", "\\1", smpls))
# simple heatmap
library(RColorBrewer)
heatmap(cor(exprs(gset)), col = rev(brewer.pal(11, "RdBu")), labCol = smpls, scale = "none", 
        margins = c(8, 8), ColSideColors = c("darkgreen", "orange", "darkviolet")[cells])

#Cluster genes using hierarchical clustering and k-means

d <- dist( t(ex) )
library(rafalib)
mypar()
hc <- hclust(d)
plot(hc,labels=colnames(ex),cex=0.5)

set.seed(1)
km <- kmeans(t(ex), centers=2)
names(km)

plot(ex, col=km$cluster)

data(iris)
# this is a little tweak so that things line up nicely later on
iris$Species <- factor(iris$Species,
                       levels = c("versicolor","virginica","setosa"))
head(iris)

round(cor(iris[,1:4]), 2)

pc <- princomp(iris[,1:4], cor=TRUE, scores=TRUE)

summary(pc)

plot(pc,type="lines")

biplot(pc)


library(rgl)
plot3d(pc$scores[,1:3], col=iris$Species)
