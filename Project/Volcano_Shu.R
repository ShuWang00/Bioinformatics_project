# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Jun 2 16:17:25 EDT 2017

################################################################
#   Differential expression analysis with limma
source("http://bioconductor.org/biocLite.R") 

library(Biobase)
library(GEOquery)
library(limma)
# load series and platform data from GEO

gset <- getGEO("GSE37648", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "000000111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
ex <- exprs(gset)
patient1 <- ex[ ,1:6]
patient2 <- ex[ ,7:12]

View(ex)
View(patient2)
boxplot(patient1)
boxplot(patient2)


# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

View(exprs(gset))

x<-exprs(gset)

View(x)

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)

options(download.file.method = "wget")
install.packages("lubridate", dependencies=TRUE, repos='http://cran.rstudio.com/')
options(install.packages.check.source = "no")
setRepositories(ind = c(1:6, 8))
ap <- available.packages()
biocLite("lmFit")

fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=54675)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
res <- write.table(tT, file=stdout(), row.names=F, sep="\t")

View(tT)

plot(tT$logFC, 1-tT$adj.P.Val, xlim=c(-6, 6), 
     main="Effect of intratracheal exposure to LPS",
     xlab="log2Ratio", ylab="1-adj.P.Val")
abline(h=0.95, col="red")



# Download the data from github (click the "raw" button, save as a text file called "results.txt").
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
head(tT)

# Make a basic volcano plot
with(tT, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2)))


# Add colored points: red if adj.P.Val<0.05, orange of log2FC>1, green if both)
with(subset(tT, -log10(P.Value)>2 & logFC>1 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(tT, -log10(P.Value)>2 & logFC< -1), points(logFC, -log10(P.Value), pch=20, col="green"))
with(subset(tT, -log10(P.Value)>2 & logFC> -1&logFC< 1), points(logFC, -log10(P.Value), pch=20, col="gray"))
abline(a = 2, b = 0, h = NULL, v = NULL, reg = NULL,
       coef = NULL, untf = FALSE)
abline(a = 0, b = 0, h = NULL, v = 1, reg = NULL,
       coef = NULL, untf = FALSE)
abline(a = 0, b = 0, h = NULL, v = -1, reg = NULL,
       coef = NULL, untf = FALSE)
grey.count = nrow(subset(tT, P.Value<.01 & abs(logFC)<1))
grey.count = paste(c("N =", grey.count), collapse = " ")
red.count = nrow(subset(tT, P.Value<.01 & logFC>1))
red.count = paste(c("N =", red.count), collapse = " ")
green.count = nrow(subset(tT, P.Value<.01 & logFC<(-1)))
green.count = paste(c("N =", green.count), collapse = " ")
left.black.count = nrow(subset(tT, P.Value>=.01 & logFC<(-1)))
left.black.count = paste(c("N =", left.black.count), collapse = " ")
right.black.count = nrow(subset(tT, P.Value>=.01 & logFC>1))
right.black.count = paste(c("N =", right.black.count), collapse = " ")

legend(x = "topleft", legend = green.count, bty = "n", text.col = "red")
legend(x = "top", legend = grey.count, bty = "n", text.col = "red")
legend(x = "topright", legend = red.count, bty = "n", text.col = "red")
legend(x = "bottomleft", legend = left.black.count, bty = "n", text.col = "red")
legend(x = "bottomright", legend = right.black.count, bty = "n", text.col = "red")

legend("topleft", legend="N = 466", cex=0.8,box.lty = 0)
legend("top", legend="N = 467", cex=0.8,box.lty = 0)
legend("topright", legend="N = 469", cex=0.8,box.lty = 0)
legend("bottomleft", legend="N = 48", cex=0.8,box.lty = 0)
legend("bottomright", legend="N = 134", cex=0.8,box.lty = 0)

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE37648", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "000000111111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Patient1","Patient2")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE37648", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

save(title,fl,ex,labels, file = 'file.Rdata')
save()

patient1 <- x[ ,1:6]

patient2 <- x[ ,7:12]

View(ex)
boxplot(patient1)
boxplot(patient2)                               
