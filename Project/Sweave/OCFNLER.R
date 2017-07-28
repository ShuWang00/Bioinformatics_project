#!/usr/bin/Rscript --vanilla
# OCFNLER current manuscript analysis
# Feb 10 2011
# Turned into Sweave document for Nat Med submission

library(limma)
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(hgu133plus2probe)
library(hgu133a.db)
library(hgu133acdf)
library(hgu133aprobe)
library(gplots)

#New libraries added
library(org.Hs.eg.db)
library(annotate)
library(survdiff)
#Unable to install
#source("http://bioconductor.org/biocLite.R")
#biocLite("survdiff")
#install.packages('matchprobes', lib='/Users/wangshu/Downloads')
#library(matchprobes)

library(ClassDiscovery)
library(affy)
library(statmod)
library(biomaRt)
library(globaltest)
#library(Design) - obsolete, replaced by rms
library(rms)
library(mclust)
library(survival)
library(survdiff)

#***************************************************************************
# Begin by loading cell line data

load("/Users/wangshu/Documents/Sweave/Data/cellLines_expressionSet.rdat")
targets.no52 <- readTargets("/Users/wangshu/Documents/Sweave/Data/Targets_no52.txt")
target.names.no52 <- paste(targets.no52[, 1], "CEL", "gz", sep=".")
eset.no52 <- exprs(cells.normalized)[, target.names.no52]
View(eset.no52)
#***************************************************************************
# IDENTIFICATION OF DIFFERENTIALLY EXPRESSED GENES USING LIMMA
# To demonstrate how I created the genelists (e.g. lm.oce.fne.1157.Feb11) above
origintrt.no52 = paste(targets.no52$FO, targets.no52$TN, sep=".")
# dm2b (limma design matrix 2b), as opposed to many other options I tried,
# cell origin comparisons (ft vs ov, within immortalized & transformed cells), no pt 52
# Aim: apply duplicate correlation function to block for patient differences
# need statmod package for correlation computation

fitdm2b.no52 <- factor(origintrt.no52, levels=c("ft.n","ov.n","ft.t","ov.t"))
designdm2b.no52 <- model.matrix(~0+fitdm2b.no52)
colnames(designdm2b.no52) <- c("ft.n","ov.n","ft.t","ov.t")

# create a vector ptdiff which indicates the blocks that are different patients
ptdiff <- c(targets.no52$patient)
# use the duplicateCorrelation function to identify the correlation btw samples derived from
# the same patient
ptdiffcorr <- duplicateCorrelation(eset.no52, design=designdm2b.no52, block=ptdiff)
#View(ptdiffcorr)
# The block function estimates the correlation between repeated observations on the 
# blocking variable. Typically the blocks are biological replicates and the repeated observations are 
# technical replicates. In either case, the correlation is estimated by fitting a mixed linear model by 
# REML individually for each gene. The function also returns a consensus correlation, which is a robust 
# average of the individual correlations, which can be used as input for functions lmFit or gls.series. 

# Value: a list with components 
# consensus.correlation the average estimated inter-duplicate correlation. The average is the trimmed 
# mean of the individual correlations on the atanh-transformed scale. 

# the dupcorfit can then be applied to the linear model as one of the coefficients
# ptdiffcorr$consensus.correlation estimates the average correlation within the blocks

fit.ptdiffcorr <- lmFit(eset.no52, designdm2b.no52, block=ptdiff, cor=ptdiffcorr$consensus.correlation)

# set up contrast matrix
matrix.dm2b.no52.dupcorr <- makeContrasts (ft.n-ov.n, ft.t-ov.t, levels=designdm2b.no52)

# apply the fit to extract the contrasts of interest
fit2dm2b.no52.dupcorr <- contrasts.fit(fit.ptdiffcorr, matrix.dm2b.no52.dupcorr)
fit2dm2b.no52.dupcorr <- eBayes(fit2dm2b.no52.dupcorr)

# create a results file
resultsdm2b.no52.dupcorr <- decideTests(fit2dm2b.no52.dupcorr,method = "separate",adjust.method = "BH",p.value=.05)
#View(resultsdm2b.no52.dupcorr)
# write the results to a file
write.fit(fit2dm2b.no52.dupcorr, results=resultsdm2b.no52.dupcorr, 
          file="/Users/wangshu/Documents/Sweave/Results/dm2b.vsn.no52_dupcorrtest_2011-06-30.txt", 
          digits=3, adjust="BH", method="sep", row.names = TRUE)
View(fit2dm2b.no52.dupcorr)
View(resultsdm2b.no52.dupcorr)
# read in again to annotate
final.result.dm2b.vsn.no52.dupcorr.test <- read.delim(file="/Users/wangshu/Documents/Sweave/Results/dm2b.vsn.no52_dupcorrtest_2011-06-30.txt",as.is=TRUE)
View(final.result.dm2b.vsn.no52.dupcorr.test)

##Copy and Paste the probe column to result matrix
x<- row.names(fit2dm2b.no52.dupcorr)

length<-dim(final.result.dm2b.vsn.no52.dupcorr.test)[2]

final.result.dm2b.vsn.no52.dupcorr.test[,length+1]<-x

names(final.result.dm2b.vsn.no52.dupcorr.test)[length+1]<-"ID"

View(x)


# Focus on OCE-FNE genelist
# Select only those probes that are significant
# > colnames(final.result.dm2b.vsn.no52.dupcorr.test)
# [1] "A"                       "Coef.ft.n...ov.n"
# [3] "Coef.ft.t...ov.t"        "t.ft.n...ov.n"
# [5] "t.ft.t...ov.t"           "p.value.ft.n...ov.n"
# [7] "p.value.ft.t...ov.t"     "p.value.adj.ft.n...ov.n"
# [9] "p.value.adj.ft.t...ov.t" "F"
# [11] "F.p.value"               "Res.ft.n...ov.n"
# [13] "Res.ft.t...ov.t"         "ID"

oce.fne.1157.dupcorr <- final.result.dm2b.vsn.no52.dupcorr.test[final.result.dm2b.vsn.no52.dupcorr.test[,12]!=0,]
# this method identifies more significant (n=1157) probes
gene.symbols.dupcorr <- unlist(mget(oce.fne.1157.dupcorr$ID, hgu133plus2SYMBOL, ifnotfound=NA))
final.result.oce.fne.1157.dupcorr <- cbind (oce.fne.1157.dupcorr,gene.symbols=gene.symbols.dupcorr)

write.table(final.result.oce.fne.1157.dupcorr,file="/Users/wangshu/Documents/Sweave/Results/final.result.oce.fne.1157_lm.fdr.P.05.sep.dupcorr_annot_2011.06.30.txt", row.names=F, sep="\t")

lm.oce.fne.1157.Feb11 <- read.delim(file="/Users/wangshu/Documents/Sweave/Results/final.result.oce.fne.1157_lm.fdr.P.05.sep.dupcorr_annot_2011.06.30.txt",as.is=TRUE)


#***************************************************************************
#***************************************************************************
# Begin Matt's analysis
# Analysis of Melissa's FTE/OSE genes using the 10-gene model (5 OCE genes and 5 FNE genes)

# Function to get glm logistic p-value
get.logistic.p <- function(model) {
  pval <- pchisq(model$null.deviance - model$deviance, 
                 model$df.null - model$df.residual, 
                 lower.tail = FALSE)
  return(pval)
}

lmp <- function(lm.object) {
  numdf <- as.numeric(summary(lm.object)$fstatistic["numdf"])
  dendf <- as.numeric(summary(lm.object)$fstatistic["dendf"])
  val <- as.numeric(summary(lm.object)$fstatistic["value"])
  pval <- pf(val, numdf, dendf, lower.tail = FALSE)
  return(pval)
}

# Function to find p-value of survdiff
# surv.diff is a survdiff object

p.value <- function(surv.diff) {
  df <- length(surv.diff$obs) - 1
  chisq <- surv.diff$chisq
  pchisq(chisq, df, lower.tail = FALSE)
}

mel.tab <- read.delim(file="/Users/wangshu/Documents/Sweave/Results/final.result.oce.fne.1157_lm.fdr.P.05.sep.dupcorr_annot_2011.06.30.txt",as.is=TRUE)

nbig <- 10
n <- 5

# Separate into OSE and FT genes
# Separate into UR and DR genes first
ose <- mel.tab[mel.tab$Coef.ft.n...ov.n < 0,] #down-regulated genes
ft <- mel.tab[mel.tab$Coef.ft.n...ov.n > 0,] #up-regulated genes

# > dim(ose)
# [1] 525  15
# > dim(ft)
# [1] 632  15

# Order the FT and OSE genes according to fdr adjusted p.value
ft.adj <- ft[order(ft$p.value.adj.ft.n...ov.n),][1:nbig,]
ose.adj <- ose[order(ose$p.value.adj.ft.n...ov.n),][1:nbig,]

# Get HGNC symbols for genes
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# > human
# Object of class 'Mart':
# Using the ensembl BioMart database
# Using the hsapiens_gene_ensembl dataset

ft.both <- getBM(c("hgnc_symbol","affy_hg_u133_plus_2"), "affy_hg_u133_plus_2", ft.adj$ID, human)
## getBM function: This function is the main biomaRt query function.  Given a set of
# filters and corresponding values, it retrieves the user specified
# attributes from the BioMart database one is connected to

# > ft.both
# hgnc_symbol affy_hg_u133_plus_2
# 1      HS6ST3          1552842_at
# 2      OSBPL3         209627_s_at
# 3        DPP6           228546_at
# 4        CD47           227259_at
# 5     NOSTRIN           226992_at
# 6        RAI2           219440_at
# 7      KCNK13           221325_at
# 8        DOK5        1554863_s_at
# 9        DOK5         214844_s_at

rownames(ft.both) <- ft.both$affy_hg_u133_plus_2
# > ft.both
# hgnc_symbol affy_hg_u133_plus_2
# 1552842_at        HS6ST3          1552842_at
# 209627_s_at       OSBPL3         209627_s_at
# 228546_at           DPP6           228546_at
# 227259_at           CD47           227259_at
# 226992_at        NOSTRIN           226992_at
# 219440_at           RAI2           219440_at
# 221325_at         KCNK13           221325_at
# 1554863_s_at        DOK5        1554863_s_at
# 214844_s_at         DOK5         214844_s_at

ft.both <- na.omit(ft.both[ft.adj$ID,])
## na.omit function: Description:
# These generic functions are useful for dealing with NAs in e.g.,
# data frames. na.omit
# returns the object with incomplete cases removed.
# reminder ft.adj is the top10 FDR adj p.values

# hgnc_symbol affy_hg_u133_plus_2
# 214844_s_at         DOK5         214844_s_at
# 1552842_at        HS6ST3          1552842_at
# 227259_at           CD47           227259_at
# 228546_at           DPP6           228546_at
# 1554863_s_at        DOK5        1554863_s_at
# 209627_s_at       OSBPL3         209627_s_at
# 219440_at           RAI2           219440_at
# 221325_at         KCNK13           221325_at
# 226992_at        NOSTRIN           226992_at

# select top probes with unique HGNC symbols
ft.genes <- NULL
ft.probes <- NULL
j <- 1
for(i in 1:nrow(ft.both)) {
  if(!(ft.both$hgnc_symbol[i] %in% ft.genes)) {
    ft.probes[j] <- ft.both$affy_hg_u133_plus_2[i]
    ft.genes[j] <- ft.both$hgnc_symbol[i]
    j <- j+1
  }
}
ft.adj <- cbind(ft.adj[ft.adj$ID %in% ft.probes[1:5],], hgnc = ft.genes[1:5])


# Do the same for the OSE genes
ose.both <- getBM(c("hgnc_symbol","affy_hg_u133_plus_2"), "affy_hg_u133_plus_2", ose.adj$ID, human)
ose.both <- ose.both[!(ose.both$hgnc_symbol == ""),]
rownames(ose.both) <- ose.both$affy_hg_u133_plus_2
ose.both <- na.omit(ose.both[ose.adj$ID,])

# select top 
ose.genes <- NULL
ose.probes <- NULL
j <- 1
for(i in 1:nrow(ose.both)) {
  if(!(ose.both$hgnc_symbol[i] %in% ose.genes)) {
    ose.probes[j] <- ose.both$affy_hg_u133_plus_2[i]
    ose.genes[j] <- ose.both$hgnc_symbol[i]
    j <- j+1
  }
}
ose.adj <- cbind(ose.adj[ose.adj$ID %in% ose.probes[1:5],], hgnc = ose.genes[1:5])

# Combine the top FT and OSE genes
mel.tab <- rbind(ft.adj, ose.adj)

# > dim(mel.tab)
# [1] 10 18


# ------------------------------ WU DATASET ---------------
# load the normalized Wu expressionset
load("/Users/wangshu/Documents/Sweave/Data/Wu.eset.vsnrma.9Nov10.Rdata")
wu.data <- Wu.eset

# Get U133a genes for the Wu dataset
probes.both <- getBM(c("affy_hg_u133a","affy_hg_u133_plus_2"), "affy_hg_u133_plus_2", mel.tab$ID, human)
for(i in 1:nrow(probes.both)) probes.both$bin[i] <- mel.tab[mel.tab$ID == probes.both$affy_hg_u133_plus_2[i], "Res.ft.n...ov.n"]

## result for probes.both
## thus 1 is FT genes, 2 is OV genes (bin grabbed "Res.ft.n...ov.n")
# > probes.both
# affy_hg_u133a affy_hg_u133_plus_2 bin
# 1                         1552842_at   1
# 2    209627_s_at         209627_s_at   1
# 3    209626_s_at         209627_s_at   1
# 4    207789_s_at           228546_at   1
# 5    202037_s_at         202035_s_at  -1
# 6    202036_s_at         202035_s_at  -1
# 7    202035_s_at         202035_s_at  -1
# 8    213857_s_at           227259_at   1
# 9      213055_at           227259_at   1
# 10   211075_s_at           227259_at   1
# 11   203439_s_at           203438_at  -1
# 12     203438_at           203438_at  -1
# 13                         229065_at  -1
# 14   214844_s_at         214844_s_at   1
# 15   220486_x_at         220486_x_at  -1
# 16   214096_s_at         214437_s_at  -1
# 17     214095_at         214437_s_at  -1
# 18   214437_s_at         214437_s_at  -1

# Combine probes for the same gene in the Wu dataset
probes.both.adj <- probes.both[!(probes.both$affy_hg_u133a == ""),]
## this removes those that are missing a value for u133a

wu.adj <- matrix(ncol = ncol(wu.data), nrow = nrow(mel.tab))
colnames(wu.adj) <- colnames(wu.data)  # these are the CEL file names
rownames(wu.adj) <- mel.tab$hgnc  # rownames are the gene symbols
for(i in 1:nrow(mel.tab)) {
  u133a <- probes.both.adj[grep(mel.tab$ID[i], probes.both.adj$affy_hg_u133_plus_2), "affy_hg_u133a"]
  if(length(u133a >= 0)) {
    curr.rows <- which(rownames(wu.data) %in% u133a)
    if(length(curr.rows) >= 0) {
      if (length(curr.rows)==1) wu.adj[i,] <- as.numeric(wu.data[curr.rows,])
      else wu.adj[i,] <- as.numeric(medianpolish(2^wu.data[curr.rows,])$exprs)
    }
  }
}

## for multiple probesets for the same gene, we apply medianpolish (robust
# way to average) used for multiple probes within a probeset
# medianpolish function is taking un-logged data (affy package)

# Get rid of genes that did not map to u133a
## use na.omit to remove all fo the NA's
wu.adj <- na.omit(wu.adj)
bin <- mel.tab$Res.ft.n...ov.n[mel.tab$hgnc %in% rownames(wu.adj)]

## this just tells again which are ft vs ov genes
# > bin <- mel.tab$Res.ft.n...ov.n[mel.tab$hgnc %in% rownames(wu.adj)]
# > bin
# [1]  1  1  1  1 -1 -1 -1 -1

# Get pdata

pdata.rough <- read.csv("/Users/wangshu/Documents/Sweave/Data/E-GEOD-6008.sdrf.csv",stringsAsFactors=FALSE)

## stringsAsFactors: logical: should character vectors be converted to factors?
## > pdata.rough$Source.Name
# [1] "GSM139463" "GSM139412", etc

rownames(pdata.rough) <- paste(pdata.rough[,1],"CEL","gz",sep=".")
## name the rows as the CEL file names
# > rownames(pdata.rough)
# [1] "GSM139463.CEL.gz" "GSM139412.CEL.gz"

pdata.rough <- pdata.rough[colnames(wu.data),]
## reminder wu.data is the ExpressionSet, so the column names are the CEL file names
## this is making sure the phenodata (rows) are compatible with column names from wu.data
# colnames(wu.data) = rownames(pdata.rough)


# Create a unified histologic subtype variable
# > colnames(pdata.rough)
# [1] "Source.Name" "Description"
Type <- NULL
for(i in 1:nrow(pdata.rough)) {
  string <- strsplit(pdata.rough[i,"Description"], split="")[[1]]
  char.id <- c("T","u","m","o","r","_","T","y","p","e",":"," ")
  semicol <- grep(";",string)
  for(j in 1:(length(string)-length(char.id))) {
    if(sum(string[j:(length(char.id)+j-1)] == char.id) == length(char.id)) string.num <- j+length(char.id)
  }
  Type[i] <- paste(string[string.num:(min(semicol[semicol > string.num])-1)],collapse="")
}
Type[Type == "N/A"] <- NA

## this gives: 
# > Type
# [1] "Clear_Cell"   "Clear_Cell"   "Clear_Cell"   "Clear_Cell"   "Clear_Cell"
# [6] "Clear_Cell"   "Clear_Cell" etc

names(Type) <- rownames(pdata.rough)
# > Type
# GSM139377.CEL.gz GSM139378.CEL.gz GSM139379.CEL.gz GSM139380.CEL.gz
# "Clear_Cell"     "Clear_Cell"     "Clear_Cell"     "Clear_Cell"
# GSM139381.CEL.gz GSM139382.CEL.gz GSM139383.CEL.gz GSM139384.CEL.gz
# "Clear_Cell"     "Clear_Cell" etc

# > table(Type)
# Type
# Clear_Cell Endometrioid     Mucinous          OSE       Serous
# 8           37           13            4           41

Stage <- NULL
for(i in 1:nrow(pdata.rough)) {
  string <- strsplit(pdata.rough[i,"Description"], split="")[[1]]
  char.id <- c("s","t","a","g","e",":"," ")
  semicol <- grep(";",string)
  for(j in 1:(length(string)-length(char.id))) {
    if(sum(string[j:(length(char.id)+j-1)] == char.id) == length(char.id)) string.num <- j + length(char.id)
  }
  Stage[i] <- paste(string[string.num:(min(semicol[semicol > string.num])-1)],collapse="")
}
Stage[Stage == "N/A"] <- NA
Stage[Stage == "I a"] <- "1"
Stage[Stage == "I c"] <- "1"
Stage[Stage == "1A"] <- "1"
Stage[Stage == "1C"] <- "1"
Stage[Stage == "1a"] <- "1"
Stage[Stage == "1c"] <- "1"
Stage[Stage == "2A"] <- "2"
Stage[Stage == "2B"] <- "2"
Stage[Stage == "2C"] <- "2"
Stage[Stage == "2a"] <- "2"
Stage[Stage == "2c"] <- "2"
Stage[Stage == "3B"] <- "3"
Stage[Stage == "3C"] <- "3"
Stage[Stage == "3c"] <- "3"
Stage[Stage == "3D"] <- "3"

# > table(Stage)
# Stage
# 1  2  3  4
# 35 11 44  9

Grade <- NULL
for(i in 1:nrow(pdata.rough)) {
  string <- strsplit(pdata.rough[i,"Description"], split="")[[1]]
  char.id <- c("g","r","a","d","e",":"," ")
  semicol <- grep(";",string)
  for(j in 1:(length(string)-length(char.id))) {
    if(sum(string[j:(length(char.id)+j-1)] == char.id) == length(char.id)) string.num <- j+length(char.id)
  }
  Grade[i] <- paste(string[string.num:(min(semicol[semicol > string.num])-1)],collapse="")
}
Grade[is.na(Type)] <- NA
Grade[Grade == "N/A"] <- NA

# > table(Grade)
# Grade
# 1   2 2-3   3
# 19  17   5  33

# Do various exploratory tests

stage <- ordered(Stage)
## gives you a stage variable with levels
# e.g. Levels: 1 < 2 < 3 < 4

gene.stage <- wu.adj[,!is.na(stage)]
## reminder, wu.adj is 8 rows (selected genes), 103 columns (CELfiles)
## tells it to only look at the CEL files where there is a value for stage
## previously noted that there were four with "NA" values
## output is
# > dim(gene.stage)
# [1]  8 99 (8 genes, 99 samples with stage info)

stage <- stage[!is.na(stage)]
## this removstagees the missings from stage variable as well
gt.stage <- gt(ordered(stage), t(gene.stage))
## ordered(stage) is the response variable
## tests the null hypothesis that the covariates in alternative (gene.stage)
# are not associated with the response

## Here is the code for the global test ("gt" is the function implemented in 2009)
# Description:
# Tests a low-dimensional null hypothesis against a potentially
# high-dimensional alternative in regression models (linear
# regression, logistic regression, poisson regression, Cox
# proportional hazards model).

# In gene set testing in
# microarray data analysis alternative may be a matrix of gene
# expression measurements, and the aim is to find which of a
# collection of predefined subsets of the genes (e.g. Gene
# Ontology terms or KEGG pathways) is most associated with the
# response

## Note these 3 are results from GlobalTests
# > gt.stage
# p-value Statistic Expected Std.dev #Cov
# 1 9.55e-12      12.3     1.03   0.527    8

grade <- ordered(Grade)
gene.grade <- wu.adj[,!is.na(grade)]
grade <- grade[!is.na(grade)]
gt.grade <- gt(ordered(grade), t(gene.grade))

# > gt.grade
# p-value Statistic Expected Std.dev #Cov
# 1 2.03e-07      9.95     1.39     0.7    8

gene.type <- wu.adj[,!is.na(Type)]
type <- Type[!is.na(Type)]
gt.type <- gt(as.factor(as.character(type)), t(gene.type))
## Reminder
## OSE is used to generate this p-value
# > table(Type)
# Type
#   Clear_Cell Endometrioid     Mucinous       Serous 
#            8           37           13           41 
# > gt.type
#    p-value Statistic Expected Std.dev #Cov
# 1 3.74e-10      10.8     1.03   0.524    8


## reminder: gene.stage is a matrix like this

# GSM139473.CEL.gz GSM139474.CEL.gz GSM139475.CEL.gz
# DOK5           10.276914         9.044567        10.149500
# CD47           10.422950         9.604965        10.051869
# DPP6            7.144967         7.470085         7.327839
# OSBPL3          7.718119         7.512260         6.417728
# STC2            7.839947         6.949597         7.038743
# SFRP1           6.876388         7.161776         6.937320
# SHMT2           7.631528         7.144393         8.147059
# TMEM164         8.115610         8.494884         9.123595

## reminder 
# > dim(gene.stage)
# [1]  8 99 (8 genes, 99 samples with stage info)

# > probes.both.adj
# affy_hg_u133a affy_hg_u133_plus_2 bin
# 2    209627_s_at         209627_s_at   1
# 3    209626_s_at         209627_s_at   1
# 4    207789_s_at           228546_at   1


## See how each gene is related to stage and grade

stage.p <- NULL
stage.coef <- NULL
for(i in 1:nrow(gene.stage)) {
  stage.model <- lrm(stage ~ as.numeric(gene.stage[i,]), na.action=na.pass)
  stage.p[i] <- stage.model$stats["P"]
  stage.coef[i] <- stage.model$coefficients["gene.stage"]
}

# > stage.p
# [1] 2.305900e-11 1.944181e-05 1.503802e-02 1.690472e-02 1.529746e-03
# [6] 3.873047e-05 2.424461e-01 4.103876e-01

# > stage.coef
# [1]  1.3252886  1.0111846  1.6303119 -0.8396347 -1.7238029 -1.4155974  0.4220139
# [8]  0.3439440

grade.p <- NULL
grade.coef <- NULL
for(i in 1:nrow(gene.grade)) {
  grade.model <- lrm(grade ~ as.numeric(gene.grade[i,]), na.action=na.pass)
  grade.p[i] <- grade.model$stats["P"]
  grade.coef[i] <- grade.model$coefficients["gene.grade"]
}

# > grade.coef
# [1]  1.1255299  1.0402291 -0.6878462 -0.2835020 -0.6887416 -1.2559906  0.8345742
# [8]  0.9199047

type.p <- NULL
type.coef <- NULL
type <- as.factor(Type)
for(i in 1:nrow(wu.adj)) {
  type.p[i] <- lmp(lm(as.numeric(wu.adj[i,]) ~ type))
  type.coef[i] <- lm(as.numeric(wu.adj[i,]) ~ as.numeric(type == "Serous"))$coefficients[2]
}

type.p.adj <- p.adjust(type.p, method = "BH")
grade.p.adj <- p.adjust(grade.p, method = "BH") 
stage.p.adj <- p.adjust(stage.p, method = "BH") 

sum(grade.p.adj < 0.05)
# [1] 3
sum(stage.p.adj < 0.05)
# [1] 6
sum(type.p.adj < 0.05)
# [1] 6

tissue <- sub(-1, "OV", bin)
## sub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE,
# fixed = FALSE, useBytes = FALSE)

tissue <- sub("1", "FT", tissue)

## this code essentially writes out the results to a table
sig <- function(val) return(signif(val, 3))
ind.genes <- data.frame(Gene = rownames(wu.adj), Tissue = tissue, Grade.pval = sig(grade.p.adj), Grade.coef = sig(grade.coef), 
                        Stage.pval = sig(stage.p.adj), Stage.coef = sig(stage.coef), Type.pval = sig(type.p.adj), Serous.coef = sig(type.coef))

## reminder: > rownames(wu.adj)
# [1] "DOK5"    "CD47"    "DPP6"    "OSBPL3"  "STC2"    "SFRP1"   "SHMT2"
# [8] "TMEM164"

# > ind.genes
# Gene Tissue Grade.pval Grade.coef Stage.pval Stage.coef Type.pval
# 1    DOK5     FT   8.87e-06      1.130   1.84e-10      1.330  5.48e-09
# 2    CD47     FT   5.62e-04      1.040   7.78e-05      1.010  4.00e-05
# 3    DPP6     FT   5.82e-01     -0.688   2.25e-02      1.630  4.23e-01
# 4  OSBPL3     FT   5.18e-01     -0.284   2.25e-02     -0.840  7.12e-02
# 5    STC2     OV   2.33e-01     -0.689   3.06e-03     -1.720  1.84e-04
# 6   SFRP1     OV   2.50e-03     -1.260   1.03e-04     -1.420  2.43e-06
# 7   SHMT2     OV   1.14e-01      0.835   2.77e-01      0.422  1.24e-02
# 8 TMEM164     OV   1.35e-01      0.920   4.10e-01      0.344  3.61e-02
# Serous.coef
# 1       1.380
# 2       0.783
# 3       0.119
# 4      -0.116
# 5      -0.108
# 6      -0.519
# 7       0.248
# 8       0.175	

ind.genes[grade.p.adj < 0.05 & stage.p.adj < 0.05 & type.p.adj < 0.05,]
## find out which genes are associated with all clinical characteristics
# e.g. driving the association
# > ind.genes[grade.p.adj < 0.05 & stage.p.adj < 0.05 & type.p.adj < 0.05,]
# Gene Tissue Grade.pval Grade.coef Stage.pval Stage.coef Type.pval
# 1  DOK5     FT   8.87e-06       1.13   1.84e-10       1.33  5.48e-09
# 2  CD47     FT   5.62e-04       1.04   7.78e-05       1.01  4.00e-05
# 6 SFRP1     OV   2.50e-03      -1.26   1.03e-04      -1.42  2.43e-06
# Serous.coef
# 1       1.380
# 2       0.783
# 6      -0.519

write.table(ind.genes, "/Users/wangshu/Documents/Sweave/Results/Wu.dataset.clinical.associations.txt", quote=F, sep="\t")

# Compute score
score.wu <- t(wu.adj) %*% bin
## Note: bin defined above as 
# bin <- mel.tab$Res.ft.n...ov.n[mel.tab$hgnc %in% rownames(wu.adj)]
# > bin
# [1]  1  1  1  1 -1 -1 -1 -1

## use the pdf for high res
pdf("/Users/wangshu/Documents/Sweave/Results/OCFNLER.genes.in.Wu-cut.both.lists.pdf")
plot(density(score.wu), xlab="Score", ylab="Density", main="FNE/OCE gene signature score in Wu dataset")
dev.off()
## classification will give you 1's and 2's and to get back to zero and 1 subtract 1
classes <- Mclust(score.wu, 2)$classification - 1
## this gives each CELfile a "0" or a "1"

# Mclust                 package:mclust                  R Documentation
# Model-Based Clustering
# Description:
# The optimal model according to BIC for EM initialized by
# hierarchical clustering for parameterized Gaussian mixture models.

# Usage:
# Mclust(data, G=NULL, modelNames=NULL, prior=NULL, control=emControl(),
# initialization=NULL, warn=FALSE, ...)

# References:

# C. Fraley and A. E. Raftery (2006, revised 2010).  MCLUST Version
# 3 for R: Normal Mixture Modeling and Model-Based Clustering,
# Technical Report no. 504, Department of Statistics, University of
# Washington.

# C. Fraley and A. E. Raftery (2002).  Model-based clustering,
# discriminant analysis, and density estimation.  _Journal of the
# American Statistical Association 97:611:631_.

# C. Fraley and A. E. Raftery (2005, revised 2009).  Bayesian
# regularization for normal mixture estimation and model-based
# clustering.  Technical Report, Department of Statistics,
# University of Washington.

# C. Fraley and A. E. Raftery (2007).  Bayesian regularization for
# normal mixture estimation and model-based clustering. _Journal of
# Classification 24:155-181_.

#save p-values for Sweave
wu.s.p <- get.logistic.p(glm(classes ~ ordered(Stage), binomial))
wu.g.p <- get.logistic.p(glm(classes ~ ordered(Grade), binomial))
wu.t.p  <- fisher.test(xtabs(~classes + Type))$p.value
wu.bipart.table <- cbind(xtabs(~classes + Grade), xtabs(~classes + Stage), xtabs(~classes + Type))
wt<-rbind(cbind(sum(wu.bipart.table[1,5:8]),0,sum(wu.bipart.table[2,5:8]),0),
          cbind(wu.bipart.table[1,12],round(wu.bipart.table[1,12]*100/sum(wu.bipart.table[1,9:12])),
                wu.bipart.table[2,12],round(wu.bipart.table[2,12]*100/sum(wu.bipart.table[2,9:12]))),
          cbind(wu.bipart.table[1,10],round(wu.bipart.table[1,10]*100/sum(wu.bipart.table[1,9:12])),
                wu.bipart.table[2,10],round(wu.bipart.table[2,10]*100/sum(wu.bipart.table[2,9:12]))),
          cbind(wu.bipart.table[1, 9],round(wu.bipart.table[1, 9]*100/sum(wu.bipart.table[1,9:12])),
                wu.bipart.table[2, 9],round(wu.bipart.table[2, 9]*100/sum(wu.bipart.table[2,9:12]))),
          cbind(wu.bipart.table[1,11],round(wu.bipart.table[1,11]*100/sum(wu.bipart.table[1,9:12])),
                wu.bipart.table[2,11],round(wu.bipart.table[2,11]*100/sum(wu.bipart.table[2,9:12]))),
          cbind(wu.bipart.table[1,5],round(wu.bipart.table[1,5]*100/sum(wu.bipart.table[1,5:8])),
                wu.bipart.table[2,5],round(wu.bipart.table[2,5]*100/sum(wu.bipart.table[2,5:8]))),
          cbind(wu.bipart.table[1,6],round(wu.bipart.table[1,6]*100/sum(wu.bipart.table[1,5:8])),
                wu.bipart.table[2,6],round(wu.bipart.table[2,6]*100/sum(wu.bipart.table[2,5:8]))),
          cbind(wu.bipart.table[1,7],round(wu.bipart.table[1,7]*100/sum(wu.bipart.table[1,5:8])),
                wu.bipart.table[2,7],round(wu.bipart.table[2,7]*100/sum(wu.bipart.table[2,5:8]))),
          cbind(wu.bipart.table[1,8],round(wu.bipart.table[1,8]*100/sum(wu.bipart.table[1,5:8])),
                wu.bipart.table[2,8],round(wu.bipart.table[2,8]*100/sum(wu.bipart.table[2,5:8]))),
          cbind(wu.bipart.table[1,1],round(wu.bipart.table[1,1]*100/sum(wu.bipart.table[1,1:4])),
                wu.bipart.table[2,1],round(wu.bipart.table[2,1]*100/sum(wu.bipart.table[2,1:4]))),
          cbind(wu.bipart.table[1,2],round(wu.bipart.table[1,2]*100/sum(wu.bipart.table[1,1:4])),
                wu.bipart.table[2,2],round(wu.bipart.table[2,2]*100/sum(wu.bipart.table[2,1:4]))),
          cbind(wu.bipart.table[1,3],round(wu.bipart.table[1,3]*100/sum(wu.bipart.table[1,1:4])),
                wu.bipart.table[2,3],round(wu.bipart.table[2,3]*100/sum(wu.bipart.table[2,1:4]))),
          cbind(wu.bipart.table[1,4],round(wu.bipart.table[1,4]*100/sum(wu.bipart.table[1,1:4])),
                wu.bipart.table[2,4],round(wu.bipart.table[2,4]*100/sum(wu.bipart.table[2,1:4])))
)
# > wu.bipart.table
#    1  2 2-3  3  1 2  3 4 Clear_Cell Endometrioid Mucinous Serous
# 0 14 13   4  9 28 8 12 2          7           25       11      7
# 1  5  4   1 24  7 3 32 7          1           12        2     34

write.table(wu.bipart.table, "/Users/wangshu/Documents/Sweave/Results/OCFNLER.genes-Wu.bipartition.clinical.table.txt", quote=F, sep="\t")

summary(lm(score.wu ~ ordered(Stage)))
# Call:
# lm(formula = score.wu ~ ordered(Stage))
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -6.0086 -0.9708 -0.0759  1.2345  3.6007 
# 
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.0428     0.2288   8.930 3.20e-14 ***
# ordered(Stage).L   2.6112     0.4760   5.486 3.39e-07 ***
# ordered(Stage).Q   0.1447     0.4575   0.316    0.752    
# ordered(Stage).C  -0.7038     0.4383  -1.606    0.112    
# ---
# Signif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01 \u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1 
# 
# Residual standard error: 1.818 on 95 degrees of freedom
#   (4 observations deleted due to missingness)
# Multiple R-squared: 0.3405,	Adjusted R-squared: 0.3196 
# F-statistic: 16.35 on 3 and 95 DF,  p-value: 1.209e-08 


# ----------------------TOTHILL DATA ------------------------

# Load Tothill data
load("/Users/wangshu/Documents/Sweave/Data/tothillExpressionData.rdat")
exprs <- tothill.e.data$expr
pdata <- tothill.e.data$pheno

## this code restricts to the malignant samples
exprs <- exprs[,pdata$Type == "MAL"]
pdata <- pdata[pdata$Type == "MAL",]

# Make survival variables
t.rfs <- pdata$Time.to.Relapse
t.os <- pdata$Time.to.Death
e.rfs <- as.numeric(pdata$Status == "R" | pdata$Status == "D")
# this creates a numeric matrix of 1's and 0's, guessing the 1's had "R" or "D"
## rfs = relapse free survival
# R = Relapse
# D = Death
# D* =
# PF = progression free

## > table(pdata$Status)
# D      D*      PF       R Unknown
# 111       2      74      77       3

e.os <- as.numeric(pdata$Status == "D")

pdata <- cbind(pdata, t.rfs, t.os, e.rfs, e.os)

# Edit clinical data
# Changing Stage
## > table(pdata$Stage)
# I      II     III      IV Unknown
# 16      14     212      21       4
pdata$Stage <- sub("Unknown", NA, pdata$Stage)
pdata$Stage <- ordered(pdata$Stage)
## > table(pdata$Stage)
# I  II III  IV
# 16  14 212  21


# Changing Grade and treatments
pdata$Grade = sub("Unknown", NA, pdata$Grade)
pdata$Grade <- ordered(pdata$Grade)
# > table(pdata$Grade)
# 
#   1   2   3 
#  11  97 155 

# Change residual disease
## pdata$res.dis ==1 (>1 cm)
## pdata$res.dis ==0 (<1 cm)

pdata$res.dis <- NA
for(i in 1:nrow(pdata)) {
  if(pdata$Residual.Disease[i] == ">1<2cm" | pdata$Residual.Disease[i] == "<2cm") {
    pdata$res.dis[i] <- 1
  }
  else if(pdata$Residual.Disease[i] == "<1cm" | pdata$Residual.Disease[i] == "NIL") {
    pdata$res.dis[i] <- 0
  }
  else {pdata$res.dis[i] <- NA}
}

resid.disease = pdata$res.dis
# > table(resid.disease)
# resid.disease
# 0   1
# 145  81

# change subtype 
# (I don't analyze adenocarcinoma, because there is only 1)
pdata$Subtype = sub("Adeno", NA, pdata$Subtype)
# > table(pdata$Subtype)
# 
# Endo  Ser 
#   20  246 

# Get score
good.exprs <- exprs[mel.tab$ID,]
bin <- mel.tab$Res.ft.n...ov.n
score <- t(as.matrix(good.exprs)) %*% as.numeric(bin)

pdf("/Users/wangshu/Documents/Sweave/Results/OCFNLER.genes.in.Tothill-cut.both.lists.pdf")
plot(density(score), xlab = "Score", ylab = "Density", main = "FNE/OCE gene signature score in Tothill dataset")
dev.off()

# Get bipartition
bipart.tot <- Mclust(score,2)$classification -1

# Compare bipartition to clinical variables
to.t.p <- fisher.test(xtabs(~pdata$Subtype + bipart.tot))$p.value
to.s.p <- get.logistic.p(glm(bipart.tot ~ Stage, data = pdata, binomial))
to.g.p <- get.logistic.p(glm(bipart.tot ~ Grade, data = pdata, binomial))

lmp(lm(as.numeric(score) ~ Subtype, data = pdata))
lmp(lm(as.numeric(score) ~ Stage, data = pdata))
lmp(lm(as.numeric(score) ~ Grade, data = pdata))

tot.bipart.table <- cbind(xtabs(~bipart.tot + pdata$Grade), xtabs(~bipart.tot + pdata$Stage), xtabs(~bipart.tot + pdata$Subtype))
# > tot.bipart.table
# 1  2   3  I II III IV Endo Ser
# 0 2 12  20  3  3  29  0    7  28
# 1 9 85 135 13 11 183 21   13 218
write.table(tot.bipart.table, "/Users/wangshu/Documents/Sweave/Results/OCFNLER.genes-Tothil.bipartition.clinical.table.txt", quote=F, sep="\t")
tt<-rbind(cbind(sum(tot.bipart.table[1,8:9]),0,sum(tot.bipart.table[2,8:9]),0),
          cbind(tot.bipart.table[1,9],round(tot.bipart.table[1,9]*100/sum(tot.bipart.table[1,8:9])),
                tot.bipart.table[2,9],round(tot.bipart.table[2,9]*100/sum(tot.bipart.table[2,8:9]))),
          cbind(tot.bipart.table[1,8],round(tot.bipart.table[1,8]*100/sum(tot.bipart.table[1,8:9])),
                tot.bipart.table[2,8],round(tot.bipart.table[2,8]*100/sum(tot.bipart.table[2,8:9]))),
          cbind(tot.bipart.table[1,4],round(tot.bipart.table[1,4]*100/sum(tot.bipart.table[1,4:7])),
                tot.bipart.table[2,4],round(tot.bipart.table[2,4]*100/sum(tot.bipart.table[2,4:7]))),
          cbind(tot.bipart.table[1,5],round(tot.bipart.table[1,5]*100/sum(tot.bipart.table[1,4:7])),
                tot.bipart.table[2,5],round(tot.bipart.table[2,5]*100/sum(tot.bipart.table[2,4:7]))),
          cbind(tot.bipart.table[1,6],round(tot.bipart.table[1,6]*100/sum(tot.bipart.table[1,4:7])),
                tot.bipart.table[2,6],round(tot.bipart.table[2,6]*100/sum(tot.bipart.table[2,4:7]))),
          cbind(tot.bipart.table[1,7],round(tot.bipart.table[1,7]*100/sum(tot.bipart.table[1,4:7])),
                tot.bipart.table[2,7],round(tot.bipart.table[2,7]*100/sum(tot.bipart.table[2,4:7]))),
          cbind(tot.bipart.table[1,1],round(tot.bipart.table[1,1]*100/sum(tot.bipart.table[1,1:3])),
                tot.bipart.table[2,1],round(tot.bipart.table[2,1]*100/sum(tot.bipart.table[2,1:3]))),
          cbind(tot.bipart.table[1,2],round(tot.bipart.table[1,2]*100/sum(tot.bipart.table[1,1:3])),
                tot.bipart.table[2,2],round(tot.bipart.table[2,2]*100/sum(tot.bipart.table[2,1:3]))),
          cbind(tot.bipart.table[1,3],round(tot.bipart.table[1,3]*100/sum(tot.bipart.table[1,1:3])),
                tot.bipart.table[2,3],round(tot.bipart.table[2,3]*100/sum(tot.bipart.table[2,1:3])))
)

# Survival analysis
library(survival)
library(survdiff)
coxph(Surv(t.rfs, e.rfs) ~ score, pdata)
coxph(Surv(t.os, e.os) ~ score, pdata)
surv.subtypes <- survfit(Surv(t.rfs, e.rfs) ~ bipart.tot, pdata)
surv.subtypes.os <- survfit(Surv(t.os, e.os) ~ bipart.tot, pdata)
surv.subtypes.diff <- survdiff(Surv(t.rfs, e.rfs) ~ bipart.tot, pdata)
surv.subtypes.diff.os <- survdiff(Surv(t.os, e.os) ~ bipart.tot, pdata)

tot.age = as.numeric(pdata$Age)

# Create KM curves
bitmap("/Users/wangshu/Documents/Sweave/Results/OCFNLER.survival.curves.inTothill.png", res = 300, height = 10, width = 20)
par(mfrow = c(1,2))

plot(surv.subtypes, main = "Disease-free survival", ylab = "Probability of survival", xlab = "Months", col=3:6,
     cex.main = 2, lwd = 3, cex.axis = 1.5, cex.lab = 1.5)
n <- surv.subtypes$n   
groups <- names(surv.subtypes$strata)
groups <- sub("0", "OV-like", groups)
groups <- sub("1", "FT-like", groups)
legend("topright",paste(groups, ", n = ", n,sep = ""), col=3:6, fill=3:6, cex = 1.5)
legend("bottomleft", paste("p-value",signif(p.value(surv.subtypes.diff), digits=3),sep=": "), cex = 1.5)

plot(surv.subtypes.os, main = "Overall survival", ylab = "Probability of survival", xlab = "Months", col=3:6,
     cex.main = 2, lwd = 3, cex.axis = 1.5, cex.lab = 1.5)
n <- surv.subtypes.os$n   
groups <- names(surv.subtypes.os$strata)
groups <- sub("0", "OV-like", groups)
groups <- sub("1", "FT-like", groups)
legend("topright",paste(groups, ", n = ", n,sep = ""), col=3:6, fill=3:6, cex = 1.5)
legend("bottomleft", paste("p-value",signif(p.value(surv.subtypes.diff.os), digits=3),sep=": "), cex = 1.5)
dev.off()

# See how each gene is related to stage and grade
gene.stage <- good.exprs[,!is.na(pdata$Stage)]
stage <- na.omit(pdata$Stage)
stage.p <- NULL
stage.coef <- NULL
for(i in 1:nrow(gene.stage)) {
  stage.model <- lrm(stage ~ as.numeric(gene.stage[i,]))
  stage.p[i] <- stage.model$stats["P"]
  stage.coef[i] <- stage.model$coefficients["gene.stage"]
}

# See how each gene is related to grade and grade
gene.grade <- good.exprs[,!is.na(pdata$Grade)]
grade <- na.omit(pdata$Grade)
grade.p <- NULL
grade.coef <- NULL
for(i in 1:nrow(gene.grade)) {
  grade.model <- lrm(grade ~ as.numeric(gene.grade[i,]))
  grade.p[i] <- grade.model$stats["P"]
  grade.coef[i] <- grade.model$coefficients["gene.grade"]
}

gene.type <- good.exprs[,pdata$Subtype != "Adeno"]
type.p <- NULL
type.coef <- NULL
type <- as.factor(pdata$Subtype[pdata$Subtype != "Adeno"])
for(i in 1:nrow(gene.type)) {
  type.p[i] <- lmp(lm(as.numeric(gene.type[i,]) ~ type))
  type.coef[i] <- lm(as.numeric(gene.type[i,]) ~ as.numeric(type == "Ser"))$coefficients[2]
}

type.p.adj <- p.adjust(type.p, method = "BH")
grade.p.adj <- p.adjust(grade.p, method = "BH") 
stage.p.adj <- p.adjust(stage.p, method = "BH") 

sum(grade.p.adj < 0.05)
sum(stage.p.adj < 0.05)
sum(type.p.adj < 0.05)

tissue <- sub(-1, "OSE", bin)
tissue <- sub("1", "FT", tissue)
sig <- function(val) return(signif(val, 3))
ind.genes <- data.frame(Gene = mel.tab$hgnc, Tissue = tissue, Grade.pval = sig(grade.p.adj), Grade.coef = sig(grade.coef), 
                        Stage.pval = sig(stage.p.adj), Stage.coef = sig(stage.coef), Type.pval = sig(type.p.adj), Serous.coef = sig(type.coef))
ind.genes[stage.p.adj < 0.05 & type.p.adj < 0.05,]

# gt.stage <- gt(stage, t(gene.stage))
# gt.grade <- gt(grade, t(gene.grade))
# gt.type <- gt(type, t(gene.type))

write.table(ind.genes, "Results/OCFNLER.genes-Tothill.dataset.clinical.associations.txt", quote=F, sep="\t")

# coxph(Surv(t.rfs, e.rfs) ~ bipart.tot + Grade + Stage + Subtype, data = pdata)
# coxph(Surv(t.os, e.os)   ~ bipart.tot + Grade + Stage + Subtype, data = pdata)
# contribution of bipart.tot

tot.coxph.rfs.p.value <- summary(coxph(Surv(t.rfs, e.rfs) ~ bipart.tot + Grade + Stage + Subtype + tot.age + resid.disease, data = pdata))$coefficients[1,5]


# > coxph(Surv(t.rfs, e.rfs) ~ bipart.tot + Grade + Stage + Subtype + tot.age + resid.disease, data = pdata)
# Call:
# coxph(formula = Surv(t.rfs, e.rfs) ~ bipart.tot + Grade + Stage +
# Subtype + tot.age + resid.disease, data = pdata)
# coef exp(coef) se(coef)       z       p
# bipart.tot     0.6844     1.983  0.27974  2.4466 0.01400
# Grade.L       -0.1361     0.873  0.37023 -0.3676 0.71000
# Grade.Q        0.0566     1.058  0.24008  0.2356 0.81000
# Stage.L        1.9937     7.342  0.53870  3.7009 0.00021
# Stage.Q       -0.4459     0.640  0.44290 -1.0067 0.31000
# Stage.C       -0.0179     0.982  0.35098 -0.0511 0.96000
# SubtypeSer     0.0763     1.079  0.43272  0.1763 0.86000
# tot.age        0.0121     1.012  0.00903  1.3358 0.18000
# resid.disease  0.3323     1.394  0.17920  1.8544 0.06400

# Likelihood ratio test=49.8  on 9 df, p=1.17e-07  n=219 (48 observations deleted due to missingness)
# Warning message:
# In model.matrix.default(Terms, m) :
# variable 'Subtype' converted to a factor

tot.coxph.os.p.value <- summary(coxph(Surv(t.os, e.os) ~ bipart.tot + Grade + Stage + Subtype + tot.age + resid.disease, data = pdata))$coefficients[1,5]
# > coxph(Surv(t.os, e.os) ~ bipart.tot + Grade + Stage + Subtype + tot.age + resid.disease, data = pdata)
# Call:
# coxph(formula = Surv(t.os, e.os) ~ bipart.tot + Grade + Stage +
# Subtype + tot.age + resid.disease, data = pdata)
# coef exp(coef) se(coef)        z      p
# bipart.tot     0.3394  1.40e+00 3.55e-01  0.95717 0.3400
# Grade.L        0.1341  1.14e+00 4.35e-01  0.30810 0.7600
# Grade.Q        0.0443  1.05e+00 2.93e-01  0.15119 0.8800
# Stage.L        1.5937  4.92e+00 5.85e-01  2.72653 0.0064
# Stage.Q        0.8279  2.29e+00 6.39e-01  1.29611 0.1900
# Stage.C       -0.6653  5.14e-01 7.05e-01 -0.94317 0.3500
# SubtypeSer    16.4430  1.38e+07 2.16e+03  0.00763 0.9900
# tot.age        0.0250  1.03e+00 1.19e-02  2.09870 0.0360
# resid.disease  0.2590  1.30e+00 2.19e-01  1.18331 0.2400

# Likelihood ratio test=36.7  on 9 df, p=2.94e-05  n=221 (46 observations deleted due to missingness)
# Warning messages:
# 1: In model.matrix.default(Terms, m) :
# variable 'Subtype' converted to a factor
# 2: In fitter(X, Y, strats, offset, init, control, weights = weights,  :
# Loglik converged before variable  7 ; beta may be infinite.

save(list=c('wt', 'wu.s.p', 'wu.g.p', 'wu.t.p', 
            'tt', 'to.s.p', 'to.g.p', 'to.t.p', 'tot.coxph.rfs.p.value', 'tot.coxph.os.p.value',
            'surv.subtypes', 'surv.subtypes.diff', 'surv.subtypes.os', 'surv.subtypes.diff.os'
), 
file="/Users/wangshu/Documents/Sweave/Results/forSweave.RData")

