Bioinformatics_project
My project in Bioimformatics class

# Introduction
1.	They hypothesized that distinct cells-of-origin may play a role in determining ovarian tumor phenotype and outcome. 
2.	they describe a new cell culture medium for in vitro culture of paired normal human ovarian (OV) and fallopian tube (FT) epithelial cells from donors without cancer. 
3.	These cells have been cultured individually for short periods of time, this is the first long-term culture of both cell types from the same donors. 
4.	Do analysis of the gene expression profiles of the cultured OV/FT cells 
5.	identified a normal cell-of-origin gene signature that classified primary ovarian cancers into OV-like and FT-like subgroups; 
6.	the normal cell-of-origin may be a source of ovarian tumor heterogeneity and the associated differences in tumor outcome.


## 1- What is the goal of the study?
objective was to identify a gene expression signature that could identify the normal cell-of-origin of ovarian carcinomas

Towards this goal, it was necessary to develop a new cell culture medium and methods to isolate and propagate normal ovarian epithelium and fallopian tube epithelium as paired cultured cells from the same individuals

## 2- What type of data the authors analyzed?
The authors used microarray RNA data for the study. The platform was Affymetrix Human X3P Array.

## 3- What analysis did they performed?
Microarray analysis, datasets generated using total unamplified RNA.

## 4- What were the conclusions of the study?
Previous tissue-based studies have reported similarities in gene expression between normal fallopian tube epithelium and papillary serous carcinoma. And we observed that >80% of serous carcinomas in two independent ovarian cancer gene expression datasets were classified as FT-like using the FNE versus OCE cell-of-origin signature. This result provides that a large proportion of high grade serous carcinomas may arise from the non-ciliated fallopian tube epithelium.

# Datasets
1.	Merritt MA: analyzed 12 samples from two donor patients and established cultures of both ovarian epithelium and fallopian tube epithelium (hTERT immortalized), each with 3 replicates (different culture passages)
2.	Wu R: 99 individual ovarian tumors (37 endometrioid, 41 serous, 13 mucinous, and 8 clear cell carcinomas) and 4 individual normal ovary samples. RNA expression was analyzed using one Affymetrix HG_U133A array per sample.
3.	Tothill R: Randomly selected 285 samples from the AOCS ( Australian Ovarian Cancer Study) was expression profiled on the affymetrix U133_plus2 platform to identify novel subtypes of ovarian tumor

# Methods
1.	Load Dataset
2.	Log2 transform
3.	Limma (Identify different expression genes)
4.	duplicateCorrelation (identify the correlation samples derived from the same patient)
5.	write.fit(write the results to a file)
6.	Three different functions to find P-Value
7.	Survival analysis

# Results




