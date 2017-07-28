Bioinformatics_project
My project in Bioimformatics class

# Normal cell-of-origin predicts ovarian tumor outcomes.


## Introduction
1.	They hypothesized that distinct cells-of-origin may play a role in determining ovarian tumor phenotype and outcome. 
2.	they describe a new cell culture medium for in vitro culture of paired normal human ovarian (OV) and fallopian tube (FT) epithelial cells from donors without cancer. 
3.	These cells have been cultured individually for short periods of time, this is the first long-term culture of both cell types from the same donors. 
4.	Do analysis of the gene expression profiles of the cultured OV/FT cells 
5.	identified a normal cell-of-origin gene signature that classified primary ovarian cancers into OV-like and FT-like subgroups; 
6.	the normal cell-of-origin may be a source of ovarian tumor heterogeneity and the associated differences in tumor outcome.

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/cells.png)

### 1- What is the goal of the study?
objective was to identify a gene expression signature that could identify the normal cell-of-origin of ovarian carcinomas

Towards this goal, it was necessary to develop a new cell culture medium and methods to isolate and propagate normal ovarian epithelium and fallopian tube epithelium as paired cultured cells from the same individuals

### 2- What type of data the authors analyzed?
The authors used microarray RNA data for the study. The platform was Affymetrix Human X3P Array.

### 3- What analysis did they performed?
Microarray analysis, datasets generated using total unamplified RNA.

### 4- What were the conclusions of the study?
Previous tissue-based studies have reported similarities in gene expression between normal fallopian tube epithelium and papillary serous carcinoma. And we observed that >80% of serous carcinomas in two independent ovarian cancer gene expression datasets were classified as FT-like using the FNE versus OCE cell-of-origin signature. This result provides that a large proportion of high grade serous carcinomas may arise from the non-ciliated fallopian tube epithelium.

## Datasets
1.	Merritt MA: analyzed 12 samples from two donor patients and established cultures of both ovarian epithelium and fallopian tube epithelium (hTERT immortalized), each with 3 replicates (different culture passages)
2.	Wu R: 99 individual ovarian tumors (37 endometrioid, 41 serous, 13 mucinous, and 8 clear cell carcinomas) and 4 individual normal ovary samples. RNA expression was analyzed using one Affymetrix HG_U133A array per sample.
3.	Tothill R: Randomly selected 285 samples from the AOCS ( Australian Ovarian Cancer Study) was expression profiled on the affymetrix U133_plus2 platform to identify novel subtypes of ovarian tumor

## Methods
1.	Load Dataset
2.	Log2 transform
3.	Limma (Identify different expression genes)
4.	duplicateCorrelation (identify the correlation samples derived from the same patient)
5.	write.fit(write the results to a file)
6.	Three different functions to find P-Value
7.	Survival analysis

## Results
1.	Normalized Expression Plot

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/Plot2.png)

2.	Volcano Plot

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/Volcano_Plot%20.png)


3.	Heatmap from Wu R datasets: 99 individual ovarian tumors.

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/HeatMap.png)

4.  Performs a principal components analysis on the given data matrix and returns the results as an object of class prcomp.

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/PCA.png)


5.	Cluster genes using hierarchical clustering and k-means

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/Cluster%20genes%20.png)

6.	Differences in the disease-free survival and overall survival between OV-like and FT-like subgroups in the Tothill dataset.

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/survival_analysis.png)

![My image](https://github.com/ShuWang00/Bioinformatics_project/blob/master/Project/Figures/SurvivalPlot.png)


## Discussion

As we can see from the plots above, The possibility that tumor behavior is influenced by the intrinsic characteristics of the normal cell-of-origin in which the cancer-promoting mutations emerge has been raised. we observed that >80% of serous carcinomas in two independent ovarian cancer gene expression datasets were classified as fallopian tube (FT)-like using the FNE versus OCE cell-of-origin signature. This result provides independent support for the hypothesis that a large proportion of high grade serous carcinomas may arise from the non-ciliated fallopian tube epithelium.

Although based on small numbers we observed that ≥85% of mucinous and clear cell tumors were classified as ovary (OV)-like, suggesting that these tumors may arise from ovarian epithelium. In contrast, endometrioid adenocarcinomas had a broader spectrum of phenotypes. Based on the cell-of-origin signature in the Wu dataset, 67% of endometrioid cancers were classified as OV-like. In contrast, 65% of endometrioid cancers were in the FT-like subgroup in the Tothill dataset. These observations suggest multiple candidate cells-of-origin for endometrioid tumors.

## Reference
1. Merritt MA, Bentink S: "Gene signatures of normal hTERT immortalized ovarian epithelium and fallopian tube epithelium", PLoS One. 2013; 8(11): e80314. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3841174/#pone.0080314.s005
2. Wu R1, Hendrix-Lucas N: "Mouse model of human ovarian endometrioid adenocarcinoma based on somatic defects in the Wnt/beta-catenin and PI3K/Pten signaling pathways", Cancer Cell. 2007 Apr;11(4):321-33.https://www.ncbi.nlm.nih.gov/pubmed/17418409
3. Tothill RW1, Tinker AV: "Novel molecular subtypes of serous and endometrioid ovarian cancer linked to clinical outcome", Clin Cancer Res. 2008 Aug 15;14(16):5198-208. doi: 10.1158/1078-0432.CCR-08-0196. https://www.ncbi.nlm.nih.gov/pubmed/18698038
4. Differential gene expression analysis I (microarray data), https://biomedizin.unibas.ch/fileadmin/DKBW/redaktion/Group_Directories/Bioinformatics/IntroBioc2016/03_MicroarrayFromCEL_html.html
5. GEO: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE37648
