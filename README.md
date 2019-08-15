# cancerCellNet (CCN)

[Shortcut to Setup CCN](#setup_ccn)

[Shortcut to Broad Training CCN](#broadTrain_ccn)

[Shortcut to Broad Validation CCN](#broadVal_ccn)

[Shortcut to Subclass Training CCN](#subTrain_ccn)

[Shortcut to Subclass Validation CCN](#subVal_ccn)

[Shortcut to Application of CCN](#app_BRCA)

[Shortcut to Application of Broad Classifier on BRCA CCLE data](#app_broad)

[Shortcut to Gene Visualization/Comparison Tool](#gene_comp)

[Shortcut to subclass application on BRCA CCLE data](#sub_app)

[Shortcut to class application on PAAD GEMM data](#GEMM_geneConvert)

### <a name="setup_ccn">Set up CCN</a>
```R
library(devtools)
install_github("pcahan1/cancerCellNet", ref="master", auth="your_token_here")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("randomForest")
install.packages("ggplot2")
```

#### Fetch the required files if you have not already donwloaded them
```R
- Ask Patrick to upload all the files that needs to be uploaded 
```

### <a name="broadTrain_ccn">Broad Class Training</a>
Load in the necessary files first. 
```{R}
library(cancerCellNet)
expGDC = utils_loadObject("Named_expGDC_20181218.rda")
stGDC = utils_loadObject("Named_stGDC_20181218.rda")

CCLE_sample = utils_loadObject("CCLE_UCEC.rda")
GEMM_sample = utils_loadObject("GEMM_UCEC.rda")
PDX_sample = utils_loadObject("PDX_UCEC.rda")
```
#### Find intersecting genes between query samples and training samples
```{R}
iGenes = Reduce(intersect, list(rownames(CCLE_sample), rownames(GEMM_sample), rownames(PDX_sample), rownames(expGDC)))
```
#### Split TCGA data into training and validation 
```{r}
stList = splitCommon_proportion(sampTab = stGDC, proportion = 0.66, dLevel = "project_id")
stTrain = stList$trainingSet
expTrain = expGDC[, stTrain$barcode]
```
#### Train the broad class classifier 
Because the training data is not balanced in this case, we would have to use stratified sampling in this case. The samplesize parameter indicates the samplesize of stratified sampling. Additionally, because the process of gene pair transform is resource intensive and time consuming, we developed a modified method to perform quick pair transform that is much quicker. 
```{R}
broad_return = broadClass_train(stTrain = stTrain, 
                                expTrain = expTrain[iGenes, ], 
                                colName_cat = "project_id", 
                                colName_samp = "barcode", 
                                nRand = 70,
                                nTopGenes = 25, 
                                nTopGenePairs = 70, 
                                nTrees = 2000, 
                                stratify=TRUE, 
                                sampsize=60, 
                                quickPairs=TRUE)

save(broad_return, file="BroadClassifier_return.rda")
```
### <a name="broadVal_ccn">Validate Broadclass classifier</a>
#### Classify the validation set
```{R}
stVal_Broad = stList$validationSet
stVal_Broad_ord = stVal_Broad[order(stVal_Broad$project_id), ] #order by broadClass
expVal_Broad = expGDC[iGenes, rownames(stVal_Broad_ord)]
cnProc_broad = broad_return$cnProc #select the cnProc from the broadclass training earlier 

classMatrix_broad = broadClass_predict(cnProc_broad, expVal_Broad, nrand = 60)
```

#### Visualize the validation classification 
```{R}
stValRand_broad = addRandToSampTab(classMatrix_broad, stVal_Broad_ord, "project_id", "barcode")
grps = as.vector(stValRand_broad$project_id)
names(grps)<-rownames(stValRand_broad)
ccn_hmClass(classMatrix_broad, grps=grps, fontsize_row=10)
```
![](md_img/TCGA_validate_heatmap.png)

#### Adding gaps between groups for more clear visualization 
```{R}
breakVector = c() # create a vector of number indicating the column at which the gap will be placed 
for (uniqueClass in unique(grps)) {
   myBreak = max(which(grps %in% uniqueClass))
   breakVector = c(breakVector, myBreak)
}
ccn_hmClass(classMatrix_broad, grps=grps, fontsize_row=10, gaps_col = breakVector) 
```
![](md_img/TCGA_validate_heatmap_split.png)


#### Classifier Assessment 
```{R}
assessmentDat = ccn_classAssess(classMatrix_broad, stValRand_broad, "project_id","barcode")
plot_class_PRs(assessmentDat)
```
![](md_img/TCGA_PR.png)








