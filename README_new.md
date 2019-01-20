# cancerCellNet (CCN)

[Shortcut to Setup CCN](#setup_ccn)
[Shortcut to Broad Training CCN](#broadTrain_ccn)
[Shortcut to Subclass Training CCN](#subTrain_ccn)

#### <a name="setup_ccn">Setting up CCN</a>
```R
library(devtools)
install_github("pcahan1/cancerCellNet", ref="master", auth="your_token_here")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("randomForest")
install.packages("ggplot2")

```

### Fetch the required files if you have not already donwloaded them
```R
download.file("https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/ccn_classifier_Jun_29_2018.rda", "ccn_classifier_Jun_29_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/expHeldOut_Jun_30_2018.rda", "expHeldOut_Jun_29_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/stHeldOut_Jun_30_2018.rda", "stHeldOut_Jun_29_2018.rda")

```

### <a name="broadTrain_ccn">Broad Class Training</a>
```R
library(cancerCellNet)
expGDC = utils_loadObject("expGDC_raw_20190118.rda")
exampleSampTab = utils_loadObject("ExampleSampTab_20190118.rda")
iGenes = utils_loadObject("intersecting_genes_20190120.rda")

exampleSampTab
```
Load in the training data, example sample table and intersecting genes. 

![](md_img/sampTab_20190120.PNG)

The sample table shown above maps each sample in the training data to a broad class and a sub class. 
```R
stList = splitCommon(exampleSampTab, ncells=40, dLevel = "BroadClass")
stTrain = stList[[1]]

expTrain = expGDC[iGenes, as.vector(stTrain$SampleBarCodes)]

returnBroad = broadClass_train(stTrain = stTrain, expTrain = expTrain, colName_cat = "BroadClass", colName_samp = "SampleBarCodes")
save(returnBroad, file="BroadClassifier_return.rda")
```
Randomly select 40 training samples in each broad category using "splitCommon" function and pass into "broadClass_train" for training a broad classifier. 

```R
names(returnBroad)
[1] "expTnorm"    "sampTab"     "cgenes_list" "cnProc"
```
In the returnBroad list, there are 4 items. "expTnorm" is the normalized expression matrix for the training samples. "sampTab" is the sample table of training samples. "cgenes_list" is a named list with all the genes belonging in each category. "cnProc" is a list that contains various components needed for prediction including the classifier object. 

### <a name="subTrain_ccn">Subclass Training</a>
To start the subclass training, you need the cnProc from the broad class training. 

```R
stList_sub = splitCommon(exampleSampTab, ncells=20, dLevel="SubClass")
stTrain_sub = stList_sub[[1]]

expTrain_sub = expGDC[iGenes, as.vector(stTrain_sub$SampleBarCodes)]
cnProc_broad = returnBroad$cnProc

returnSubClass = subClass_train(cnProc_broad = cnProc_broad, stTrain = stTrain_sub, expTrain = expTrain_sub, colName_broadCat = "BroadClass", colName_subClass = "SubClass", name_broadCat = "TCGA-BRCA", colName_samp="SampleBarCodes")
save(returnSubClass, file="SubClassifier_return.rda")
```





