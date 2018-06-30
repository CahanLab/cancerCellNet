# cancerCellNet (cCN)


This is a cursory walk-thru for using cancerCellNet. In addition to the R package, you will also need a few files:

- (classifier list)[https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/ccn_classifier_Jun_29_2018.rda]
- (held out valdiation expression data)[https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/expHeldOut_Jun_30_2018.rda]
- (held out valdiation meta data table)[https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/stHeldOut_Jun_30_2018.rda]


#### Setup
```R
library(devtools)
install_github("pcahan1/cancerCellNet", ref="master", auth="your_token_here")
library(cancerCellNet)

library(RColorBrewer)
library(pheatmap)
library(randomForest)
library(ggplot2)
```

### Fetch the required files if you have not already donwloaded them
#### Fetch the data (optional if you have alread done this)
```R
download.file("https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/ccn_classifier_Jun_29_2018.rda", "ccn_classifier_Jun_29_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/expHeldOut_Jun_30_2018.rda", "expHeldOut_Jun_29_2018.rda")

download.file("https://s3.amazonaws.com/cnobjects/cancerCellNet/resources/stHeldOut_Jun_30_2018.rda", "stHeldOut_Jun_29_2018.rda")

```

#### Load Classifier list and held-out validation data
```R
mydate<-utils_myDate()
ccnList<-utils_loadObject("ccn_classifier_Jun_29_2018.rda")
rf_tsp<-ccnList[['classifier']]
cgenes<-ccnList[['cgenes']]
xpairs<-ccnList[['xpairs']]

expVal<-utils_loadObject("expHeldOut_Jun_30_2018.rda")
stVal<-utils_loadObject("stHeldOut_Jun_30_2018.rda")
```

#### Transform the query/validation data
```R
expValTrans<-query_transform(expVal[cgenes,], xpairs)
```

#### Classify the query/validation data, and add some randomized profiles, too
```R
nrand<-50
classRes_val<-rf_classPredict(rf_tsp, expValTrans, numRand=nrand)
```

The results of the ^^ are different than the traditional CellNet in that they only (for now) return a classification matrix.

#### Plot the classification results
```R
stValRand<-addRandToSampTab(classRes_val, stVal, "description2", "sample_name")
grps<-as.vector(stValRand$description2)
names(grps)<-rownames(stValRand)
ccn_hmClass(classRes_val, grps=grps)
####

![](md_img/hmClass_val_Jun_29_2018.png)


Use cnn_classAssess() and plot_class_PRs() to assess the performance of this classifier.

For the CCLE analysis, write a function to reorder the classRes columns in increasing order of target tumor category.











