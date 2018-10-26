#' Plot results of cn_classAssess
#'
#' Plot one precision recall curve per CT
#' @param assessed result of runnung cn_classAssess
#'
#' @return ggplot pbject
#'
#' @examples
#' testAssTues<-cn_splitMakeAssess(stTrain, expTrain, ctGRNs, prop=.5)
#' plot_class_PRs(testAssTues$ROCs)
#'
#' @export
plot_class_PRs<-function
(assessed
){
  ctts<-names(assessed)
  df<-data.frame()
  for(ctt in ctts){
    tmp<-assessed[[ctt]]
    tmp<-cbind(tmp, ctype=ctt)
    df<-rbind(df, tmp)
  }

  ggplot(data=df, aes(x=Sens, y=Prec)) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
    theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
    theme(axis.text = element_text(size=5)) + ggtitle("Classifier performance")
}
