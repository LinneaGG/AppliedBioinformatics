# change file for other performance measures and softwares
concordance <- read.csv("all_values_cross_entropy_phaser.csv", header = FALSE)

colnames(concordance) <- c('uncertainty', 'concordance')

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     quant25 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.25)),
                     quant75 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.75)),
                     quant10 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.10)),
                     quant90 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.90)),
                     quant00 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.00)),
                     quant01 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.01)),
                     quant99 = quantile (xx[[col]], na.rm=na.rm, probs=c(0.99)),
                     quant100 = quantile (xx[[col]], na.rm=na.rm, probs=c(1.00))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac <- rename(datac, c("quant25.25%" = "quant25"))
  datac <- rename(datac, c("quant10.10%" = "quant10"))
  datac <- rename(datac, c("quant75.75%" = "quant75"))
  datac <- rename(datac, c("quant90.90%" = "quant90"))
  datac <- rename(datac, c("quant00.0%" = "quant00"))
  datac <- rename(datac, c("quant01.1%" = "quant01"))
  datac <- rename(datac, c("quant99.99%" = "quant99"))
  datac <- rename(datac, c("quant100.100%" = "quant100"))
  
  return(datac)
}

performance <- summarySE(concordance, measurevar="concordance", groupvars=c("uncertainty"))

library(ggplot2)
library(ggthemes)
ggplot(performance, aes(x=uncertainty, y=concordance)) +
  geom_line(aes(x=uncertainty, y=quant25, color='25% to 75% percentile'), size=1, alpha=0.8)+
  geom_line(aes(x=uncertainty, y=quant01, color="1% to 99% percentile"), size=1, alpha=0.8)+
  geom_line(aes(x=uncertainty, y=quant75, color="25% to 75% percentile"), size=1, alpha=0.8)+
  geom_line(aes(x=uncertainty, y=quant99, color="1% to 99% percentile"), size=1, alpha=0.8)+
  geom_line(aes(x=uncertainty, y=quant00, color="0% to 100% percentile"), size=1, alpha=0.8)+
  geom_line(aes(x=uncertainty, y=quant100, color="0% to 100% percentile"), size=1, alpha=0.8)+
  geom_ribbon(aes(ymin=quant25,ymax=quant75), fill="orange", alpha=0.7) +
  geom_ribbon(aes(ymin=quant01,ymax=quant99), fill="orange", alpha=0.45) +
  geom_ribbon(aes(ymin=quant00,ymax=quant100), fill="orange", alpha=0.2) +
  geom_line(aes(x=uncertainty, y=concordance, color="Mean concordance"), size=1.5, alpha=0.8)+
  labs(color="")+
  xlab('Genotype likelihood')+
  ylab('Cross entropy')+
  theme_minimal()+
  theme(axis.title=element_text(), text=element_text(size=20))+
  ggtitle("Cross entropy across genotype likelihoods")+
  labs(subtitle="Samples imputed with Phaser")+
  scale_color_brewer(palette="Oranges")
