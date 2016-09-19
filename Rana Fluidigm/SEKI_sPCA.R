## sPCA for Haplotypes ##
## 09/19/16 Andy Rothstein ##
setwd("~/SEKI_YOSE/cohortvcf/")
require(adegenet)
require(ade4)
require(ggplot2)
require(ggmap)
rm(list=ls())
## bring in hap data
datSeki = read.table("~/SEKI_YOSE/cohortvcf/data.seki.pops.txt", header = T, stringsAsFactors = F, sep="\t")
## determine % missing data
datSeki$PercMissing = apply(datSeki[,-c(1:9)], 1, function(x) length(which(is.na(x)))/(ncol(datSeki)-9))
datSeki = datSeki[,c(1:9,47,10:46)]
hist(datSeki$PercMissing, breaks=50)
# filter based on % missing
SEKI = datSeki[which(datSeki$PercMissing < 0.5),]

# use adegenet
library(adegenet)
# library(wordcloud)
gen = df2genind(as.matrix(SEKI[,-c(1:10)]), sep = ",", ploidy = 2, pop = as.factor(SEKI$pops))
gen$pop

# adegenet PCA
X <- tab(gen, NA.method="mean")