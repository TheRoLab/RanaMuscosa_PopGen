### PCA and IBD for Haplotypes and SNPs ###
### Fall 2016 Andy Rothstein ###

## Read in haplotype data## 
rm(list=ls())
datSeki = read.table("~/SEKI_YOSE/cohortvcf/data.seki.pops.txt", header = T, stringsAsFactors = F, sep="\t")
## determine % missing data
datSeki$PercMissing = apply(datSeki[,-c(1:9)], 1, function(x) length(which(is.na(x)))/(ncol(datSeki)-9))
datSeki = datSeki[,c(1:9,47,10:46)]
hist(datSeki$PercMissing, breaks=50)
# filter based on % missing
SEKI = datSeki[which(datSeki$PercMissing < 0.5),]

#write.table(SEKI, file = "./haps.strucutre.txt", sep = "\t", col.names = T, row.names = F)

#plot lat/lon map of populations
library(ggmap)
dev.off()
# mean lat long for orientation
location_SEKI <- c(lon=-118.5086, lat=36.85772)
map = get_map(location =location_SEKI, source = "google", maptype = "hybrid", zoom = 9)

# map plot for grid
map.plot = ggmap(map) + geom_jitter(aes(x = SEKI$Longitude, y = SEKI$Latitude, 
                                        color=as.factor(SEKI$pops)), 
                                    data = SEKI , size = 2) + guides(color=guide_legend(title="Population"))
## Start PCA Haplotypes ##
# use adegenet
library(adegenet)
# library(wordcloud)
gen = df2genind(as.matrix(SEKI[,-c(1:10)]), sep = ",", ploidy = 2, pop = as.factor(SEKI$pops))
gen$pop

# adegenet PCA
X <- tab(gen, NA.method="mean")

## make PCA
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)

## basic plot
plot(pca1$li, cex=1, pch=16)
# plot(pca1$li, col=myCol, cex=1, pch=16, xlim=c(-8,8),ylim=c(-6,6))
## use wordcloud for non-overlapping labels
# textplot(pca1$li[,1], pca1$li[,2], words=pops, cex=1, new=FALSE)
## legend the axes by adding loadings
abline(h=0,v=0,col="grey",lty=2)
#dev.off()

# eigenvalues
pca1$eig[1:3]
# eigenvalues - proportion of summed eigenvalues. a measure of how much variation explained by the axes. 
#         this is typically reported in the axis labels
round(pca1$eig[1:3] / sum(pca1$eig),3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

# pca with pop clusters
s.class(pca1$li, pop(gen))
add.scatter.eig(pca1$eig[1:20], 3,1,2, posi="bottomright")

#######
## plot PCA with Populations
pca1Plot = pca1$li
pca1Plot = cbind(pca1Plot, SEKI[,1:10])
colnames(pca1Plot)

library(ggplot2)
library(ggrepel)

## by amphib species 
gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=pca1Plot$amphib_species))
gp + theme_bw() + geom_jitter(color="grey80", size = 1) + geom_text(size=2) +
  xlab("PC1") + 
  ylab("PC2") 

## by populations
gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=pops, color=factor(pca1Plot$pops)))
pca.plot = gp + theme_bw() + geom_jitter(color="grey80", size = 1) + geom_text() +
  xlab("PC1") + 
  ylab("PC2") + theme(legend.position="none")

plot_grid(pca.plot, map.plot)

#library(ggplot2)
#library(ggrepel)
#gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=pop_id))
#gp + theme_bw() + geom_point(color="grey80", size = 1) + geom_text(size=2) +
#  xlab("PC1") + 
#  ylab("PC2") 

#gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=drainage_id, color=factor(drainage_id)))
#gp + theme_bw() + geom_point(color="grey80", size = 1) + geom_text() +
#  xlab("PC1") + 
#  ylab("PC2") 

# gp + theme_bw() + geom_point() +
#   xlab("PC1") + 
#   ylab("PC2") +
#   geom_text_repel(
#     size = 3, color="grey60",
#     box.padding = unit(0.35, "lines"),
#     point.padding = unit(0.3, "lines"))


# get centroid coords for pop factor
#pop_centroids = aggregate(cbind(Axis1,Axis2)~amphib_species,pca1Plot,mean)
#gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=drainage_name2, color=as.character(drainage_name2)))

#pop_centroids = aggregate(cbind(Axis1,Axis2)~pops,pca1Plot,mean)
#gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=drainage_name3, color=drainage_name3))
#gp + theme_bw() + geom_jitter(color="grey80") + geom_text(show.legend = F, size=3) +
#  xlab("PC1") + 
#  ylab("PC2") +
#  geom_point(data=pop_centroids, size = 2) + guides(color=guide_legend(title="Drainage Forks"))


# Check where low quality samples are located
gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=lake_id, color=PercMissing))
gp + theme_bw() + geom_jitter() +
  xlab("PC1") + 
  ylab("PC2") 

# Subset by PercMissing?
gp = ggplot(pca1Plot[which(pca1Plot$PercMissing < .8),], aes(Axis1, Axis2, label=pca1Plot$pops, color=PercMissing))
gp + theme_bw() + geom_point() +
  xlab("PC1") + 
  ylab("PC2") 

#### IBD Haplotypes ####
rm(list = ls())
setwd("~/SEKI_YOSE/cohortvcf/")
load("vcf_structure/data.seki.pops.robj")
#load("seki.data.plink.Robj")

## plot lat lon
col.scheme = as.factor(data.seki.pops$pops)
plot(jitter(data.seki.pops$Longitude,amount=0.01), xlab="Longitude",
     jitter(data.seki.pops$Latitude,amount=0.01), ylab="Latitude",
     pch=19,col=col.scheme)
# bring in required packages and build genid/genpop for haplotypes
require(adegenet)
require(PopGenReport)
require(dismo)
gen = df2genind(as.matrix(data.seki.pops[, -c(1:9)]), sep = ",", pop=as.factor(data.seki.pops$pops), ploidy = 2)
latlong = data.seki.pops[,c(5,7,8)]
gen$other$longlat = latlong[, 2:1]
gen$other$xy = Mercator(gen$other$longlat)
gen.pop = genind2genpop(gen)

# calculate genetic distance var
Dgen = dist.genpop(gen.pop, method=1)

## read in lat long for pop ##
pop.latlong = aggregate(latlong, list(latlong$pops), mean)
pop.latlong = pop.latlong[,c(2,3,4)]
gen.pop$other$longlat = pop.latlong[, 3:2]
gen.pop$other$xy = Mercator(gen.pop$other$longlat)

# calculate Euclidean distance var
Dgeo = dist(gen.pop$other$xy, method = "euclidean")

# calculate IBD
ibd = mantel.randtest(Dgen, Dgeo)
# plot IBD
plot(ibd)

# Geo/Gen plot
col.scheme = as.factor(data.seki.pops$pops)
plot(Dgeo, Dgen, pch=19, xlab="Geographic Distance (Euclidean)", ylab="Genetic Distance(Nei's)")
abline(lm(Dgen~Dgeo), col="blue", lty=2)

# calculate pairwise Fst
library(hierfstat)
fst.seki = pairwise.fst(gen, res.type = c("dist", "matrix"))
seki.fst.matrix = as.matrix(fst.seki)

# plot Geo/Fst
plot(Dgeo, fst.seki)

geoDist.seki <- fields::rdist.earth(cbind(pop.latlong$Longitude, pop.latlong$Latitude),miles=FALSE)
index.matrix <- upper.tri(geoDist.seki,diag=TRUE)
plot(geoDist.seki[index.matrix],seki.fst.matrix[index.matrix], xlab="Geo Distance (km)", ylab ="Fst")
abline(lm(seki.fst.matrix~geoDist.seki), col="blue", lty=2)

################################################################################################
## VCF PCAs and IBD ## 
## 09/07/16 Andy Rothstein ##
require(adegenet)
require(VariantAnnotation)
require(SNPRelate)
require(pegas)
setwd("~/SEKI_YOSE/cohortvcf/")
rm(list = ls())

vcf= readVcf("cohort.snp.filtered.subset.vcf","MYLFTargets")
mat <- genotypeToSnpMatrix(vcf)
mat = mat$genotypes
meta = read.table("seki_meta.txt", header = T, sep = "\t")

mat.summ = col.summary(mat)
hist(mat.summ$MAF, breaks=500, xlim=c(0,.02))
hist(mat.summ$Call.rate, breaks=500)
table(mat.summ$MAF >= .006); table(mat.summ$MAF >= .007)
table(mat.summ$Calls == 0)
plot(mat.summ$Call.rate, mat.summ$MAF)

vcf = vcf[-which(mat.summ$Calls == 0),]
dim(vcf)

mat <- genotypeToSnpMatrix(vcf)
mat = mat$genotypes
mat.summ = col.summary(mat)
hist(1/(mat.summ$Calls*2))
hist(mat.summ$MAF, breaks=10000, xlim=c(0.002,.004))
hist(mat.summ$MAF)
dim(vcf)

mat <- genotypeToSnpMatrix(vcf)
mat = mat$genotypes
mat.rowsumm = row.summary(mat)
hist(mat.rowsumm$Call.rate, breaks=20)
mean(mat.rowsumm$Call.rate);median(mat.rowsumm$Call.rate); 
vcf = vcf[,which(mat.rowsumm$Call.rate > .6)]
dim(vcf)

#cutSamples = c("RKS12265", "RKS11739")
#if(!all(is.na(match(cutSamples,colnames(vcf))))){ 
#  vcf = vcf[,-na.omit(match(cutSamples,colnames(vcf)))] 
#}               
#dim(vcf)

meta = meta[match(colnames(vcf),meta$swab_id),]
coord = meta[match(colnames(vcf),meta$sample_ID),c(12,13)]
pop = meta$pop_id[match(colnames(vcf),meta$swab_id)]

## setup vcf matrix for genid
snp = t(geno(vcf)$GT)
dim(snp)
snp = gsub("\\/","",snp)
snp[snp=="."] = "NA"
rownames(snp) = rownames(t(geno(vcf)$GT))
snp[2,2]
dim(snp)

## make genid
snptogenind = snp
snptogenind = as.matrix(snptogenind)
colnames(snptogenind) = 1:ncol(snptogenind)
gen <- df2genind(snptogenind, ploidy=2, pop=meta$pop_id,sep="/",NA.char="NA")
gen.pop = genind2genpop(gen)

# incorporate lat/lon meta
meta.latlong = meta[,c(2,12,13,15)]

# aggregate avg lat lon for each pop (drainage)
drainage.latlong = meta.latlong[,2:4]
drainage.latlong = aggregate(drainage.latlong, list(drainage.latlong$pop), mean)
drainage.latlong = drainage.latlong[,c(2:4)]
data.seki.geoavg = merge(meta, drainage.latlong, by="pop_id")
data.seki.geoavg = data.seki.geoavg[,c(1,3,16,17)]
data.seki.geoavg = rename(data.seki.geoavg, c("Latitude.y"="lat", "Longitude.y"="long"))

# add lat lon to genid
latlong.pop = drainage.latlong[,c(2,1)]
latlong = data.seki.geoavg[,c(4,3)]
gen.pop$other$longlat = latlong.pop
gen.pop$other$xy = Mercator(gen.pop$other$longlat)

gen$other$longlat = latlong
gen$other$xy = Mercator(gen$other$longlat)

# pca of genid
X <- tab(gen, NA.method="mean")
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
plot(pca1$li, cex=1, pch=16)
abline(h=0,v=0,col="grey",lty=2)

pca1$eig[1:3]
round(pca1$eig[1:3] / sum(pca1$eig),3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

s.class(pca1$li, pop(gen))
add.scatter.eig(pca1$eig[1:20], 3,1,2, posi="bottomright")

# bind PCA val with meta
pca1Plot = pca1$li
pca1Plot = cbind(pca1Plot, meta)
colnames(pca1Plot)

require(ggplot2)
require(ggrepel)
require(cowplot)
## by pop
gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=pca1Plot$pop_id))
gp + theme_bw() + geom_jitter(color="grey80", size = 1) + geom_text(size=8) +
  xlab("PC1") + 
  ylab("PC2") 

gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=pca1Plot$pop_id, color=factor(pca1Plot$pop_id)))
gp + theme_bw() + geom_jitter(color="grey80", size = 4) + geom_text(size=8) +
  xlab("PC1") + 
  ylab("PC2") + guides(color=guide_legend(title="Population"))

library(ggmap)
# mean lat long for each location
location_SEKI <- c(lon=-118.5086, lat=36.85772)
map = get_map(location =location_SEKI, source = "google", maptype = "hybrid", zoom = 9)

map.plot = ggmap(map) + geom_jitter(aes(x = meta$Longitude, y = meta$Latitude, 
                                        color=as.factor(pca1Plot$pop_id)), 
                                    data = meta, size = 1) + guides(color=guide_legend(title="Population"))

gp = ggplot(pca1Plot, aes(Axis1, Axis2, label=pca1Plot$pop_id, color=factor(pca1Plot$pop_id)))
pca.plot = gp + theme_bw() + geom_jitter(color="grey80", size = 1) + geom_text() + 
  xlab("PC1") + ylab("PC2") + theme(legend.position="none")
plot_grid(pca.plot, map.plot)

################
