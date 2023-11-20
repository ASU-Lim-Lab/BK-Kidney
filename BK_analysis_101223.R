setwd("/Users/rmaqsood/Documents/BK/virome")
library(nlme)
library(lme4)
require(car)
require(MASS)
library(lmerTest)
library(RColorBrewer)
library(gplots)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)
library(vegan)
library(tidyverse)
library(tidyr)
library(nlme)
library(phyloseq)
library(ape)
library(corrplot)
library(taxonomizr)
#library(VirFinder)
library("qvalue")
library(glmnet)
library(Rcpp)
library(decontam)

####Taxonomy_blastn
library(taxonomizr)
setwd('/Users/rmaqsood/Documents/')
getNamesAndNodes()
read.names.sql('names.dmp','accessionTaxa.sql', overwrite = TRUE)
read.nodes.sql('nodes.dmp','accessionTaxa.sql',overwrite = TRUE)
blastResults<-read.table("/Users/rmaqsood/Documents/BK/virome/finalViralContigsCandidate_viralNRBlastx.out",header=FALSE,stringsAsFactors=FALSE)
summary(blastResults)
accessions<-blastResults$V2
taxaId<-accessionToTaxa(accessions,"accessionTaxa.sql")
taxonomy<-getTaxonomy(taxaId,'accessionTaxa.sql')
write.csv(taxonomy,"/Users/rmaqsood/Documents/BK/virome/finalViralContigsCandidate_viralNRBlastx_taxonomy.csv")

set.seed(100)
data <- read.delim("viromeRawCountsTable_090123.txt", row.names=1)
data<-as.matrix(t(data))
metadata<-read.table("metadata_preDecontam.txt", sep = "\t", header = T, row.names = 1, stringsAsFactors = TRUE, comment.char = "")
contam<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.1)
contam.intermediate<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.25)
contam.strict<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.5)
write.csv(contam, "LK_SpeciesDefault_090123.csv")
write.csv(contam.intermediate,"LK_SpeciesIntermediate_090123.csv")
write.csv(contam.strict,"LK_SpeciesStrict_090123.csv")

d.1<-read.delim("BK_cleanNormalized_counts.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
data_vegan<-diversity(dataTransposed, index = 'shannon')
write.table(data_vegan,"BK_ViromeAlphaDiversity_092623.txt", sep = '\t')

#gender differences
setwd("/Users/rmaqsood/Documents/BK/virome/")
m.1 <- read.delim("Analysis_BK_092623.txt", header=TRUE, row.names = 1)
m.1$Patient.ID<-as.character(m.1$Patient.ID)
#Gender significance
set.seed(100)
allrichness<-lme(Richness~ Gender*Time + Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ Gender*Time + Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)

set.seed(100)
allrichness<-lme(Richness~Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)

#difference by BKV status
set.seed(100)
allrichness<-lme(Richness~Status*Time+Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~Status*Time+Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)

#split BKV
m.1 <- read.delim("Analysis_BK_092623.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$Patient.ID)
summary(m.1)
BKV.virome<-m.1[ which(m.1$Status=='BK+'),]
Control.virome<-m.1[ which(m.1$Status=='Control'),]

set.seed(100)
allrichness<-lme(Richness~Gender+Time, random=~1|Patient.ID, data=BKV.virome)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~Gender + Time, random=~1|Patient.ID, data=BKV.virome)
summary(allalphaDiversity)

set.seed(100)
allrichness<-lme(Richness~Gender+Time, random=~1|Patient.ID, data=Control.virome)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~Gender + Time, random=~1|Patient.ID, data=Control.virome)
summary(allalphaDiversity)

loessPlotRichness <- ggplot(m.1, aes(x = Time, y = Richness, color = Status, fill = Status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Time), max(m.1$Time), by =1 ),1)) + labs (x= "Post-transplant week",y="Virome Richness", title = "Females vs males" )
loessPlotRichness

loessPlotShannon <- ggplot(m.1, aes(x = Time, y = AlphaDiversity, color = Status, fill = Status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Time), max(m.1$Time), by =1 ),1)) + labs (x= "Post-transplant week",y="Virome Alpha Diversity", title = "Females vs males" )
loessPlotShannon

m.1 <- read.delim("Analysis_BK_092623.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$Patient.ID)
d.1<-read.delim("BK_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
beta_dist <- vegdist(t(d.1),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetavirome_BKFemaleMale_092623.csv")

res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC3 <- as.matrix(res$vectors[,3])
PC1
PC2
PC3
data<-read.delim("PC_femaleMale.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
data$PC3 <- PC3

#personCode
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Gender))) + 
  geom_point() + xlim(-.6, .8) + ylim(-.6, .8) +
  stat_ellipse() +geom_point(size=2.5)+
  scale_color_manual(breaks = c("female","male"),
                     values = c("salmon", "darkgreen")) +
  theme(aspect.ratio=1)+labs(title="All samples by gender")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#BKV
p <- ggplot(data, aes(x=PC1, y=PC3, color=as.factor(Status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +geom_point(size=2.5)+
  scale_color_manual(breaks = c("BK+","Control"),
                     values = c("#b2182b", "#045a8d")) +
  theme(aspect.ratio=1)+labs(title="All samples by BKV status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p <- ggplot(data, aes(x=PC1, y=PC3, color=(Time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All sample over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# 
sub_info <- cbind( smp = m.1$Sample, # sample ID
                   id = m.1$Patient.ID, # subject ID 
                   tp = m.1$Time,
                   status = m.1$Status,
                   gender=m.1$Gender)
sub_info<-as.data.frame(sub_info)

WUDM<-read.csv("weightedBetavirome_BKFemaleMale_092623.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA

sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)

ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)

set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))#, plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   gender + tp +status , data = sub_clust, permutations = perm, by = "margin")

set.seed(100)
adonis2(wudm ~   gender*tp + gender + tp +status , data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   gender*status + gender + tp +status , data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   status*tp + gender + tp +status , data = sub_clust, permutations = perm, by = "margin")

library(reshape2)
df <- melt(as.matrix(beta_dist))
p    <- t(apply(df[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])
p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))
df   <- df[-c(rmv1,rmv2),] 
write.csv(df,"Pair_WeightedBeta_v2.csv")


##################Maaslin2##################
setwd('/Users/rmaqsood/Documents/BK/virome/')
m.1 <- read.delim("Analysis_BK_092623.txt", header=TRUE)
d.1<-read.delim("BK_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
# m.2<-m.1
# #female Virome
# female.virome<-m.2[ which(m.2$Gender=='female'),]
# d.2<-d.1[,colnames(d.1) %in% female.virome$Sample]
# write.table(female.virome,"female.virome_metadata.txt")
# write.table(d.2,"female.virome_Relabun.txt")
# #male
# male.virome<-m.2[ which(m.2$Gender=='male'),]
# d.2<-d.1[,colnames(d.1) %in% male.virome$Sample]
# write.table(male.virome,"male.virome_metadata.txt")
# write.table(d.2,"male.virome_Relabun.txt")

library(Maaslin2)
input_data <- read.delim("BK_cleanNormalized_relativeAbundance.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("Analysis_BK_092623.txt", row.names = 1)
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_all_101223', transform = "NONE",
  fixed_effects = c('Time','Status','Gender'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)

library(Maaslin2)
input_data <- read.delim("male.virome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("male.virome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_Male_092823', transform = "NONE",
  fixed_effects = c('Time','Status'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)

library(Maaslin2)
input_data <- read.delim("female.virome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("female.virome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_female_092823', transform = "NONE",
  fixed_effects = c('Time','Status'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)


##################k-means##################
library(factoextra)
library(cluster)
setwd('/Users/rmaqsood/Documents/BK/virome/')
m.1 <- read.delim("Analysis_BK_092623.txt", header=TRUE)
# m.1$Patient.ID<-as.character(m.1$Patient.ID)
# d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
# write.csv(d.1,"LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance_033023.txt")
# dataTransposed<-t(d.1)
# m.2<-m.1[!m.1$Richness=="#N/A",]
# d.2<-d.1[,colnames(d.1) %in% m.1$Sample]
# beta_dist <- vegdist(t(d.2),index = "bray")
# distMatrix <- as.matrix(beta_dist)
# write.csv(distMatrix,"weightedBetaBacteriomeMotherInfants_033023.csv")

#calculate optimal number of clusters (sum of squares)
WUDM<-read.csv("weightedBetavirome_BKFemaleMale_092623.csv",row.names=1,header = T)
a<-fviz_nbclust(WUDM, pam, method = "wss")
#plot optimal number of clusters (sum of squares)
a
b<-fviz_nbclust(WUDM, pam, method = "silhouette")
b
#plot number of clusters vs. gap statistic
gap_stat <- clusGap(WUDM,
                    FUN = pam,
                    K.max = 20, #max clusters to consider
                    B = 50) #total bootstrapped iterations
fviz_gap_stat(gap_stat)


#set.seed(123)
km.res<-kmeans(WUDM, 2, iter.max = 10, nstart = 25)
km.res
km.res1<-kmeans(WUDM, 2, iter.max = 10, nstart = 25)
km.res1
set.seed(123)
km.res2<-kmeans(WUDM, 4, iter.max = 10, nstart = 25)
km.res2
set.seed(123)
km.res3<-kmeans(WUDM, 5, iter.max = 10, nstart = 25)
km.res3
set.seed(123)
km.res4<-kmeans(WUDM, 6, iter.max = 10, nstart = 25)
km.res4



##################k-meansNoPolyNoBKDom##################
library(factoextra)
library(cluster)
library(vegan)
setwd('/Users/rmaqsood/Documents/BK/virome/')
m.1 <- read.delim("subsetAnalysis_BK_092623.txt", header=TRUE)
m.1$Patient.ID<-as.character(m.1$Patient.ID)
d.1<-read.delim("BK_cleanNormalized_relativeAbundance_noPolyomaNoBKDom.txt", row.names = 1, header = TRUE)
#write.csv(d.1,"LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance_033023.txt")
dataTransposed<-t(d.1)
#m.2<-m.1[!m.1$Richness=="#N/A",]
#d.2<-d.1[,colnames(d.1) %in% m.1$Sample]
beta_dist <- vegdist(t(d.1),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaNoPolyomaNoBKDom_100523.csv")

d.1<-read.delim("subsetNoBK_cleanNormalized_counts_.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
data_vegan<-diversity(dataTransposed, index = 'shannon')
write.table(data_vegan,"subsetNoBK_ViromeAlphaDiversity_100623.txt", sep = '\t')

#gender differences
setwd("/Users/rmaqsood/Documents/BK/virome/")
m.1 <- read.delim("subsetAnalysis_BK_092623.txt", header=TRUE, row.names = 1)
m.1$Patient.ID<-as.character(m.1$Patient.ID)
#Gender significance
set.seed(100)
allrichness<-lme(Richness~ Gender*Time + Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ Gender*Time + Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)

set.seed(100)
allrichness<-lme(Richness~ Status*Time + Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ Status*Time + Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)

set.seed(100)
allrichness<-lme(Richness~Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)


#BKVsubset
m.1 <- read.delim("subsetAnalysis_BK_092623.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$Patient.ID)
summary(m.1)
BKV.virome<-m.1[ which(m.1$Status=='BK+'),]
Control.virome<-m.1[ which(m.1$Status=='Control'),]

set.seed(123)
allalphaDiversity<-lme(AlphaDiversity~Gender + Time, random=~1|Patient.ID, data=BKV.virome)
summary(allalphaDiversity)
set.seed(123)
allalphaDiversity<-lme(AlphaDiversity~Gender + Time, random=~1|Patient.ID, data=Control.virome)
summary(allalphaDiversity)

loessPlotRichness <- ggplot(m.1, aes(x = Time, y = Richness, color = Status, fill = Status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Time), max(m.1$Time), by =1 ),1)) + labs (x= "Post-transplant week",y="Virome Richness", title = " BKV status without BK dominant" )
loessPlotRichness

loessPlotShannon <- ggplot(m.1, aes(x = Time, y = AlphaDiversity, color = Status, fill = Status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Time), max(m.1$Time), by =1 ),1)) + labs (x= "Post-transplant week",y="Virome Alpha Diversity", title = "BKV status without BK dominant" )
loessPlotShannon


#difference by BKV status
set.seed(100)
allrichness<-lme(Richness~Status*Time+Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~Status*Time+Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)

#beta
m.1 <- read.delim("subsetAnalysis_BK_092623.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$Patient.ID)
d.1<-read.delim("weightedBetaNoPolyomaNoBKDom_100523.csv", row.names = 1, header = TRUE, sep = ",")
distMatrix <- as.matrix(d.1)

res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC3 <- as.matrix(res$vectors[,3])
PC1
PC2
PC3
data<-read.delim("PC_NoPolyNoBKDom.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
data$PC3 <- PC3

#KT BKV
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Status))) +
  geom_point() + xlim(-.6, .6) + ylim(-.6, .6) +
  stat_ellipse() +geom_point(size=2.5)+
  scale_color_manual(breaks = c("BK+","Control"),
                     values = c("#b2182b", "#045a8d")) +
  theme(aspect.ratio=1)+labs(title="Females by BKV status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Gender))) +
  geom_point() +xlim(-.6, .6) + ylim(-.6, .6) +
  stat_ellipse() +geom_point(size=2.5)+
  scale_color_manual(breaks = c("female","male"),
                     values = c("salmon", "darkgreen")) +
  theme(aspect.ratio=1)+labs(title="All samples by gender")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

data<-read.delim("PC_NoPolyNoBKDom.txt", header=TRUE)
sub_info <- cbind( smp = data$Sample, # sample ID
                   id = data$PatientID, # subject ID
                   tp = data$Time,
                   status = data$Status,
                   gender=data$Gender)
sub_info<-as.data.frame(sub_info)

WUDM<-read.csv("weightedBetaNoPolyomaNoBKDom_100523.csv", header=T, row.names = 1, sep = ",")
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA

sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)

ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)

set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))#, plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   status*tp+gender + tp +status , data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   gender*tp+gender + tp +status , data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   status*gender+gender + tp +status , data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   gender + tp +status , data = sub_clust, permutations = perm, by = "margin")


#difference by BKV status and Gender
m.1 <- read.delim("subsetAnalysis_BK_092623.txt", header=TRUE)
m.1$Patient.ID<-as.character(m.1$Patient.ID)
d.1<-read.delim("BK_cleanNormalized_relativeAbundance_noPolyomaNoBKDom.txt", row.names = 1, header = TRUE)
m.2<-m.1
#BKV+
BKV.virome<-m.2[ which(m.2$Status=='BK+'),]
d.2<-d.1[,colnames(d.1) %in% BKV.virome$Sample]
write.table(BKV.virome,"BKVNOBKDom.virome_metadata.txt")
write.table(d.2,"BKVNOBKDom.virome_Relabun.txt")
data<-read.delim("BKVNOBKDom.virome_Relabun.txt", header = TRUE,row.names = 1, sep = " ")
dataTransposed<-t(data)
beta_dist <- vegdist(t(data),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBKVBetaNoPolyomaNoBKDom_100523.csv")

data<-read.delim("BKVNOBKDom.virome_metadata.txt", header=TRUE,row.names = 1, sep = " ")
sub_info <- cbind( smp = data$Sample, # sample ID
                   id = data$Patient.ID, # subject ID
                   tp = data$Time,
                   status = data$Status,
                   gender=data$Gender)
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBKVBetaNoPolyomaNoBKDom_100523.csv", header=T, row.names = 1, sep = ",")
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))#, plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   tp+gender , data = sub_clust, permutations = perm, by = "margin")

#BKV-
Control.virome<-m.2[ which(m.2$Status=='Control'),]
d.2<-d.1[,colnames(d.1) %in% Control.virome$Sample]
write.table(Control.virome,"ControlNOBKDom.virome_metadata.txt")
write.table(d.2,"ControlNOBKDom.virome_Relabun.txt")
data<-read.delim("ControlNOBKDom.virome_Relabun.txt", header = TRUE,row.names = 1, sep = " ")
dataTransposed<-t(data)
beta_dist <- vegdist(t(data),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedControlBetaNoPolyomaNoBKDom_100523.csv")

data<-read.delim("ControlNOBKDom.virome_metadata.txt", header=TRUE,row.names = 1, sep = " ")
sub_info <- cbind( smp = data$Sample, # sample ID
                   id = data$Patient.ID, # subject ID
                   tp = data$Time,
                   status = data$Status,
                   gender=data$Gender)
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedControlBetaNoPolyomaNoBKDom_100523.csv", header=T, row.names = 1, sep = ",")
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))#, plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   tp+gender , data = sub_clust, permutations = perm, by = "margin")

#female
female.virome<-m.2[ which(m.2$Gender=='female'),]
d.2<-d.1[,colnames(d.1) %in% female.virome$Sample]
write.table(female.virome,"femaleNOBKDom.virome_metadata.txt")
write.table(d.2,"femaleNOBKDom.virome_Relabun.txt")
data<-read.delim("femaleNOBKDom.virome_Relabun.txt", header = TRUE,row.names = 1, sep = " ")
dataTransposed<-t(data)
beta_dist <- vegdist(t(data),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedFemaleBetaNoPolyomaNoBKDom_100523.csv")

data<-read.delim("femaleNOBKDom.virome_metadata.txt", header=TRUE,row.names = 1, sep = " ")
sub_info <- cbind( smp = data$Sample, # sample ID
                   id = data$Patient.ID, # subject ID
                   tp = data$Time,
                   status = data$Status,
                   gender=data$Gender)
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedFemaleBetaNoPolyomaNoBKDom_100523.csv", header=T, row.names = 1, sep = ",")
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))#, plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   tp+status , data = sub_clust, permutations = perm, by = "margin")

#male
male.virome<-m.2[ which(m.2$Gender=='male'),]
d.2<-d.1[,colnames(d.1) %in% male.virome$Sample]
write.table(male.virome,"maleNOBKDom.virome_metadata.txt")
write.table(d.2,"maleNOBKDom.virome_Relabun.txt")
data<-read.delim("maleNOBKDom.virome_Relabun.txt", header = TRUE,row.names = 1, sep = " ")
dataTransposed<-t(data)
beta_dist <- vegdist(t(data),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedMaleBetaNoPolyomaNoBKDom_100523.csv")

data<-read.delim("maleNOBKDom.virome_metadata.txt", header=TRUE,row.names = 1, sep = " ")
sub_info <- cbind( smp = data$Sample, # sample ID
                   id = data$Patient.ID, # subject ID
                   tp = data$Time,
                   status = data$Status,
                   gender=data$Gender)
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedMaleBetaNoPolyomaNoBKDom_100523.csv", header=T, row.names = 1, sep = ",")
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))#, plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   tp+status , data = sub_clust, permutations = perm, by = "margin")

m.1 <- read.delim("maleNOBKDom.virome_metadata.txt", header=TRUE,row.names = 1, sep = " ")
m.1$ptnum<-as.character(m.1$Patient.ID)
d.1<-read.delim("weightedMaleBetaNoPolyomaNoBKDom_100523.csv", row.names = 1, header = TRUE, sep = ",")
distMatrix <- as.matrix(d.1)
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC3 <- as.matrix(res$vectors[,3])
PC1
PC2
PC3
data<-read.delim("PC_NoPolyNoBKDom_Male.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
data$PC3 <- PC3
#KT BKV
p <- ggplot(data, aes(x=PC1, y=PC3, color=as.factor(Status))) +
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +geom_point(size=2.5)+
  scale_color_manual(breaks = c("BK+","Control"),
                     values = c("#b2182b", "#045a8d")) +
  theme(aspect.ratio=1)+labs(title="No BK dominate males by BKV status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#calculate optimal number of clusters (sum of squares)
WUDM<-read.csv("weightedBetaNoPolyomaNoBKDom_100523.csv",row.names=1,header = T)
a<-fviz_nbclust(WUDM, pam, method = "wss")
#plot optimal number of clusters (sum of squares)
a
b<-fviz_nbclust(WUDM, pam, method = "silhouette")
b
#plot number of clusters vs. gap statistic
gap_stat <- clusGap(WUDM,
                    FUN = pam,
                    K.max = 20, #max clusters to consider
                    B = 50) #total bootstrapped iterations
fviz_gap_stat(gap_stat)
set.seed(123)
km.res1<-kmeans(WUDM, 6, iter.max = 10, nstart = 25)
km.res1
set.seed(123)
km.res2<-kmeans(WUDM, 7, iter.max = 10, nstart = 25)
km.res2
set.seed(123)
km.res3<-kmeans(WUDM, 8, iter.max = 10, nstart = 25)
km.res3
set.seed(123)
km.res4<-kmeans(WUDM, 9, iter.max = 10, nstart = 25)
km.res4
set.seed(123)
km.res5<-kmeans(WUDM, 10, iter.max = 10, nstart = 25)
km.res5

library(mclogit)
m <- read.delim("Analysis_BK_092623.txt", header=TRUE)
m.1<-m[!m$k.means.noBKdom=="#N/A",]
m.1$Patient.ID<-as.factor(m.1$Patient.ID)
m.1$Status<-as.factor(m.1$Status)
m.1$Gender<-as.factor(m.1$Gender)
m.1$k.means.noBKdom<-as.factor(m.1$k.means.noBKdom)

Subset1<-subset(m.1,k.means.noBKdom =='2'|k.means.noBKdom =='3'|k.means.noBKdom =='4'|k.means.noBKdom =='5'|k.means.noBKdom =='6')
Subset2<-subset(m.1,k.means.noBKdom =='3'|k.means.noBKdom =='4'|k.means.noBKdom =='5'|k.means.noBKdom =='6')
Subset3<-subset(m.1,k.means.noBKdom =='4'|k.means.noBKdom =='5'|k.means.noBKdom =='6')
Subset4<-subset(m.1,k.means.noBKdom =='5'|k.means.noBKdom =='6')

set.seed(123)
(comm.mblogit <- mblogit(k.means.noBKdom~Status+Gender+Time,  data = m.1, random=~1|Patient.ID))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(k.means.noBKdom~Status+Gender+Time,  data = Subset1, random=~1|Patient.ID))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(k.means.noBKdom~Status+Gender+Time,  data = Subset2, random=~1|Patient.ID))
summary(comm.mblogit.subset2)
set.seed(123)
(comm.mblogit.subset3 <- mblogit(k.means.noBKdom~Status+Gender+Time,  data = Subset3, random=~1|Patient.ID))
summary(comm.mblogit.subset3)
set.seed(123)
(comm.mblogit.subset4 <- mblogit(k.means.noBKdom~Status+Gender+Time,  data = Subset4, random=~1|Patient.ID))
summary(comm.mblogit.subset4)


library(Maaslin2)
input_data <- read.delim("BK_cleanNormalized_relativeAbundance_noPolyomaNoBKDom.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("subsetAnalysis_BK_092623.txt", row.names = 1)
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_NOBKDom_101323', transform = "NONE",
  fixed_effects = c('Time','Status','Gender'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)


setwd("/Users/rmaqsood/Documents/BK/virome/")
d.1<-read.delim("subsetNoBK_cleanNormalized_counts_.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
data_vegan<-diversity(dataTransposed, index = 'shannon')
write.table(data_vegan,"subsetNoBK_ViromeAlphaDiversity_100623.txt", sep = '\t')

m.1 <- read.delim("subsetAnalysis_BK_092623.txt", header=TRUE, row.names = 1)
m.1$Patient.ID<-as.character(m.1$Patient.ID)

set.seed(100)
allrichness<-lme(Richness~Gender+Time+Status, random=~1|Patient.ID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~Gender + Time +Status, random=~1|Patient.ID, data=m.1)
summary(allalphaDiversity)


loessPlotRichness <- ggplot(m.1, aes(x = Time, y = Richness, color = Gender, fill = Gender)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Time), max(m.1$Time), by =1 ),1)) + labs (x= "Post-transplant week",y="Virome Richness", title = "Females vs males" )
loessPlotRichness

loessPlotShannon <- ggplot(m.1, aes(x = Time, y = AlphaDiversity, color = Gender, fill = Gender)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Time), max(m.1$Time), by =1 ),1)) + labs (x= "Post-transplant week",y="Virome Alpha Diversity", title = "Females vs males" )
loessPlotShannon


setwd('/Users/rmaqsood/Documents/BK/virome/')
library(Maaslin2)
input_data <- read.delim("BK_cleanNormalized_relativeAbundance_noPolyomaNoBKDom.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("subsetAnalysis_BK_092623.txt", row.names = 1)
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_oPolyomaNoBKDom_100623', transform = "NONE",
  fixed_effects = c('Time','Status','Gender'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)


##################k-meansBKDom##################
library(factoextra)
library(cluster)
library(vegan)
setwd('/Users/rmaqsood/Documents/BK/virome/')
m.1 <- read.delim("Analysis_BK_092623.txt", header=TRUE)
m.1$Patient.ID<-as.character(m.1$Patient.ID)
d.1<-read.delim("dominantBK_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
#write.csv(d.1,"LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance_033023.txt")
dataTransposed<-t(d.1)
#m.2<-m.1[!m.1$Richness=="#N/A",]
#d.2<-d.1[,colnames(d.1) %in% m.1$Sample]
beta_dist <- vegdist(t(d.1),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaBKDom_100623.csv")

#calculate optimal number of clusters (sum of squares)
WUDM<-read.csv("weightedBetaBKDom_100623.csv",row.names=1,header = T)
a<-fviz_nbclust(WUDM, pam, method = "wss")
#plot optimal number of clusters (sum of squares)
a
b<-fviz_nbclust(WUDM, pam, method = "silhouette")
b
#plot number of clusters vs. gap statistic
gap_stat <- clusGap(WUDM,
                    FUN = pam,
                    K.max = 20, #max clusters to consider
                    B = 50) #total bootstrapped iterations
fviz_gap_stat(gap_stat)
set.seed(123)
km.res<-kmeans(WUDM, 3, iter.max = 10, nstart = 25)
km.res


library(mclogit)
m <- read.delim("Analysis_BK_092623.txt", header=TRUE)
m.1<-m[!m$k.means.BKdom=="#N/A",]
m.1$Patient.ID<-as.factor(m.1$Patient.ID)
m.1$Status<-as.factor(m.1$Status)
m.1$Gender<-as.factor(m.1$Gender)
m.1$k.means.BKdom<-as.factor(m.1$k.means.BKdom)

Subset1<-subset(m.1,k.means.BKdom =='2'|k.means.BKdom =='3')

set.seed(123)
(comm.mblogit <- mblogit(k.means.BKdom~Status+Gender+Time,  data = m.1, random=~1|Patient.ID))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(k.means.BKdom~Gender+Time,  data = Subset1, random=~1|Patient.ID))
summary(comm.mblogit.subset1)


# m <- read.delim("Analysis_BK_092623.txt", header=TRUE)
# m.1<-m[!m$k.means.BKdom=="#N/A",]
# m.1$Patient.ID<-as.factor(m.1$Patient.ID)
# m.1$Status<-as.factor(m.1$Status)
# m.1$Gender<-as.factor(m.1$Gender)
# m.1$k.means.BKdom<-as.factor(m.1$k.means.BKdom)
# write.table(m.1, "bkDom_metadata.txt")
library(Maaslin2)
input_data <- read.delim("dominantBK_cleanNormalized_relativeAbundance.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("bkDom_metadata.txt", row.names = 1)
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_BKDom_101323', transform = "NONE",
  fixed_effects = c('Time','Status','Gender'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)

library(Maaslin2)
input_data <- read.delim("BKVNOBKDom.virome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("BKVNOBKDom.virome_metadata.txt", row.names = 1,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_BKNoDom_110823', transform = "NONE",
  fixed_effects = c('Time','Status','Gender'),
  random_effects = c('Patient.ID'),
  normalization = 'NONE',
  standardize = TRUE)

setwd('/Users/rmaqsood/Documents/BK/virome/')
library(ggplot2)
library(RColorBrewer)
data<-read.delim("maaslin2_female_BKV.txt")
ggplot(data, aes(x = Status, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

setwd('/Users/rmaqsood/Documents/BK/virome/')
library(ggplot2)
library(RColorBrewer)
data<-read.delim("maaslin2_female_time.txt")
ggplot(data, aes(x = Time, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

setwd('/Users/rmaqsood/Documents/BK/virome/')
library(ggplot2)
library(RColorBrewer)
data<-read.delim("maaslin2_allsamples_status.txt")
ggplot(data, aes(x = Status, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

setwd('/Users/rmaqsood/Documents/BK/virome/')
library(ggplot2)
library(RColorBrewer)
data<-read.delim("maaslin2_allsamples_time.txt")
ggplot(data, aes(x = Time, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()



