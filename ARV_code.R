setwd("/Users/rabiamaqsood/Documents/ARV/16s_arv/silvaQiiime/silva_012122/")
library(nlme)
library(vegan)
library(ggplot2)
library(ape)
library(decontam)
library(ggplot2)
library(RColorBrewer)
library(Maaslin2)
library(mclogit)
library(factoextra)
library(cluster)
library(dendextend)
library(mclogit)

data<-read.delim("Analysis_012422.txt", row.names = 1) #Use metadata tab in supplementary excel file
set.seed(100)
m1 <- lme(shannon ~ Week + treated.2 + antibiotics.original , random=~1|Subject, data=data)
summary(m1)

loessPlotRichness <- ggplot(data, aes(x = Week, y = Richness, color = antibiotics.original, fill = antibiotics.original)) +
  geom_point() +
  stat_smooth(method="loess", se=TRUE, span=0.5, level=0.95)+
  theme_bw()
loessPlotRichness

set.seed(100)
m2 <- lme(Richness ~ Week + treated.2 + antibiotics.original , random=~1|Subject, data=data)
summary(m2)

loessPlotShannon <- ggplot(data, aes(x = Week, y = shannon, color = antibiotics.original, fill = antibiotics.original)) +
  geom_point() +
  stat_smooth(method="loess", se=TRUE, span=0.5, level=0.95)+
  theme_bw()
loessPlotShannon
################betaDiversity########################
input <- read.csv("weighted_unifrac_distance_matrix.tsv", header=TRUE, row.names = 1, sep="\t") #Use data matrix in betaDiversity.2 tab in supplementary excel file
metadata <- read.delim("ArvMetadata.txt") #Use metadata tab in supplementary excel file
#data <- subset(metadata, metadata[1,] %in% rownames(input))
#data <- data[order(as.numeric(data$X.SampleID)),]
distMatrix <- as.matrix(input)

res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_Data.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2

p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Antibiotics))) + 
  geom_point() + xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Yes", "No"),
                     values = c("#ee3a90", "#c2c1c0")) +
  theme(aspect.ratio=1)
p
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

data<-read.delim("PC_Data.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file

q <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(treated))) + 
  geom_point() + xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  stat_ellipse() +
  scale_color_manual(breaks = c("trt", "untrt"),
                     values = c("#039254", "#941c52")) +
  theme(aspect.ratio=1)
q
q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

data$Week.2<-as.factor(data$Week)
r <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Week.2))) + 
  geom_point() + xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  stat_ellipse() +
  scale_color_manual(breaks = c("1", "2","3","4"),
                     values = c("#F6eb16", "#2baee4","#fdc010","#703a96")) +
  theme(aspect.ratio=1)
r
r + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

setwd("/Users/rabiamaqsood/Documents/ARV/16s_arv/silvaQiiime/silva_012122/")

####Full
MetaDat <- read.csv("Analysis_012422.txt", header=TRUE, sep='\t',check.names = FALSE)
# subject information, 180 samples
sub_info <- cbind( smp = MetaDat$SampleID, # sample ID
                   id = MetaDat$Subject, # subject ID 
                   grp = MetaDat$treated, # trt, unt
                   tp = MetaDat$Week, 
                   ab = MetaDat$Antibiotics) # time point
sub_info<-as.data.frame(sub_info)

WUDM<-read.table("weighted_unifrac_distance_matrix.tsv", header=T, row.names = 1) #Use weighted unifrac table from betaDiversity.2 tab
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA

sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)

ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)

perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~   grp +tp + ab, data = sub_clust, permutations = perm, by = "margin")

pair_dat <- function(dist_mat, info, t1, t2) {
  ind <- info$tp == t1 | info$tp == t2
  pair_info <- info[ind, ]
  pair_mat <- as.dist(dist_mat[ind, ind])
  return(list(dm = pair_mat, pi = pair_info))
}

pair_tp <- combn(unique(sub_clust$tp), 2)
pval <- c()
for (i in 1:ncol(pair_tp)) {
  x <- pair_dat(WUDM, sub_clust, pair_tp[1,i], pair_tp[2,i])
  perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = x$pi$id))
  res <- adonis2(x$dm ~  + tp, data = x$pi, permutations = perm, by = "margin")
  p <- res$`Pr(>F)`[2]
  pval <- c(pval, p)
}
pval

padj <- p.adjust(pval, method = "BH")
rbind(pair_tp, padj)

########################
AbxNo
MetaDat <- read.csv("Analysis_012422.txt", header=TRUE, sep='\t',check.names = FALSE)
# subject information, 180 samples
sub_info.1 <- cbind( smp = MetaDat$SampleID, # sample ID
                   id = MetaDat$Subject, # subject ID 
                   grp = MetaDat$treated, # trt, unt
                   tp = MetaDat$Week, 
                   ab = MetaDat$Antibiotics) # time point
sub_info.1<-as.data.frame(sub_info.1)

WUDM.1<-read.table("weighted_unifrac_distance_matrix_AbxNo.tsv", header=T, row.names = 1) #Use weighted unifrac table for ABX No from betaDiversity.3 tab
x.1 <- setdiff(sub_info.1$smp, rownames(WUDM.1))
sub_info.1[sub_info.1$smp %in% x.1, 'smp'] <- NA

sub_wide.1 <- reshape(sub_info.1, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust.1 <- sub_info.1[which(sub_info.1$smp %in% rownames(WUDM.1)), ]
sub_clust.1$id <- as.factor(sub_clust.1$id)
sub_clust.1$tp <- as.factor(sub_clust.1$tp)

ord.1 <- match(sub_clust.1$smp, rownames(WUDM.1))
WUDM.1 <- WUDM.1[ord.1, ord.1]
wudm.1 <- as.dist(WUDM.1)

perm.1 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust.1$id))
adonis.1<-adonis2(wudm.1 ~ tp , data = sub_clust.1, permutations = perm.1, by = "margin")
adonis.1

########################
MetaDat <- read.csv("Analysis_012422.txt", header=TRUE, sep='\t',check.names = FALSE)
# subject information, 180 samples
sub_info.2 <- cbind( smp = MetaDat$SampleID, # sample ID
                   id = MetaDat$Subject, # subject ID 
                   grp = MetaDat$treated, # trt, unt
                   tp = MetaDat$Week, 
                   ab = MetaDat$Antibiotics) # time point
sub_info.2<-as.data.frame(sub_info.2)

WUDM.2<-read.table("weighted_unifrac_distance_matrix_AbxYes.tsv", header=T, row.names = 1) #Use weighted unifrac table for ABX Yes from betaDiversity.3 tab
x.2 <- setdiff(sub_info.2$smp, rownames(WUDM.2))
sub_info.2[sub_info.2$smp %in% x.2, 'smp'] <- NA

sub_wide.2 <- reshape(sub_info.2, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust.2 <- sub_info.2[which(sub_info.2$smp %in% rownames(WUDM.2)), ]
sub_clust.2$id <- as.factor(sub_clust.2$id)
sub_clust.2$tp <- as.factor(sub_clust.2$tp)

ord.2 <- match(sub_clust.2$smp, rownames(WUDM.2))
WUDM.2 <- WUDM.2[ord.2, ord.2]
wudm.2 <- as.dist(WUDM.2)

perm.2 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust.2$id))
adonis.2<-adonis2(wudm.2 ~ tp, data = sub_clust.2, permutations = perm.2, by = "margin")
adonis.2

########################
MetaDat <- read.csv("Analysis_012422.txt", header=TRUE, sep='\t',check.names = FALSE)
# subject information, 180 samples
sub_info.3 <- cbind( smp = MetaDat$SampleID, # sample ID
                     id = MetaDat$Subject, # subject ID 
                     grp = MetaDat$treated, # trt, unt
                     tp = MetaDat$Week, 
                     ab = MetaDat$Antibiotics) # time point
sub_info.3<-as.data.frame(sub_info.3)

WUDM.3<-read.table("weighted_unifrac_distance_matrix_cARTexposed.tsv", header=T, row.names = 1)#Use weighted unifrac table for exposed from betaDiversity.3 tab
x.3 <- setdiff(sub_info.3$smp, rownames(WUDM.3))
sub_info.3[sub_info.3$smp %in% x.3, 'smp'] <- NA

sub_wide.3 <- reshape(sub_info.3, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust.3 <- sub_info.3[which(sub_info.3$smp %in% rownames(WUDM.3)), ]
sub_clust.3$id <- as.factor(sub_clust.3$id)
sub_clust.3$tp <- as.factor(sub_clust.3$tp)

ord.3 <- match(sub_clust.3$smp, rownames(WUDM.3))
WUDM.3 <- WUDM.3[ord.3, ord.3]
wudm.3 <- as.dist(WUDM.3)

perm.3 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust.3$id))
adonis.3<-adonis2(wudm.3 ~ tp, data = sub_clust.3, permutations = perm.3, by = "margin")
adonis.3

########################
MetaDat <- read.csv("Analysis_012422.txt", header=TRUE, sep='\t',check.names = FALSE)
# subject information, 180 samples
sub_info.4 <- cbind( smp = MetaDat$SampleID, # sample ID
                     id = MetaDat$Subject, # subject ID 
                     grp = MetaDat$treated, # trt, unt
                     tp = MetaDat$Week, 
                     ab = MetaDat$Antibiotics) # time point
sub_info.4<-as.data.frame(sub_info.4)

WUDM.4<-read.table("weighted_unifrac_distance_matrix_cARTunexposed.tsv", header=T, row.names = 1)#Use weighted unifrac table for unexposed from betaDiversity.3 tab
x.4 <- setdiff(sub_info.4$smp, rownames(WUDM.4))
sub_info.4[sub_info.4$smp %in% x.4, 'smp'] <- NA

sub_wide.4 <- reshape(sub_info.4, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust.4 <- sub_info.4[which(sub_info.4$smp %in% rownames(WUDM.4)), ]
sub_clust.4$id <- as.factor(sub_clust.4$id)
sub_clust.4$tp <- as.factor(sub_clust.4$tp)

ord.4 <- match(sub_clust.4$smp, rownames(WUDM.4))
WUDM.4 <- WUDM.4[ord.4, ord.4]
wudm.4 <- as.dist(WUDM.4)

perm.4 <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust.4$id))
adonis.4<-adonis2(wudm.4 ~ tp, data = sub_clust.4, permutations = perm.4, by = "margin")
adonis.4



########################Maaslin########################
input_data <- read.delim("feature-table-arv-8k-rarefied.txt", row.names = 1) #Use data in masslin tab in supplementary excel file
input_data2<-t(input_data)
input_metadata <- read.delim("ArvMetadata.txt", row.names = 1) #Use metadata in masslin tab in supplementary excel file
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_021522', transform = "NONE",
  fixed_effects = c('treated','Week','Antibiotics'),
  random_effects = c('idNum'),
  normalization = 'NONE',
  standardize = FALSE)

data<-read.delim("Significant_Results_PrevalanceRelAbun.txt") #Use data in masslin tab in supplementary excel file
ggplot(data, aes(x = Time, y = ASV, size = Prevalance, color = MeanAbundance)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 0.25, 0.5, 0.75, 1),
             labels = c("0", "0.25", "0.5", "0.75", "1"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

########################Multinomial########################
data<-read.delim2("Analysis_012422.txt", row.names = 1) #Use data in metadata tab in supplementary excel file
data$Subject<-as.factor(data$Subject)
data$treated<-as.factor(data$treated)
data$Antibiotics<-as.factor(data$Antibiotics)
data$community.State<-as.factor(data$community.State)

Subset1<-subset(data,community.State =='2'|community.State =='3'|community.State =='4'|community.State =='5'|community.State =='6'|community.State =='7'|community.State =='8')
Subset2<-subset(data,community.State =='3'|community.State =='4'|community.State =='5'|community.State =='6'|community.State =='7'|community.State =='8')
Subset3<-subset(data,community.State =='4'|community.State =='5'|community.State =='6'|community.State =='7'|community.State =='8')
Subset4<-subset(data,community.State =='5'|community.State =='6'|community.State =='7'|community.State =='8')  
Subset5<-subset(data,community.State =='6'|community.State =='7'|community.State =='8')  
Subset6<-subset(data,community.State =='7'|community.State =='8')  


(comm.mblogit <- mblogit(community.State~treated+Week+Antibiotics,  data = data, random=~1|Subject))
summary(comm.mblogit)
(comm.mblogit.subset1 <- mblogit(community.State~treated+Week+Antibiotics,  data = Subset1, random=~1|Subject))
summary(comm.mblogit.subset1)
(comm.mblogit.subset2 <- mblogit(community.State~treated+Week+Antibiotics,  data = Subset2, random=~1|Subject))
summary(comm.mblogit.subset2)
(comm.mblogit.subset3 <- mblogit(community.State~treated+Week+Antibiotics,  data = Subset3, random=~1|Subject))
summary(comm.mblogit.subset3)
(comm.mblogit.subset4 <- mblogit(community.State~treated+Week+Antibiotics,  data = Subset4, random=~1|Subject))
summary(comm.mblogit.subset4)
(comm.mblogit.subset5 <- mblogit(community.State~treated+Week+Antibiotics,  data = Subset5, random=~1|Subject))
summary(comm.mblogit.subset5)
(comm.mblogit.subset6 <- mblogit(community.State~treated+Week+Antibiotics,  data = Subset6, random=~1|Subject))
summary(comm.mblogit.subset6)


########################Decontam########################
set.seed(100)
data <- read.delim("table-arv-v1.txt", row.names=1) #use data in rawdata tab in supplementary excel file
data<-as.matrix(t(data))
metadata<-read.delim("ArvMetadata.txt")#use data in Metadata tab in excel in supplementary excel fil
contam<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.1)
contam.intermediate<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.25)
contam.strict<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.5)
write.csv(contam, "16s_ARV_Default.csv")
write.csv(contam.intermediate,"16s_ARV_Intermediate.csv")
write.csv(contam.strict,"16s_ARV_Strict.csv")



