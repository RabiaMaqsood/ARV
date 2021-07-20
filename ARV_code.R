library(nlme)
library(vegan)
library(ggplot2)
library(ape)
library(decontam)
library(ggplot2)
library(RColorBrewer)
library(Maaslin2)
library(mclogit)


########################LME########################
data<-read.delim("arv_analysis.txt", row.names = 1) #Use metadata tab in supplementary excel file
m1 <- lme(shannon ~ Week + treated+antibiotics , random=~1|Subject, data=data)
summary(m1)

loessPlotRichness <- ggplot(data, aes(x = Week, y = Richness, color = treated, fill = treated)) +
  geom_point() +
  stat_smooth(method="loess", se=TRUE, span=0.5, level=0.95)+
  theme_bw()
loessPlotRichness

m2 <- lme(Richness ~ Week + treated+antibiotics , random=~1|Subject, data=data)
summary(m2)

loessPlotShannon <- ggplot(data, aes(x = Week, y = shannon, color = treated, fill = treated)) +
  geom_point() +
  stat_smooth(method="loess", se=TRUE, span=0.5, level=0.95)+
  theme_bw()
loessPlotShannon

################betaDiversity########################
input <- read.csv("weighted-unifrac-distance-matrix.tsv", header=TRUE, row.names = 1, sep="\t") #Use data matrix in betaDiversity.2 tab in supplementary excel file
metadata <- read.delim("ArvMetadata.txt", header=TRUE, row.names = 1) #Use metadata tab in supplementary excel file
data <- subset(metadata, metadata[1,] %in% rownames(input))
data <- data[order(as.numeric(data$X.SampleID)),]
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

data<-read.delim("PC_TimeAbx.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file

q <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(treated))) + 
  geom_point() + xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  stat_ellipse() +
  scale_color_manual(breaks = c("trt", "untrt"),
                     values = c("#039254", "#941c52")) +
  theme(aspect.ratio=1)
q
q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

########################Adonis########################
input <- read.csv("weighted-unifrac-distance-matrix.tsv", header=TRUE, row.names = 1, sep="\t") #Use data matrix in betaDiversity.2 tab in supplementary excel file
metadata <- read.csv("ARVMetadata.txt", header=TRUE, row.names=1, sep='\t')
metasubset <- subset(metadata, rownames(metadata) %in% colnames(input))

folate <- metasubset$Folate
iron <- metasubset$Iron
antibiotics <- metasubset$Antibiotics
Week<-metasubset$Week
Treatment<-metasubset$treated
trtTime<-metasubset$TimeWeek
patients<-metasubset$idNum
abxTime<-metasubset$TimeAbx

adon.results<-adonis(t(input) ~ folate,perm=10000)
print(adon.results)

adon.results<-adonis(t(input) ~ iron,perm=10000)
print(adon.results)

adon.results<-adonis(t(input) ~ antibiotics,perm=10000)
print(adon.results)

adon.results<-adonis(t(input) ~ Week,perm=10000)
print(adon.results)

adon.results<-adonis(t(input) ~ Treatment,perm=10000)
print(adon.results)

adon.results<-adonis(t(input) ~ Week,perm=10000, strata = patients)
print(adon.results)

adon.results<-adonis(t(input) ~ trtTime,perm=10000,strata = patients)
print(adon.results)

adon.results<-adonis(t(input) ~ abxTime,perm=10000,strata = patients)
print(adon.results)

########################Maaslin########################
input_data <- read.delim("feature-table-arv-8k-rarefied.txt", row.names = 1) #Use data in masslin tab in supplementary excel file
input_data2<-t(input_data)
input_metadata <- read.delim("ArvMetadata.txt", row.names = 1) #Use data in metadata tab in supplementary excel file
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_092220', transform = "NONE",
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
data<-read.delim2("arv_analysisCommState.txt", row.names = 1) #Use data in metadata tab in supplementary excel file
data$Subject<-as.factor(data$Subject)
data$treated<-as.factor(data$treated)
data$Antibiotics<-as.factor(data$Antibiotics)
data$community.State<-as.factor(data$community.State)

Subset1<-subset(data,community.State =='grp2'|community.State =='grp3'|community.State =='grp4'|community.State =='grp5'|community.State =='grp6'|community.State =='grp7')
Subset2<-subset(data,community.State =='grp3'|community.State =='grp4'|community.State =='grp5'|community.State =='grp6'|community.State =='grp7')
Subset3<-subset(data,community.State =='grp4'|community.State =='grp5'|community.State =='grp6'|community.State =='grp7')
Subset4<-subset(data,community.State =='grp5'|community.State =='grp6'|community.State =='grp7')  
Subset5<-subset(data,community.State =='grp6'|community.State =='grp7')  



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

shannonModel <- lme(shannon ~ Week + Antibiotics + treated.2 , random=~1|Subject, data = data)
summary(shannonModel)


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



