#analysis of tadpole body sizes relative to mesocosm treatments
library(stats)
library(ggplot2)
library(reshape)
library(reshape2)
tadsize<-read.csv("C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2010 Mesocosms//Tadpole.csv", header=T)
str(tadsize)
#1 $ TankNumber : int  10 10 10 10 10 10 10 10 10 10 ...
#2 MayflyTx   : int  0 0 0 0 0 0 0 0 0 0 ...
#3 $ BdTx       : Factor w/ 2 levels "Bd-","Bd+": 2 2 2 2 2 2 2 2 2 2 ...
#4 $ SampleDate : Factor w/ 5 levels "7/26/2010","8/11/2010",..: 1 1 1 1 5 5 5 5 2 2 ...
#5 $ AnimalID   : int  1 2 3 4 1 2 3 4 1 2 ...
#6 $ soakID     : int  16 19 23 26 70 69 71 72 -32 -33 ...
#7 $ TadpoleMark: Factor w/ 5 levels "OO","OY","unmarked",..: 5 4 1 2 1 5 4 2 4 5 ...
#8 $ GosnerStage: int  36 37 36 36 37 37 38 37 41 40 ...
#9 $ UBPct      : int  100 100 100 100 100 100 100 100 100 100 ...
#10 $ LBPct      : int  100 100 100 100 100 100 100 100 100 100 ...
#11 $ Trpig      : Factor w/ 3 levels "","Toothrows fully pigmented",..: 2 2 2 2 2 2 2 2 2 2 ...
#12 $ TadWt      : num  2.7 3.9 3.4 2.8 3.7 2.9 4.9 3.1 4.4 3.2 ...
#13 $ BL         : num  23.4 25.9 25.5 23.9 26.8 26.3 29.8 23.9 28 25.7 ...
#14 $ TAL        : num  48.5 54.6 49.5 47.5 47.2 47.9 52.2 46.7 54.7 48.6 ...
#15 $ TMW        : num  5.9 6.7 6.2 4.9 5.7 5.7 6.3 5.6 6.6 6.2 ...
#16 $ TMH        : num  7.8 8.8 8.9 8.4 7.4 8.1 7.4 7.8 7 6.8 ...


tadsize<-tadsize[,-c(6,7,9,10,11)]
tadsize<-tadsize[-c(168, 157,221),]
levels(tadsize$SampleDate)
#[1] "7/26/2010" "8/11/2010" "8/15/2010" "8/16/2010" "8/2/2010" 
tadsize<-tadsize[-c(which(tadsize$SampleDate=="7/26/2010" | tadsize$SampleDate=="8/2/2010"| tadsize$SampleDate=="8/11/2010")),]
tadsize<-tadsize[-c(which(tadsize$GosnerStage>=42)),]
summary(tadsize[which(tadsize$BdTx=="Bd-"),])
summary(tadsize[which(tadsize$BdTx=="Bd+"),])
plot(tadsize)
aggregate(tadsize$AnimalID, by=list(as.factor=(tadsize$TankNumber)), FUN=max)

cor.test(tadsize$BL, tadsize$TadWt, method="sp")#r = 0.57, p<<0.0001
cor.test(tadsize$BL, tadsize$TadWt)#r = 0.57, p<<0.0001
hist(tadsize$BL[-which(tadsize$BL>35)])
hist(tadsize$TadWt)

#body length
aggregate(tadsize$BL, by=list(as.factor(tadsize$TankNumber)), FUN=mean)
SummaryMeans<-aggregate(tadsize$BL, by=list(tadsize$TankNumber, tadsize$BdTx, as.factor(tadsize$MayflyTx)), FUN=mean)
anova(aov(x~Group.2, data=SummaryMeans))
TukeyHSD(aov(x~Group.2, data=SummaryMeans))
anova(aov(x~Group.3, data=SummaryMeans))
TukeyHSD(aov(x~Group.3, data=SummaryMeans))
anova(aov(x~as.factor+Group.2, data=SummaryMeans))
anova(aov(x~as.factor*Group.2, data=SummaryMeans))
ggplot(SummaryMeans, aes(y=x, x=Group.2))+
	geom_boxplot()
ggplot(SummaryMeans, aes(y=x, x=Group.3))+
	geom_boxplot()

#Stage
SummaryMeans<-aggregate(tadsize$GosnerStage, by=list(tadsize$TankNumber, tadsize$BdTx, as.factor(tadsize$MayflyTx)), FUN=mean)
anova(aov(x~Group.2, data=SummaryMeans))
TukeyHSD(aov(x~Group.2, data=SummaryMeans))
anova(aov(x~as.factor, data=SummaryMeans))
anova(aov(x~as.factor+Group.2, data=SummaryMeans))
anova(aov(x~as.factor*Group.2, data=SummaryMeans))
ggplot(SummaryMeans, aes(y=x, x=Group.2))+
	geom_boxplot()

#Weight
SummaryMeans<-aggregate(tadsize$TadWt, by=list(tadsize$TankNumber, tadsize$BdTx, as.factor=(tadsize$MayflyTx)), FUN=mean)
anova(aov(x~Group.2, data=SummaryMeans))
TukeyHSD(aov(x~Group.2, data=SummaryMeans))
anova(aov(x~as.factor, data=SummaryMeans))
anova(aov(x~as.factor+Group.2, data=SummaryMeans))
anova(aov(x~as.factor*Group.2, data=SummaryMeans))
ggplot(SummaryMeans, aes(y=x, x=Group.2))+
	geom_boxplot()

#Tail length
SummaryMeans<-aggregate(tadsize$TAL, by=list(tadsize$TankNumber, tadsize$BdTx, as.factor=(tadsize$MayflyTx)), FUN=mean)
anova(aov(x~Group.2, data=SummaryMeans))
TukeyHSD(aov(x~Group.2, data=SummaryMeans))
anova(aov(x~as.factor, data=SummaryMeans))
anova(aov(x~as.factor+Group.2, data=SummaryMeans))
anova(aov(x~as.factor*Group.2, data=SummaryMeans))
ggplot(SummaryMeans, aes(y=x, x=Group.2))+
	geom_boxplot()

#TMW
SummaryMeans<-aggregate(tadsize$TMW, by=list(tadsize$TankNumber, tadsize$BdTx, as.factor=(tadsize$MayflyTx)), FUN=mean)
anova(aov(x~Group.2, data=SummaryMeans))
TukeyHSD(aov(x~Group.2, data=SummaryMeans))
anova(aov(x~as.factor, data=SummaryMeans))
anova(aov(x~as.factor+Group.2, data=SummaryMeans))
anova(aov(x~as.factor*Group.2, data=SummaryMeans))
ggplot(SummaryMeans, aes(y=x, x=Group.2))+
	geom_boxplot()

#TMH
SummaryMeans<-aggregate(tadsize$TMH, by=list(tadsize$TankNumber, tadsize$BdTx, as.factor=(tadsize$MayflyTx)), FUN=mean)
anova(aov(x~Group.2, data=SummaryMeans))
TukeyHSD(aov(x~Group.2, data=SummaryMeans))
anova(aov(x~as.factor, data=SummaryMeans))
anova(aov(x~as.factor+Group.2, data=SummaryMeans))
anova(aov(x~as.factor*Group.2, data=SummaryMeans))
ggplot(SummaryMeans, aes(y=x, x=Group.2))+
	geom_boxplot()





















boxplot(tadsize$BL~as.factor(tadsize$TankNumber))

summary(tadsize$BdTx)


summary(tadsize[which(tadsize$BdTx=="Bd-"),])
summary(tadsize[which(tadsize$BdTx=="Bd+"),])

anova(aov(GosnerStage~MayflyDensity*BdTx, data=tadsize))
anova(aov(GosnerStage~BdTx, data=tadsize))
anova(aov(GosnerStage~MayflyDensity, data=tadsize))

anova(aov(BL~MayflyDensity*BdTx, data=tadsize))
anova(aov(BL~BdTx, data=tadsize))
anova(aov(BL~MayflyDensity, data=tadsize))
ggplot(tadsize, aes(x=BdTx, y=BL))+
	geom_boxplot(aes(fill=BdTx))
	
anova(aov(TadWt~MayflyDensity*BdTx, data=tadsize))
anova(aov(TadWt~MayflyDensity, data=tadsize))
anova(aov(TadWt~BdTx, data=tadsize))
ggplot(tadsize, aes(x=BdTx, y=TadWt))+
	geom_boxplot(aes(fill=BdTx))

anova(aov(TMW~MayflyDensity*BdTx, data=tadsize))
anova(aov(TMW~MayflyDensity, data=tadsize))
anova(aov(TMW~BdTx, data=tadsize))
ggplot(tadsize, aes(x=BdTx, y=TMW))+
	geom_boxplot(aes(fill=BdTx))

anova(aov(TMH~MayflyDensity*BdTx, data=tadsize))
anova(aov(TMH~MayflyDensity, data=tadsize))
anova(aov(TMH~BdTx, data=tadsize))
ggplot(tadsize, aes(x=BdTx, y=TMH))+
	geom_boxplot(aes(fill=BdTx))

tadmelt<-melt(tadsize[,c(4,5,6,7,8,9)], id=c("BdTx"))
p<-ggplot(tadmelt, aes(x=BdTx, y=value))+
	geom_boxplot(aes(fill=BdTx))+
	facet_grid(.~variable)

f_labels <- data.frame(variable = c("GosnerStage", "TadWt", "BL", "TMW", "TMH"), 
				label = c("p=0.28", "p=0.014", "p<<0.01", "p=0.36", "p<<0.01"))
p + geom_text(x=1.5, y=15, aes(label=label), data=f_labels)


#-----------------------
#since stg 42 tadpoles shirnk, remove stage 42, rerun:
tadsize<-tadsize[-c(which(tadsize$GosnerStage>=42)),]
summary(tadsize[which(tadsize$BdTx=="Bd-"),])
summary(tadsize[which(tadsize$BdTx=="Bd+"),])
plot(tadsize)

anova(aov(BL~SampleDate, data=tadsize))# itested sample date in all models below, it gets dropped as not significant
anova(aov(GosnerStage~BdTx, data=tadsize))
anova(aov(TadWt~BdTx, data=tadsize))
anova(aov(BL~BdTx, data=tadsize))
anova(aov(TMW~BdTx, data=tadsize))
anova(aov(TMH~BdTx, data=tadsize))

tadmelt<-melt(tadsize[,c(1,4,5,6,7,8,9)], id=c("SampleDate", "BdTx"))
p<-ggplot(tadmelt, aes(x=BdTx, y=value))+
	geom_boxplot(aes(fill=BdTx))+
	facet_grid(SampleDate~variable)
f_labels <- data.frame(variable = c("GosnerStage", "TadWt", "BL", "TMW", "TMH"), 
				label = c("p=0.9", "p=0.08", "p=0.0008", "p=0.2", "p<<0.01"))
p + geom_text(x=1.5, y=15, aes(label=label), data=f_labels)

#-----------------------
#differences in starting tadpole sizes?
startsize<-read.csv("C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2010 Mesocosms//MesocosmTreatements_TadpoleBodySizeMeasurements.csv", header=T)
startsize<-startsize[-c(33, 111, 216, 223),]
levels(startsize$SampleDate)
#[1] "7/26/2010" "8/11/2010" "8/15/2010" "8/16/2010" "8/2/2010" 
startsize<-startsize[c(which(startsize$SampleDate=="7/26/2010")),]

anova(aov(GosnerStage~BdTx, data=startsize))
anova(aov(TadWt~BdTx, data=startsize))
anova(aov(BL~BdTx, data=startsize))
anova(aov(TMW~BdTx, data=startsize))
anova(aov(TMH~BdTx, data=startsize))
#no differences in tadpoles in treatment categories at start.

