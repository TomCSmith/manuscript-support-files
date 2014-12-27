##analysis of 2010 mesocosm data
#############
detach()
rm(list=ls())
#############
#load packages:
library(sqldf)
library(nlme)
library(lme4)
library(ggplot2)
library(plyr)
library(ggthemes)
library(reshape2)
library(reshape)
library(scales)
#############
#load and organize data:
whole<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Summer 2010/2010MesocosmCSVdata/ExperimentalAlgae.csv", header=T)
str(whole)
names(whole)
#exclude Shelf Tiles:
whole<-whole[which(whole$SubstrateType=="BottomTile"),]
#Sample 2 includes all non-independent samples, 
#exclude Sample2:
whole<-whole[-which(whole$SampleNumber=="S2"),]
#sample 3 includes SOME non-independent samples:
s1<-which(whole$SampleNumber=="S1")
whole[s1, c(1,5)]

sqldf()

#############
PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	qqnorm(resid(x))
	plot(resid(x)~fitted(x), data=y)
	plot(resid(x)~y[,3], data=y, xlab=colnames(y[3]))
	plot(resid(x)~y[,13], data=y, xlab=colnames(y[13]))
	boxplot(resid(x)~y[,15], data=y, xlab=colnames(y[15]))
	boxplot(resid(x)~y[,16], data=y, xlab=colnames(y[16]))
	plot(resid(x)~y[,17], data=y, xlab=colnames(y[17]))
	boxplot(resid(x)~y[,18], data=y, xlab=colnames(y[18]))
	}
#PlotRes(model, dataset, plot dimensions)

##############
#plotting:
#Figure 4: 2010 mesocosm algal abudnance:
AFDMsummary<-ddply(whole, "Treatment",  summarize,
	meanAFDM=mean(AFDMadjPerMperD), 
	seAFDM=sd(AFDMadjPerMperD)/sqrt(length(AFDMadjPerMperD)))
q<-ggplot(whole, aes(x=factor(Treatment, levels=c("none", "mayflies", "tadpoles", "both")), y=log(AFDMadjPerMperD+0.01)))+
	geom_boxplot(fill="white", color="black", size=1.5, outlier.size=2.5)+
	ylab(expression("log AFDM m"^-2))+
	stat_summary(fun.y="mean", geom="point", shape=23, size=5, color="black", fill="grey")+
	theme(panel.background=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.line=element_line(color="black", size=1.5),
		axis.text.x = element_text(face="bold",colour="black", size=14), 
		axis.text.y = element_text(face="bold",colour="black", size=14), 
		axis.title.x=element_blank(), 
		axis.title.y=element_text(face="bold", color="black", size=14)
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2010_AlgalAbundanceBoxplot.png", 
	res=300, width=5, height=5, units="in", restoreConsole=F)
q
dev.off()

qqnorm(log10(whole$AFDMadjPerMperD+0.0000001))
bartlett.test((AFDMadjPerMperD+0.0000001)^(1/10)~Treatment, data=whole)
#outliers are zeroes:
whole<-whole[-c(which(whole$AFDMadjPerMperD<=0.0001)),]
qqnorm(log10(whole$AFDMadjPerMperD))
bartlett.test(log10(AFDMadjPerMperD)~Treatment, data=whole)
bartlett.test(log10(AFDMadjPerMperD)~TadpolePDensity, data=whole)
bartlett.test(log10(AFDMadjPerMperD)~MFDensity, data=whole)


##############
#models:

model1<-aov(log10(AFDMadjPerMperD)~Treatment, data=whole)
anova(model1)
TukeyHSD(model1)
PlotRes(model1, whole, c(4,5))
summary(model1)
AIC(model1)

#test for interaction
model2<-lm(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity, data=whole)
PlotRes(model2, whole, c(4,5))
summary(model2)
AIC(model2)

model3<-lm(log10(AFDMadjPerMperD)~TadpolePDensity+MFDensity, data=whole)
PlotRes(model3, whole, c(4,5))
summary(model3)
AIC(model3)

model4<-gls(log10(AFDMadjPerMperD)~as.factor(TadpolePDensity)*as.factor(MFDensity), data=whole, method="ML")
PlotRes(model4, whole, c(3,4))
summary(model4)
AIC(model4)	#375.7

TadVar<-varIdent(form=~1|TadpolePDensity)
model5<-gls(log10(AFDMadjPerMperD)~
	as.factor(TadpolePDensity)*as.factor(MFDensity), 
	weights=TadVar,
	method="ML",
	data=whole)
summary(model5)
PlotRes(model5, whole, c(3,4))
AIC(model5)	#377

MFVar<-varIdent(form=~1|MFDensity)
model6<-gls(log10(AFDMadjPerMperD)~
	as.factor(TadpolePDensity)*as.factor(MFDensity), 
	weights=MFVar,
	method="ML",
	data=whole)
summary(model6)
PlotRes(model6, whole, c(3,4))
AIC(model6)	#350.99

#other work, not saved, showed that model 6 is the best for the whole dataset.
model7<-gls(log10(AFDMadjPerMperD)~
	as.factor(TadpolePDensity)*as.factor(MFDensity), 
	weights=MFVar,
	method="REML",
	data=whole)
summary(model7)
anova(model7)

##############
#for bottom tiles only:
bottom<-whole[which(whole$SubstrateType=="BottomTile"),]

qqnorm(log10(bottom$AFDMadjPerMperD))
bartlett.test(log10(AFDMadjPerMperD)~Treatment, data=bottom)
bartlett.test(log10(AFDMadjPerMperD)~TadpolePDensity, data=bottom)
bartlett.test(log10(AFDMadjPerMperD)~MFDensity, data=bottom)

ggplot(bottom, aes(x=Treatment, y=log10(AFDMadjPerMperD)))+
	geom_point(aes(color=SampleDate))
ggplot(bottom, aes(x=Treatment, y=log10(AFDMadjPerMperD)))+
	geom_boxplot()
ggplot(bottom, aes(x=as.factor(TankNumber), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(bottom, aes(x=as.factor(SampleDate), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(bottom, aes(x=as.factor(TadpolePDensity), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(bottom, aes(x=as.factor(MFDensity), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()

bottom1<-lm(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity, data=bottom)
summary(bottom1)
PlotRes(bottom1, bottom, c(3,4) )
AIC(bottom1)	#321

bottom2<-lm(log10(AFDMadjPerMperD)~TadpolePDensity+MFDensity, data=bottom)
summary(bottom2)
PlotRes(bottom2, bottom, c(3,4) )
AIC(bottom2)	#335

bottom3<-gls(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity,method="ML", data=bottom)
summary(bottom3)
PlotRes(bottom3, bottom, c(3,4) )
AIC(bottom3)	#321

bottom4<-gls(log10(AFDMadjPerMperD)~TadpolePDensity+MFDensity,method="ML",data=bottom)
summary(bottom4)
PlotRes(bottom4, bottom, c(3,4) )
AIC(bottom4)	#335.9
#in this case, interaction is not supported

TadVar<-varIdent(form=~1|TadpolePDensity)
bottom5<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+MFDensity, 
	weights=TadVar,
	method="ML",
	data=bottom)
summary(bottom5)
PlotRes(bottom5, bottom, c(3,4) )
AIC(bottom5)	#335

MF<-varIdent(form=~1|MFDensity)
bottom6<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+MFDensity, 
	weights=MF,
	method="ML",
	data=bottom)
summary(bottom6)
PlotRes(bottom6, bottom, c(3,4) )
AIC(bottom6)	#312

Combvar<-varComb(TadVar, MF)
bottom7<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+MFDensity, 
	weights=Combvar,
	method="ML",
	data=bottom)
summary(bottom7)
PlotRes(bottom7, bottom, c(3,4) )
AIC(bottom7)	#314
#no reason to include tad variance

#refit model 6
bottom8<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+MFDensity, 
	weights=MF,
	method="REML",
	data=bottom)
summary(bottom8)

##############
#reporting:
#differences between treatments:
mean(log10(bottom$AFDMadjPerMperD)[])	#-1.25




##############
#for shelf tiles only
shelf<-whole[which(whole$SubstrateType=="ShelfTile"),]

qqnorm(log10(shelf$AFDMadjPerMperD))
bartlett.test(log10(AFDMadjPerMperD)~Treatment, data=shelf)
bartlett.test(log10(AFDMadjPerMperD)~TadpolePDensity, data=shelf)
bartlett.test(log10(AFDMadjPerMperD)~MFDensity, data=shelf)

ggplot(shelf, aes(x=Treatment, y=log10(AFDMadjPerMperD)))+
	geom_point(aes(color=SampleDate))
ggplot(shelf, aes(x=Treatment, y=log10(AFDMadjPerMperD)))+
	geom_boxplot()
ggplot(shelf, aes(x=as.factor(TankNumber), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(shelf, aes(x=as.factor(SampleDate), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(shelf, aes(x=as.factor(TadpolePDensity), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(shelf, aes(x=as.factor(MFDensity), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()

shelfAOV<-aov(log10(AFDMadjPerMperD)~Treatment, data=shelf)
anova(shelfAOV)
TukeyHSD(shelfAOV)#no differences among groups

shelf1<-lm(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity, data=shelf)
summary(shelf1)
PlotRes(shelf1, shelf, c(3,4) )
AIC(shelf1)	#-1.23
#confirms no differences.

shelf2<-lm(log10(AFDMadjPerMperD)~TadpolePDensity+MFDensity, data=shelf)
summary(shelf2)
PlotRes(shelf2, shelf, c(3,4) )
AIC(shelf2)	#-1.033

shelf3<-gls(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity, method="ML", data=shelf)
summary(shelf3)
PlotRes(shelf3, shelf, c(3,4) )
AIC(shelf3)	#-1.23

shelf4<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+MFDensity, 
	method="ML", 
	data=shelf)
summary(shelf4)
PlotRes(shelf4, shelf, c(3,4) )
AIC(shelf4)	#-1.033

############

remodel1<-lm(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity+SubstrateType, data=whole)
PlotRes(remodel1, whole, c(3,4))
summary(remodel1)

remodel2<-gls(log10(AFDMadjPerMperD)~TadpolePDensity*MFDensity+SubstrateType, method="ML", data=whole)
PlotRes(remodel2, whole, c(3,4))
summary(remodel2)
AIC(remodel2)	#349.3

remodel3<-gls(log10(AFDMadjPerMperD)~TadpolePDensity+MFDensity+SubstrateType, method="ML", data=whole)
PlotRes(remodel3, whole, c(3,4))
summary(remodel3)
AIC(remodel3)	#365.5; keep interaction

remodel4<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	method="ML", data=whole)
summary(remodel4)
PlotRes(remodel4, whole, c(3,4))
AIC(remodel4)	#357

remodel5<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	weights=varIdent(form=~1|SubstrateType),
	method="ML", data=whole)
summary(remodel5)
PlotRes(remodel5, whole, c(3,4))
AIC(remodel5)	#346.73

remodel6<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	weights=varIdent(form=~1|TadpolePDensity),
	method="ML", data=whole)
summary(remodel6)
PlotRes(remodel6, whole, c(3,4))
AIC(remodel6)	#355

remodel7<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	weights=varIdent(form=~1|MFDensity),
	method="ML", data=whole)
summary(remodel7)
PlotRes(remodel7, whole, c(3,4))
AIC(remodel7)	#331.99

MFVar<-varIdent(form=~1|MFDensity)
TadVar<-varIdent(form=~1|TadpolePDensity)
SubstrateVar<-varIdent(form=~1|SubstrateType)
CombVar<-varComb(MFVar, TadVar, SubstrateVar)
remodel8<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	weights=CombVar,
	method="ML", 
	data=whole)
summary(remodel8)
PlotRes(remodel8, whole, c(3,4))
AIC(remodel8)	#324.5201

CombVar<-varComb(MFVar, SubstrateVar)
remodel9<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	weights=CombVar,
	method="ML", 
	data=whole)
summary(remodel9)
PlotRes(remodel9, whole, c(3,4))
AIC(remodel9)	#323.2
#bestbest model

remodel9<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*MFDensity,
	random=~1|SubstrateType, 
	weights=CombVar,
	method="REML", 
	data=whole)
summary(remodel9)
anova(remodel9)



##############
#Results reported in Manuscript:
#Linear mixed-effects model fit by REML
# Data: whole 
#377.1346 415.1994 -180.5673
#
#Random effects:
# Formula: ~1 | SubstrateType
#        (Intercept)  Residual
#StdDev:  0.08813041 0.2672232
#
#Combination of variance functions: 
# Structure: Different standard deviations per stratum
# Formula: ~1 | MFDensity 
# Parameter estimates:
#       0      250 
#1.000000 1.277082 
# Structure: Different standard deviations per stratum
# Formula: ~1 | SubstrateType 
# Parameter estimates:
#BottomTile  ShelfTile 
# 1.0000000  0.8145705 
#Fixed effects: log10(AFDMadjPerMperD) ~ TadpolePDensity * MFDensity 
#                               Value  Std.Error  DF    t-value p-value
#(Intercept)               -1.2324850 0.06470813 860 -19.046833  0.0000
#TadpolePDensity           -0.0149701 0.00152219 860  -9.834561  0.0000
#MFDensity                 -0.0000418 0.00010889 860  -0.384255  0.7009
#TadpolePDensity:MFDensity  0.0000412 0.00000990 860   4.164528  0.0000
# Correlation: 
#                          (Intr) TdplPD MFDnst
#TadpolePDensity           -0.182              
#MFDensity                 -0.163  0.431       
#TadpolePDensity:MFDensity  0.112 -0.615 -0.687
#
#Standardized Within-Group Residuals:
#      Min         Q1        Med         Q3        Max 
#-2.6350907 -0.6614767 -0.0759585  0.6443546  3.8049975 
#
#Number of Observations: 865
#Number of Groups: 2 

TadpoleOnly<-which(whole$TadpolePDensity==16 & whole$MFDensity==0)
mean(log10(whole$AFDMadjPerMperD[TadpoleOnly]))#-1.44
mean((whole$AFDMadjPerMperD[TadpoleOnly]))#0.0411
mean(10^resid(remodel9)[TadpoleOnly])#1.13

MFOnly<-which(whole$TadpolePDensity==0 & whole$MFDensity==250)
mean(log10(whole$AFDMadjPerMperD[MFOnly]))#-1.2
mean((whole$AFDMadjPerMperD[MFOnly]))#0.0813
mean(10^resid(remodel9)[MFOnly])#1.29

Both<-which(whole$TadpolePDensity==16 & whole$MFDensity==250)
mean(log10(whole$AFDMadjPerMperD[Both]))#-1.29
mean((whole$AFDMadjPerMperD[Both]))#0.0797
mean(10^resid(remodel9)[Both])#1.512

None<-which(whole$TadpolePDensity==0 & whole$MFDensity==0)
mean(log10(whole$AFDMadjPerMperD[None]))#-1.189
mean((whole$AFDMadjPerMperD[None]))#0.0817
mean(10^resid(remodel9)[None])#1.262

#raw abundances relative to one another
0.0411/0.0817	#Algae about 50% as abundance in present of tadpoles relative to with no consumers
0.0411/0.0813	#51% tadpole/mayfly
0.0411/0.0797	#52% tadpole/both

0.0797/0.0817	
0.0813/0.0817

############
#rerun with actual ending mayfly abundance and ending algal measures. add 'sample number' field
SampleThree<-whole[which(whole$SampleNumber=="S3"),]

qqnorm(log10(SampleThree$AFDMadjPerMperD))
bartlett.test(log10(AFDMadjPerMperD)~Treatment, data=SampleThree)
bartlett.test(log10(AFDMadjPerMperD)~TadpolePDensity, data=SampleThree)
bartlett.test(log10(AFDMadjPerMperD)~MFDensity, data=SampleThree)

ggplot(whole, aes(x=SampleNumber, y=log10(AFDMadjPerMperD)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=Treatment, y=log10(AFDMadjPerMperD)))+
	geom_point(aes(color=SampleDate))
ggplot(SampleThree, aes(x=Treatment, y=log10(AFDMadjPerMperD)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=as.factor(TankNumber), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=as.factor(SampleDate), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=as.factor(TadpolePDensity), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=as.factor(MFDensity), y=log(AFDMadjPerMperD+0.0001)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=as.factor(FinalMFCount), y=log(AFDMadjPerMperD)))+
	geom_boxplot()
ggplot(SampleThree, aes(x=SubstrateType, y=log(AFDMadjPerMperD)))+
	geom_boxplot()

MFCountModel1<-lm(log10(AFDMadjPerMperD)~TadpolePDensity*FinalMFCount+SubstrateType, data=SampleThree)
PlotRes(MFCountModel1, SampleThree, c(3,4))
summary(MFCountModel1)
AIC(MFCountModel1)	#313

MFCountModel2<-lm(log10(AFDMadjPerMperD)~TadpolePDensity+FinalMFCount+SubstrateType, data=SampleThree)
PlotRes(MFCountModel2, SampleThree, c(3,4))
summary(MFCountModel2)
AIC(MFCountModel2)	#314

MFCountModel3<-gls(log10(AFDMadjPerMperD)~TadpolePDensity*FinalMFCount+SubstrateType, data=SampleThree, method="ML")
PlotRes(MFCountModel3, SampleThree, c(3,4))
summary(MFCountModel3)
AIC(MFCountModel3)	#313.77

MFCountModel4<-lme(log10(AFDMadjPerMperD)~
	TadpolePDensity*FinalMFCount, 
	random=~1|SubstrateType,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel4, SampleThree, c(3,4))
summary(MFCountModel4)
AIC(MFCountModel4)	#317; random intercept not that helpful; but could drop interaction and substrate

MFVar<-varIdent(form=~1|MFDensity)
TadVar<-varIdent(form=~1|TadpolePDensity)
SubstrateVar<-varIdent(form=~1|SubstrateType)
CombVar<-varComb(MFVar, TadVar, SubstrateVar)
MFCountModel5<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity*FinalMFCount+SubstrateType,
	weights=MFVar,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel5, SampleThree, c(3,4))
summary(MFCountModel5)
AIC(MFCountModel5)	#300.0903

MFCountModel6<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity*FinalMFCount+SubstrateType,
	weights=TadVar,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel6, SampleThree, c(3,4))
summary(MFCountModel6)
AIC(MFCountModel6)	#303.89

MFCountModel7<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity*FinalMFCount+SubstrateType,
	weights=SubstrateVar,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel7, SampleThree, c(3,4))
summary(MFCountModel7)
AIC(MFCountModel7)	#309.0672

MFCountModel8<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity*FinalMFCount+SubstrateType,
	weights=CombVar,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel8, SampleThree, c(3,4))
summary(MFCountModel8)
AIC(MFCountModel8)	#291.49

#step down model 8, drop interaction:
MFCountModel9<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+FinalMFCount+SubstrateType,
	weights=CombVar,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel9, SampleThree, c(3,4))
summary(MFCountModel9)
AIC(MFCountModel9)	#290.81; drop MFCount
anova(MFCountModel9)

CombVar<-varComb(TadVar, SubstrateVar)
MFCountModel10<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+SubstrateType,
	weights=CombVar,
	data=SampleThree, 
	method="ML")
PlotRes(MFCountModel10, SampleThree, c(3,4))
summary(MFCountModel10)
AIC(MFCountModel10)	#297.7; dropping MFCount from fixed and variance structure leads to a poorer fit.
anova(MFCountModel10)

#refit modle 9:
MFCountModel9<-gls(log10(AFDMadjPerMperD)~
	TadpolePDensity+FinalMFCount+SubstrateType,
	weights=CombVar,
	data=SampleThree, 
	method="REML")
summary(MFCountModel9)
anova(MFCountModel9)

AFDMsummary<-ddply(SampleThree, "SubstrateType",  summarize,
	meanAFDM=mean(AFDMadjPerMperD), 
	seAFDM=sd(AFDMadjPerMperD)/sqrt(length(AFDMadjPerMperD)))

#-----------------------------------------------------------------------------------------------------------------------------
#patterns across time?
str(whole)
plot(whole$GrowthPeriod)
quantile(whole$AFDM, p=0.0075)
whole[which(whole$AFDM<0.0005),]


whole2<-whole[-which(whole$AFDM<0.0005),]
whole2$GrowthPeriodCat[whole2$GrowthPeriod>19 & whole2$GrowthPeriod<24]<-"22"
whole2$GrowthPeriodCat[whole2$GrowthPeriod>11 & whole2$GrowthPeriod<17]<-"14"
whole2$GrowthPeriodCat[whole2$GrowthPeriod>7 & whole2$GrowthPeriod<10]<-"9"
whole2$GrowthPeriodCat[whole2$GrowthPeriod>2 & whole2$GrowthPeriod<8]<-"7"
boxplot(whole$AFDM~whole2$GrowthPeriodCat)
ggplot(whole2, aes(x=as.numeric(GrowthPeriodCat), y=log(AFDM+0.001)))+
	geom_boxplot(aes(GrowthPeriodCat))
ggplot(whole2, aes(x=GrowthPeriodCat, y=log(AFDM+0.001)))+
	geom_boxplot(aes(color=Treatment))
ggplot(whole2, aes(x=as.numeric(GrowthPeriod), y=log(AFDM)))+
	geom_point()
ggplot(whole2, aes(x=as.numeric(GrowthPeriod), y=AFDM))+
	geom_point(shape=1)

#lines?
tapply(whole$AFDM, whole$GrowthPeriod, length)
whole2$TxOrder<-factor(whole2$Treatment, levels = c("none", "mayflies", "tadpoles", "both"))
ggplot(whole2, aes(x=GrowthPeriod, y=log(AFDM+0.00001)))+
	geom_point()+
	facet_grid(.~TxOrder)+
	geom_smooth(method=lm)

llcoolJ<-ggplot(whole2, aes(x=GrowthPeriod, y=log(AFDM+0.00001), color=TxOrder))+
	geom_smooth(method=lm, alpha=0.25, size=2)+	#alternatively use method=lm
	xlab("Days")+ylab("log(AFDM)")+
	theme(panel.background=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=14, face="bold"),	
		axis.title=element_text(color="black", size=14, face="bold")
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2010_AlgalGrowthRates_byTreatment.png", 
	res=300, width=5, height=5, units="in", restoreConsole=F)
llcoolJ
dev.off()


	

hist(log(whole2$AFDM+0.0001))
shapiro.test(log(whole2$AFDM+0.0001))
qqnorm(log(whole2$AFDM+0.0001))
lm1<-lm(log(AFDM+0.00001)~GrowthPeriod*Treatment,data=whole2)
summary(lm1)
PlotRes(lm1, whole, c(3,3))




#------------------------------------------------------------------------------------------------------------
#############
#Analyze growth at time S3 with respect to treatment, time, and initial algal abundance from S0, as a fixed continuous effect
#create the dataset to use:

#reload "whole":
whole<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Summer 2010/2010MesocosmCSVdata/ExperimentalAlgae_Sample3.csv", header=T)
str(whole)
start<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Summer 2010/2010MesocosmCSVdata/StartingPtAlgaeSamples.csv", header=T)
str(start)
whole2<-rbind(whole, start)
str(whole)
#############
PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	qqnorm(resid(x))
	plot(resid(x)~fitted(x), data=y)
	plot(resid(x)~y[,1], data=y, xlab=colnames(y[1]))
	plot(resid(x)~y[,2], data=y, xlab=colnames(y[2]))
	plot(resid(x)~y[,3], data=y, xlab=colnames(y[3]))
	plot(resid(x)~y[,5], data=y, xlab=colnames(y[5]))
	boxplot(resid(x)~y[,6], data=y, xlab=colnames(y[6]))
	boxplot(resid(x)~y[,7], data=y, xlab=colnames(y[7]))
	plot(resid(x)~y[,8], data=y, xlab=colnames(y[8]))
	}
#PlotRes(model, dataset, plot dimensions)
##############
#remove outliers?
hist(log10(whole$AFDM))
#none seen!
str(whole2)

#create dataframe AlgaeVsT0 that has mean values for each tank at each sample period,
#with respect to the mean at T0, treatment, and number of days:
whole3<-whole2[,c(1,4,12,18)]
meltwhole<-melt(whole3, id=c("TankNumber","SampleNumber", "Treatment"))
TreatMeans<-cast(meltwhole, TankNumber+Treatment~SampleNumber, max)

whole4<-whole2[,c(1,4,18,13)]
meltwhole4<-melt(whole4, id=c("TankNumber","SampleNumber", "Treatment"))
TankSamplePeriods<-cast(meltwhole4, TankNumber~SampleNumber, mean)

whole5<-whole2[,c(1,4,18, 15)]
meltwhole5<-melt(whole5, id=c("TankNumber","SampleNumber", "Treatment"))
TankTadDens<-cast(meltwhole5, TankNumber~SampleNumber, mean)

whole6<-whole2[,c(1,4,18,16)]
meltwhole6<-melt(whole6, id=c("TankNumber","SampleNumber", "Treatment"))
TankMFDens<-cast(meltwhole6, TankNumber~SampleNumber, mean)

whole7<-whole2[,c(1,4,18,17)]
meltwhole7<-melt(whole7, id=c("TankNumber","SampleNumber", "Treatment"))
TankMFFinal<-cast(meltwhole7, TankNumber~SampleNumber, mean)

AlgaeVsT0<-data.frame(rep(TreatMeans$TankNumber, 1), 
		rep(TreatMeans$Treatment, 1), 
		rep(TreatMeans$S0, 1),
		c(TreatMeans$S3),
		c(TankSamplePeriods$S3),
		rep(TankTadDens$S0, 1),
		rep(TankMFDens$S0, 1),
		rep(TankMFFinal$S0, 1))
names(AlgaeVsT0)
names(AlgaeVsT0)<-c("TankNumber", "Treatment", "AlgaeT0", "TreatmentMeans","GrowthPeriod", "TadDens", "MFDens", "MFFinal")
str(AlgaeVsT0)

#analyse AlgaeVsT0
lm<-lm(TreatmentMeans~TadDens*MFFinal+GrowthPeriod+log(AlgaeT0), data=AlgaeVsT0)
PlotRes(lm, AlgaeVsT0, c(2,5))
summary(lm)

lm<-lm(log(TreatmentMeans)~TadDens*MFFinal+GrowthPeriod+log(AlgaeT0), data=AlgaeVsT0)
PlotRes(lm, AlgaeVsT0, c(2,5))
summary(lm)

lm2<-gls(log(TreatmentMeans)~
	TadDens*MFFinal+GrowthPeriod+log(AlgaeT0),
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm2, AlgaeVsT0, c(2,5))
summary(lm2)

TadVar<-varIdent(form=~1|TadDens)
MFVar<-varIdent(form=~1|MFDens)
lm3<-gls(log(TreatmentMeans)~
	TadDens*MFDens+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm3, AlgaeVsT0, c(2,5))
summary(lm3)
anova(lm2, lm3)	#that improves model... 

lm4<-gls(log(TreatmentMeans)~
	TadDens*MFDens+GrowthPeriod+log(AlgaeT0),
	weights=MFVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm4, AlgaeVsT0, c(2,5))
summary(lm4)
anova(lm2, lm4)	#that does not imporve the model

#try lm3 sans interaction:
lm5<-gls(log(TreatmentMeans)~
	TadDens+MFDens+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm5, AlgaeVsT0, c(2,5))
summary(lm5)
anova(lm3, lm5)	#not better, not worse, drop interaction, then drop...MFDens

lm6<-gls(log(TreatmentMeans)~
	TadDens+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm6, AlgaeVsT0, c(2,5))
summary(lm6)
anova(lm5, lm6)	#stop here.

lm6.reml<-gls(log(TreatmentMeans)~
	TadDens+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="REML",
	data=AlgaeVsT0)
PlotRes(lm6.reml, AlgaeVsT0, c(2,5))
summary(lm6.reml)
anova(lm6.reml)
#plot: residuals against tadpole density?
lm6.reml.noTad<-gls(log(TreatmentMeans)~
	GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="REML",
	data=AlgaeVsT0)
PlotRes(lm6.reml.noTad, AlgaeVsT0, c(2,5))
summary(lm6.reml.noTad)

boxplot(resid(lm6.reml.noTad)~AlgaeVsT0$TadDens)


#analyse AlgaeVsT0 using MFFinal abundance
lm<-lm(log(TreatmentMeans)~TadDens*MFFinal+GrowthPeriod+log(AlgaeT0), data=AlgaeVsT0)
PlotRes(lm, AlgaeVsT0, c(2,5))
summary(lm)

lm2<-gls(log(TreatmentMeans)~
	TadDens*MFFinal+GrowthPeriod+log(AlgaeT0),
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm2, AlgaeVsT0, c(2,5))
summary(lm2)

TadVar<-varIdent(form=~1|TadDens)
MFVar<-varIdent(form=~1|MFDens)
lm3<-gls(log(TreatmentMeans)~
	TadDens*MFFinal+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm3, AlgaeVsT0, c(2,5))
summary(lm3)
anova(lm2, lm3)	#that improves model... 

lm4<-gls(log(TreatmentMeans)~
	TadDens*MFFinal+GrowthPeriod+log(AlgaeT0),
	weights=MFVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm4, AlgaeVsT0, c(2,5))
summary(lm4)
anova(lm2, lm4)	#that does not improve model... 

lm5<-gls(log(TreatmentMeans)~
	TadDens+MFFinal+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="ML",
	data=AlgaeVsT0)
PlotRes(lm5, AlgaeVsT0, c(2,5))
summary(lm5)
anova(lm3, lm5)	#significantly worse model.  so stick with lm4, cannot drop anything else.

lm4.reml<-gls(log(TreatmentMeans)~
	TadDens*MFFinal+GrowthPeriod+log(AlgaeT0),
	weights=TadVar,
	method="REML",
	data=AlgaeVsT0)
PlotRes(lm4.reml, AlgaeVsT0, c(2,5))
summary(lm4.reml)

#plotting: - with goal of displaying the results of the interaction bw tadpole density and mf abundance
so_rad<-ggplot(AlgaeVsT0, aes(x=MFFinal, y=log(TreatmentMeans)))+
	geom_smooth(method=lm, aes(linetype=as.factor(TadDens), color=as.factor(TadDens)), size=2, fill="light grey", fullrange=TRUE)+
	geom_point(aes(color=as.factor(TadDens), shape=as.factor(TadDens)), size=5)+
	scale_linetype_manual(values=c("dotted", "solid"), labels=c(" no tadpoles ", " tadpoles "))+
	scale_color_manual(values=c("black","black"), labels=c(" no tadpoles ", " tadpoles "))+
	scale_shape_manual(values=c(24,16), labels=c(" no tadpoles ", " tadpoles "))+
	guides(linetype=guide_legend(keywidth=10, keyheight=0.5, label.position="bottom", label.hjust = 0.5))+
	xlab("Final Mayfly Abundance")+ylab("ln (mg AFDM)")+
	theme(legend.position=c("bottom"),
		legend.title=element_blank(),
		legend.text=element_text(size=14, face="bold"),
		panel.background=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=14, face="bold"),	
		axis.title=element_text(color="black", size=14, face="bold")
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2010_logAlgalGrowth_byMFAbundance_byTadpolePresence.png", 
	res=300, width=5, height=5, units="in", restoreConsole=F)
so_rad
dev.off()




