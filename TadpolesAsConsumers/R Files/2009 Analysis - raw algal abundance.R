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
library(scales)
library(gridExtra)
#############
#load and organize data:
whole<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\ConsumerResourceAnalysis_for_analysis_withTadpoleAFDM.csv", header=T)
str(whole)
#############

PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	qqnorm(resid(x))
	plot(resid(x)~fitted(x), data=y)
	plot(resid(x)~y[,1], data=y, xlab=colnames(y[1]))
	plot(resid(x)~y[,3], data=y, xlab=colnames(y[3]))
	plot(resid(x)~y[,4], data=y, xlab=colnames(y[4]))
	boxplot(resid(x)~y[,9], data=y, xlab=colnames(y[9]))
	boxplot(resid(x)~y[,12], data=y, xlab=colnames(y[12]))
	plot(resid(x)~y[,16], data=y, xlab=colnames(y[16]))
	plot(resid(x)~y[,17], data=y, xlab=colnames(y[17]))
	}
#PlotRes(model, dataset, plot dimensions)

#############
#calculations:
	#1. Experimental AFDM per m^2
	#2. Control AFDM per m^2 per day
	#3. Differences between Experimental and Control Abundance
	#4. Experimental AFDM per m^2 per day
	#5. Control AFDM per m^2 per day
	#6. Differences between Experimental and Control Grow Rate
	#7. Change in Tadpole Wet Weight
	#8. Change in MF Density (Density-Count)
	#9. Change in MF Biomass

#1. Experimental AFDM per m^2 
	#each tile is 2.4 * 2.4 cm = 0.024 * 0.024 m
	#there are 24 experimental tiles
ExperimentalTileArea<-(0.024^2)*24	#0.013824
ExperimentalAFDMpM<-whole$ExperimentalAFDM/ExperimentalTileArea

#2. Control AFDM per m^2 
	#each tile is 2.4 * 2.4 cm = 0.024 * 0.024 m
	#there are 12 Control tiles
ControlTileArea<-(0.024^2)*12		#0.006912
ControlAFDMpM<-whole$ControlAFDM/ControlTileArea

#3. Differences between Experimental and Control Abundance
ExperimentalControlDifference_Abundance<-ExperimentalAFDMpM-ControlAFDMpM

#4. Experimental AFDM per m^2 per day
ExperimentalAFDMpMpD<-whole$ExperimentalAFDM/ExperimentalTileArea/whole$SamplePeriod

#5. Control AFDM per m^2 per day
ControlAFDMpMpD<-whole$ControlAFDM/ControlTileArea/whole$SamplePeriod

#6. Differences between Experimental and Control Growth Rate
ExperimentalControlDifference_Growth<-ExperimentalAFDMpMpD-ControlAFDMpMpD

#7. Change in Tadpole Wet Weight
TadpoleWeightChange<-whole$TadWtStart-whole$TadWtEnd

#8. Change in MF Density (Density-Count)
MFDensChange<-whole$MFDens-whole$MFCount

#9. Change in MF Biomass
MFBiomassChange<-whole$MFBiomassStart-whole$MFBiomassEnd

#lakeControl tiles
LakeControl<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\LakeControlTilesAFDM.csv", header=T)
#2. Control AFDM per m^2 
	#each tile is 2.4 * 2.4 cm = 0.024 * 0.024 m
	#there are 12 Control tiles
ControlTileArea<-(0.024^2)*12		#0.006912
LakeAFDMpM<-LakeControl$AFDM/ControlTileArea
#compare experimental, control, and lake-control
df<-data.frame(c(ControlAFDMpM, ExperimentalAFDMpM, LakeAFDMpM), c(rep("ControlAFDMpM", 102),rep("ExperimentalAFDMpM", 102), rep("Lake",8)))
str(df)
df<-rename(df, c("c.ControlAFDMpM..ExperimentalAFDMpM..LakeAFDMpM."="Value" , "c.rep..ControlAFDMpM...102...rep..ExperimentalAFDMpM...102..."="Type"))
ggplot(df, aes(y=log(Value), x=Type))+
	geom_boxplot()
TukeyHSD(aov(Value~Type, data=df))
l<-aggregate(df, by=list(df$Type), mean)
m<-aggregate(df, by=list(df$Type), sd)
m$Value/c(102,102,8)
ll<-aggregate(LakeControl, by=list(LakeControl$LakeID), mean)
ml<-aggregate(LakeControl, by=list(LakeControl$LakeID), sd)
ml$AFDM/c(1,4,3)


#############
#basic plotting
#groups:
leconte<-which(whole$LakeName=="LeConte")
spur<-which(whole$LakeName=="Spur")

plot(ExperimentalAFDMpMpD, ControlAFDMpMpD)

ggplot(whole, aes(x=LakeName, y=ExperimentalAFDMpMpD))+
  geom_boxplot()
ggplot(whole, aes(ExperimentalAFDMpMpD))+
  geom_histogram()
ggplot(whole[leconte,], aes(ExperimentalAFDMpMpD))+
	geom_histogram()
hist(ExperimentalAFDMpMpD[leconte])
hist(ExperimentalAFDMpMpD[spur])
mean(ExperimentalAFDMpMpD)
sd(ExperimentalAFDMpMpD)

ggplot(whole, aes(factor(LakeName), log(ExperimentalAFDMpMpD)))+
  geom_boxplot()

ggplot(whole, aes(x=as.factor(SamplePeriod), y=log(ExperimentalAFDMpM)))+
	geom_boxplot()

#############
#summaries:
#groups:
leconte<-which(whole$LakeName=="LeConte")
spur<-which(whole$LakeName=="Spur")

#between lake comparison of algal abundance
mean(ExperimentalAFDMpMpD[leconte])#
sd(ExperimentalAFDMpMpD[leconte])#
mean(ExperimentalAFDMpMpD[spur])#
sd(ExperimentalAFDMpMpD[spur])#
(mean(ExperimentalAFDMpMpD[leconte]))/(mean(ExperimentalAFDMpMpD[spur]))#

anova(aov(ExperimentalAFDMpMpD~whole$LakeName))#
wholesummarized<-summarySE(whole, measurevar="ExperimentalAFDM", groupvars=c("LakeName","SampleNumber"))
ggplot(wholesummarized, aes(y=ExperimentalAFDM, x=LakeName, color=SampleNumber))+
	geom_errorbar(aes(ymin=ExperimentalAFDM-se, ymax=ExperimentalAFDM+se), width=.1, position=position_dodge(0.1))
#requires running/calling summarySE from another file.

summary(aov(log(ExperimentalAFDMpM)~as.factor(SamplePeriod), data=whole))
TukeyHSD(aov(log(ExperimentalAFDMpM)~as.factor(SamplePeriod), data=whole))

summary(aov(log(ExperimentalAFDMpM)~as.factor(SampleNumber), data=whole))
TukeyHSD(aov(log(ExperimentalAFDMpM)~as.factor(SampleNumber), data=whole))

ggplot(whole, aes(x=SamplePeriod, y=log(ExperimentalAFDMpM)))+
	geom_boxplot(aes(color=SampleNumber))+
	facet_grid(.~LakeName)
	

############
#models of Raw algal abundance (not subtracted from control)
ggplot(whole, aes(x=CageNumber, y=ExperimentalAFDMpM))+
	geom_point(aes(color=SampleNumber))+
	facet_grid(.~LakeName)
hist(ExperimentalAFDMpM)
hist(log(ExperimentalAFDMpM))
shapiro.test(log(ExperimentalAFDMpM))
shapiro.test((ExperimentalAFDMpM)^(1/4))

hist(log(ExperimentalAFDMpM[leconte]))
shapiro.test(log(ExperimentalAFDMpM[leconte]))
hist((ExperimentalAFDMpM[leconte])^(1/2))
shapiro.test((ExperimentalAFDMpM[leconte])^(1/2))	#sqrt transform normalizes leconte
hist(log(ExperimentalAFDMpM[spur]))
shapiro.test(log(ExperimentalAFDMpM[spur]))
hist((ExperimentalAFDMpM[spur]))
boxplot((ExperimentalAFDMpM[spur])^(1/10))
shapiro.test((ExperimentalAFDMpM[spur])^(1/10))	#

lm1<-lm(ExperimentalAFDMpM~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,4))
shapiro.test(resid(lm1))
AIC(lm1)

lm2<-lm(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm2)
PlotRes(lm2, whole, c(3,4))
shapiro.test(resid(lm2))	#p=0.08
AIC(lm2)	#327.75

#the log transform of the raw algal abundance achieves normality of the residuals
#random effects for 

lm3<-gls(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation+SampleNumber, 
	method="ML",
	data=whole)
summary(lm3)
PlotRes(lm3, whole, c(3,4))
shapiro.test(resid(lm3))	#p=0.09
AIC(lm3)	#327.75

lm4<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	method="ML",
	data=whole)
summary(lm4)
PlotRes(lm4, whole, c(3,4))
shapiro.test(resid(lm4))	#p=0.015
AIC(lm4)	#329.5
anova(lm3, lm4)	#p=0.05, so, not really worse, visually better fit of fitted to residuals, logically better

lm5<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+Silt+SamplePeriod+Radiation+SampleNumber, 
	random=~1|LakeName,
	method="ML",
	data=whole)
summary(lm5)
PlotRes(lm5, whole, c(3,4))
shapiro.test(resid(lm5))	#0.032
AIC(lm5)	#334

lm6<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+SamplePeriod+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	method="ML",
	data=whole)
summary(lm6)
PlotRes(lm6, whole, c(3,4))
shapiro.test(resid(lm6))	#.005
AIC(lm6)	#335.5

#so use lm4, add heterogeneity
LakeVar<-varIdent(form=~1|LakeName)
lm7<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	weights=LakeVar,
	method="ML",
	data=whole)
summary(lm7)
PlotRes(lm7, whole, c(3,4))
shapiro.test(resid(lm7))	#p=0.003
AIC(lm7)	#300.4
anova(lm4, lm7)	#p<0.0001

SampleVar<-varIdent(form=~1|SampleNumber)
lm8<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	weights=SampleVar,
	method="ML",
	data=whole)
summary(lm8)
PlotRes(lm8, whole, c(3,4))
shapiro.test(resid(lm8))	#p=0.005
AIC(lm8)	#327
anova(lm4, lm8)	#p=0.052, only a marginal improvement

PeriodVar<-varIdent(form=~1|SamplePeriod)
lm9<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	weights=PeriodVar,
	method="ML",
	data=whole)
summary(lm9)
PlotRes(lm9, whole, c(3,4))
shapiro.test(resid(lm9))	#p=0.00013
AIC(lm9)	#300.91
anova(lm4, lm9)	#p<0.0001; but residuals dont look very good.

TadVar<-varIdent(form=~1|TadDens)
lm10<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	weights=TadVar,
	method="ML",
	data=whole)
summary(lm10)
PlotRes(lm10, whole, c(3,4))
shapiro.test(resid(lm10))	#p=0.006
AIC(lm10)	#331
anova(lm4, lm10)	#p=0.23, not really an improvement

MFVar<-varIdent(form=~1|MFDens)
lm11<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	weights=MFVar,
	method="ML",
	data=whole)
summary(lm11)
PlotRes(lm11, whole, c(3,4))
shapiro.test(resid(lm11))	#p=0.0022
AIC(lm11)	#330.4
anova(lm4, lm11)	#p=0.17

#improvements using Period Var, Lake Var, and SampleVar
CombVar<-varComb(SampleVar, LakeVar)
lm12<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	weights=CombVar,
	method="ML",
	data=whole)
summary(lm12)
PlotRes(lm12, whole, c(3,4))
shapiro.test(resid(lm12))	#p=0.0001
AIC(lm12)	#303.6
anova(lm4, lm12)	#p<0.0001; 
#all combinations of variance structure offer improvement in AIC, 
#but with a cost to normality of residuals
#the best is the model without variance structure... lm4
summary(lm4)

#stepdown, dropping consumer interaction first:
lm13<-lme(log(ExperimentalAFDMpM)~
	TadDens+MFDens+LakeName+SamplePeriod+Silt+Radiation, 
	random=~1|SampleNumber,
	method="ML",
	data=whole)
summary(lm13)
PlotRes(lm13, whole, c(3,4))
shapiro.test(resid(lm13))	#p=0.025
AIC(lm13)	#328
anova(lm4, lm13)	#p=0.25 no difference without interaction
#try keeping interaction, drop SamplePeriod

lm14<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1|SampleNumber,
	method="ML",
	data=whole)
summary(lm14)
PlotRes(lm14, whole, c(3,4))
shapiro.test(resid(lm14))	#p=0.02
AIC(lm14)	#326
anova(lm13, lm14)	#p=0.78
#drop interaction

lm15<-lme(log(ExperimentalAFDMpM)~
	TadDens+MFDens+LakeName+Silt+Radiation, 
	random=~1|SampleNumber,
	method="ML",
	data=whole)
summary(lm15)
PlotRes(lm15, whole, c(3,4))
shapiro.test(resid(lm15))	#p=0.18
AIC(lm15)	#326.8
anova(lm14, lm15)	#p=0.5
#drop silt

lm16<-lme(log(ExperimentalAFDMpM)~
	TadDens+MFDens+LakeName+Radiation, 
	random=~1|SampleNumber,
	method="ML",
	data=whole)
summary(lm16)
PlotRes(lm16, whole, c(3,4))
shapiro.test(resid(lm16))	#p=0.012
AIC(lm16)	#326
anova(lm15, lm16)	#p=0.3 no improvement
#drop radiation

lm17<-lme(log(ExperimentalAFDMpM)~
	TadDens+MFDens+LakeName, 
	random=~1|SampleNumber,
	method="ML",
	data=whole)
summary(lm17)
PlotRes(lm17, whole, c(3,4))
shapiro.test(resid(lm17))	#p=0.007
AIC(lm17)	#324.9
anova(lm16, lm17)	#p=0.34

#refit with REML:
lm17.reml<-lme(log(ExperimentalAFDMpM)~
	TadDens+MFDens+LakeName, 
	random=~1|SampleNumber,
	method="REML",
	data=whole)
summary(lm17.reml)
PlotRes(lm17.reml, whole, c(3,4))
shapiro.test(resid(lm17.reml))	#p=0.011
intervals(lm17.reml)

############
#PLOTTING:
#
#plotting residuals:
#run sans MF, plot resid against MF
lm17.noMF<-lme(log(ExperimentalAFDMpM)~
	TadDens+LakeName, 
	random=~1|SampleNumber,
	method="REML",
	data=whole)
#run sans Tad, plot resid against Tad
lm17.noTad<-lme(log(ExperimentalAFDMpM)~
	MFDens+LakeName, 
	random=~1|SampleNumber,
	method="REML",
	data=whole)
lm17.noLake<-lme(log(ExperimentalAFDMpM)~
	TadDens+MFDens, 
	random=~1|SampleNumber,
	method="REML",
	data=whole)

par(mfrow=c(1,3))
#################
#Plots of original raw algal abundances
library(gridExtra)
#rename S2... etc to Block 1
whole2<-whole
levels(whole2$SampleNumber)[levels(whole2$SampleNumber)=="S2"]<-"Early August"
levels(whole2$SampleNumber)[levels(whole2$SampleNumber)=="S3"]<-"Late August"
levels(whole2$SampleNumber)[levels(whole2$SampleNumber)=="S4"]<-"Mid September"

#boxplots:
l<-ggplot(whole2, aes(y=log(ExperimentalAFDMpM+0.00001), x=as.factor(TadDens)))+
	geom_boxplot(outlier.shape = 16,outlier.size = 4, size=1.5)+
	facet_grid(LakeName~., scales="free_y")+
	xlab("Tadpole Abundance")+ylab("ln (mg algal AFDM)")+
	guides(fill=F)+
	theme(panel.background = element_rect(color="grey", fill="white"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.title=element_text(color="black", face="bold", size=18),
		axis.text.y=element_text(color="black", face="bold", size=14),
		axis.text.x=element_text(color="black", face="bold", size=14), 
		axis.title.x=element_text(color="black", face="bold", size=18),
		axis.title.y=element_text(color="black", face="bold", size=18), 
		strip.background=element_blank(),
		strip.text=element_blank(),
		legend.text=element_text(face="bold", size=14),
		legend.title=element_text(face="bold", size=14), 
		axis.line = element_line(colour="black", size=1.5)
		)
m<-ggplot(whole2, aes(y=log(ExperimentalAFDMpM+0.00001), x=as.factor(MFDens)))+
	geom_boxplot(outlier.shape = 16,outlier.size = 4, size=1.5)+
	facet_grid(LakeName~., scales="free_y")+
	xlab("Mayfly Abundance")+ylab("")+
	guides(fill=F)+
	theme(panel.background = element_rect(color="grey", fill="white"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.title=element_text(color="black", face="bold", size=18),
		axis.text.y=element_blank(),
		axis.text.x=element_text(color="black", face="bold", size=14), 
		axis.title.x=element_text(color="black", face="bold", size=18),
		strip.background=element_rect(fill="light grey"),
		strip.text=element_text(face="bold", size=14), 
		axis.line = element_line(colour="black", size=1.5),
		axis.line.y = element_blank(),
		axis.line.x = element_line(colour="black", size=1.5)
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_rawalgalabundace_byconsumer_bylake.png", 
	res=300, width=6, height=9, units="in", restoreConsole=F)
grid.arrange(l,m, ncol=2)
dev.off()


#heat map tile plot:
one<-ggplot(whole2, aes(x=as.factor(MFDens), y=as.factor(TadDens)), group=LakeName)+
	geom_tile(aes(fill=log(ExperimentalAFDMpM)))+
	scale_fill_gradient2(midpoint=0, low="white", mid=muted("green"), high="brown", limits=c(-4.5,2.5), guide="legend")+
	xlab("Mayfly Abundance")+ylab("Tadpole Abundance")+guides(fill=guide_legend(title=NULL))+
	theme(legend.position="right", 
		panel.background = element_rect(fill="white"),
		panel.border=element_rect(fill=NA, color="light grey", size=20),
		axis.text.x=element_text(color="black", face="bold", size=18), 
		axis.text.y=element_text(color="black", face="bold", size=18),
		axis.title.x=element_text(color="black", face="bold", size=20),
		axis.title.y=element_text(color="black", face="bold", size=20), 
		strip.background=element_rect(fill="light grey"),
		strip.text=element_text(face="bold", size=20), 
		legend.text=element_text(face="bold", size=18),
		axis.ticks=element_blank()
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_logalgalabundace_AllOne_heatmap.png", 
	res=300, width=9, height=6, units="in", restoreConsole=F)
one
dev.off()



#both lakes together, all samples summarized as one...
both<-ggplot(whole2, aes(x=as.factor(MFDens), y=as.factor(TadDens)), group=LakeName)+
	geom_tile(aes(fill=log(ExperimentalAFDMpM)))+
	facet_grid(.~LakeName)+
	scale_fill_gradient2(midpoint=0, low="white", mid=muted("green"), high="brown", limits=c(-4.5,2.5), guide="legend")+
#	scale_fill_gradient2(midpoint=-1.25, low=muted("blue"), mid="white", high=muted("green"), limits=c(-4.5,2.5), guide="legend")+
	xlab("Mayfly Abundance")+ylab("Tadpole Abundance")+guides(fill=guide_legend(title=NULL))+
	theme(legend.position="right", 
		panel.background = element_rect(fill="white"),
		panel.border=element_rect(fill=NA, color="light grey", size=10),
		axis.text.x=element_text(color="black", face="bold", size=18), 
		axis.text.y=element_text(color="black", face="bold", size=18),
		axis.title.x=element_text(color="black", face="bold", size=20),
		axis.title.y=element_text(color="black", face="bold", size=20), 
		strip.background=element_rect(fill="light grey"),
		strip.text=element_text(face="bold", size=20), 
		legend.text=element_text(face="bold", size=18),
		axis.ticks=element_blank()
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_logalgalabundace_bylakeheatmap.png", 
	res=300, width=9, height=6, units="in", restoreConsole=F)
both
dev.off()

#lakes and blocks as facets:
HeatMap<-ggplot(whole2, aes(x=as.factor(MFDens), y=as.factor(TadDens) ))+
	geom_tile(aes(fill=log(ExperimentalAFDMpM)))+
	scale_fill_gradient2(midpoint=0, low="white", mid=muted("green"), high="brown", limits=c(-4.5,2.5), guide="legend")+
#	scale_fill_gradient2(midpoint=-1.25, low=muted("blue"), mid="white", high=muted("green"), limits=c(-4.5,2.5), guide="legend")+
	facet_grid(LakeName~SampleNumber, scales=("free_y"))+
	xlab("Mayfly Abundance")+ylab("Tadpole Abundance")+guides(fill=guide_legend(title=NULL))+
	theme(legend.position="right", 
		panel.background = element_rect(fill="white"),
		panel.border=element_rect(fill=NA, color="light grey", size=7),
		axis.text.x=element_text(color="black", face="bold", size=18), 
		axis.text.y=element_text(color="black", face="bold", size=18),
		axis.title.x=element_text(color="black", face="bold", size=20),
		axis.title.y=element_text(color="black", face="bold", size=20), 
		strip.background=element_rect(fill="light grey"),
		strip.text=element_text(face="bold", size=20), 
		legend.text=element_text(face="bold", size=18),
		axis.ticks=element_blank()
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_logalgalabundance_bylakebyblock_heatmap.png", 
	res=300, width=9, height=6, units="in", restoreConsole=F)
HeatMap
dev.off()






































################
#allow different slopes?
ggplot(whole, aes(y=resid(lm1), x=as.factor(MFDens)))+
	geom_boxplot(aes(color=MFDens))+
	facet_grid(.~LakeName)
dev.new()
ggplot(whole, aes(y=resid(lm1), x=as.factor(TadDens)))+
	geom_boxplot(aes(color=TadDens))+
	facet_grid(.~LakeName)
#the two consumers certainly appear to have different effects in each lake, though that may largely be result of the different variances

lmRS1<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+SamplePeriod+Silt+Radiation+SampleNumber, 
	random=~1+TadDens|LakeName, 
	method="ML",
	data=whole)
summary(lmRS1)
PlotRes(lmRS1, whole, c(3,4))
shapiro.test(resid(lmRS1))	#
AIC(lmRS1)	#337 with nice looking residuals

lmRS2<-lme(log(ExperimentalAFDMpM)~
	TadDens*MFDens+SamplePeriod+Silt+Radiation+SampleNumber, 
	random=~1+MFDens|LakeName, 
	method="ML",
	data=whole)
summary(lmRS2)
PlotRes(lmRS2, whole, c(3,4))
shapiro.test(resid(lmRS2))	#0.005
AIC(lmRS2)	#333.1; the residuals are uglier

#random slope does not help.


#with biomass?
























