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
	plot(resid(x)~y[,9], data=y, xlab=colnames(y[9]))
	plot(resid(x)~y[,12], data=y, xlab=colnames(y[12]))
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

#############
#basic plotting
#groups:
leconte<-which(whole$LakeName=="LeConte")
spur<-which(whole$LakeName=="Spur")

plot(ExperimentalAFDMpM, ControlAFDMpM)

ggplot(whole, aes(factor(LakeName), ExperimentalControlDifference_Abundance))+
  geom_point()
ggplot(whole, aes(ExperimentalControlDifference_Abundance))+
  geom_histogram()
ggplot(whole[leconte,], aes(ExperimentalControlDifference_Abundance))+
	geom_histogram()
hist(ExperimentalControlDifference_Abundance[leconte])
hist(ExperimentalControlDifference_Abundance[spur])
mean(ExperimentalControlDifference_Abundance)
sd(ExperimentalControlDifference_Abundance)


ggplot(whole, aes(factor(LakeName), log10(ExperimentalAFDMpM)))+
  geom_boxplot()

#############
#summaries:
#groups:
leconte<-which(whole$LakeName=="LeConte")
spur<-which(whole$LakeName=="Spur")

#between lake comparison of algal abundance
mean(ExperimentalAFDMpM[leconte])#
sd(ExperimentalAFDMpM[leconte])#
mean(ExperimentalAFDMpM[spur])#
sd(ExperimentalAFDMpM[spur])#
mean(mean(ExperimentalAFDMpM[leconte]))/(mean(ExperimentalAFDMpM[spur]))#

anova(aov(ExperimentalAFDMpM~whole$LakeName))#
wholesummarized<-summarySE(whole, measurevar="ExperimentalAFDM", groupvars=c("LakeName","SampleNumber"))
ggplot(wholesummarized, aes(y=ExperimentalAFDM, x=LakeName, color=SampleNumber))+
	geom_errorbar(aes(ymin=ExperimentalAFDM-se, ymax=ExperimentalAFDM+se), width=.1, position=position_dodge(0.1))
#requires running/calling summarySE from another file.

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






############
#models of Algal Abundance:
lm1<-lm(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,4))
AIC(lm1)	#362.03

#outliers?
which(whole$ExperimentalAFDM<=0.001)
which(resid(lm1)<=-5)	#72
which(resid(lm1)>=2.5)	#73, 80, 88
whole<-whole[-72,]
ExperimentalControlDifference_Abundance<-ExperimentalControlDifference_Abundance[-72]
hist(ExperimentalControlDifference_Abundance[leconte])
shapiro.test(ExperimentalControlDifference_Abundance[leconte])#not normal
hist(ExperimentalControlDifference_Abundance[spur])
shapiro.test(ExperimentalControlDifference_Abundance[spur])#not normal

lm1<-lm(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,4))
AIC(lm1)	#362.03

lm1.1<-gls(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	method="ML",
	data=whole)
summary(lm1.1)
PlotRes(lm1.1, whole, c(3,4))
AIC(lm1.1)	#362.03

lm2<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1|SampleNumber,
	data=whole,
	method="ML"
	)
summary(lm2)
PlotRes(lm2, whole, c(3,4))
AIC(lm2)	#367.23

lm3<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation+SampleNumber, 
	random=~1|LakeName,
	data=whole,
	method="ML"
	)
summary(lm3)
PlotRes(lm3, whole, c(3,4))
AIC(lm3)	#362.07

lm4<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	data=whole,
	method="ML"
	)
summary(lm4)
PlotRes(lm4, whole, c(3,4))
AIC(lm4)	#367.1

LakeVar<-varIdent(form=~1|LakeName)
SampleVar<-varIdent(form=~1|SampleNumber)
lm5<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=LakeVar,
	data=whole,
	method="ML"
	)
summary(lm5)
PlotRes(lm5, whole, c(3,4))
AIC(lm5)	#237.95

lm6<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=SampleVar,
	data=whole,
	method="ML"
	)
summary(lm6)
PlotRes(lm6, whole, c(3,3))
AIC(lm6)	#367.4

lm7<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm7)
PlotRes(lm7, whole, c(3,3))
AIC(lm7)	#232.54

lm8<-gls(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation, 
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm8)
PlotRes(lm8, whole, c(3,3))
AIC(lm8)	#233.8

lm9<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm9)
PlotRes(lm9, whole, c(3,3))
AIC(lm9)	#230.95; drop radiation

lm10<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm10)
PlotRes(lm10, whole, c(3,3))
AIC(lm10)	#228.9

lm11<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm11)
PlotRes(lm11, whole, c(3,3))
AIC(lm11)	#229.5172

lm12<-lme(ExperimentalControlDifference_Abundance~
	TadDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm12)
PlotRes(lm12, whole, c(3,3))
AIC(lm12)	#238.8

lm13<-lme(ExperimentalControlDifference_Abundance~
	MFDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm13)
PlotRes(lm13, whole, c(3,3))
AIC(lm13)	#229.72

lm10.reml<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="REML"
	)
summary(lm10.reml)
PlotRes(lm10.reml, whole, c(3,3))
AIC(lm10.reml)	#


##############
#plotting:
whole<-data.frame(whole, ExperimentalControlDifference_Abundance)
LeConte<-which(whole$LakeName=="LeConte")
Spur<-which(whole$LakeName=="Spur")
p<-ggplot(whole[LeConte,], aes(x=as.factor(TadDens), y=ExperimentalControlDifference_Abundance, color=SampleNumber))
p+geom_boxplot()
p<-ggplot(whole[LeConte,], aes(x=as.factor(MFDens), y=ExperimentalControlDifference_Abundance, color=SampleNumber))
p+geom_boxplot()
p<-ggplot(whole[Spur,], aes(x=as.factor(TadDens), y=ExperimentalControlDifference_Abundance, color=SampleNumber))
p+geom_boxplot()
p<-ggplot(whole[Spur,], aes(x=as.factor(MFDens), y=ExperimentalControlDifference_Abundance, color=SampleNumber))
p+geom_boxplot()

q<-ggplot(whole[LeConte,], aes(x=as.factor(MFDens), y=as.factor(TadDens) ))
	q+
	geom_tile(aes(fill=ExperimentalControlDifference_Abundance))+
	scale_fill_gradient2(midpoint=0, low="black", mid="grey", high="white", limits=c(-2,2))+
	facet_grid(.~SampleNumber)
z<-ggplot(whole[Spur,], aes(x=as.factor(MFDens), y=as.factor(TadDens) ))
	z+
	geom_tile(aes(fill=ExperimentalControlDifference_Abundance))+
	scale_fill_gradient2(midpoint=0, low="black", mid="grey", high="white", limits=c(-11,11))+
	facet_grid(.~SampleNumber)
#using residuals does not give a different or stronger pattern

##############
#models using log-modulus transformed 

#transforms?
#add a constant and log transform
shiftlog<-log(ExperimentalControlDifference_Abundance+10)
hist(shiftlog)
qqnorm(shiftlog)
hist(shiftlog[leconte])
qqnorm(shiftlog[leconte])
shapiro.test(shiftlog[leconte])#
hist(shiftlog[spur])
qqnorm(shiftlog[spur])
shapiro.test(shiftlog[spur])#

#log-modulus, as in http://www.statsblogs.com/2014/07/14/a-log-transformation-of-positive-and-negative-values/
#J.A.John and N.R.Draper 1980 J. Roy. Stat. Soc. Series C
#"An alternative family of transformations"

logmod<-sign(ExperimentalControlDifference_Abundance)*log(abs(ExperimentalControlDifference_Abundance)+1)
hist(logmod)
qqnorm(logmod)
hist(logmod[leconte])
qqnorm(logmod[leconte])
shapiro.test(logmod[leconte])#nearly normal - certainly an improvement
hist(logmod[spur])
qqnorm(logmod[spur])
shapiro.test(logmod[spur])#nearl normal

lm1<-lm(logmod~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,4))
AIC(lm1)	#171.63

lm2<-gls(logmod~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	method="ML",
	data=whole)
summary(lm2)
PlotRes(lm2, whole, c(3,4))
AIC(lm2)	#171.63

lm3<-lme(logmod~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1|SampleNumber,
	data=whole,
	method="ML"
	)
summary(lm3)
PlotRes(lm3, whole, c(3,4))
AIC(lm3)	#178.2; no help

lm4<-lme(logmod~
	TadDens*MFDens+SampleNumber+Silt+Radiation, 
	random=~1|LakeName,
	data=whole,
	method="ML"
	)
summary(lm4)
PlotRes(lm4, whole, c(3,4))
AIC(lm4)	#172.22; no big help

lm5<-lme(logmod~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	data=whole,
	method="ML"
	)
summary(lm5)
PlotRes(lm5, whole, c(3,4))
AIC(lm5)	#179.5; maybe a little more normal

LakeVar<-varIdent(form=~1|LakeName)
SampleVar<-varIdent(form=~1|SampleNumber)
lm6<-gls(logmod~
	TadDens*MFDens+Silt+Radiation+LakeName+SampleNumber, 
	weights=LakeVar,
	data=whole,
	method="ML"
	)
summary(lm6)
PlotRes(lm6, whole, c(3,4))
AIC(lm6)	#98.74

lm7<-gls(logmod~
	TadDens*MFDens+Silt+Radiation+LakeName+SampleNumber, 
	weights=SampleVar,
	data=whole,
	method="ML"
	)
summary(lm7)
PlotRes(lm7, whole, c(3,4))
AIC(lm7)	#175.03

lm8<-gls(logmod~
	TadDens*MFDens+Silt+Radiation+LakeName+SampleNumber, 
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm8)
PlotRes(lm8, whole, c(3,4))
AIC(lm8)	#95.35

lm9<-lme(logmod~
	TadDens*MFDens+Silt+Radiation+LakeName, 
	random=~1|SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm9)
PlotRes(lm9, whole, c(3,4))
AIC(lm9)	#103.599

lm10<-lme(logmod~
	TadDens*MFDens+Silt+Radiation+SampleNumber, 
	random=~1|LakeName,
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm10)
PlotRes(lm10, whole, c(3,4))
AIC(lm10)	#95.36

lm11<-lme(logmod~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm11)
PlotRes(lm11, whole, c(3,4))
AIC(lm11)	#106.36

#lm8 is best.... start dropping terms:
#test Sample Number
lm8.noSampleNumber<-gls(logmod~
	TadDens*MFDens+Silt+Radiation+LakeName, 
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm8.noSampleNumber)
PlotRes(lm8.noSampleNumber, whole, c(3,4))
anova(lm8, lm8.noSampleNumber)	#aic=110, p=0.0001; dont drop block.

lm8.noLake<-gls(logmod~
	TadDens*MFDens+Silt+Radiation+SampleNumber, 
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm8.noLake)
PlotRes(lm8.noLake, whole, c(3,4))
anova(lm8, lm8.noLake)	#aic=, p=0.927, drop lake, or put it in random component; then drop radiation

lm12<-gls(logmod~
	TadDens*MFDens+Silt+SampleNumber, 
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm12)
PlotRes(lm12, whole, c(3,4))
anova(lm8, lm12)	#aic=91.5, p=0.94; drop interaction

lm13<-gls(logmod~
	TadDens+MFDens+Silt+SampleNumber, 
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="ML"
	)
summary(lm13)
PlotRes(lm13, whole, c(3,4))
anova(lm8, lm13)	#aic=89.9, p=0.88; don't drop anything else. Reift with REML:

lm13.reml<-gls(logmod~
	TadDens+MFDens+Silt+SampleNumber, 
	weights=varComb(LakeVar, SampleVar),
	data=whole,
	method="REML"
	)
summary(lm13.reml)
lm13.reml.nosamplenumber<-gls(logmod~TadDens+MFDens+Silt, weights=varComb(LakeVar, SampleVar),data=whole,method="REML")
anova(lm13.reml, lm13.reml.nosamplenumber)	#for SampleNumber, p=0.0045, L.Ratio=10.822

#Combination of variance functions: 
# Structure: Different standard deviations per stratum
# Formula: ~1 | LakeName 
# Parameter estimates:
# LeConte     Spur 
#1.000000 4.572348 
# Structure: Different standard deviations per stratum
# Formula: ~1 | SampleNumber 
# Parameter estimates:
#       S2        S3        S4 
#1.0000000 0.7357974 0.6165303 
#
#Coefficients:
#                     Value  Std.Error   t-value p-value
#(Intercept)    -0.08135782 0.06302189 -1.290945  0.1999
#TadDens        -0.00630614 0.00302051 -2.087770  0.0395
#MFDens         -0.00089943 0.00022423 -4.011140  0.0001
#Silt           -0.00138773 0.00071189 -1.949348  0.0542
#SampleNumberS3  0.29735400 0.06371543  4.666907  0.0000
#SampleNumberS4  0.23362048 0.06026127  3.876793  0.0002
#
# Correlation: 
#               (Intr) TadDns MFDens Silt   SmpNS3
#TadDens        -0.427                            
#MFDens         -0.347 -0.016                     
#Silt           -0.357  0.380  0.094              
#SampleNumberS3 -0.656  0.001  0.000  0.000       
#SampleNumberS4 -0.693  0.000  0.000  0.000  0.685
#
#Standardized residuals:
#        Min          Q1         Med          Q3         Max 
#-2.54653303 -0.54805284 -0.01252839  0.51667328  3.53840506 
#
#Residual standard error: 0.216497 
#Degrees of freedom: 101 total; 95 residual

#non parametric methods?  
#laplace distribution?


##############
#models using estimated total biomass of consumers

#redefine PlotRes for biomass:
PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	plot(resid(x)~fitted(x), data=y)
	plot(resid(x)~y[,1], data=y, xlab=colnames(y[1]))
	plot(resid(x)~y[,3], data=y, xlab=colnames(y[3]))
	plot(resid(x)~y[,4], data=y, xlab=colnames(y[4]))
	plot(resid(x)~y[,11], data=y, xlab=colnames(y[11]))
	plot(resid(x)~y[,12], data=y, xlab=colnames(y[12]))
	plot(resid(x)~y[,13], data=y, xlab=colnames(y[13]))
	plot(resid(x)~y[,14], data=y, xlab=colnames(y[14]))
	plot(resid(x)~y[,15], data=y, xlab=colnames(y[15]))
	plot(resid(x)~y[,16], data=y, xlab=colnames(y[16]))
	plot(resid(x)~y[,17], data=y, xlab=colnames(y[17]))
	plot(resid(x)~y[,18], data=y, xlab=colnames(y[18]))
}
#PlotRes(model, dataset, plot dimensions)
###############

#LeConte Mayflies
q<-ggplot(whole[LeConte,], aes(x=MFBiomassEnd, y=ExperimentalControlDifference_Abundance))
	q+geom_point(aes(color=SampleNumber, position="dodge"), shape=16, size=5)+
	facet_grid(.~SampleNumber)
#Spur mayflies
q<-ggplot(whole[Spur,], aes(x=MFBiomassEnd, y=ExperimentalControlDifference_Abundance))
	q+geom_point(aes(color=SampleNumber, position="dodge"), shape=16, size=5)+
	facet_grid(.~SampleNumber)

q<-ggplot(whole[Spur,], aes(x=CageNumber, y=ExperimentalControlDifference_Abundance))
	q+geom_point(aes(color=TadpoleDensity, position="dodge"), shape=16, size=5)+
	facet_grid(.~SampleNumber)

################
#Tadpoles
TadAFDM<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Data Files/2009 Enclosures/TadpoleEndweightAFDM.csv", header=T)
SumsTadAFDM<-sqldf("select LakeID, SampleNumber, CageNumber, sum(EstimatedTadpoleAFDM) from TadAFDM group by LakeID, SampleNumber, CageNumber")

ggplot(TadAFDM, aes(x=as.factor(TadpoleDensity), y=EstimatedTadpoleAFDM))+
	geom_boxplot(aes(fill=SampleNumber), color="black")+
	facet_grid(.~LakeID)
ggplot(TadAFDM, aes(x=as.factor(TadpoleDensity), y=EstimatedTadpoleAFDM))+
	geom_boxplot(color="black")+
	facet_grid(.~LakeID)
ggplot(TadAFDM, aes(x=as.factor(TadpoleDensity), y=EstimatedTadpoleAFDM))+
	geom_point(aes(color=SampleNumber, position="dodge"))+
	facet_grid(.~LakeID)
ggplot(TadAFDM, aes(x=as.factor(MayflyDensity), y=EstimatedTadpoleAFDM))+
	geom_boxplot(aes(fill=SampleNumber), color="black")+
	facet_grid(.~LakeID)
ggplot(TadAFDM, aes(x=as.factor(MayflyDensity), y=EstimatedTadpoleAFDM))+
	geom_boxplot(color="black")+
	facet_grid(.~LakeID)


#
whole<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Data Files/2009 Enclosures/ConsumerResourceAnalysis_for_analysis_withTadpoleAFDM.csv", header=T)
str(whole)

lm1<-lm(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM+MFBiomassEnd+LakeName+SamplePeriod+Silt+Radiation+SampleNumber+CageNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,4))
AIC(lm1)	#415

#low-end outliers: 
which(resid(lm1)<=-5)	#72
#remove value 72, which is essentially zero and may be an error in measurement
whole<-whole[-72,]
ExperimentalControlDifference_Abundance<-ExperimentalControlDifference_Abundance[-72]
spur<-which(whole$LakeName=="Spur")
#high-end outliers:
which(resid(lm5)>=0.005)
AIC(lm1)	#415

lm1int<-lm(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+LakeName+SamplePeriod+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1int)
PlotRes(lm1int, whole, c(3,5))
AIC(lm1int)	#361; interaction is better model

lm2<-gls(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+LakeName+SamplePeriod+Silt+Radiation+SampleNumber, 
	method="ML",
	data=whole)
summary(lm2)
PlotRes(lm2, whole, c(3,5))
AIC(lm2)	#361.3

lm3<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+SamplePeriod+Silt+Radiation+SampleNumber, 
	random=~1|LakeName,
	method="ML",
	data=whole)
summary(lm3)
PlotRes(lm3, whole, c(3,5))
AIC(lm3)	#361.3

lm4<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+SamplePeriod+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	method="ML",
	data=whole)
summary(lm4)
PlotRes(lm4, whole, c(3,5))
AIC(lm4)	#370

LakeVar<-varIdent(form=~1|LakeName)
SampleVar<-varIdent(form=~1|SampleNumber)
lm6<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+SamplePeriod+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=LakeVar, 
	method="ML",
	data=whole)
summary(lm6)
PlotRes(lm6, whole, c(3,5))
AIC(lm6)	#233.84

lm7<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+SamplePeriod+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=SampleVar, 
	method="ML",
	data=whole)
summary(lm7)
PlotRes(lm7, whole, c(3,5))
AIC(lm7)	#371, not helpful, combo var?

lm8<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM*MFBiomassEnd+SamplePeriod+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar), 
	method="ML",
	data=whole)
summary(lm8)
PlotRes(lm8, whole, c(3,5))
AIC(lm8)	#226.091

#fixed variances, decreasing with increasing consumer biomass, did not function; 

#step down, drop interaction:
lm9<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM+MFBiomassEnd+SamplePeriod+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar), 
	method="ML",
	data=whole)
summary(lm9)
PlotRes(lm9, whole, c(3,5))
AIC(lm9)	#224, drop radiation

lm10<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM+MFBiomassEnd+SamplePeriod+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar), 
	method="ML",
	data=whole)
summary(lm10)
PlotRes(lm10, whole, c(3,5))
AIC(lm10)	#222.23, drop silt

lm11<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM+MFBiomassEnd+SamplePeriod, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar), 
	method="ML",
	data=whole)
summary(lm11)
PlotRes(lm11, whole, c(3,5))
AIC(lm11)	#222.9, drop tadpole

lm12<-lme(ExperimentalControlDifference_Abundance~
	MFBiomassEnd+SamplePeriod, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar), 
	method="ML",
	data=whole)
summary(lm12)
PlotRes(lm12, whole, c(3,5))
AIC(lm12)	#222.6, drop tadpole
anova(lm11, lm12)#leave tadpole in, refit lm11


lm11<-lme(ExperimentalControlDifference_Abundance~
	SumOfEstimatedTadpoleAFDM+MFBiomassEnd+SamplePeriod, 
	random=~1|LakeName/SampleNumber,
	weights=varComb(LakeVar, SampleVar), 
	method="REML",
	data=whole)
summary(lm11)
PlotRes(lm11, whole, c(3,5))
AIC(lm11)	#255
anova(lm11)


######################################
#rerun as algal growth rate, to eliminate the sample period factor
#but, do run with tadpole and mayfly "growth" rates?
#are they reall growth rates? NO, I don't think so, esp. because they are estimated
#but could test growth rate against densities

##############
#models using Density of consumers and algal growth rate

#redefine PlotRes for biomass:
PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	plot(resid(x)~fitted(x), data=y)
	plot(resid(x)~y[,1], data=y, xlab=colnames(y[1]))
	plot(resid(x)~y[,3], data=y, xlab=colnames(y[3]))
	plot(resid(x)~y[,4], data=y, xlab=colnames(y[4]))
	boxplot(resid(x)~y[,9], data=y, xlab=colnames(y[9]))
	plot(resid(x)~y[,11], data=y, xlab=colnames(y[11]))
	boxplot(resid(x)~y[,12], data=y, xlab=colnames(y[12]))
	plot(resid(x)~y[,13], data=y, xlab=colnames(y[13]))
	plot(resid(x)~y[,14], data=y, xlab=colnames(y[14]))
	plot(resid(x)~y[,15], data=y, xlab=colnames(y[15]))
	plot(resid(x)~y[,16], data=y, xlab=colnames(y[16]))
	plot(resid(x)~y[,17], data=y, xlab=colnames(y[17]))
	plot(resid(x)~y[,18], data=y, xlab=colnames(y[18]))
}
#PlotRes(model, dataset, plot dimensions)
###############

lm1<-lm(ExperimentalControlDifference_Growth~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,5))
AIC(lm1)	#same pattern of outliers as in other models:
which(resid(lm1)<=-0.2|resid(lm1)>=0.2)#72, 73, 80, 88
whole.not<-whole[-which(resid(lm1)<=-0.2),]#|resid(lm1)>=0.2),]#72, 73, 80, 88
ExperimentalControlDifference_Growth.whole.not<-ExperimentalControlDifference_Growth[-which(resid(lm1)<=-0.2)]#|resid(lm1)>=0.2)]

lm1nooutlier<-lm(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	data=whole.not)
summary(lm1nooutlier)
PlotRes(lm1nooutlier, whole.not, c(3,5))
AIC(lm1nooutlier)	#-229

lm2<-gls(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	method="ML",
	data=whole.not)
summary(lm2)
PlotRes(lm2, whole.not, c(3,5))
AIC(lm2)	#-230.5

lm3<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation+SampleNumber, 
	random=~1|LakeName,
	method="ML",
	data=whole.not)
summary(lm3)
PlotRes(lm3, whole.not, c(3,5))
AIC(lm3)	#-230.5

lm4<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	method="ML",
	data=whole.not)
summary(lm4)
PlotRes(lm4, whole.not, c(3,5))
AIC(lm4)	#-225.12

#variance sturctures?
boxplot(resid(lm4)~whole.not$TadDens)
boxplot(resid(lm4)~whole.not$MFDens)

LakeVar<-varIdent(form=~1|LakeName)
SampleVar<-varIdent(form=~1|SampleNumber)
MFVar<-varIdent(form=~1|MFDens)
TadVar<-varIdent(form=~1|TadDens)

lm5<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=LakeVar,
	method="ML",
	data=whole.not)
summary(lm5)
PlotRes(lm5, whole.not, c(3,5))
AIC(lm5)	#-340.9; big improvement

lm6<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=SampleVar,
	method="ML",
	data=whole.not)
summary(lm6)
PlotRes(lm6, whole.not, c(3,5))
AIC(lm6)	#-229; small improvement 

lm7<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=MFVar,
	method="ML",
	data=whole.not)
summary(lm7)
PlotRes(lm7, whole.not, c(3,5))
AIC(lm7)	#-259; improvemt

lm8<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=TadVar,
	method="ML",
	data=whole.not)
summary(lm8)
PlotRes(lm8, whole.not, c(3,5))
AIC(lm8)	#-258 improvement

ConsumerVar<-varComb(TadVar, MFVar)
lm9<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=ConsumerVar,
	method="ML",
	data=whole.not)
summary(lm9)
PlotRes(lm9, whole.not, c(3,5))
AIC(lm9)	#-302, combined is improvment

LakeConsumerVar<-varComb(TadVar, MFVar, LakeVar)
lm10<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=LakeConsumerVar,
	method="ML",
	data=whole.not)
summary(lm10)
PlotRes(lm10, whole.not, c(3,5))
AIC(lm10)	#-365

SampleLakeConsumerVar<-varComb(TadVar, MFVar, LakeVar, SampleVar)
lm11<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt+Radiation, 
	random=~1|LakeName/SampleNumber,
	weights=SampleLakeConsumerVar,
	method="ML",
	control=lmeControl(maxIter = 1000, msMaxIter=100),
	data=whole.not)
summary(lm11)
PlotRes(lm11, whole.not, c(3,5))
AIC(lm11)	#-369.07

#step down and summarize model 11; drop Radiation
lm12<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens*MFDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=SampleLakeConsumerVar,
	method="ML",
	control=lmeControl(maxIter = 1000, msMaxIter=100),
	data=whole.not)
summary(lm12)
PlotRes(lm12, whole.not, c(3,5))
AIC(lm12)	#-370; drop interaction

lm13<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens+MFDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=SampleLakeConsumerVar,
	method="ML",
	control=lmeControl(maxIter = 1000, msMaxIter=100),
	data=whole.not)
summary(lm13)
PlotRes(lm13, whole.not, c(3,5))
AIC(lm13)	#-371.14; now drop silt

lm14<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens+MFDens, 
	random=~1|LakeName/SampleNumber,
	weights=SampleLakeConsumerVar,
	method="ML",
	control=lmeControl(maxIter = 1000, msMaxIter=100),
	data=whole.not)
summary(lm14)
PlotRes(lm14, whole.not, c(3,5))
AIC(lm14)	#-371.64
anova(lm13, lm14)	# 14 is not an improvement, leave silt in, refit lm13

lm13<-lme(ExperimentalControlDifference_Growth.whole.not~
	TadDens+MFDens+Silt, 
	random=~1|LakeName/SampleNumber,
	weights=SampleLakeConsumerVar,
	method="REML",
	control=lmeControl(maxIter = 1000, msMaxIter=100),
	data=whole.not)
summary(lm13)
PlotRes(lm13, whole.not, c(3,5))
AIC(lm13)	#-308

#summarize and report in table in text!


###
BUT:!
Lake as fixed?
Is including the Block Duration not necessary or redundant IF you are including block as a random effect?
It seems redundant to me.

#refer to zuur about random effects interpretation:  
	# the variance of the normally distributed intercepts (a^2); 
	# the variance for the normally distributed residuals (b^2); 
	# and the correlation among observations within a random effect class as (a^2) / (a^2 + b^2)
	# but what is the log-ratio test that shows significance of this correlation
	# the log ratio test of the model without the random effect and the one with?  







