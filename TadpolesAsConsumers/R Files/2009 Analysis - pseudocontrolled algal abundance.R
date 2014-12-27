#analysis of the pseudo-controlled algal abundance in enclosures
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
whole_original<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\ConsumerResourceAnalysis_for_analysis_withTadpoleAFDM.csv", header=T)
str(whole_original)
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
	plot(resid(x)~y[,4], data=y, xlab=colnames(y[4]))
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
ExperimentalAFDMpM<-whole_original$ExperimentalAFDM/ExperimentalTileArea

#2. Control AFDM per m^2 
	#each tile is 2.4 * 2.4 cm = 0.024 * 0.024 m
	#there are 12 Control tiles
ControlTileArea<-(0.024^2)*12		#0.006912
ControlAFDMpM<-whole_original$ControlAFDM/ControlTileArea

#3. Differences between Experimental and Control Abundance
ExperimentalControlDifference_Abundance<-ExperimentalAFDMpM-ControlAFDMpM
#so, a negative value means less algae in the experiment

#4. Experimental AFDM per m^2 per day
ExperimentalAFDMpMpD<-whole_original$ExperimentalAFDM/ExperimentalTileArea/whole_original$SamplePeriod

#5. Control AFDM per m^2 per day
ControlAFDMpMpD<-whole_original$ControlAFDM/ControlTileArea/whole_original$SamplePeriod

#6. Differences between Experimental and Control Growth Rate
ExperimentalControlDifference_Growth<-ExperimentalAFDMpMpD-ControlAFDMpMpD

#7. Change in Tadpole Wet Weight
TadpoleWeightChange<-whole_original$TadWtStart-whole_original$TadWtEnd

#8. Change in MF Density (Density-Count)
MFDensChange<-whole_original$MFDens-whole_original$MFCount

#9. Change in MF Biomass
MFBiomassChange<-whole_original$MFBiomassStart-whole_original$MFBiomassEnd

#############
#Define groups
#by lake:
leconte<-which(whole_original$LakeName=="LeConte")
spur<-which(whole_original$LakeName=="Spur")

#############
#basic plotting
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
#between lake comparison of algal abundance
mean(ExperimentalAFDMpM[leconte])#
sd(ExperimentalAFDMpM[leconte])#
mean(ExperimentalAFDMpM[spur])#
sd(ExperimentalAFDMpM[spur])#
mean(mean(ExperimentalAFDMpM[leconte]))/(mean(ExperimentalAFDMpM[spur]))#

anova(aov(ExperimentalAFDMpM~whole_original$LakeName))#
wholesummarized<-summarySE(whole_original, measurevar="ExperimentalAFDM", groupvars=c("LakeName","SampleNumber"))
ggplot(wholesummarized, aes(y=ExperimentalAFDM, x=LakeName, color=SampleNumber))+
	geom_errorbar(aes(ymin=ExperimentalAFDM-se, ymax=ExperimentalAFDM+se), width=.1, position=position_dodge(0.1))
#requires running/calling summarySE from another file.

############
#Set a working data set from teh original:
whole<-whole_original
############
#models of Algal Abundance:

lm1<-lm(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation+SampleNumber, 
	data=whole)
summary(lm1)
PlotRes(lm1, whole, c(3,4))
shapiro.test(resid(lm1))
AIC(lm1)	#362.03

#outliers?
quantile(whole$ExperimentalAFDM, p=c(0.05, 0.95))
quantile(resid(lm1), p=c(0.01, 0.05, 0.95, 0.99))
#      1%        5%       95%       99% 
#-3.505851 -1.726993  1.153512  6.304898 
which(resid(lm1)<=-3.5)	#54, 72
which(resid(lm1)>=6.3)	#80, 88
#despite these quantile based selections, I can see that value #73 is still outside normal, and closer to 80 and 88;
#meanwhile, value #54 sits quite close to a lot of other observations and far from #72
#I select outliers visually, as those that were most strongly influencing normality (QQplot), as:
#low: 72
#high: 73,80,88

#remove low end, "systematic" outliers (where controlled algal abundance was extremely negative b/c of methods)
#run model with this removed a priori, because there is a good reason to exclude it
whole<-whole[-72,]
ExperimentalControlDifference_Abundance<-ExperimentalControlDifference_Abundance[-72]

#rerun models with these removed after wards, because their removal is exploratory to see their influence, and compare to initial model
#remove high end "statistical" outliers
#whole2<-whole[-c(73,80,88)]
#ExperimentalControlDifference_Abundance2<-ExperimentalControlDifference_Abundance[-c(73,80,88)]

#run with low outliers removed:
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

lmc<-lmeControl(niterEM = 5000, msMaxIter = 1000)
lm2<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod, 
	random=~1|SampleNumber,
	data=whole,
	method="ML", 
	control=lmc
	)
summary(lm2)
PlotRes(lm2, whole, c(3,4))
AIC(lm2)	#365.7


lm4<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=c(~1+SamplePeriod|SampleNumber),
	data=whole,
	method="ML",
	control=lmc
	)
summary(lm4)
PlotRes(lm4, whole, c(3,4))
AIC(lm4)	#370.26
anova(lm2, lm4)# so the random slope for days within block does not help.
ggplot(whole, aes(y=resid(lm2), x=SamplePeriod))+
	geom_boxplot(aes(color=LakeName))+
	facet_grid(.~SampleNumber)
#and there is no real pattern to support that it would work...

LakeVar<-varIdent(form=~1|LakeName)
SampleVar<-varIdent(form=~1|SampleNumber)
lm5<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation+SamplePeriod+LakeName, 
	random=~1|SampleNumber,
	weights=LakeVar,
	data=whole,
	method="ML"
	)
summary(lm5)
PlotRes(lm5, whole, c(3,4))
AIC(lm5)	#237.95
anova(lm2, lm5)
0.5 * (1 - pchisq(131.06, 1))

lm6<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation+SamplePeriod+LakeName, 
	random=~1|SampleNumber,
	weights=SampleVar,
	data=whole,
	method="ML"
	)
summary(lm6)
PlotRes(lm6, whole, c(3,4))
AIC(lm6)	#366
anova(lm6, lm2)#no difference, not an improvement

lm7<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+Radiation+SamplePeriod+LakeName, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML"
	)
summary(lm7)
PlotRes(lm7, whole, c(3,4))
AIC(lm7)	#
anova(lm2, lm5, lm7)	#this is an improvement

#step down lm7, drop radiation first
lm8<-lme(ExperimentalControlDifference_Abundance~
	TadDens*MFDens+Silt+SamplePeriod+LakeName, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML",
	control=lmc
	)
summary(lm8)
PlotRes(lm8, whole, c(3,4))
AIC(lm8)	#227
anova(lm7, lm8)	#slight improvement, drop interaction

lm9<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens+Silt+SamplePeriod+LakeName, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML",
	control=lmc
	)
summary(lm9)
PlotRes(lm9, whole, c(3,4))
AIC(lm9)	#
anova(lm8, lm9)	#225; drop lake name

lm10<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens+Silt+SamplePeriod, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML",
	control=lmc
	)
summary(lm10)
PlotRes(lm10, whole, c(3,4))
AIC(lm10)	#
anova(lm9, lm10)	#225, no different, drop silt

lm11<-lme(ExperimentalControlDifference_Abundance~
	TadDens+MFDens+SamplePeriod, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML",
	control=lmc
	)
summary(lm11)
PlotRes(lm11, whole, c(3,4))
AIC(lm11)	#
anova(lm10, lm11)	#225.9, no different but slightly higher, try drop tadpoles

lm12<-lme(ExperimentalControlDifference_Abundance~
	MFDens+SamplePeriod, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="ML",
	control=lmc
	)
summary(lm12)
PlotRes(lm12, whole, c(3,4))
AIC(lm12)	#
anova(lm11, lm12)	#225.1

#refit with REML:
lm12.reml<-lme(ExperimentalControlDifference_Abundance~
	MFDens+SamplePeriod, 
	random=~1|SampleNumber,
	weights=varComb(SampleVar,LakeVar),
	data=whole,
	method="REML",
	control=lmc
	)
summary(lm12.reml)
PlotRes(lm12.reml, whole, c(3,4))
intervals(lm12.reml)

#Linear mixed-effects model fit by REML
# Data: whole 
#       AIC      BIC    logLik
#  249.2583 269.9381 -116.6292
#Random effects:
# Formula: ~1 | SampleNumber
#       (Intercept)  Residual
#StdDev:   0.4836506 0.3257609
#
#Combination of variance functions: 
# Structure: Different standard deviations per stratum
# Formula: ~1 | SampleNumber 
# Parameter estimates:
#       S2        S3        S4 
#1.0000000 0.5765392 0.4645842 
# Structure: Different standard deviations per stratum
# Formula: ~1 | LakeName 
# Parameter estimates:
#LeConte    Spur 
#  1.000  10.175 
#Fixed effects: ExperimentalControlDifference_Abundance ~ MFDens + SamplePeriod 
#                  Value Std.Error DF   t-value p-value
#(Intercept)   1.4398106 1.1768280 96  1.223467  0.2241
#MFDens       -0.0010252 0.0002720 96 -3.769186  0.0003
#SamplePeriod -0.0735639 0.0590301 96 -1.246210  0.2157
# Correlation: 
#             (Intr) MFDens
#MFDens       -0.020       
#SamplePeriod -0.971 -0.001
#
#Standardized Within-Group Residuals:
#        Min          Q1         Med          Q3         Max 
#-2.87957421 -0.38191853 -0.01766917  0.42277204  4.57597301 
#
#Number of Observations: 101
#Number of Groups: 3 

#Approximate 95% confidence intervals
#
# Fixed effects:
#                    lower         est.         upper
#(Intercept)  -0.896174414  1.439810611  3.7757956365
#MFDens       -0.001565081 -0.001025184 -0.0004852863
#SamplePeriod -0.190737781 -0.073563917  0.0436099460
#attr(,"label")
#[1] "Fixed effects:"
#
# Random Effects:
#  Level: SampleNumber 
#                     lower      est.    upper
#sd((Intercept)) 0.09283645 0.4836506 2.519677
#
# Variance function:
#           lower       est.      upper
#A.S3   0.3871702  0.5765392  0.8585308
#A.S4   0.3041154  0.4645842  0.7097257
#B.Spur 7.1309480 10.1750016 14.5184986
#attr(,"label")
#[1] "Variance function:"
#
# Within-group standard error:
#    lower      est.     upper 
#0.2551285 0.3257609 0.4159479 



#############
#Rerun with out high-end outliers
#rerun models with these removed, because their removal is exploratory to see their influence, and compare to initial model
#remove high end "statistical" outliers
whole2<-whole[-c(73,80,88),]
ExperimentalControlDifference_Abundance2<-ExperimentalControlDifference_Abundance[-c(73,80,88)]

lm1<-lm(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod+SampleNumber, 
	data=whole2)
summary(lm1)
PlotRes(lm1, whole2, c(3,4))
AIC(lm1)	#362.03

#...now have more outliers!
which(resid(lm1)>=2)# values in slots 57,72,79,87
whole2<-whole2[-c(57,72,79,87),]
ExperimentalControlDifference_Abundance2<-ExperimentalControlDifference_Abundance2[-c(57,72,79,87)]

lm1<-lm(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod+SampleNumber, 
	data=whole2)
summary(lm1)
PlotRes(lm1, whole2, c(3,4))
AIC(lm1)	#
#...now have more outliers!
which(resid(lm1)>=4)# values in slots 76, 82
whole2<-whole2[-c(76,82),]
ExperimentalControlDifference_Abundance2<-ExperimentalControlDifference_Abundance2[-c(76,82)]

lm1<-lm(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod+SampleNumber, 
	data=whole2)
summary(lm1)
PlotRes(lm1, whole2, c(3,4))
AIC(lm1)	#
#...now have more outliers!
which(resid(lm1)<=-1.5)# values in slots 54 55 58 81
whole2<-whole2[-c(54,55,58,81),]
ExperimentalControlDifference_Abundance2<-ExperimentalControlDifference_Abundance2[-c(54,55,58,81)]

lm1<-lm(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod+SampleNumber, 
	data=whole2)
summary(lm1)
PlotRes(lm1, whole2, c(3,4))
AIC(lm1)	#133
#OK, that is about as normal and symmetrical as anything... proceed with model selection!

lm2<-gls(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod+SampleNumber, 
	method="ML",
	data=whole2)
summary(lm2)
PlotRes(lm2, whole2, c(3,4))
AIC(lm2)	#131.8

lmc<-lmeControl(niterEM = 5000, msMaxIter = 1000)
lm3<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation+SamplePeriod, 
	random=~1|SampleNumber,
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm3)
PlotRes(lm3, whole2, c(3,4))
AIC(lm3)	#140.42, no change in residuals

lm4<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1+SamplePeriod|SampleNumber,
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm4)
PlotRes(lm4, whole2, c(3,4))
AIC(lm4)	#140.5, maybe barely better. certainly not worse

LakeVar<-varIdent(form=~1|LakeName)
SampleVar<-varIdent(form=~1|SampleNumber)
lm5<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1+SamplePeriod|SampleNumber,
	weights=LakeVar,
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm5)
PlotRes(lm5, whole2, c(3,4))
AIC(lm5)	#110.45 but reduces normality of residuals

lm6<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1+SamplePeriod|SampleNumber,
	weights=SampleVar,
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm6)
PlotRes(lm6, whole2, c(3,4))
AIC(lm6)	#142

lm7<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt+Radiation, 
	random=~1+SamplePeriod|SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm7)
PlotRes(lm7, whole2, c(3,4))
AIC(lm7)	#98
#ultimately, removing outliers does not really change the pattern of the residuals, you still have long tails no matter how many values removed...

#step down lm7, drop radiation
lm8<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+LakeName+Silt, 
	random=~1+SamplePeriod|SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm8)
PlotRes(lm8, whole2, c(3,4))
AIC(lm8)	#96.5, drop lakename

lm9<-lme(ExperimentalControlDifference_Abundance2~
	TadDens*MFDens+Silt, 
	random=~1+SamplePeriod|SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm9)
PlotRes(lm9, whole2, c(3,4))
AIC(lm9)	#96.7, drop interaction

lm10<-lme(ExperimentalControlDifference_Abundance2~
	TadDens+MFDens+Silt, 
	random=~1+SamplePeriod|SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm10)
PlotRes(lm10, whole2, c(3,4))
AIC(lm10)	#95.5, drop silt

lm11<-lme(ExperimentalControlDifference_Abundance2~
	TadDens+MFDens, 
	random=~1+SamplePeriod|SampleNumber,
	weights=varComb(LakeVar, SampleVar),
	data=whole2,
	method="ML", 
	control=lmc
	)
summary(lm11)
PlotRes(lm11, whole2, c(3,4))
AIC(lm11)	#96.06
anova(lm10, lm11)#not different
#using the more restricted data set, with about 10 outliers removed, does not change the story, or ultimately remove the non-normal distribution of the residuals


##############
#plotting:
#heat map plotting of pseudo controlled algal abundance:
whole<-data.frame(whole, ExperimentalControlDifference_Abundance)
leconte<-which(whole$LakeName=="LeConte")
spur<-which(whole$LakeName=="Spur")
#redefine data:
whole_original<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\ConsumerResourceAnalysis_for_analysis_withTadpoleAFDM.csv", header=T)
wholeplot<-whole_original
levels(wholeplot$SampleNumber)[levels(wholeplot$SampleNumber)=="S2"]<-"Early August"
levels(wholeplot$SampleNumber)[levels(wholeplot$SampleNumber)=="S3"]<-"Late August"
levels(wholeplot$SampleNumber)[levels(wholeplot$SampleNumber)=="S4"]<-"Mid September"
#Recalculate Differences between Experimental and Control Abundance
ExperimentalControlDifference_Abundance<-ExperimentalAFDMpM-ControlAFDMpM
#here, log modulus transform used only for display purposes, with the original low outliers removed:
logmod<-sign(ExperimentalControlDifference_Abundance[-72])*log(abs(ExperimentalControlDifference_Abundance[-72])+1)

#Summarized into AllOne
HeatMap<-ggplot(wholeplot[-72,], aes(x=as.factor(MFDens), y=as.factor(TadDens) ))+
	geom_tile(aes(fill=logmod))+
	scale_fill_gradient2(midpoint=0, low=muted("blue"), mid="white", high=muted("green"), limits=c(-2,2.1), guide="legend")+
	xlab("Mayfly Abundance")+ylab("Tadpole Abundance")+guides(fill=guide_legend(title=NULL))+
	theme(legend.position="right", 
		panel.background = element_rect(fill="white"),
		panel.border=element_rect(fill=NA, color="light grey", size=15),
		axis.text.x=element_text(color="black", face="bold", size=18), 
		axis.text.y=element_text(color="black", face="bold", size=18),
		axis.title.x=element_text(color="black", face="bold", size=20),
		axis.title.y=element_text(color="black", face="bold", size=20), 
		strip.background=element_rect(fill="light grey"),
		strip.text=element_text(face="bold", size=20), 
		legend.text=element_text(face="bold", size=18),
		axis.ticks=element_blank()
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_pseudocontrolledalgalabundance_AllOne_heatmap.png", 
	res=300, width=9, height=6, units="in", restoreConsole=F)
HeatMap
dev.off()

#lakes as facets:
HeatMap<-ggplot(wholeplot[-72,], aes(x=as.factor(MFDens), y=as.factor(TadDens) ))+
	geom_tile(aes(fill=logmod))+
#	scale_fill_gradient2(midpoint=-0.06, low="black", mid="grey", high="white", limits=c(-2,2.1), guide="legend")+
	scale_fill_gradient2(midpoint=0, low=muted("blue"), mid="white", high=muted("green"), limits=c(-2,2.1), guide="legend")+
	facet_grid(.~LakeName)+
	xlab("Mayfly Abundance")+ylab("Tadpole Abundance")+guides(fill=guide_legend(title=NULL))+
	theme(legend.position="right", 
		panel.background = element_rect(fill="white"),
		panel.border=element_rect(fill=NA, color="light grey", size=15),
		axis.text.x=element_text(color="black", face="bold", size=18), 
		axis.text.y=element_text(color="black", face="bold", size=18),
		axis.title.x=element_text(color="black", face="bold", size=20),
		axis.title.y=element_text(color="black", face="bold", size=20), 
		strip.background=element_rect(fill="light grey"),
		strip.text=element_text(face="bold", size=20), 
		legend.text=element_text(face="bold", size=18),
		axis.ticks=element_blank()
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_pseudocontrolledalgalabundance_bylake_heatmap.png", 
	res=300, width=9, height=6, units="in", restoreConsole=F)
HeatMap
dev.off()

#lakes and blocks as facets:
HeatMap<-ggplot(wholeplot[-72,], aes(x=as.factor(MFDens), y=as.factor(TadDens) ))+
	geom_tile(aes(fill=logmod))+
#	scale_fill_gradient2(midpoint=-0.06, low="black", mid="grey", high="white", limits=c(-2,2.1), guide="legend")+
	scale_fill_gradient2(midpoint=0, low=muted("blue"), mid="white", high=muted("green"), limits=c(-2,2.1), guide="legend")+
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
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_Enclosures_pseudocontrolledalgalabundance_bylakebyblock_heatmap.png", 
	res=300, width=9, height=6, units="in", restoreConsole=F)
HeatMap
dev.off()







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

logmod<-sign(ExperimentalControlDifference_Abundance[-72])*log(abs(ExperimentalControlDifference_Abundance[-72])+1)
hist(logmod)
qqnorm(logmod)
hist(logmod[leconte])
qqnorm(logmod[leconte])
shapiro.test(logmod[leconte])#nearly normal - certainly an improvement
hist(logmod[spur])
qqnorm(logmod[spur])
shapiro.test(logmod[spur])#nearl normal

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







