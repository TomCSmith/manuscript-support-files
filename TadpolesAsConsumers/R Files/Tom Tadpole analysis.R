#analysis of tadpole endweights relative to consumer densities.
##############
detach()
rm(list=ls())
##############
#load required packages
library(stats)
library(reshape2)
library(reshape)
library(plyr)
library(ggplot2)
library(nlme)
require(gridExtra)
##############
#fit observed AFDM to observed stages for marmot tadpoles
Marmot<-read.csv("C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2009 Enclosures//MarmotLakeTadpole_AFDMlengthstages.csv", header=T)
str(Marmot)
plot(Marmot$AFDM~Marmot$GosnerStage)
nControl<-nls.control(maxiter = 100, tol = 1e-05, minFactor = 1/1024,
            printEval = FALSE, warnOnly = FALSE)
n1<-nls(AFDM~a*GosnerStage^b, start=list(a=1, b=1), data=Marmot, control=nControl)
summary(n1)
AFDMest2<-0.0001028*Marmot$GosnerStage^3.6974552

ggplot(Marmot, aes(GosnerStage, AFDM))+
	geom_point()+
	scale_y_log10()+scale_x_log10()+
	stat_smooth()


###############
#use estimated AFDM for Experimental tadpoles based on observed gosner stages and regression:
#Load data:
#endweights <- read.csv("\\\\BRIGGS-NAS1\\Thomas Smith\\Tadpole Data\\Endweights.csv", header = T)
endweights <- read.csv("C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2009 Enclosures//Endweights.csv", header = T)
# $ LakeID         : int  10102 10102 10102 10102 10102 10102 10102 10102 10102 10102 ...
# $ SampleDate     : Factor w/ 8 levels "1-Aug-09","12-Aug-09",..: 4 4 4 4 4 4 4 4 4 4 ...
# $ CageNumber     : int  1 1 1 1 1 1 1 1 1 1 ...
# $ AnimalNumber   : int  1 2 3 4 5 6 7 8 9 10 ...
# $ GosnerStage    : int  37 38 31 35 31 31 34 35 34 35 ...
# $ Predicted.AFDM : num  55.1 63.4 21.7 41.1 21.7 ...
# $ Tadpole.Density: int  20 20 20 20 20 20 20 20 20 20 ...
# $ MayflyDensity  : int  125 125 125 125 125 125 125 125 125 125 ...
str(endweights)

##############
#data formatting:
levels(endweights$SampleDate)
endweights$SampleNumber<-endweights$SampleDate
endweights$SampleNumber<-revalue(endweights$SampleNumber, c("15-Jul-09"="July","1-Aug-09"="early August","23-Aug-09"="late August","14-Sep-09"="September","21-Jul-09"="July","12-Aug-09"="early August",	"31-Aug-09"="late August",	"21-Sep-09"="September"))
endweights$LakeName<-as.factor(endweights$LakeID)
endweights$LakeName<-revalue(endweights$LakeName, c("10475"="Spur", "10102"="LeConte"))

#summarize tad wts
endweightsTadWt<-endweights[,c(9,10,3,6)]
endweightsmelt<-melt(endweightsTadWt, id=c("LakeName", "SampleNumber","CageNumber"))
TadWt<-cast(endweightsmelt, LakeName+SampleNumber+CageNumber~variable, mean)

#Summarize Gosner Stage
endweightsGos<-endweights[,c(9,10,3,5)]
Gosmelt<-melt(endweightsGos, id=c("LakeName", "SampleNumber","CageNumber"))
Gos<-cast(Gosmelt, LakeName+SampleNumber+CageNumber~variable, mean)[4]

#summarize Tadpole.Density
TadDens<-endweights[,c(9,10,3,7)]
TadDensmelt<-melt(TadDens, id=c("LakeName", "SampleNumber","CageNumber"))
TadpoleDensity<-cast(TadDensmelt, LakeName+SampleNumber+CageNumber~variable, mean)[4]

#Summarize MayflyDensity
mf<-endweights[,c(9,10,3,8)]
mfmelt<-melt(mf, id=c("LakeName", "SampleNumber","CageNumber"))
MFDens<-cast(mfmelt, LakeName+SampleNumber+CageNumber~variable, mean)[4]

#define df with recast summarized data:
TadWtDF<-data.frame(TadWt, Gos, TadpoleDensity, MFDens)
str(TadWtDF)
#'data.frame':   96 obs. of  7 variables:
# $ LakeName       : Factor w/ 2 levels "LeConte","Spur": 1 1 1 1 1 1 1 1 1 1 ...
# $ SampleNumber   : Factor w/ 4 levels "2","4","1","3": 1 1 1 1 1 1 1 1 1 1 ...
# $ CageNumber     : int  1 2 3 6 7 8 9 10 12 13 ...
# $ Predicted.AFDM : num  49.1 35.6 38.4 44.5 38.2 ...
# $ GosnerStage    : num  35.7 34 34 35.1 34.5 ...
# $ Tadpole.Density: num  20 2 2 20 2 10 20 20 10 2 ...
# $ MayflyDensity  : num  125 125 0 25 25 12

TadWtDF$SampleNumber<-factor(TadWtDF$SampleNumber, levels=c("July","early August","late August","September"))

####################
#analysis:

#load function for linear model diagnostic plots
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
	}
#PlotRes(model, dataset, plot dimensions)
##################
#linear models:
lm1<-lm(Predicted.AFDM~Tadpole.Density+MayflyDensity+SampleNumber+LakeName, data=TadWtDF)
PlotRes(lm1, TadWtDF, c(2,5))
min(TadWtDF$Predicted.AFDM)
which(TadWtDF$Predicted.AFDM<=22)
lm1<-lm(Predicted.AFDM~Tadpole.Density+MayflyDensity+SampleNumber+LakeName, data=TadWtDF[-55,])
PlotRes(lm1, TadWtDF[-55,], c(2,5))
#may need variance structures for SampleNumbers and for TadpoleDensity and MFDensity
#could use random intercepts for SampleNumbers
summary(lm1)
lm2<-lm(Predicted.AFDM~Tadpole.Density+MayflyDensity+LakeName, 
	data=TadWtDF)
summary(lm2)
PlotRes(lm2, TadWtDF, c(2,5))
anova(lm1, lm2)	#<0.0001
#so def. include samplenumber as random.
#include different slopes for different lakes - as an interaction

lm3<-gls(Predicted.AFDM~Tadpole.Density+MayflyDensity+LakeName+SampleNumber, 
	method="ML",
	data=TadWtDF)
lm4<-gls(Predicted.AFDM~Tadpole.Density*LakeName+MayflyDensity+SampleNumber, 
	method="ML",
	data=TadWtDF)
summary(lm4)
anova(lm3, lm4)# lm4 is a little better

#add random intercepts and random slopes for sample NUmber
lm5<-lme(Predicted.AFDM~Tadpole.Density*LakeName+MayflyDensity, 
	random=~1|SampleNumber,			
	method="ML",
	data=TadWtDF)
summary(lm5)
PlotRes(lm5, TadWtDF, c(2,5))
AIC(lm4, lm5)

lm6<-lme(Predicted.AFDM~Tadpole.Density*LakeName+MayflyDensity, 
	random=~1+SampleNumber|SampleNumber,			
	method="ML",
	data=TadWtDF)
summary(lm6)
PlotRes(lm6, TadWtDF, c(2,5))
anova(lm5, lm6)
#random slopes not necessary

#add variance structure:
#among sample dates:
SampleVar<-varIdent(form=~1|SampleNumber)
lm7<-lme(Predicted.AFDM~Tadpole.Density*LakeName+MayflyDensity, 
	random=~1|SampleNumber,			
	weights=SampleVar,
	method="ML",
	data=TadWtDF)
summary(lm7)
PlotRes(lm7, TadWtDF, c(2,5))
anova(lm5, lm7)	#that helps

#among tadpole groups
TadVar<-varFixed(~Tadpole.Density)
TadVar<-varIdent(form=~1|Tadpole.Density)
lm8<-lme(Predicted.AFDM~Tadpole.Density*LakeName+MayflyDensity, 
	random=~1|SampleNumber,			
	weights=TadVar,
	method="ML",
	data=TadWtDF)
summary(lm8)
PlotRes(lm8, TadWtDF, c(2,5))
anova(lm5, lm7, lm8)	#neither formulation helps, varFixed is much worse, no need for variance with tad density

#among mayfly
MFVar<-varFixed(~MayflyDensity)	#this returns an error...
MFVar<-varIdent(form=~1|MayflyDensity)
lm9<-lme(Predicted.AFDM~Tadpole.Density*LakeName+MayflyDensity, 
	random=~1|SampleNumber,			
	weights=MFVar,
	method="ML",
	data=TadWtDF)
summary(lm9)
PlotRes(lm9, TadWtDF, c(2,5))
anova(lm5, lm7, lm8, lm9)	#not an improvement

#step down lm7, try drop lake x tadpole interaction, then MF density
lm10<-lme(Predicted.AFDM~Tadpole.Density+LakeName+MayflyDensity, 
	random=~1|SampleNumber,			
	weights=SampleVar,
	method="ML",
	data=TadWtDF)
summary(lm10)
PlotRes(lm10, TadWtDF, c(2,5))
anova(lm7, lm10)	#worse model by AIC

lm11<-lme(Predicted.AFDM~Tadpole.Density*LakeName, 
	random=~1|SampleNumber,			
	weights=SampleVar,
	method="ML",
	data=TadWtDF)
summary(lm11)
PlotRes(lm1, TadWtDF, c(2,5))
anova(lm7, lm11)	#

#refit with REML
lm11<-lme(Predicted.AFDM~Tadpole.Density*LakeName, 
	random=~1|SampleNumber,			
	weights=SampleVar,
	method="REML",
	data=TadWtDF)
summary(lm11)


#interaction appears to be that there is a negative response in Spur, but a positive response in LeConte
#test for each lake, to help interpret the interaction:

lm11.L<-lme(Predicted.AFDM~Tadpole.Density, 
	random=~1|SampleNumber,			
	weights=SampleVar,
	method="REML",
	data=TadWtDF[TadWtDF$LakeName=="LeConte",])
summary(lm11.L)

lm11.S<-lme(Predicted.AFDM~Tadpole.Density, 
	random=~1|SampleNumber,			
	weights=SampleVar,
	method="REML",
	data=TadWtDF[TadWtDF$LakeName=="Spur",])
summary(lm11.S)

#in Spur, there is no effect of tadpoles (but maybe a negative trend...) - 
#where there is large variance and large differences in variance among blocks
#In LeConte, there is a POSITIVE effect of tadpole density on tadpole AFDM; 
#there is little difference in variance among blocks



#---------------------
#plotting:
#################
#basic plots
TadWtDFConsCat<-TadWtDF
attach(TadWtDFConsCat)
TadWtDFConsCat$MayflyDensity[MayflyDensity==25]<-"low"
TadWtDFConsCat$MayflyDensity[MayflyDensity==125]<-"medium"
TadWtDFConsCat$MayflyDensity[MayflyDensity==250]<-"high"
TadWtDFConsCat$MayflyDensity[MayflyDensity==0]<-"zero"
TadWtDFConsCat$Tadpole.Density[Tadpole.Density==2]<-"low"
TadWtDFConsCat$Tadpole.Density[Tadpole.Density==10]<-"medium"
TadWtDFConsCat$Tadpole.Density[Tadpole.Density==20]<-"high"
detach(TadWtDFConsCat)
str(TadWtDFConsCat)
TadWtDFMelt<-melt(TadWtDFConsCat[,c(1,2,3,4,6,7)], id=c("LakeName", "CageNumber", "SampleNumber", "Predicted.AFDM"))
TadWtDFMelt<-rename(TadWtDFMelt, c(variable="consumer"))
attach(TadWtDFMelt)
TadWtDFMelt$consumertype[consumer=="Tadpole.Density"]<-"Tadpole"
TadWtDFMelt$consumertype[consumer=="MayflyDensity"]<-"Mayfly"
detach(TadWtDFMelt)
str(TadWtDFMelt)
plot3<-ggplot(TadWtDFMelt, aes(x=value, y=Predicted.AFDM))+
	geom_boxplot(size=1.5, color="black")+
	facet_grid(LakeName~consumertype, scales="free")+
#	scale_fill_manual(values=c("white", "light grey", "grey", "dark grey"))+
	xlab("\nConsumer Density")+ylab("Tadpole AFDM (mg)\n")+
	scale_x_discrete(limits=c("zero","low","medium","high"))+
	theme(panel.background=element_blank(),
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank(),
		panel.border=element_rect(fill=NA),
		strip.background=element_rect(color="black", fill=NA),
		strip.text=element_text(color="black", size=14, face="bold"),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=10),	
		axis.title=element_text(color="black", size=14, face="bold")
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_TadpoleAFDM_boxplots_byconsumer_lake.png", 
	res=300, width=5, height=4, units="in", restoreConsole=F)
plot3
dev.off()






legend.position=c(0,0.5),
		legend.key=element_rect(fill=NA, color=NA),
		legend.key.size=unit(1.25, "cm"),
		legend.justification=c(0,0),
		legend.title=element_blank(),
		legend.background=element_rect(fill=NA, color="black"),
		legend.text=element_text(size=14, face="bold"),








#-------------------------------------------
#maintenance of treatments
str(endweights)
over39<-function(x){
	length(which(x>=39))}
a<-aggregate(data=endweights[,c(1,2,3,4,5)], GosnerStage~CageNumber+SampleDate+LakeID, over39)
b<-aggregate(data=endweights[,c(1,2,3,4,5)], AnimalNumber~CageNumber+SampleDate+LakeID, max)
c<-aggregate(data=endweights[,c(1,2,3,4,7)], Tadpole.Density~CageNumber+SampleDate+LakeID, max)

d<-cbind(a, c$Tadpole.Density)
#write.table(d, "C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2009 Enclosures//TadpoleDensityDynamics.csv", sep="\t") 
d<-read.csv("C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2009 Enclosures//TadpoleDensityDynamics.csv", header=T)
str(d)
e<-aggregate(data=endweights[,c(1,2,3,4,5)], GosnerStage~SampleDate+LakeID, mean)


#recode variables:
d$Block[d$SampleDate=="15-Jul-09" | d$SampleDate=="21-Jul-09"]<-0
d$Block[d$SampleDate=="1-Aug-09" | d$SampleDate=="12-Aug-09"]<-1
d$Block[d$SampleDate=="23-Aug-09" | d$SampleDate=="31-Aug-09"]<-2
d$Block[d$SampleDate=="14-Sep-09" | d$SampleDate=="21-Sep-09"]<-3
#calculate percent change:
d$PpnChange<-d$GosnerStage/d$TadpoleTx
#write.csv(d, "C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2009 Enclosures//TadpoleDensityDynamics.csv", sep=",") 
d$BlockName[d$Block=="1"]<-"Jul-Aug"
d$BlockName[d$Block=="2"]<-"Mid Aug"
d$BlockName[d$Block=="3"]<-"Aug-Sep"
d$LakeName[d$LakeID=="10102"]<-"LeConte"
d$LakeName[d$LakeID=="10475"]<-"Spur"
d$BlockName<-factor(d$BlockName, levels=c("Jul-Aug","Mid Aug","Aug-Sep"))
d$Days[d$Block=="1" & d$LakeID=="10102"]<-16
d$Days[d$Block=="2" & d$LakeID=="10102"]<-23
d$Days[d$Block=="3" & d$LakeID=="10102"]<-19
d$Days[d$Block=="1" & d$LakeID=="10475"]<-23
d$Days[d$Block=="2" & d$LakeID=="10475"]<-18
d$Days[d$Block=="3" & d$LakeID=="10475"]<-18
str(d)
notZero<-which(d$Block>0)
Block1<-which(d$Block==1)
q<-which(d$Block==3 & d$LakeName=="Spur")
v<-aov((PpnChange/Days)~as.factor(TadpoleTx), data=d[q,])
summary(v)
drop1(v,~.,test="F")
TukeyHSD(v)

png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_TadpoleDensityDynamics_byBlockLake.png", 
	res=300, width=8, height=6, units="in", restoreConsole=F)
ggplot(d[notZero,], aes(x=as.factor(TadpoleTx), y=(((PpnChange+0.0001)*100)/Days)))+
	geom_boxplot(aes(fill=as.factor(TadpoleTx)), size=2)+
	facet_grid(LakeName~BlockName)+
	scale_fill_manual(values=c("white", "white", "white"))+
	guides(fill=FALSE)+
	scale_y_continuous(limits=c(0,3))+#breaks=c(10,20,30,40,50), 
	xlab("\nTadpole Treatment")+ylab("% tadpoles > stage 38\n")+
	theme(panel.background=element_blank(),
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank(),
		panel.border=element_rect(fill=NA),
		strip.background=element_rect(color="black", fill=NA),
		strip.text=element_text(color="black", size=20, face="bold"),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=18, face="bold"),	
		axis.title=element_text(color="black", size=20, face="bold")
		)
dev.off()
