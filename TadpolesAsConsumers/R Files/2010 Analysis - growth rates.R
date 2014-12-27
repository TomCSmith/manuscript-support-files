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
library(RColorBrewer)
#############
#load and organize data:
whole<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Summer 2010/2010MesocosmCSVdata/ExperimentalAlgae.csv", header=T)
start<-read.csv("C:/Users/thsmith/Desktop/Consumer Resource Experiment/Summer 2010/2010MesocosmCSVdata/StartingPtAlgaeSamples.csv", header=T)
str(whole)
str(start)
names(whole)
levels(whole$Treatment)
levels(start$Treatment)
whole<-rbind(whole, start)

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
#############
#############
#restrict data set for growth rate analysis
bottom<-whole[which(whole$SubstrateType=="BottomTile"),]
whole2<-bottom[-which(whole$AFDM<0.0005),]
names(whole2)
# [1] "TankNumber"      "SampleDate"      "SampleNumber"    "SubstrateType"   "TileNumber"      "FilterNumber"    "FilterWt"       
# [8] "DryWt"           "AshedWt"         "AFDM"            "FilterVol"       "AFDMadjbyVol"    "GrowthPeriod"    "AFDMadjPerMperD"
#[15] "TadpolePDensity" "MFDensity"       "FinalMFCount"    "Treatment" 

#manipulate data to calculate mean abundances in each tank/treatment combination, across sample dates.
S0all<-which(whole2$SampleNumber=="S0")
S1all<-which(whole2$SampleNumber=="S1")
S3over20<-which(whole2$SampleNumber=="S3" & whole2$GrowthPeriod>=20)
whole3<-whole2[c(S0all,S1all,S3over20),c(1,3,12,18)]
meltwhole<-melt(whole3, id=c("TankNumber","SampleNumber", "Treatment"))
TreatMeans<-cast(meltwhole, TankNumber+Treatment~SampleNumber, max)
#calculate growth rates:
S0S1<-log(TreatMeans$S1/TreatMeans$S0)
S0S3<-log(TreatMeans$S3/TreatMeans$S0)
whole4<-whole2[c(S0all,S1all,S3over20),c(1,3,18,13)]
meltwhole4<-melt(whole4, id=c("TankNumber","SampleNumber", "Treatment"))
TankSamplePeriods<-cast(meltwhole4, TankNumber~SampleNumber, mean)

Growth<-data.frame(rep(TreatMeans$TankNumber, 2), 
			rep(TreatMeans$Treatment, 2), 
			c(S0S1, S0S3), c(rep("S0S1", 16),rep("S0S3", 16)), 
			c(TankSamplePeriods$S1,TankSamplePeriods$S3))
names(Growth)
#"rep.TreatMeans.TankNumber..3."                       
#"rep.TreatMeans.Treatmen..3."                         
#"c.S0S1..S0S2..S0S3."                                 
#"c.rep..S0S1...16...rep..S0S2...16...rep..S0S3...16.."
#"c.rep..S0S1...16...rep..S0S3...16.." 
names(Growth)<-c("TankNumber", "Treatment", "GrowthRate", "SamplePeriod","GrowthPeriod")
str(Growth)
Growth$RealGrowthRate<--Growth$GrowthRate/Growth$GrowthPeriod
ggplot(Growth[which(Growth$SamplePeriod=="S0S3"),], aes(x=Treatment, y=RealGrowthRate))+
	geom_boxplot()+
	facet_grid(.~SamplePeriod)
ggplot(Growth[which(Growth$SamplePeriod=="S0S1"),], aes(x=Treatment, y=RealGrowthRate))+
	geom_boxplot()+
	facet_grid(.~SamplePeriod)
ggplot(Growth, aes(x=Treatment, y=RealGrowthRate))+
	geom_boxplot()
anova(aov(GrowthRate~Treatment, data=Growth))
TukeyHSD(aov(GrowthRate~Treatment, data=Growth[which(Growth$SamplePeriod=="S0S1"),]))



ggplot(whole2[whole2$SampleNumber=="S0",], aes(Treatment, AFDM))+
	geom_boxplot()
anova(aov(log(AFDM)~Treatment, data=whole2[whole2$SampleNumber=="S0",]))
TukeyHSD(aov(log(AFDM)~Treatment, data=whole2[whole2$SampleNumber=="S0",]))
#there is no difference among these, but the trend follows the trend that we end up with,
#so, definitely do the growth rates of the S1, S2, S3 abundance relative to the S0




#################
#PLOTTING:
#
whole3$TxOrder<-factor(whole3$Treatment, levels=c("none", "mayflies", "tadpoles", "both"))
llcoolJ<-ggplot(whole3[whole3$SampleNumber!="S0",], aes(x=GrowthPeriod, y=log(AFDMadjbyVol+0.00001), group=TxOrder))+
	geom_smooth(method=lm, aes(linetype=TxOrder, color=TxOrder), size=2, fill="light grey", fullrange=TRUE)+	#alternatively use method=lm
	scale_linetype_manual(values=c("dotted", "longdash", "solid", "dotdash"))+
	scale_color_manual(values=c("black", "black", "black", "black"))+
	guides(linetype=guide_legend(keywidth=5, keyheight=0.5, label.position="bottom", label.hjust = 0.5))+
	xlab("Days")+ylab("ln (mg AFDM)")+
	theme(legend.position=c("bottom"),
		legend.title=element_blank(),
		legend.text=element_text(size=14, face="bold"),
		panel.background=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=14, face="bold"),	
		axis.title=element_text(color="black", size=14, face="bold")
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2010_logAlgalGrowth_byTreatment.png", 
	res=300, width=5, height=5, units="in", restoreConsole=F)
llcoolJ
dev.off()


ggplot(whole3[whole3$SampleNumber!="S0",], aes(x=as.factor(GrowthPeriod), y=log(AFDMadjbyVol+0.00001), group=TxOrder))+
	geom_boxplot(aes(color=TxOrder), size=2)+	#alternatively use method=lm
	xlab("Days")+ylab("ln (mg AFDM)")+
	theme(legend.position=c("bottom"),
		legend.title=element_blank(),
		legend.text=element_text(size=14, face="bold"),
		panel.background=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=14, face="bold"),	
		axis.title=element_text(color="black", size=14, face="bold")
		)


ggplot(whole3[whole3$SampleNumber!="S0",], aes(x=Treatment, y=log(AFDMadjbyVol+0.00001)))+
	geom_boxplot(aes(color=SampleNumber))










