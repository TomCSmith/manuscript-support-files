qdetach()
rm(list=ls())
#############
#load packages:
library(nlme)
library(lme4)
library(ggplot2)
library(plyr)
library(ggthemes)
library(reshape2)
library(reshape)
library(scales)
library(gridExtra)
#############
PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	qqnorm(resid(x))
	plot(resid(x)~fitted(x), data=y)
	plot(resid(x)~y[,1], data=y, xlab=colnames(y[1]))
	plot(resid(x)~y[,3], data=y, xlab=colnames(y[3]))
	plot(resid(x)~log(y[,6]), data=y, xlab=colnames(y[6]))
	boxplot(resid(x)~y[,7], data=y, xlab=colnames(y[7]))
	boxplot(resid(x)~y[,8], data=y, xlab=colnames(y[8]))
	plot(resid(x)~y[,9], data=y, xlab=colnames(y[9]))
	plot(resid(x)~y[,10], data=y, xlab=colnames(y[10]))
	}
#PlotRes(model, dataset, plot dimensions)
#####################

mayfly<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\Mayfly data\\YISHEN\\mayfly.csv", header=T)
mayfly<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\Mayfly data\\MFLength_Genus_FinalSample_Densities_EnclosureCharacters.csv", header=T)
str(mayfly)
# [1] "LakeName"        
# [2] "SampleNumber"    
# [3] "CageNumber"      
# [4] "AlgMeans"
# [5] "TadDens"         
# [6] "MFDens"          
# [7] "Silt"            
# [8] "Radiation"       
# [9] "GenusSpp"        
# [10] "MFLength"        
#get mayfly length means and SEs
df1<-mayfly[,c(1,3,10,9)]
df1m<-melt(df1, id=c("LakeName", "CageNumber", "GenusSpp"))
str(df1m)
MFlengthMeans<-cast(df1m, LakeName+CageNumber+GenusSpp~variable, mean)
StanErr<-function(x){
	sd(x)/sqrt(length(x))}
MFlengthSE<-cast(df1m, LakeName+CageNumber+GenusSpp~variable, StanErr)[,4]

#get Algal AFDM
dfAlg<-mayfly[,c(1,3,9,4)]
dfAlgm<-melt(dfAlg, id=c("LakeName", "CageNumber", "GenusSpp"))
AlgMeans<-cast(dfAlgm, LakeName+CageNumber+GenusSpp~variable, mean)[,4]

#get MF Dens
dfMF<-mayfly[,c(1,3,9,6)]
dfMFm<-melt(dfMF, id=c("LakeName", "CageNumber", "GenusSpp"))
MF<-cast(dfMFm, LakeName+CageNumber+GenusSpp~variable, mean)[,4]

#get tadpole density
dfTad<-mayfly[,c(1,3,9,5)]
dfTadm<-melt(dfTad, id=c("LakeName", "CageNumber", "GenusSpp"))
Tad<-cast(dfTadm, LakeName+CageNumber+GenusSpp~variable, mean)[,4]

#Get radition
dfRad<-mayfly[,c(1,3,9,8)]
dfRadm<-melt(dfRad, id=c("LakeName", "CageNumber", "GenusSpp"))
Rad<-cast(dfRadm, LakeName+CageNumber+GenusSpp~variable, mean)[,4]

#get silt
dfSad<-mayfly[,c(1,3,9,7)]
dfSadm<-melt(dfSad, id=c("LakeName", "CageNumber", "GenusSpp"))
Silt<-cast(dfSadm, LakeName+CageNumber+GenusSpp~variable, mean)[,4]

#put all together in a useful data frame:
MFSummary<-data.frame(MFlengthMeans, MFlengthSE, AlgMeans, MF, Tad, Rad, Silt)
str(MFSummary)
#'data.frame':   35 obs. of  10 variables:
# $ LakeName  : Factor w/ 2 levels "LeConte","Spur": 1 1 1 1 1 1 1 1 1 1 ...
# $ CageNumber: int  1 1 2 4 5 6 7 8 10 13 ...
# $ GenusSpp  : Factor w/ 2 levels "Ameletus spp.",..: 1 2 1 1 1 1 1 1 1 1 ...
# $ MFLength  : num  8.31 5.64 10.31 9.08 9.12 ...
# $ MFlengthSE: num  0.347 0.237 0.142 0.219 0.201 ...
# $ AlgMeans  : num  0.004 0.004 0.0045 0.0036 0.002 0.0047 0.0043 0.0039 0.0016 0.0007 ...
# $ MF        : num  125 125 125 125 250 25 25 125 250 250 ...
# $ Tad       : num  20 20 2 0 0 20 2 10 20 2 ...
# $ Rad       : num  1950 1950 1950 1950 1950 1950 1970 1970 1970 1970 ...
# $ Silt      : num  0 0 100 62.5 

################
#basic plotting
plot(MFSummary[,c(1,4,5,6,7,8,9,10)])
ggplot(MFSummary, aes(x=as.factor(Tad), y=MFLength))+
	geom_boxplot()+
	facet_grid(.~LakeName)
ggplot(MFSummary, aes(x=as.factor(MF), y=MFLength))+
	geom_boxplot()+
	facet_grid(.~LakeName)
windows()
ggplot(MFSummary, aes(x=as.factor(Tad), y=MFLength))+
	geom_boxplot(aes(fill=GenusSpp), color="black")+
	facet_grid(.~LakeName)+
	theme(legend.justification=c(1,1), legend.position=c(1,1))
windows()
ggplot(MFSummary, aes(x=as.factor(MF), y=MFLength))+
	geom_boxplot(aes(fill=GenusSpp), color="black")+
	facet_grid(.~LakeName)+
	theme(legend.justification=c(1,1), legend.position=c(1,1))

ggplot(MFSummary, aes(x=MFLength))+
	geom_histogram(binwidth=1)+
	facet_grid(LakeName~GenusSpp)

###############
#models:
lm1<-lm(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+MF+GenusSpp*LakeName, 
	data=MFSummary)
summary(lm1)
PlotRes(lm1, MFSummary, c(2,5))
#there def. seems to be a downward trend in the fitted v. residuals, or something, maybe additive, though not much...
#add variance structure on Lake, maybe tadpole density, spp and MF

lm2<-gls(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+MF+GenusSpp*LakeName, 
	method="ML",
	data=MFSummary)
PlotRes(lm2, MFSummary, c(2,5))

LakeVar<-varIdent(form=~1|LakeName)
#TadVar<-varFixed(~Tad)	#should be variance across a continuous variable but that causes memory problems
TadVar<-varIdent(form=~1|Tad)	#
MFVar<-varIdent(form=~1|MF)
SppVar<-varIdent(form=~1|GenusSpp)
GlsCtrl<-glsControl(maxIter=999, msMaxIter=999, tolerance=0.001, msTol=0.001)
lm3<-gls(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+MF+GenusSpp*LakeName, 
	weights=LakeVar,
	method="ML",
	data=MFSummary)
summary(lm3)
PlotRes(lm3, MFSummary, c(2,5))
anova(lm2, lm3)	# better

lm4<-gls(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+MF+GenusSpp*LakeName, 
	weights=TadVar,
	method="ML",
	data=MFSummary, 
	control=GlsCtrl)
summary(lm4)
PlotRes(lm4, MFSummary, c(2,5))	# maybe less normal, but not much 
anova(lm2, lm4)

lm5<-gls(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+MF+GenusSpp*LakeName, 
	weights=MFVar,
	method="ML",
	data=MFSummary)
summary(lm5)
PlotRes(lm5, MFSummary, c(2,5))	#
anova(lm2, lm5)	#betterish not drmatic or significant 

lm6<-gls(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+MF+GenusSpp*LakeName, 
	weights=SppVar,
	method="ML",
	data=MFSummary)
summary(lm6)
PlotRes(lm6, MFSummary, c(2,5))	#
anova(lm2, lm6)	#big improvement

#so which do we include in varComb
lm9<-gls(MFLength~log(AlgMeans)+Tad*MF+Silt+Rad+GenusSpp*LakeName, 
	weights=varComb(SppVar, TadVar),
	method="ML",
	data=MFSummary)
summary(lm9)
PlotRes(lm9, MFSummary, c(2,5))	#
anova(lm2, lm4, lm5, lm6, lm7, lm8, lm9)	#

#additional combinations of variance structure not helpful
#step down lm9, drop each interaction in turn and compare to full model:
#dropping the Genus x Lake interaction is a much worse model; 
#dropping the Tad x MF interaction at this stage is no different, keep both interactions
#drop log(AlgMeans)
lm20<-gls(MFLength~Tad*MF+Silt+Rad+GenusSpp*LakeName, 
	weights=varComb(SppVar, TadVar),
	method="ML",
	data=MFSummary)
summary(lm20)
PlotRes(lm20, MFSummary, c(2,5))	#
anova(lm9, lm20)	#slight impvt. drop Tad x MF interaction

lm21<-gls(MFLength~Tad+MF+Silt+Rad+GenusSpp*LakeName, 
	weights=varComb(SppVar, TadVar),
	method="ML",
	data=MFSummary)
summary(lm21)
PlotRes(lm21, MFSummary, c(2,5))	#
anova(lm20, lm21)	#barely differnt, drop nothing else. refit

lm22<-gls(MFLength~Tad+MF+Rad+GenusSpp*LakeName, 
	weights=varComb(SppVar, TadVar),
	method="ML",
	data=MFSummary)
summary(lm22)
PlotRes(lm22, MFSummary, c(2,5))	#
anova(lm21, lm22)	#barely differnt, drop nothing else. refit

lm23<-gls(MFLength~Tad+MF+GenusSpp*LakeName, 
	weights=varComb(SppVar, TadVar),
	method="ML",
	data=MFSummary,
	control=GlsCtrl)
summary(lm23)
PlotRes(lm23, MFSummary, c(2,5))	#
anova(lm22, lm23)	#barely differnt, drop nothing else. refit

lm21.reml<-gls(MFLength~Tad+MF+GenusSpp*LakeName, 
	weights=varComb(SppVar, TadVar),
	method="REML",
	data=MFSummary)
summary(lm21.reml)
PlotRes(lm21.reml, MFSummary, c(2,5))	#

#report and interpret...
#percent decreases:
#effect of tads on MFlenght in LeConte
#mf length, leconte, t=0
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$Tad==0)])	#9.47
#mf length, leconte, t=2
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$Tad==2)]) #8.89
#mf length, leconte, t=10
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$Tad==10)]) #7.83
#mf length, leconte, t=20
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$Tad==20)]) #8.45
8.89/9.47# 6%
7.83/9.47# 17%

#effect of tads on MFlenght in Spur; no reduction
anova(aov(MFLength~as.factor(Tad), data=MFSummary[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Ameletus spp."),]))
anova(aov(MFLength~as.factor(Tad), data=MFSummary[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Callibaetis spp."),]))

#effect of MF on MFlenght in LeConte
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$MF==25)])	#9.65
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$MF==125)]) #8.88
mean(MFSummary$MFLength[which(MFSummary$LakeName=="LeConte" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$MF==250)]) #7.45
8.88/9.65# 8%
7.45/9.65# 23%

#effect of MF on MFlenght in Spur
anova(aov(MFLength~as.factor(MF), data=MFSummary[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Ameletus spp."),]))		#no difference
anova(aov(MFLength~as.factor(MF), data=MFSummary[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Callibaetis spp."),]))	#
TukeyHSD(aov(MFLength~as.factor(MF), data=MFSummary[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Callibaetis spp."),]))	#
#125-25  -0.46922987 -0.9895217 0.051061968 0.0755389
#250-25  -0.55266309 -1.1146425 0.009316321 0.0536064
#250-125 -0.08343322 -0.6454126 0.478546189 0.9066801

#AMEL not stat. different, but
mean(MFSummary$MFLength[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$MF==25)])	#
mean(MFSummary$MFLength[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$MF==125)]) #
mean(MFSummary$MFLength[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Ameletus spp." & MFSummary$MF==250)]) #
5.17/5.81	#11%
4.83/5.81	#17%

#effect of MF on MFlenght in Spur - CALB
mean(MFSummary$MFLength[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Callibaetis spp." & MFSummary$MF==25)])	#
mean(MFSummary$MFLength[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Callibaetis spp." & MFSummary$MF==125)]) #
mean(MFSummary$MFLength[which(MFSummary$LakeName=="Spur" & MFSummary$GenusSpp=="Callibaetis spp." & MFSummary$MF==250)]) #
5.48/5.95	#8%
5.3994/5.95	#10%



#boxplots:
MFSummaryConsCat<-MFSummary
MFSummaryConsCat$MF[MF==25]<-"low"
MFSummaryConsCat$MF[MF==125]<-"medium"
MFSummaryConsCat$MF[MF==250]<-"high"
MFSummaryConsCat$Tad[Tad==0]<-"zero"
MFSummaryConsCat$Tad[Tad==2]<-"low"
MFSummaryConsCat$Tad[Tad==10]<-"medium"
MFSummaryConsCat$Tad[Tad==20]<-"high"
MFSummaryMelt<-melt(MFSummaryConsCat[,c(1,2,3,4,7,8)], id=c("LakeName", "CageNumber", "GenusSpp", "MFLength"))
MFSummaryMelt<-rename(MFSummaryMelt, c(variable="consumer"))
attach(MFSummaryMelt)
MFSummaryMelt$consumertype[consumer=="Tad"]<-"Tadpole"
MFSummaryMelt$consumertype[consumer=="MF"]<-"Mayfly"
detach(MFSummaryMelt)
plot3<-ggplot(MFSummaryMelt, aes(x=as.factor(value), y=MFLength, fill=GenusSpp))+
	geom_boxplot(size=1.5)+
	facet_grid(LakeName~consumertype, scales="free")+
	scale_fill_manual(values=c("white", "dark grey"))+
	xlab("\nConsumer Density")+ylab("\nMayfly length (mm)")+
	scale_x_discrete(limits=c("zero","low","medium","high"))+
	theme(#legend.position=c(0,0.5),
		legend.position="bottom",
		legend.key=element_rect(fill=NA, color=NA),
		legend.key.size=unit(0.25, "in"),
		legend.justification=c(0,0),
		legend.title=element_blank(),
		legend.background=element_rect(fill=NA, color="black"),
		legend.text=element_text(size=10, face="bold"),
		panel.background=element_blank(),
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank(),
		panel.border=element_rect(fill=NA),
		strip.background=element_rect(color="black", fill=NA),
		strip.text=element_text(color="black", size=14, face="bold"),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=10, face="bold"),	
		axis.title=element_text(color="black", size=14, face="bold")
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_MayflyLength_boxplots_byconsumer_lake.png", 
	res=300, width=5, height=6, units="in", restoreConsole=F)
plot3
dev.off()





plot4<-ggplot(MFSummaryMelt, aes(x=as.factor(value), y=MFLength))+
	geom_boxplot(size=1.5, color="black")+
	facet_grid(LakeName~consumertype, scales="free")+
#	scale_fill_manual(values=c("white", "dark grey"))+
	xlab("\nConsumer Density")+ylab("\nMayfly length (mm)")+
	scale_x_discrete(limits=c("zero","low","medium","high"))+
	theme(panel.background=element_blank(),
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank(),
		panel.border=element_rect(fill=NA),
		strip.background=element_rect(color="black", fill=NA),
		strip.text=element_text(color="black", size=14, face="bold"),
		axis.line=element_line(color="black", size=2), 
		axis.text=element_text(color="black", size=10, face="bold"),	
		axis.title=element_text(color="black", size=14, face="bold")
		)
png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_MayflyLength_boxplots_byconsumer_lake_nospecies.png", 
	res=300, width=5, height=4, units="in", restoreConsole=F)
plot4
dev.off()

#-------------------------------------------------------------------
#Emergence
emerge<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\MF_Emergence_byBlockTreatment.csv", header=T)
str(emerge)
emerge$LakeID<-as.factor(emerge$LakeID)
emerge$Block<-as.factor(emerge$Block)
plot(emerge)

ggplot(emerge, aes(y=EmergedPct, x=Block))+
	geom_boxplot(size=3, alpha=0.5)+
	facet_grid(.~LakeID)
lm1<-lm(EmergedPct~TadDens*MFDens+Block+LakeID, data=emerge)
summary(lm1)
PlotRes<-function(x,y,z){
	par(mfrow=z)
	hist(resid(x))
	qqnorm(resid(x))
	plot(resid(x)~fitted(x), data=y)
	boxplot(resid(x)~y[,1], data=y, xlab=colnames(y[1]))
	boxplot(resid(x)~y[,2], data=y, xlab=colnames(y[2]))
	plot(resid(x)~log(y[,4]), data=y, xlab=colnames(y[4]))
	plot(resid(x)~y[,5], data=y, xlab=colnames(y[5]))
	}
#PlotRes(model, dataset, plot dimensions)
PlotRes(lm1, emerge, c(2,4))
lm2<-gls(EmergedPct~
	TadDens*MFDens*Block*LakeID,
	method="ML",
	data=emerge)
summary(lm2)
lm3<-lme(EmergedPct~
	TadDens*MFDens+Block,
	random=~1|LakeID,
	method="ML",
	data=emerge)
anova(lm2, lm3)
PlotRes(lm3, emerge, c(2,4))
lm4<-lme(EmergedPct~
	TadDens*MFDens,
	random=~1|LakeID/Block,
	method="ML",
	data=emerge)
anova(lm2, lm3, lm4)
PlotRes(lm3, emerge, c(2,4))
summary(lm4)

lm5<-lm(EmergedPct~TadDens, data=emerge, subset=which(emerge$LakeID=="10102" & emerge$Block=="1"))
summary(lm5)

#no effect of consumer density on emergence rates.

anova(aov(EmergedPct~Block*LakeID, data=emerge))
TukeyHSD(aov(EmergedPct~Block*LakeID, data=emerge))




#--------------------------------------
#mayfly density dynamics
MFDens<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\MayflyTreatmentDynamics.csv", header=T)
str(MFDens)
MFDens$LakeID<-as.factor(MFDens$LakeID)
MFDens$CageNumber<-as.factor(MFDens$CageNumber)

cage1<-which(MFDens$CageNumber=="1")
LC<-which(MFDens$LakeName=="LeConte")
v<-which(MFDens$LakeName=="LeConte" & MFDens$CageNumber=="1")
Zero<-which(MFDens$MayflyTx==0)
Start<-which(MFDens$StartEnd=="Start")

MFDens$BlockName[MFDens$SampleNumber=="S2"]<-"Jul-Aug"
MFDens$BlockName[MFDens$SampleNumber=="S3"]<-"Mid Aug"
MFDens$BlockName[MFDens$SampleNumber=="S4"]<-"Aug-Sep"
MFDens$BlockName<-factor(MFDens$BlockName, levels=c("Jul-Aug","Mid Aug","Aug-Sep"))
MFDens$LakeName[MFDens$LakeID=="10102"]<-"LeConte"
MFDens$LakeName[MFDens$LakeID=="10475"]<-"Spur"
MFDens$PpnEmerged<-1-(MFDens$MFCount/MFDens$MayflyTx)
MFDens$Days[MFDens$SampleNumberInt=="2" & MFDens$LakeID=="10102"]<-16
MFDens$Days[MFDens$SampleNumberInt=="3" & MFDens$LakeID=="10102"]<-23
MFDens$Days[MFDens$SampleNumberInt=="4" & MFDens$LakeID=="10102"]<-19
MFDens$Days[MFDens$SampleNumberInt=="2" & MFDens$LakeID=="10475"]<-23
MFDens$Days[MFDens$SampleNumberInt=="3" & MFDens$LakeID=="10475"]<-18
MFDens$Days[MFDens$SampleNumberInt=="4" & MFDens$LakeID=="10475"]<-18


png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2009_MayflyDensityChange_byLake_byTx.png", 
	res=300, width=8, height=6, units="in", restoreConsole=F)
ggplot(MFDens[-Zero,][-Start,], aes(x=as.factor(MayflyTx),y=(100*PpnEmerged)))+
	geom_boxplot(aes(fill=as.factor(MayflyTx)), size=2)+
	facet_grid(LakeName~BlockName)+
	scale_fill_manual(values=c("white","white","white"))+
	guides(fill=FALSE)+
	scale_y_continuous(breaks=c(25,50,75,100), limits=c(0,100))+
	xlab("\nMayfly Treatment")+ylab("% lost\n")+
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

v<-lm(PpnEmerged~LakeName*BlockName, data=MFDens[-Zero,][-Start,])
summary(v)
drop1(v,~.,test="F")
TukeyHSD(v)

q<-which(MFDens$SampleNumber!="S0" & MFDens$MayflyTx!=0 & MFDens$StartEnd=="End")
v<-aov(PpnEmerged~LakeName*BlockName, data=MFDens[q,])
summary(v)
drop1(v,~.,test="F")
TukeyHSD(v)


#-------------------------------------
#Cage light comparisons:
light<-read.csv("C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Data Files\\2009 Enclosures\\CageLight.csv", header=T)
str(light)
anova(aov(Light~LtLoc, data=light))
TukeyHSD(aov(Light~LtLoc, data=light))
sd(light$Light[which(light$LtLoc=="ExtLt")])/sqrt(length(light$Light[which(light$LtLoc=="ExtLt")]))
sd(light$Light[which(light$LtLoc=="IntLt")])/sqrt(length(light$Light[which(light$LtLoc=="IntLt")]))
