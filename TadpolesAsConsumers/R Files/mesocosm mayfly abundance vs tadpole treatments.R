library(reshape2)
library(reshape)
library(ggplot2)
#

TadDens<-rep(c(16, 16, 0,0,16,0,16,0),2)
MayflyStart<-rep(250,8)
MayflyEnd<-c(7,48,52,0,18,121,100,76)

Tank<-as.factor(seq(1:8))
MayflyPpnStart<-rep(1,8)
MayflyPpnofStart<-MayflyEnd/250
StartEnd<-factor(c(rep("Start",8),rep("End", 8)), levels=c("Start", "End"))
MFSample<-c(rep(1,8),rep(2, 8))
MFDens<-c(MayflyPpnStart, MayflyPpnofStart)
MFdfSE<-data.frame(Tank,StartEnd, MFSample, MFDens, TadDens)
MFdf<-MFdfSE[which(MFdfSE$StartEnd=="End"),]
anova(aov(MFDens~StartEnd, MFdfSE))

png(filename="C:\\Users\\thsmith\\Desktop\\Consumer Resource Experiment\\Figures\\2010_MayflyDynamics.png", 
	res=300, width=4, height=3, units="in", restoreConsole=F)
ggplot(MFdf, aes(x=as.factor(TadDens), y=MFDens*100))+
	geom_point(size=6, shape=16, fill="black")+
#	geom_point(aes(fill=as.factor(TadDens)), size=5, shape=21)+
	scale_fill_manual(values=c("black","black"), limits=c(0,0.5))+#, labels=c("No Tadpoles", "Tadpoles\nPresent"))+
#	guides(fill=guide_legend(title=NULL,keywidth = 1, keyheight = 2, label.position="left", label.hjust=0.5))+
	guides(fill=FALSE)+
#	scale_y_continuous(breaks=c(10,20,30,40,50), limits=c(0,50))+
	xlab("Tadpoles")+ylab("% remaining\n")+
#	annotate("text", x=0.85, y=0.85, label="a", fontface="bold.italic", size=7)+
#	annotate("text", x=1.2, y=0.85, label="a", fontface="bold.italic", size=7)+
#	annotate("text", x=1.8, y=0.55, label="b", fontface="bold.italic", size=7)+
#	annotate("text", x=2.2, y=0.55, label="b", fontface="bold.italic", size=7)+
	theme(
		#legend.position=c(1.025,1.1),
		#legend.background=element_blank(),
		#legend.justification=c(1,1),
		#legend.text=element_text(face="bold",size=14),
		#legend.key=element_blank(),
		panel.background=element_blank(),
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

#----------------------------------
lengths<-read.csv("C://Users//thsmith//Desktop//Consumer Resource Experiment//Data Files//2010 Mesocosms//MesocosmMayflyLengthsByTankandTreatment_YishenMiaoNov2014.csv", header=TRUE)
str(lengths)
summary(lengths)
LengthSummary<-aggregate(lengths$MayflyLength, by=list(lengths$TankNumber, lengths$TadpoleDensity, lengths$MayflyDensity), mean)
LengthSummary<-rename(LengthSummary, replace=c(Group.1="Tank", Group.2="TadDens", Group.3="MFDens", "x"="meanLength"))
anova(aov(meanLength~TadDens, data=LengthSummary))
boxplot(meanLength~TadDens, data=LengthSummary)