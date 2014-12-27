library(nlme)
mfdata<-read.csv(file.choose(), header=T)
str(mfdata)
mfdata$Lake<-as.factor(mfdata$Lake)
mfdata$TadpoleDensity<-as.factor(mfdata$TadpoleDensity)
boxplot(Length~TadpoleDensity, data=mfdata[mfdata$Lake=="10102",])

#anova of all data:
summary(aov(Length~TadpoleDensity+Lake, data=mfdata))

#anova of 10475:
anova(aov(Length~TadpoleDensity, data=mfdata[mfdata$Lake=="10475",]))
#anova of 10102
anova(aov(Length~TadpoleDensity, data=mfdata[mfdata$Lake=="10102",]))

#plot means for MF lenght 10475
meanLength10475T0<-mean(mfdata$Length[mfdata$Lake=="10475" & mfdata$TadpoleDensity==0])
seLength10475T0<-(sd(mfdata$Length[mfdata$Lake=="10475" & mfdata$TadpoleDensity==0]))/(sqrt(length(mfdata$Length[mfdata$Lake=="10475" & mfdata$TadpoleDensity==0])))
meanLength10475T20<-mean(mfdata$Length[mfdata$Lake=="10475" & mfdata$TadpoleDensity==20])
seLength10475T20<-(sd(mfdata$Length[mfdata$Lake=="10475" & mfdata$TadpoleDensity==20]))/(sqrt(length(mfdata$Length[mfdata$Lake=="10475" & mfdata$TadpoleDensity==20])))
means10475<-c(meanLength10475T0, meanLength10475T20)
se10475<-c(seLength10475T0, seLength10475T20)

par(mar=c(3,5,1,0), bg="transparent", fg="light grey", cex.lab=1.5, cex.axis=1.5, font=2, lwd=2, cex=1.5)
barplot(means10475, 
	width=1,
	space=c(0,0),
	xlim=c(0,2),
	ylim=c(0,10),
	col=c("transparent", "olivedrab4"),
	col.lab="light grey", col.axis="light grey", font.axis=2, cex.lab=2)
mtext("mean length +/- s.e. (mm)", side=2, line=3, at=5, font=2, col="light grey", adj=0.5, cex=2)
mtext("0", side=1, line=0.5, at=0.5, font=2, col="light grey", adj=0.5, cex=2)
mtext("20", side=1, line=0.5, at=1.5, font=2, col="light grey", adj=0.5, cex=2)
mtext("Tadpole density", side=1, line=2, at=1, font=2, col="light grey", adj=0.5, cex=2)
lines(c(0.5,0.5),c(means10475[1], means10475[1]+se10475[1]), lwd=5)
lines(c(0.5,0.5),c(means10475[1], means10475[1]-se10475[1]), lwd=5)
lines(c(1.5,1.5),c(means10475[2], means10475[2]+se10475[2]), lwd=5)
lines(c(1.5,1.5),c(means10475[2], means10475[2]-se10475[2]), lwd=5)
text("p = 0.11", x=1, y=9.5, cex=2, font=2)

#MF length means 10102
meanLength10102T0<-mean(mfdata$Length[mfdata$Lake=="10102" & mfdata$TadpoleDensity==0])
seLength10102T0<-(sd(mfdata$Length[mfdata$Lake=="10102" & mfdata$TadpoleDensity==0]))/(sqrt(length(mfdata$Length[mfdata$Lake=="10102" & mfdata$TadpoleDensity==0])))
meanLength10102T20<-mean(mfdata$Length[mfdata$Lake=="10102" & mfdata$TadpoleDensity==20])
seLength10102T20<-(sd(mfdata$Length[mfdata$Lake=="10102" & mfdata$TadpoleDensity==20]))/(sqrt(length(mfdata$Length[mfdata$Lake=="10102" & mfdata$TadpoleDensity==20])))
means10102<-c(meanLength10102T0, meanLength10102T20)
se10102<-c(seLength10102T0, seLength10102T20)

par(mar=c(3,5,1,0), bg="transparent", fg="light grey", cex.lab=1.5, cex.axis=1.5, font=2, lwd=2, cex=1.5)
barplot(means10102, 
	width=1,
	space=c(0,0),
	xlim=c(0,2),
	ylim=c(0,10),
	col=c("transparent", "olivedrab4"),
	col.lab="light grey", col.axis="light grey", font.axis=2, cex.lab=2)
mtext("mean length +/- s.e. (mm)", side=2, line=3, at=5, font=2, col="light grey", adj=0.5, cex=2)
mtext("0", side=1, line=0.5, at=0.5, font=2, col="light grey", adj=0.5, cex=2)
mtext("20", side=1, line=0.5, at=1.5, font=2, col="light grey", adj=0.5, cex=2)
mtext("Tadpole density", side=1, line=2, at=1, font=2, col="light grey", adj=0.5, cex=2)
lines(c(0.5,0.5),c(means10102[1], means10102[1]+se10102[1]), lwd=5)
lines(c(0.5,0.5),c(means10102[1], means10102[1]-se10102[1]), lwd=5)
lines(c(1.5,1.5),c(means10102[2], means10102[2]+se10102[2]), lwd=5)
lines(c(1.5,1.5),c(means10102[2], means10102[2]-se10102[2]), lwd=5)
text("p = 0.022", x=1, y=9.5, cex=2, font=2)


