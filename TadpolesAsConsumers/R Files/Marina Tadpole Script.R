#Marina Bozinovic, 7/23/14, This is an analysis of tadpole weights relative to tadpole density and mayfly density.
detach()
rm(list=ls())
library(stats)
library(reshape)
################
#load data:
taddata <- read.csv("\\\\BRIGGS-NAS1\\Thomas Smith\\Tadpole Data\\SummarizedTadpoleData.csv", header = T)
str(taddata)
taddata$LakeID <- as.factor(taddata$LakeID)
taddata$CageNumber <- as.factor(taddata$CageNumber)
attach(taddata)
taddata[TadpoleDensity>20,]
taddatacage <- taddata[taddata$CageNumber!="0",]
plot(taddatacage)
str(taddatacage)
attach(taddatacage)

#Categorical Barplots
Cage.freq <- table(taddatacage$CageNumber)
barplot(Cage.freq)
barplot(Cage.freq[order(Cage.freq)],horiz=T)
str(Cage.freq)
cagecolor<-c(rep("gray",5),
             rgb(59,89,152,maxColorValue=255))
barplot(Cage.freq[order(Cage.freq)],
        horiz = T,col = cagecolor)
barplot(Cage.freq[order(Cage.freq)],
        border = NA,
        main = "Frequency of Cages\nin experiment",
        col = cagecolor,
        xlim = c(0,18),
        xlab = "Frequency of Cages") 

lake1 <- taddatacage[taddatacage$LakeID=="10102",]
cage1 <- lake1[lake1$CageNumber=="1",]
cage2 <- lake1[lake1$CageNumber=="2",]
cage3 <- lake1[lake1$CageNumber=="3",]
cage4 <- lake1[lake1$CageNumber=="4",]
cage5 <- lake1[lake1$CageNumber=="5",]
cage6 <- lake1[lake1$CageNumber=="6",]
cage7 <- lake1[lake1$CageNumber=="7",]
cage8 <- lake1[lake1$CageNumber=="8",]
cage9 <- lake1[lake1$CageNumber=="9",]
cage10 <- lake1[lake1$CageNumber=="10",]
cage12 <- lake1[lake1$CageNumber=="12",]
cage13 <- lake1[lake1$CageNumber=="13",]
cage14 <- lake1[lake1$CageNumber=="14",]
cage17 <- lake1[lake1$CageNumber=="17",]

rm(cage4)
avgl1c1 <- mean(cage1$EndAvg)
avgl1c2 <- mean(cage2$EndAvg)
avgl1c3 <- mean(cage3$EndAvg)
avgl1c6 <- mean(cage6$EndAvg)
avgl1c7 <- mean(cage7$EndAvg)
avgl1c8 <- mean(cage8$EndAvg)
avgl1c9 <- mean(cage9$EndAvg)
avgl1c10 <- mean(cage10$EndAvg)
avgl1c12 <- mean(cage12$EndAvg)
avgl1c13 <- mean(cage13$EndAvg)
avgl1c14 <- mean(cage14$EndAvg)
avgl1c17 <- mean(cage17$EndAvg)

lake1.avg <- matrix(c(avgl1c1,avgl1c2,avgl1c3,avgl1c6,avgl1c7,avgl1c8,avgl1c9,avgl1c10,avgl1c12,avgl1c13,avgl1c14,avgl1c17))
rownames (lake1.avg) <- c("Cage1","Cage2","Cage3","Cage6","Cage7","Cage8","Cage9","Cage10","Cage12","Cage13","Cage14","Cage17")
colnames (lake1.avg) <- c("Average.weights")
barplot(lake1.avg, beside = T, ylim = c(0,50), names = rownames(lake1.avg))

#Creating Error Bars

#arrows ((1:12)+0.5, lake1.avg+1, (1:12)+0.5, lake1.avg-1, angle=90, code=3, length=.1)
s <- tapply (lake1$EndStndDev,lake1$CageNumber, sd)
n <- tapply (lake1$EndStndDev,lake1$CageNumber, length)
sem <- s/sqrt(n)
sem2 <- sem[c(1,2,3,4,5,6,7,8,9,10,11,12)]
UpSem2 <-data.frame(data.frame(lake1.avg)+data.frame(sem2))
LowSem2<-data.frame(data.frame(lake1.avg)-data.frame(sem2))
arrows ((1:12)+0.5, UpSem2$Average.weights, (1:12)+0.5, LowSem2$Average.weights, angle=90, code=3, length=.1)

z<-c(1,2,3)
x<-c("a","b","c")
zx<-data.frame(z,x)
zx$y<-c(25,125,250)

dens2 <- lake1[lake1$TadpoleDensity<4,]
dens10 <- lake1[lake1$TadpoleDensity >9 & lake1$TadpoleDensity <12,]
dens20 <- lake1[lake1$TadpoleDensity>18,]
#Subset mayfly densities of 0,25,125,250
dens2may0 <- dens2[dens2$MayflyDensity=="0",]
dens2may25 <- dens2[dens2$MayflyDensity=="25",]
dens2may125 <- dens2[dens2$MayflyDensity=="125",]
dens2may250 <- dens2[dens2$MayflyDensity=="250",]
dens10may0 <- dens10[dens10$MayflyDensity=="0",]
dens10may25 <- dens10[dens10$MayflyDensity=="25",]
dens10may125 <- dens10[dens10$MayflyDensity=="125",]
dens10may250 <- dens10[dens10$MayflyDensity=="250",]
dens20may0 <- dens20[dens20$MayflyDensity=="0",]
dens20may25 <- dens20[dens20$MayflyDensity=="25",]
dens20may125 <- dens20[dens20$MayflyDensity=="125",]
dens20may250 <- dens20[dens20$MayflyDensity=="250",]

d2m0avg <- mean(dens2may0$EndAvg)
d2m25avg <- mean(dens2may25$EndAvg)
d2m125avg <- mean(dens2may125$EndAvg)
d2m250avg <- mean(dens2may250$EndAvg)
d10m0avg <- mean(dens10may0$EndAvg)
d10m25avg <- mean(dens10may25$EndAvg)
d10m125avg <- mean(dens10may125$EndAvg)
d10m250avg <- mean(dens10may250$EndAvg)
d20m0avg <- mean(dens20may0$EndAvg)
d20m25avg <- mean(dens20may25$EndAvg)
d20m125avg <- mean(dens20may125$EndAvg)
d20m250avg <- mean(dens20may250$EndAvg)

#Creating Mayfly/Tad density table and barplot
maytad.density <- matrix(c(d2m0avg,d2m25avg,d2m125avg,d2m250avg,d10m0avg,d10m25avg,d10m125avg,d10m250avg,d20m0avg,d20m25avg,d20m125avg,d20m250avg),nrow=3,byrow=T)
colnames(maytad.density) <- c("0","25","125","250")
rownames(maytad.density) <- c("2","10","20")
names(dimnames(maytad.density)) <- c("Tadpole density","Mayfly density")
as.data.frame(as.table(maytad.density))
barplot(t(maytad.density),beside=T, ylim = c(0,60), main = "Lake 1",col=c("white","grey80","grey40","grey9"),ylab="Average Weights",xlab="Tadpole Density")
legend("topright",c("0","25","125","250"),cex=.6,bty="n",fill=c("white","grey80","grey40","black"),title="Mayfly Density",horiz=T)
arrows (1.5, UpSem2$Average.weights[3], 1.5, LowSem2$Average.weights[3], angle=90, code=3, length=.1)
arrows (2.5, UpSem2$Average.weights[5], 2.5, LowSem2$Average.weights[5], angle=90, code=3, length=.1)
arrows (3.5, UpSem2$Average.weights[2], 3.5, LowSem2$Average.weights[2], angle=90, code=3, length=.1)
arrows (4.5, UpSem2$Average.weights[10], 4.5, d2m250avg,col="black", angle=90, code=3, length=.1)
arrows (4.5, LowSem2$Average.weights[10], 4.5, d2m250avg,col="white", angle=90, code=1, length=.1)
arrows (6.5, UpSem2$Average.weights[9], 6.5, LowSem2$Average.weights[9], angle=90, code=3, length=.1)
arrows (7.5, UpSem2$Average.weights[12], 7.5, LowSem2$Average.weights[12], angle=90, code=3, length=.1)
arrows (8.5, UpSem2$Average.weights[6], 8.5, LowSem2$Average.weights[6], angle=90, code=3, length=.1) 
arrows (9.5, UpSem2$Average.weights[11], 9.5, d10m250avg, col="black",angle=90, code=3, length=.1)
arrows (9.5, LowSem2$Average.weights[11], 9.5, d10m250avg, col="white",angle=90, code=1, length=.1)
arrows (11.5, UpSem2$Average.weights[7], 11.5, LowSem2$Average.weights[7], angle=90, code=3, length=.1)
arrows (12.5, UpSem2$Average.weights[4], 12.5, LowSem2$Average.weights[4], angle=90, code=3, length=.1)
arrows (13.5, UpSem2$Average.weights[1], 13.5, LowSem2$Average.weights[1], angle=90, code=3, length=.1)
arrows (14.5, UpSem2$Average.weights[8], 14.5, d20m250avg, col="black",angle=90, code=3, length=.1)
arrows (14.5, LowSem2$Average.weights[8], 14.5, d20m250avg, col="white",angle=90, code=1, length=.1)

abline(0,0)
lines(c(1,5), c(mean(dens2$EndAvg),mean(dens2$EndAvg)), col="red")
abline(mean(taddata$EndAvg),0, col="brown",lty=2, lwd=2)

d2average <- mean(d2m0avg,d2m25avg,d2m125avg,d2m250avg)
d10average <- mean(d10m0avg,d10m25avg,d10m125avg,d10m250avg)
d20average <- mean(d20m0avg,d20m25avg,d20m125avg,d20m250avg)
d2stndvalues <- mean(l1dens2stnddev0,l1dens2stnddev25,l1dens2stnddev125,l1dens2stnddev250)
dens2stnddev0 <- dens2may0$EndStndDev
dens2stnddev25 <- dens2may25$EndStndDev
dens2stnddev125 <- dens2may125$EndStndDev
dens2stnddev250 <- dens2may250$EndStndDev

9.61147718
23.61253603
15.58584598
12.39606982

lake1densavg <- matrix(c(d2average,d10average,d20average))
rownames(lake1densavg) <- c("Tad Density 2","Tad Density 10", "Tad Density 20")
colnames(lake1densavg) <- c("Average weights")
barplot(lake1densavg, beside = T, names = rownames(lake1densavg), ylab = "Average weights", main = "Lake 1 Average",ylim = c(0,50), xlim = c(0,3), width = 0.75, space = 0.25)

s <- tapply (lake1densavg,, sd) #Help
n <- tapply (lake1$EndStndDev,lake1$CageNumber, length)
sem10 <- s/sqrt(no)
sem20 <- sem[c(1,2,3,4,5,6,7,8,9,10,11,12)]
UpSem2 <-data.frame(data.frame(lake1.avg)+data.frame(sem2))
LowSem2<-data.frame(data.frame(lake1.avg)-data.frame(sem2))
arrows ((1:12)+0.5, UpSem2$Average.weights, (1:12)+0.5, LowSem2$Average.weights, angle=90, code=3, length=.1)

#Lake 2
taddata <- read.csv("\\\\BRIGGS-NAS1\\Thomas Smith\\Tadpole Data\\SummarizedTadpoleData.csv", header = T)
taddatacage <- taddata[taddata$CageNumber!="0",]
lake2 <- taddatacage[taddatacage$LakeID=="10475",]
cage1l2 <- lake2[lake2$CageNumber=="1",]
cage2l2 <- lake2[lake2$CageNumber=="2",]
cage3l2 <- lake2[lake2$CageNumber=="3",]
cage4l2 <- lake2[lake2$CageNumber=="4",]
cage5l2 <- lake2[lake2$CageNumber=="5",]
cage7l2 <- lake2[lake2$CageNumber=="7",]
cage10l2 <- lake2[lake2$CageNumber=="10",]
cage11l2 <- lake2[lake2$CageNumber=="11",]
cage13l2 <- lake2[lake2$CageNumber=="13",]
cage15l2 <- lake2[lake2$CageNumber=="15",]
cage16l2 <- lake2[lake2$CageNumber=="16",]
cage17l2 <- lake2[lake2$CageNumber=="17",]

avgl2c1 <- mean(cage1l2$EndAvg)
avgl2c2 <- mean(cage2l2$EndAvg)
avgl2c3 <- mean(cage3l2$EndAvg)
avgl2c4 <- mean(cage4l2$EndAvg)
avgl2c5 <- mean(cage5l2$EndAvg)
avgl2c7 <- mean(cage7l2$EndAvg)
avgl2c10 <- mean(cage10l2$EndAvg)
avgl2c11 <- mean(cage11l2$EndAvg)
avgl2c13 <- mean(cage13l2$EndAvg)
avgl2c15 <- mean(cage15l2$EndAvg)
avgl2c16 <- mean(cage16l2$EndAvg)
avgl2c17 <- mean(cage17l2$EndAvg)

lake2.avg <- matrix(c(avgl2c1,avgl2c2,avgl2c3,avgl2c4,avgl2c5,avgl2c7,avgl2c10,avgl2c11,avgl2c13,avgl2c15,avgl2c16,avgl1c17))
rownames (lake2.avg) <- c("Cage1","Cage2","Cage3","Cage4","Cage5","Cage7","Cage10","Cage11","Cage13","Cage15","Cage16","Cage17")
colnames (lake2.avg) <- c("Average.weights")
barplot(lake2.avg, beside = T, ylim = c(0,50), names = rownames(lake2.avg))

dens2l2 <- lake2[lake2$TadpoleDensity<4,]
dens10l2 <- lake2[lake2$TadpoleDensity >7 & lake2$TadpoleDensity <11,]
dens20l2 <- lake2[lake2$TadpoleDensity>16,]

dens2may0l2 <- dens2l2[dens2l2$MayflyDensity=="0",]
dens2may25l2 <- dens2l2[dens2l2$MayflyDensity=="25",]
dens2may125l2 <- dens2l2[dens2l2$MayflyDensity=="125",]
dens2may250l2 <- dens2l2[dens2l2$MayflyDensity=="250",]
dens10may0l2 <- dens10l2[dens10l2$MayflyDensity=="0",]
dens10may25l2 <- dens10l2[dens10l2$MayflyDensity=="25",]
dens10may125l2 <- dens10l2[dens10l2$MayflyDensity=="125",]
dens10may250l2 <- dens10l2[dens10l2$MayflyDensity=="250",]
dens20may0l2 <- dens20l2[dens20l2$MayflyDensity=="0",]
dens20may25l2 <- dens20l2[dens20l2$MayflyDensity=="25",]
dens20may125l2 <- dens20l2[dens20l2$MayflyDensity=="125",]
dens20may250l2 <- dens20l2[dens20l2$MayflyDensity=="250",]

d2m0avgl2 <- mean(dens2may0l2$EndAvg)
d2m25avgl2 <- mean(dens2may25l2$EndAvg)
d2m125avgl2 <- mean(dens2may125l2$EndAvg)
d2m250avgl2 <- mean(dens2may250l2$EndAvg)
d10m0avgl2 <- mean(dens10may0l2$EndAvg)
d10m25avgl2 <- mean(dens10may25l2$EndAvg)
d10m125avgl2 <- mean(dens10may125l2$EndAvg)
d10m250avgl2 <- mean(dens10may250l2$EndAvg)
d20m0avgl2 <- mean(dens20may0l2$EndAvg)
d20m25avgl2 <- mean(dens20may25l2$EndAvg)
d20m125avgl2 <- mean(dens20may125l2$EndAvg)
d20m250avgl2 <- mean(dens20may250l2$EndAvg)

#Lake 2 Barplot

maytad.density.l2 <- matrix(c(d2m0avgl2,d2m25avgl2,d2m125avgl2,d2m250avgl2,d10m0avgl2,d10m25avgl2,d10m125avgl2,d10m250avgl2,d20m0avgl2,d20m25avgl2,d20m125avgl2,d20m250avgl2),nrow=3,byrow=T)
colnames(maytad.density.l2) <- c("0","25","125","250")
rownames(maytad.density.l2) <- c("2","10","20")
names(dimnames(maytad.density.l2)) <- c("Tadpole density","Mayfly density")
as.data.frame(as.table(maytad.density.l2))
barplot(t(maytad.density.l2),beside=T, ylim = c(0,60), main = "Lake 2",col=c("white","grey80","grey40","grey9"),ylab="Average Weights",xlab="Tadpole Density")
legend("topright",c("0","25","125","250"),cex=.6,bty="n",fill=c("white","grey80","grey40","black"),title="Mayfly Density",horiz=T)

lake2.avg.Vec <- c(avgl2c10, avgl2c15,avgl2c3, avgl2c4,avgl2c13,avgl2c7,avgl2c11,avgl2c16,avgl1c17, avgl2c2,avgl2c5, avgl2c1)

barplot(lake2.avg.Vec,beside=T, space=c(0,0,0,0,1,0,0,0,1,0,0,0), ylim = c(0,60), main = "Lake 2",col=c("white","grey80","grey40","grey9"),ylab="Average Weights",xlab="Tadpole Density")
legend("topright",c("0","25","125","250"),cex=.6,bty="n",fill=c("white","grey80","grey40","black"),title="Mayfly Density",horiz=T)


q <- tapply (lake2$EndStndDev,lake2$CageNumber, sd)
m <- tapply (lake2$EndStndDev,lake2$CageNumber, length)
sem3 <- q/sqrt(m)
sem4 <- sem3[c(2,3,4,5,6,8,12,14,16,17,18)]
UpSem4 <-data.frame(data.frame(lake2.avg[-7]) + data.frame(sem4))
LowSem4 <-data.frame(data.frame(lake2.avg[-7])-data.frame(sem4))

install.packages("reshape")
library(reshape)
UpSem4 <- rename(UpSem4, c(lake2.avg..7.="Average.weights"))
LowSem4 <- rename(LowSem4, c(lake2.avg..7.="Average.weights"))

arrows (2.5, UpSem4$Average.weights[9], 2.5, LowSem4$Average.weights[9], angle=90, code=3, length=.1)
arrows (3.5, UpSem4$Average.weights[3], 3.5, LowSem4$Average.weights[3], angle=90, code=3, length=.1)
arrows (4.5, UpSem4$Average.weights[4], 4.5, d2m250avgl2,col="black", angle=90, code=3, length=.1)
arrows (4.5, LowSem4$Average.weights[4], 4.5, d2m250avgl2,col="white", angle=90, code=1, length=.1)
arrows (6.5, UpSem4$Average.weights[8], 6.5, LowSem4$Average.weights[8], angle=90, code=3, length=.1)
arrows (7.5, UpSem4$Average.weights[6], 7.5, LowSem4$Average.weights[6], angle=90, code=3, length=.1)
# Problem: arrows (8.5, UpSem4$Average.weights[7], 8.5, LowSem4$Average.weights[7], angle=90, code=3, length=.1) 
arrows (9.5, UpSem4$Average.weights[10], 9.5, d10m250avgl2, col="black",angle=90, code=3, length=.1)
arrows (9.5, LowSem4$Average.weights[10], 9.5, d10m250avgl2, col="white",angle=90, code=1, length=.1)
arrows (11.5, UpSem4$Average.weights[11], 11.5, LowSem4$Average.weights[11], angle=90, code=3, length=.1)
arrows (12.5, UpSem4$Average.weights[2], 12.5, LowSem4$Average.weights[2], angle=90, code=3, length=.1)
arrows (13.5, UpSem4$Average.weights[5], 13.5, LowSem4$Average.weights[5], angle=90, code=3, length=.1)
arrows (14.5, UpSem4$Average.weights[1], 14.5, d20m250avgl2, col="black",angle=90, code=3, length=.1)
arrows (14.5, LowSem4$Average.weights[1], 14.5, d20m250avgl2, col="white",angle=90, code=1, length=.1)

d2averagel2 <- mean(d2m0avgl2,d2m25avgl2,d2m125avgl2,d2m250avgl2)
d10averagel2 <- mean(d10m0avgl2,d10m25avgl2,d10m125avgl2,d10m250avgl2)
d20averagel2 <- mean(d20m0avgl2,d20m25avgl2,d20m125avgl2,d20m250avgl2)

lake2densavg <- matrix(c(d2averagel2,d10averagel2,d20averagel2))
apply(lake2densavg, 1, sd)
rownames(lake2densavg) <- c("Tad Density 2","Tad Density 10", "Tad Density 20")
colnames(lake2densavg) <- c("Average weights")
barplot(lake2densavg, beside = T, names = rownames(lake2densavg), ylab = "Average weights", main = "Lake 2 Average", xlim = c(0,3), width = 0.75, space = 0.25, ylim = c(0,40))


#Summary Statistics
summary(maytad.density)
shapiro.test(lake1$EndAvg)
result <- shapiro.test(lake1$EndAvg)
result$p.value
qqline(lake1$EndAvg)
qqnorm(lake1$EndAvg)
t.test (lake1$EndAvg)
t.test(lake1$EndAvg,lake2$EndAvg)

kruskal.test(EndAvg ~ SampleDate, data = taddatacage)
#nonidentical populations
kruskal.test(TadpoleDensity ~ MayflyDensity, data = lake1)
oneway.test (EndAvg ~ SampleDate, data = lake1)
bartlett.test (TadpoleDensity ~ SampleDate, data = lake1)
summary(lm(EndAvg ~ SampleDate, data = lake1))
anova(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity, data = endweights))
#looks at whether the variation in the groups occurs due to random error, tests 3 or more variables for significance
#f-value less than 1 indicates differences among group means were actually smaller than expected by chance,
#effect of interest is very large
anova(lm(EndAvg ~ SampleDate, data = lake2))
#f-value greater than 1 indicates null hypothesis is false and not statistically significant, likely to have occurred
#by chance
anova(lm(dens2may0$EndAvg ~ dens2may250$SampleDate, data = lake1))
#Good measurement?



#ANOVA Stats
endweights <- read.csv("\\\\BRIGGS-NAS1\\Thomas Smith\\Tadpole Data\\Endweights.csv", header = T)
anova(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity, data = endweights))
anova(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity + as.factor(LakeID), data = endweights))
summary(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity, data = endweights))

summary(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity + as.factor(LakeID), data = endweights))
summary(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity + as.factor(SampleDate), data = endweights))
one <-lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity + as.factor(SampleDate), data = endweights)
two <- lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity, data = endweights)
anova(one,two)

summary(lm(Predicted.AFDM ~ Tadpole.Density * MayflyDensity + as.factor(LakeID), data = endweights))



str(lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity + LakeID, data = endweights))

lmTadMF<-lm(Predicted.AFDM ~ Tadpole.Density + MayflyDensity + as.factor(LakeID), data = endweights)

lmTadMFl1 <- lm(Predicted.AFDM ~ Tadpole.Density * MayflyDensity, data = endweights[endweights$LakeID=="10102",])
lmTadMFl2 <- lm(Predicted.AFDM ~ Tadpole.Density * MayflyDensity, data = endweights[endweights$LakeID=="10475",])

# Important for written analysis!
summary(lmTadMFl1)
summary(lmTadMFl2)
summary(lm(Predicted.AFDM ~ Tadpole.Density * MayflyDensity + as.factor(LakeID), data = endweights))

#Graphing Residuals
plot(resid(lmTadMF)~fitted(lmTadMF))
hist(resid(lmTadMF))

boxplot(resid(lmTadMFl1)~endweights$Tadpole.Density[endweights$LakeID=="10102"])
boxplot(resid(lmTadMFl2)~endweights$Tadpole.Density[endweights$LakeID=="10475"])

boxplot(resid(lmTadMF)~endweights$MayflyDensity)

boxplot(resid(lmTadMF)~factor(endweights$SampleDate,levels = c("15-Jul-09","1-Aug-09","23-Aug-09","14-Sep-09","21-Jul-09","12-Aug-09","31-Aug-09","21-Sep-09")))
boxplot(resid(lmTadMF)~endweights$LakeID)


plot(resid()~fitted(lmTadMFl1))

shapiro.test(lmTadMFl1)
