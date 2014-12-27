library(ggplot2)
#load TadBL18Aug: "...//Consumer Resource Experiment//Summer2010//TadBL18Aug.csv"
TadSize<-read.csv(file.choose(),header=T)
str(TadSize)
plot(TadSize)
summary(aov(BL~BdExposure, data=TadSize))
ggplot(TadSize, aes(x=BdExposure, y=BL))+
	geom_boxplot(col="black", aes(fill=BdExposure))
#no effect of exposure on Tadpole Size
