#bootstrap of pseudocontrolled algal abundance
#as in John Fox 2000
#Bootstrapping Regression Models
#Appendix to An R and S-PLUS Companion to Applied Regression

#random-X Resampling rather than Fixed-X resampling
#despite having fixed predictor values in the experimental design, Fox suggests using the random-X resampling method 
#b/c it does not make an assumption of identical distribution of errors, which would be violated by non-constatn error varince
#which we have...

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
#############
#calculations:
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

#############
#basic plotting
#groups:
#############
library(boot)
#statistic function for boot:
boot.algae<-function(data,indices,maxit){
	data<-data[indices,]	#select obs. in bootstrap sample
	mod<-lm(ExperimentalControlDifference_Abundance~TadDens+MFDens+LakeName+SampleNumber+SamplePeriod+Silt+Radiation, data=data, maxit=maxit)
	coefficients(mod)	#return coefficient vector
	}
booty.boot<-boot(whole, boot.algae, 1999, maxit=100)
boot.ci(booty.boot, index=1, type=c("norm", "perc", "bca"))	#TadDens
boot.ci(booty.boot, index=2, type=c("norm", "perc", "bca"))	#MFDens
boot.ci(booty.boot, index=3, type=c("norm", "perc", "bca"))	#LakeName
boot.ci(booty.boot, index=4, type=c("norm", "perc", "bca"))	#SampleNumber
boot.ci(booty.boot, index=5, type=c("norm", "perc", "bca"))	#SamplePeriod
boot.ci(booty.boot, index=6, type=c("norm", "perc", "bca"))	#Silt
boot.ci(booty.boot, index=7, type=c("norm", "perc", "bca"))	#Radiation




