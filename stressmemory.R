## USAGE: "Thermal Priming Costs Outweigh Benefits in the Staghorn Coral Acropora cervicornis (Lamarck 1816)

# Stress Memory Experiment 
# loads in SM data & eventually, all Fv/Fm data (all timepoints)
# creates boxplots for symbiont density, total chl, chl per cell, total algal protein and protein per cell

# load libraries
library(rcompanion)
library(ggplot2)
library(easyGgplot2)
library(devtools) 
library(lmmfit)
library(lme4)
library(labdsv)
library(vegan)
library(plotrix)
library(pgirmess)
library(gridExtra)
library(pbkrtest)
library(RVAideMemoire)
library(car)
library(gtools)
library(plyr)
library(quantreg)
library(calibrate)
library(MASS)
library(e1071)
library(FSA)
library(DescTools)
library(Hmisc)
library(dplyr)
library(PMCMRplus)
library(multcomp)
library(AICcmodavg)
library(nlme)
library(exactRankTests)
library(MCMCglmm)
library(psych)
source("~/Desktop/ODU/BARSHIS_LAB_ALL/Chapter1_Recent Thermal History/Data/RsquaredGLMM.R")

## set the working directory for all analyses
setwd("~/Desktop/Dissertation/Chapter3/analyses/July2019/")

####################################
## load the chapter 3 data
	data=read.csv("~/Desktop/Dissertation/Chapter3/analyses/July2019/SM2018_allSamples.csv", header=TRUE, stringsAsFactors=TRUE)
# with Chl per cell outliers gone!
	datanoout=read.csv("~/Desktop/Dissertation/Chapter3/analyses/June2019/SM2018_allSamples_outliersRemoved.csv", header=TRUE, stringsAsFactors=TRUE)

## reorder factors to appear as desired

data$genet=factor(data$genet,levels=c("A","B","C","D","E","F","G","H","I","J"))
data$trt=factor(data$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
data$tp=factor(data$tp, levels=c("1","2","3"), labels=c("1","2","3"))

datanoout$genet=factor(datanoout$genet,levels=c("A","B","C","D","E","F","G","H","I","J"))
datanoout$trt=factor(datanoout$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
datanoout$tp=factor(datanoout$tp, levels=c("1","2","3"), labels=c("1","2","3"))

## examine the complete dataset
names(data) # what variables are there
str(data) # look at the data structure
headTail(data) # look at the first few rows
#tail(data) # look at the last few rows
summary(data) # look for possible missing values or unbalanced design
# no NAs!!!

names(datanoout) # what variables are there
str(datanoout) # look at the data structure
headTail(datanoout) # look at the first few rows
summary(datanoout) # look for possible missing values or unbalanced design

library(plyr)

#chlpcell
chlpcellStats <- ddply(datanoout, c("tp","trt"), summarise,
			n = length(chlpcell),
			min = min(chlpcell),
			max = max(chlpcell),
			mean = mean(chlpcell),
			median = median(chlpcell),
			iqr = IQR(chlpcell, na.rm=TRUE),
			sd = sd(chlpcell),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
chlpcellStats

protpcellStats <- ddply(datanoout, c("tp","trt"), summarise,
			n = length(protpcell),
			min = min(protpcell),
			max = max(protpcell),
			mean = mean(protpcell),
			median = median(protpcell),
			iqr = IQR(protpcell, na.rm=TRUE),
			sd = sd(protpcell),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
protpcellStats

# This is to identify the outliers via Mahalanobis distance
 library(stats)
 library(reshape)
 library(psych)

data_wide<-reshape(data, idvar = "trt_genet", timevar = "tp", direction = "wide")
headTail(data_wide)
chlpcell<- select(data_wide, chlpcell.1, chlpcell.2, chlpcell.3); dim(chlpcell); summary(chlpcell)
mahal_chlpcell = mahalanobis(chlpcell, colMeans(chlpcell, na.rm=TRUE), cov(chlpcell,use="pairwise.complete.obs"));
mahal_chlpcell
## determine a cutoff score
cutoff=qchisq(1-.001, ncol(chlpcell))
cutoff
ncol(chlpcell) # df
summary(mahal_chlpcell < cutoff)  ## 2 outliers
noout_chlpcell = subset(chlpcell, mahal_chlpcell < cutoff) # remove samples with mahal values < the cutoff (> are outlying)
noout_chlpcell ## outliers #45, 46 were removed from all timepoints

#chl
chlStats <- ddply(data, c("tp","trt"), summarise,
			n = length(chl),
			min = min(chl),
			max = max(chl),
			mean = mean(chl),
			median = median(chl),
			iqr = IQR(chl, na.rm=TRUE),
			sd = sd(chl),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
chlStats

#cells
cellStats <- ddply(data, c("tp","trt"), summarise,
			n = length(cells),
			min = min(cells),
			max = max(cells),
			mean = mean(cells),
			median = median(cells),
			iqr = IQR(cells, na.rm=TRUE),
			sd = sd(cells),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
cellStats

protStats <- ddply(data, c("tp","trt"), summarise,
			n = length(prot),
			min = min(prot),
			max = max(prot),
			mean = mean(prot),
			median = median(prot),
			iqr = IQR(prot, na.rm=TRUE),
			sd = sd(prot),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
protStats

fvfmStats <- ddply(data, c("tp","trt"), summarise,
			n = length(fvfm),
			min = min(fvfm),
			max = max(fvfm),
			mean = mean(fvfm),
			median = median(fvfm),
			iqr = IQR(fvfm, na.rm=TRUE),
			sd = sd(fvfm),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
fvfmStats

## parse timepoints without outliers first

time1noout<-datanoout[datanoout$tp=="1",]; dim(time1noout)
time2noout<-datanoout[datanoout$tp=="2",]; dim(time2noout)
time3noout<-datanoout[datanoout$tp=="3",]; dim(time3noout)

## Chlorophyll per Cell
# time 1
#pdf(file="chlpcell_v_trt_boxplot_time1.pdf")
		black.bold.text<-element_text(family="Times New Roman",face="bold",color="black", size=24)
		black.italic.text<-element_text(family="Times New Roman",face="bold.italic",color="black", size=28)
		italic.text<-element_text(family="Times New Roman",face="italic",color="black", size=24)
		plain.text<-element_text(family="Times New Roman",face="plain",color="black", size=24)
	
	p<-ggplot(time1noout, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,30) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Chlorophyll per Cell (pg Chl cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()
# time 2
#pdf(file="chlpcell_v_trt_boxplot_time2.pdf")
	p<-ggplot(time2noout, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,30) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Chlorophyll per Cell (pg Chl cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()
# time 3
#pdf(file="chlpcell_v_trt_boxplot_time3.pdf")
	p<-ggplot(time3noout, aes(x=trt, y=chlpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,30) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Chlorophyll per Cell (pg Chl cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

### Protein per cell
# time 1
#pdf(file="protpcell_v_trt_boxplot_time1.pdf")
	p<-ggplot(time1noout, aes(x=trt, y=protpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.5) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Protein per Cell (ng cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()
# time 2
#pdf(file="protpcell_v_trt_boxplot_time2.pdf")
	p<-ggplot(time2noout, aes(x=trt, y=protpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.5) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Protein per Cell (ng cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()
# time 3
#pdf(file="protpcell_v_trt_boxplot_time3.pdf")
	p<-ggplot(time3noout, aes(x=trt, y=protpcell)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.5) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Protein per Cell (ng cell'^-1*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

### Now parse data with all the values for cells and chl

time1<-data[data$tp=="1",]; dim(time1)
time2<-data[data$tp=="2",]; dim(time2)
time3<-data[data$tp=="3",]; dim(time3)

### Chlorophyll

#pdf(file="totalchl_v_trt_boxplot_time1.pdf")
		black.bold.text<-element_text(family="Times New Roman",face="bold",color="black", size=24)
		black.italic.text<-element_text(family="Times New Roman",face="bold.italic",color="black", size=28)
		italic.text<-element_text(family="Times New Roman",face="italic",color="black", size=24)
		plain.text<-element_text(family="Times New Roman",face="plain",color="black", size=24)
	p<-ggplot(time1, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		#ylim(0,8) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Total Chlorophyll (' ~mu *'g Chl cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="totalchl_v_trt_boxplot_time2.pdf")
	p<-ggplot(time2, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		#ylim(0,8) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Total Chlorophyll (' ~mu *'g Chl cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="totalchl_v_trt_boxplot_time3.pdf")
	p<-ggplot(time3, aes(x=trt, y=chl)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		#ylim(0,8) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Total Chlorophyll (' ~ mu *'g Chl cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

### PROTEIN
#pdf(file="totalprot_v_trt_boxplot_time1.pdf")
	p<-ggplot(time1, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Algal Protein (' ~mu *'g cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="totalchl_v_trt_boxplot_time2.pdf")
	p<-ggplot(time2, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Algal Protein (' ~mu *'g cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="totalprot_v_trt_boxplot_time3.pdf")
	p<-ggplot(time3, aes(x=trt, y=prot)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3, shape=21) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Algal Protein (' ~ mu *'g cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

### CELLS
#pdf(file="cells_v_trt_boxplot_time1.pdf")
	p<-ggplot(time1, aes(x=trt, y=cells)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="cells_v_trt_boxplot_time2.pdf")
	p<-ggplot(time2, aes(x=trt, y= cells)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(alpha=0.3, size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()


### Time3 cell stats
time3stats <- ddply(time3, "trt", summarise,
			n = length(cells),
			min = min(cells),
			max = max(cells),
			mean = mean(cells),
			median = median(cells),
			iqr = IQR(cells, na.rm=TRUE),
			sd = sd(cells),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
time3stats

#pdf(file="cells_v_trt_barplot_time3.pdf")
	p<-ggplot(time3stats, aes(x=trt, y= median)) + 
		#geom_boxplot(width=.5, lwd=1) +
		geom_bar(stat="identity") +
		geom_errorbar(aes(ymin=median-se,ymax=median+se), width=.2) +
		#geom_point(alpha=0.3, size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="cells_v_trt_boxplot_time3.pdf")
	p<-ggplot(time3, aes(x=trt, y= cells)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(alpha=0.3, size=3) +
		ylim(0,5050000) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bold(paste('Symbiont Density (cells cm'^-2*") ")))) 
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

### FV/FM
#pdf(file="FvFm_v_trt_boxplot_time1.pdf")
	p<-ggplot(time1, aes(x=trt, y=fvfm)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bolditalic(paste('Photochemical Efficiency (F'[V]*'/F'[M]*") " )))  )
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="FvFm_v_trt_boxplot_time2.pdf")
	p<-ggplot(time2, aes(x=trt, y=fvfm)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bolditalic(paste('Photochemical Efficiency (F'[V]*'/F'[M]*") " )))  )
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

#pdf(file="FvFm_v_trt_boxplot_time3.pdf")
	p<-ggplot(time3, aes(x=trt, y=fvfm)) + 
		geom_boxplot(width=.5, lwd=1) +
		#geom_point(size=3) +
		ylim(0,0.7) +
		theme_bw() +
		theme(legend.position="none")
	p.labs <- p + labs(x="Treatment", y=expression(bolditalic(paste('Photochemical Efficiency (F'[V]*'/F'[M]*") " )))  )
	plot(p.labs + theme(title= black.bold.text, axis.title=black.bold.text, axis.text=black.bold.text))
#dev.off()

### Assumption Checking

#parse each response variable into trts

#Chl per cell & Prot per cell
	c <-time1noout[time1noout$trt=="C", ]; dim(c)
	n <-time1noout[time1noout$trt=="N", ]; dim(n)
	ll <-time1noout[time1noout$trt=="LL", ]; dim(ll)
	lh <-time1noout[time1noout$trt=="LL", ]; dim(ll)
	hl <-time1noout[time1noout$trt=="HL", ]; dim(hl)
	hh <-time1noout[time1noout$trt=="HH", ]; dim(hh)

# chl p cell
par(mfrow=c(2,3))
	qqnorm(c$chlpcell, main="control")
	qqline(c$chlpcell)
	qqnorm(n$chlpcell, main="naive")
	qqline(n$chlpcell)
	qqnorm(ll$chlpcell, main="LL")
	qqline(ll$chlpcell)
	qqnorm(lh$chlpcell, main="LH")
	qqline(lh$chlpcell)
	qqnorm(hl$chlpcell, main="HL")
	qqline(hl$chlpcell)
	qqnorm(hh$chlpcell, main="HH")
	qqline(hh$chlpcell)
## Chl per cell data appear to be normal in TP1, use mean

levene<-leveneTest(chlpcell ~ interaction(trt,tp), data=time1noout); levene
			if(levene[1,3] <= 0.05) {
				cat("No reason to reject null hypothesis, p < 0.05
				These data appear to be of equal variances")
			} else {
				cat("These data are not of equal variances, p > 0.05")
			}
## Data do not meet assumptions! Transform or use KW

# prot p cell
par(mfrow=c(2,3))
	qqnorm(c$protpcell, main="control")
	qqline(c$protpcell)
	qqnorm(n$protpcell, main="naive")
	qqline(n$protpcell)
	qqnorm(ll$protpcell, main="LL")
	qqline(ll$protpcell)
	qqnorm(lh$protpcell, main="LH")
	qqline(lh$protpcell)
	qqnorm(hl$protpcell, main="HL")
	qqline(hl$protpcell)
	qqnorm(hh$protpcell, main="HH")
	qqline(hh$protpcell)
## Prot per cell data do not appear to be normal in TP1, transform or use KW test
levene<-leveneTest(protpcell ~ interaction(trt,tp), data=time1noout); levene
			if(levene[1,3] <= 0.05) {
				cat("No reason to reject null hypothesis, p < 0.05
				These data appear to be of equal variances")
			} else {
				cat("These data are not of equal variances, p > 0.05")
			}
## Data do not meet assumptions, transform or use KW

## Parse for other response vars
	c <-time1[time1$trt=="C", ]; dim(c)
	n <-time1[time1$trt=="N", ]; dim(n)
	ll <-time1[time1$trt=="LL", ]; dim(ll)
	lh <-time1[time1$trt=="LL", ]; dim(ll)
	hl <-time1[time1$trt=="HL", ]; dim(hl)
	hh <-time1[time1$trt=="HH", ]; dim(hh)

# cells
par(mfrow=c(2,3))
	qqnorm(c$cells, main="control")
	qqline(c$cells)
	qqnorm(n$cells, main="naive")
	qqline(n$cells)
	qqnorm(ll$cells, main="LL")
	qqline(ll$cells)
	qqnorm(lh$cells, main="LH")
	qqline(lh$cells)
	qqnorm(hl$cells, main="HL")
	qqline(hl$cells)
	qqnorm(hh$cells, main="HH")
	qqline(hh$cells)
##  cells data do not appear to be normal in TP1, transform or use KW test
levene<-leveneTest(cells ~ interaction(trt,tp), data=time1); levene
			if(levene[1,3] <= 0.05) {
				cat("No reason to reject null hypothesis, p < 0.05
				These data appear to be of equal variances")
			} else {
				cat("These data are not of equal variances, p > 0.05")
			}
## Data do not meet assumptions, transform or use KW

# chl
par(mfrow=c(2,3))
	qqnorm(c$chl, main="control")
	qqline(c$chl)
	qqnorm(n$chl, main="naive")
	qqline(n$chl)
	qqnorm(ll$chl, main="LL")
	qqline(ll$chl)
	qqnorm(lh$chl, main="LH")
	qqline(lh$chl)
	qqnorm(hl$chl, main="HL")
	qqline(hl$chl)
	qqnorm(hh$chl, main="HH")
	qqline(hh$chl)
## chl data do not appear to be normal in TP1, transform or use KW test
levene<-leveneTest(chl ~ interaction(trt,tp), data=time1); levene
			if(levene[1,3] <= 0.05) {
				cat("No reason to reject null hypothesis, p < 0.05
				These data appear to be of equal variances")
			} else {
				cat("These data are not of equal variances, p > 0.05")
			}
## Data do not meet assumptions, transform or use KW

# Protein
par(mfrow=c(2,3))
	qqnorm(c$prot, main="control")
	qqline(c$prot)
	qqnorm(n$prot, main="naive")
	qqline(n$prot)
	qqnorm(ll$prot, main="LL")
	qqline(ll$prot)
	qqnorm(lh$prot, main="LH")
	qqline(lh$prot)
	qqnorm(hl$prot, main="HL")
	qqline(hl$prot)
	qqnorm(hh$prot, main="HH")
	qqline(hh$prot)
## prot data do not appear to be normal in TP1, transform or use KW test
levene<-leveneTest(prot ~ interaction(trt,tp), data=time1); levene
			if(levene[1,3] <= 0.05) {
				cat("No reason to reject null hypothesis, p < 0.05
				These data appear to be of equal variances")
			} else {
				cat("These data are not of equal variances, p > 0.05")
			}
## Data do not meet assumptions, transform or use KW

# FvFm
par(mfrow=c(2,3))
	qqnorm(c$fvfm, main="control")
	qqline(c$fvfm)
	qqnorm(n$fvfm, main="naive")
	qqline(n$fvfm)
	qqnorm(ll$fvfm, main="LL")
	qqline(ll$fvfm)
	qqnorm(lh$fvfm, main="LH")
	qqline(lh$fvfm)
	qqnorm(hl$fvfm, main="HL")
	qqline(hl$fvfm)
	qqnorm(hh$fvfm, main="HH")
	qqline(hh$fvfm)
## fvfm data do not appear to be normal in TP1, transform or use KW test
levene<-leveneTest(fvfm ~ interaction(trt,tp), data=time1); levene
			if(levene[1,3] <= 0.05) {
				cat("No reason to reject null hypothesis, p < 0.05
				These data appear to be of equal variances")
			} else {
				cat("These data are not of equal variances, p > 0.05")
			}
## Data do not meet assumptions, transform or use KW

### KW tests
library(FSA)
library(pgirmess)
library(multcompView)

### chl per cell
		#inter<-interaction(time1$trt,time1$tp) # specify interaction term

	k<-kruskal.test(chlpcell~interaction(trt,tp), datanoout) 
			if(kruskal.test(chlpcell ~ interaction(trt,tp), datanoout)$p.value <= 0.05) {
						cat("There is an effect of trt*tp on chl per cell, p < 0.05")
						tuk<-kruskalmc(chlpcell ~ interaction(trt,tp), datanoout) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsChlpcell_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsChlpcell_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsChlpcell_TrtxTime.txt"))
				} else {
						cat("There are no significant differences in chl per cell across trt or tp, p > 0.05")
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsChlpcell_TrtxTime.txt"))
				}
				
	# # # if(kruskal.test(chlpcell ~ trt, time1noout)$p.value <= 0.05) {
				# # # dunit<-DunnettTest(time1noout$chlpcell, time1noout$trt, control="C", conf.level=0.95)
						# # # print(k)
						# # # print(dunit)
						# # # cat(capture.output(print(k), file="kruskalResultsChlpcell_x_Trt_Time1_postPriming.txt"))
						# # # cat(capture.output(print(dunit), file="dunnettPosthocTest_Chlpcell_x_Trt_Time1_postPriming.txt"))
					# # # plot(time1noout$chlpcell ~ time1noout$trt)
				# # # } else {
						# # # cat("No posthoc test required")
						# # # cat(capture.output(print(k), file="kruskalResultsChlpcell_x_Trt_Time1_postPriming.txt"))
				# # # }

### Protein per Cell	
	k<-kruskal.test(protpcell~interaction(trt,tp), datanoout) 
			if(kruskal.test(protpcell ~ interaction(trt,tp), datanoout)$p.value <= 0.05) {
						cat("There is an effect of trt*tp on prot per cell, p < 0.05")
						tuk<-kruskalmc(protpcell ~ interaction(trt,tp), datanoout) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsProtpcell_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsProtpcell_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsProtpcell_TrtxTime.txt"))
				} else {
						cat("There are no significant differences in prot per cell across trt or tp, p > 0.05")
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsProtpcell_TrtxTime.txt"))
				}
				

	# # # ### time 1 after priming
	# # # k<-kruskal.test(protpcell~trt, time1noout) 
			# # # if(kruskal.test(protpcell ~ trt, time1noout)$p.value <= 0.05) {
						# # # cat("There is an effect of trt on prot per cell, p < 0.05")
				# # # } else {
						# # # cat("There are no significant differences in prot per cell across trt, p > 0.05")
				# # # }
				
	# # # if(kruskal.test(protpcell ~ trt, time1noout)$p.value <= 0.05) {
				# # # dunit<-DunnettTest(time1noout$protpcell, time1noout$trt, control="C", conf.level=0.95)
						# # # print(k)
						# # # print(dunit)
						# # # cat(capture.output(print(k), file="kruskalResultsProtpcell_x_Trt_Time1_postPriming.txt"))
						# # # cat(capture.output(print(dunit), file="dunnettPosthocTest_Protpcell_x_Trt_Time1_postPriming.txt"))
					# # # plot(time1noout$protpcell ~ time1noout$trt)
				# # # } else {
						# # # cat("No posthoc test required")
						# # # cat(capture.output(print(k), file="kruskalResultsProtpcell_x_Trt_Time1_postPriming.txt"))
				# # # }

### total chl
	k<-kruskal.test(chl~interaction(trt,tp), data) 
			if(kruskal.test(chl ~ interaction(trt,tp), data)$p.value <= 0.05) {
						cat("There is an effect of trt*tp on chl , p < 0.05")
						tuk<-kruskalmc(chl ~ interaction(trt,tp), data) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsChl_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsChl_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsChl_TrtxTime.txt"))
				} else {
						cat("There are no significant differences in chl across trt or tp, p > 0.05")
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsChl_TrtxTime.txt"))
				}


	### time 1 after priming
	# # # k<-kruskal.test(chl ~ trt, time1) 
			# # # if(kruskal.test(chl ~ trt, time1)$p.value <= 0.05) {
						# # # cat("There is an effect of trt on chl, p < 0.05")
				# # # } else {
						# # # cat("There are no significant differences in chl across trt, p > 0.05")
				# # # }
	# # # if(kruskal.test(chl ~ trt, time1)$p.value <= 0.05) {
						# # # dunit<-DunnettTest(time1$chl, time1$trt, control="C",conf.level=0.95)
						# # # print(k)
						# # # print(dunit)
						# # # cat(capture.output(print(k), file="kruskalResultsTotalChl_x_Trt_Time1_postPriming.txt"))
						# # # cat(capture.output(print(dunit), file="dunnettPosthocTest_TotalChl_x_Trt_Time1_postPriming.txt"))
					# # # plot(time1$chl ~ time1$trt)
				# # # } else {
						# # # cat("No posthoc test required")
						# # # cat(capture.output(print(k), file="kruskalResultsTotalChl_x_Trt_Time1_postPriming.txt"))

				# # # }

### total prot

	k<-kruskal.test(prot~interaction(trt,tp), data) 
			if(kruskal.test(prot ~ interaction(trt,tp), data)$p.value <= 0.05) {
						cat("There is an effect of trt*tp on prot , p < 0.05")
						tuk<-kruskalmc(prot ~ interaction(trt,tp), data) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsProt_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsProt_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsProt_TrtxTime.txt"))
				} else {
						cat("There are no significant differences in prot across trt or tp, p > 0.05")
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsProt_TrtxTime.txt"))
				}

	### time 1 after priming
	# # # k<-kruskal.test(prot ~ trt, time1) 
			# # # if(kruskal.test(prot ~ trt, time1)$p.value <= 0.05) {
						# # # cat("There is an effect of prot on chl, p < 0.05")
				# # # } else {
						# # # cat("There are no significant differences in prot across trt, p > 0.05")
				# # # }
	# # # if(kruskal.test(prot ~ trt, time1)$p.value <= 0.05) {
						# # # dunit<-DunnettTest(time1$prot, time1$trt, control="C",conf.level=0.95)
						# # # print(k)
						# # # print(dunit)
						# # # cat(capture.output(print(k), file="kruskalResultsTotalProt_x_Trt_Time1_postPriming.txt"))
						# # # cat(capture.output(print(dunit), file="dunnettPosthocTest_TotalProt_x_Trt_Time1_postPriming.txt"))
					# # # plot(time1$prot ~ time1$trt)
				# # # } else {
						# # # cat("No posthoc test required")
						# # # cat(capture.output(print(k), file="kruskalResultsTotalProt_x_Trt_Time1_postPriming.txt"))
				# # # }

### symbiont density
k<-kruskal.test(cells~interaction(trt,tp), data) 
			if(kruskal.test(cells ~ interaction(trt,tp), data)$p.value <= 0.05) {
						cat("There is an effect of trt*tp on cells , p < 0.05")
						tuk<-kruskalmc(cells ~ interaction(trt,tp), data) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsCells_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsCells_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsCells_TrtxTime.txt"))
				} else {
						cat("There are no significant differences in Cells across trt or tp, p > 0.05")
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsCells_TrtxTime.txt"))
				}
				
				
	### time 1 after priming
	# # # k<-kruskal.test(cells~trt, time1) 
			# # # if(kruskal.test(cells~trt, time1)$p.value <= 0.05) {
						# # # cat("There is an effect of trt on sym, p < 0.05")
				# # # } else {
						# # # cat("There are no significant differences in sym across trt, p > 0.05")
				# # # }
				
	# # # if(kruskal.test(cells ~ trt, time1)$p.value <= 0.05) {
				# # # dunit<-DunnettTest(time1$cells, time1$trt, control="C", conf.level=0.95)
						# # # print(k)
						# # # print(dunit)
						# # # cat(capture.output(print(k), file="kruskalResultsCells_x_Trt_Time1_postPriming.txt"))
						# # # cat(capture.output(print(dunit), file="dunnettPosthocTest_Cells_x_Trt_Time1_postPriming.txt"))
					# # # plot(time1$cells ~time1$trt)
				# # # } else {
						# # # cat("No posthoc test required")
						# # # cat(capture.output(print(k), file="kruskalResultsCells_x_Trt_Time1_postPriming.txt"))
				# # # }

### Fv/Fm
k<-kruskal.test(fvfm~interaction(trt,tp), data) 
			if(kruskal.test(fvfm ~ interaction(trt,tp), data)$p.value <= 0.05) {
						cat("There is an effect of trt*tp on fvfm , p < 0.05")
						tuk<-kruskalmc(fvfm ~ interaction(trt,tp), data) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsFvFm_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsFvFm_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsFvFm_TrtxTime.txt"))
				} else {
						cat("There are no significant differences in FvFm across trt or tp, p > 0.05")
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsFvFm_TrtxTime.txt"))
				}

	### time 1 after priming
	k<-kruskal.test(fvfm ~ trt, time1) 
			if(kruskal.test(fvfm ~ trt, time1)$p.value <= 0.05) {
						cat("There is an effect of trt on Fv/Fm, p < 0.05")
				} else {
						cat("There are no significant differences in Fv/Fm across trt, p > 0.05")
				}
	if(kruskal.test(fvfm ~ trt, time1)$p.value <= 0.05) {
						dunit<-DunnettTest(time1$fvfm, time1$trt, control="C",conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsFvFm_x_Trt_Time1_postPriming.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_FvFm_x_Trt_Time1_postPriming.txt"))
					plot(time1$prot ~ time1$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsFvFm_x_Trt_Time1_postPriming.txt"))
				}

			
#######   Timepoint 2   ########				

### chl per cell
	### time 2 after priming
	k<-kruskal.test(chlpcell ~ trt, time2noout) 
			if(kruskal.test(chlpcell ~ trt, time2noout)$p.value <= 0.05) {
						cat("There is an effect of trt on chl per cell, p < 0.05")
				} else {
						cat("There are no significant differences in chl per cell across trt, p > 0.05")
				}
	if(kruskal.test(chlpcell ~ trt, time2noout)$p.value <= 0.05) {
						dunit<-DunnettTest(time2noout$chlpcell, time2noout$trt, control="C",conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsChlpcell_x_Trt_Time2_postRecovery.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_Chlpcell_x_Trt_Time2_postRecovery.txt"))
					plot(time2noout$chlpcell~time2noout$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsChlpcell_x_Trt_Time2_postRecovery.txt"))
				}
	
### Protein per Cell	
	### time 2 after priming
	k<-kruskal.test(protpcell~trt, time2noout) 
			if(kruskal.test(protpcell ~ trt, time2noout)$p.value <= 0.05) {
						cat("There is an effect of trt on prot per cell, p < 0.05")
				} else {
						cat("There are no significant differences in prot per cell across trt, p > 0.05")
				}
				
	if(kruskal.test(protpcell ~ trt, time2noout)$p.value <= 0.05) {
				dunit<-DunnettTest(time2noout$protpcell, time2noout$trt, control="C", conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsProtpcell_x_Trt_Time2_postPriming.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_Protpcell_x_Trt_Time2_postPriming.txt"))
					plot(time2noout$protpcell ~ time2noout$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsProtpcell_x_Trt_Time2_postPriming.txt"))
				}

				
### total chl
	### time 2 after recovery
	k<-kruskal.test(chl ~ trt, time2) 
			if(kruskal.test(chl ~ trt, time2)$p.value <= 0.05) {
						cat("There is an effect of trt on chl, p < 0.05")
				} else {
						cat("There are no significant differences in chl across trt, p > 0.05")
				}

	if(kruskal.test(chl ~ trt, time2)$p.value <= 0.05) {
						dunit<-DunnettTest(time2$chl, time2$trt, control="C",conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsTotalChl_x_Trt_Time2_postRecovery.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_TotalChl_x_Trt_Time2_postRecovery.txt"))
					plot(time2$chl ~ time2$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsTotalChl_x_Trt_Time2_postRecovery.txt"))
				}
				
### total prot
	### time 2 after priming
	k<-kruskal.test(prot ~ trt, time2) 
			if(kruskal.test(prot ~ trt, time2)$p.value <= 0.05) {
						cat("There is an effect of prot on chl, p < 0.05")
				} else {
						cat("There are no significant differences in prot across trt, p > 0.05")
				}
	if(kruskal.test(prot ~ trt, time2)$p.value <= 0.05) {
						dunit<-DunnettTest(time2$prot, time2$trt, control="C",conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsTotalProt_x_Trt_Time2_postPriming.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_TotalProt_x_Trt_Time2_postPriming.txt"))
					plot(time2$prot ~ time2$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsTotalProt_x_Trt_Time2_postPriming.txt"))
				}

### symbiont density
	### time 2 after priming
	k<-kruskal.test(cells ~ trt, time2) 
			if(kruskal.test(cells ~ trt, time2)$p.value <= 0.05) {
						cat("There is an effect of trt on sym, p < 0.05")
				} else {
						cat("There are no significant differences in sym across trt, p > 0.05")
				}
	if(kruskal.test(cells ~ trt, time2)$p.value <= 0.05) {
						dunit<-DunnettTest(time2$cells, time2$trt, control="C",conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsCells_x_Trt_Time2_postRecovery.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_Cells_x_Trt_Time2_postRecovery.txt"))
					plot(time2$cells~time2$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsCells_x_Trt_Time2_postRecovery.txt"))
				}

### Fv/Fm
	### time 2 after priming
	k<-kruskal.test(fvfm ~ trt, time2) 
			if(kruskal.test(fvfm ~ trt, time2)$p.value <= 0.05) {
						cat("There is an effect of trt on Fv/Fm, p < 0.05")
				} else {
						cat("There are no significant differences in Fv/Fm across trt, p > 0.05")
				}
	if(kruskal.test(fvfm ~ trt, time2)$p.value <= 0.05) {
						dunit<-DunnettTest(time2$fvfm, time2$trt, control="C",conf.level=0.95)
						print(k)
						print(dunit)
						cat(capture.output(print(k), file="kruskalResultsFvFm_x_Trt_Time2_postPriming.txt"))
						cat(capture.output(print(dunit), file="dunnettPosthocTest_FvFm_x_Trt_Time2_postPriming.txt"))
					plot(time2$fvfm ~ time2$trt)
				} else {
						cat("No posthoc test required")
						cat(capture.output(print(k), file="kruskalResultsFvFm_x_Trt_Time2_postPriming.txt"))
				}

				
#######   Timepoint 3   ########				

# use Dunn's Test
# use Letters for multiple comparisons

### Chl per cell
	### time 3 after bleaching
	k<-kruskal.test(chlpcell ~ trt, time3noout) 
			if(kruskal.test(chlpcell ~ trt, time3noout)$p.value <= 0.05) {
						cat("There is an effect of trt on chl per cell, p < 0.05")
				} else {
						cat("There are no significant differences in chl per cell across trt, p > 0.05")
				}
	if(kruskal.test(chlpcell ~ trt, time3noout)$p.value <= 0.05) {
						tuk<-kruskalmc(fvfm ~ interaction(trt,tp), data) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						cat(capture.output(print(k), file="kruskalResultsFvFm_TrtxTime.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsFvFm_TrtxTime.txt"))
						cat(capture.output(print(let), file="lettersResultsFvFm_TrtxTime.txt"))
						
						duns<-dunnTest(chlpcell ~ trt, data=time3noout, method="bh")
						print(k)
						print(duns)
						cat(capture.output(print(k), file="kruskalResultsChlpcell_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(duns), file="dunnsPosthocTest_Chlpcell_x_Trt_Time3_postBleaching.txt"))
					plot(time3noout$chlpcell ~ time3noout$trt)
				} else {
						cat("No posthoc test required")
						print(k)
						cat(capture.output(print(k), file="kruskalResultsChlpcell_x_Trt_Time3_postBleaching.txt"))
				}

### Prot per cell
	### time 3 after bleaching
	k<-kruskal.test(protpcell ~ trt, time3noout) 
			if(kruskal.test(protpcell ~ trt, time3noout)$p.value <= 0.05) {
						cat("There is an effect of trt on prot per cell, p < 0.05")
				} else {
						cat("There are no significant differences in prot per cell across trt, p > 0.05")
				}
	if(kruskal.test(protpcell ~ trt, time3noout)$p.value <= 0.05) {
						duns<-dunnTest(protpcell ~ trt, data=time3noout, method="bh")
						print(k)
						print(duns)
						cat(capture.output(print(k), file="kruskalResultsProtpcell_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(duns), file="dunnsPosthocTest_Protpcell_x_Trt_Time3_postBleaching.txt"))
					plot(time3noout$protpcell ~ time3noout$trt)
				} else {
						cat("No posthoc test required")
						print(k)
						cat(capture.output(print(k), file="kruskalResultsProtpcell_x_Trt_Time3_postBleaching.txt"))
				}


### total chl
	### time 3 after bleaching
	k<-kruskal.test(chl ~ trt, time3) 
			if(kruskal.test(chl ~ trt, time3)$p.value <= 0.05) {
						cat("There is an effect of trt on chl, p < 0.05")
				} else {
						cat("There are no significant differences in chl across trt, p > 0.05")
				}
	if(kruskal.test(chl ~ trt, time3)$p.value <= 0.05) {
						duns<-dunnTest(chl ~ trt, data=time3, method="bh")
						print(k)
						print(duns)
						cat(capture.output(print(k), file="kruskalResultsTotalChl_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(duns), file="dunnsPosthocTest_TotalChl_x_Trt_Time3_postBleaching.txt"))
					plot(time3$chl ~ time3$trt)
				} else {
						cat("No posthoc test required")
						print(k)
						cat(capture.output(print(k), file="kruskalResultsTotalChl_x_Trt_Time3_postBleaching.txt"))
				}

### total prot
	### time 3 after bleaching
	k<-kruskal.test(prot ~ trt, time3) 
			if(kruskal.test(prot ~ trt, time3)$p.value <= 0.05) {
						cat("There is an effect of trt on prot, p < 0.05")
				} else {
						cat("There are no significant differences in prot across trt, p > 0.05")
				}
	if(kruskal.test(prot ~ trt, time3)$p.value <= 0.05) {
						duns<-dunnTest(prot ~ trt, data=time3, method="bh")
						print(k)
						print(duns)
						cat(capture.output(print(k), file="kruskalResultsTotalProt_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(duns), file="dunnsPosthocTest_TotalProt_x_Trt_Time3_postBleaching.txt"))
					plot(time3$prot ~ time3$trt)
				} else {
						cat("No posthoc test required")
						print(k)
						cat(capture.output(print(k), file="kruskalResultsTotalProt_x_Trt_Time3_postBleaching.txt"))
				}

### symbiont density
	### time 3 after bleaching
	k<-kruskal.test(cells ~ trt, time3) 
			if(kruskal.test(cells ~ trt, time3)$p.value <= 0.05) {
						cat("There is an effect of trt on sym, p < 0.05")
				} else {
						cat("There are no significant differences in sym across trt, p > 0.05")
				}
	if(kruskal.test(cells ~ trt, time3)$p.value <= 0.05) {
						tuk<-kruskalmc(cells ~ trt, time3) # multiple-comparison test
						test<- tuk$dif.com$difference # select logical vector
						names(test)<-row.names(tuk$dif.com) # add comparison names
						let<- multcompLetters(test,compare="<",threshold=0.05,Letters=c(letters,LETTERS,"."),reversed = FALSE); 
						print(k)
						print(tuk)
						print(let)
						duns<-dunnTest(cells ~ trt, data= time3, method="bh")
						print(k)
						print(duns)
						cat(capture.output(print(k), file="kruskalResultsCells_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(duns), file="dunnsPosthocTest_Cells_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(tuk), file="tukeysResultsCells_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(let), file="lettersResultsCells_x_Trt_Time3_postBleaching.txt"))
						
					plot(time3$cells ~ time3$trt)
				} else {
						cat("No posthoc test required")
						print(k)
						cat(capture.output(print(k), file="kruskalResultsCells_x_Trt_Time3_postBleaching.txt"))
				}

### Fv/Fm
	### time 3 after bleaching
	k<-kruskal.test(fvfm ~ trt, time3) 
			if(kruskal.test(fvfm ~ trt, time3)$p.value <= 0.05) {
						cat("There is an effect of trt on Fv/Fm, p < 0.05")
				} else {
						cat("There are no significant differences in Fv/Fm across trt, p > 0.05")
				}
	if(kruskal.test(fvfm ~ trt, time3)$p.value <= 0.05) {
						duns<-dunnTest(fvfm ~ trt, data= time3, method="bh")
						print(k)
						print(duns)
						cat(capture.output(print(k), file="kruskalResultsFvFm_x_Trt_Time3_postBleaching.txt"))
						cat(capture.output(print(duns), file="dunnsPosthocTest_FvFm_x_Trt_Time3_postBleaching.txt"))
					plot(time3$fvfm ~ time3$trt)
				} else {
						cat("No posthoc test required")
						print(k)
						cat(capture.output(print(k), file="kruskalResultsFvFm_x_Trt_Time3_postBleaching.txt"))
				}
				
###################################

### Some summary stats

chlpcellSM<- ddply(datanoout, c("day","trt"), summarise,
			n = length(chlpcell),
			mean = mean(chlpcell),
			median = median(chlpcell),
			iqr = IQR(chlpcell, na.rm=TRUE),
			sd = sd(chlpcell),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
chlpcellSM
write.csv(chlpcellSM,file="summaryChlpCell_x_Trt_Time.csv")

protpcellSM<- ddply(datanoout, c("day","trt"), summarise,
			n = length(protpcell),
			mean = mean(protpcell),
			median = median(protpcell),
			iqr = IQR(protpcell, na.rm=TRUE),
			sd = sd(protpcell),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
protpcellSM
write.csv(protpcellSM,file="summaryProtpCell_x_Trt_Time.csv")

totalchlSM<- ddply(data, c("day","trt"), summarise,
			n = length(chl),
			mean = mean(chl),
			median = median(chl),
			iqr = IQR(chl, na.rm=TRUE),
			sd = sd(chl),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
totalchlSM
write.csv(totalchlSM,file="summaryTotalChl_x_Trt_Time.csv")

totalprotSM<- ddply(data, c("day","trt"), summarise,
			n = length(prot),
			mean = mean(prot),
			median = median(prot),
			iqr = IQR(prot, na.rm=TRUE),
			sd = sd(prot),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
totalprotSM
write.csv(totalprotSM,file="summaryTotalProt_x_Trt_Time.csv")

cellsSM<- ddply(data, c("day","trt"), summarise,
			n = length(cells),
			mean = mean(cells),
			median = median(cells),
			min = min(cells),
			max = max(cells),
			iqr = IQR(cells, na.rm=TRUE),
			sd = sd(cells),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
cellsSM
write.csv(cellsSM,file="summaryCells_x_Trt_Time.csv")

fvfmSM<- ddply(data, c("day","trt"), summarise,
			n = length(fvfm),
			mean = mean(fvfm),
			median = median(fvfm),
			iqr = IQR(fvfm, na.rm=TRUE),
			sd = sd(fvfm),
			se = sd/ sqrt(n),
			ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
		)
fvfmSM
write.csv(fvfmSM,file="summaryFvFm_x_Trt_Time.csv")


### Cells over time

# first, add a row of control values for time zero: the median Control values
# dummy_data <- data.frame(day= c(0,0,0,0,0,0), trt= c('C','N','LL','LH','HL','HH'), n= c(8,8,8,8,8,8), mean= c(3504030,3504030,3504030,3504030,3504030,3504030) , median= c(3504030,3504030,3504030,3504030,3504030,3504030), iqr= c(730466.8,730466.8,730466.8,730466.8,730466.8,730466.8), sd=c(482316.5, 482316.5, 482316.5, 482316.5, 482316.5, 482316.5), se=c(170524.64, 170524.64, 170524.64, 170524.64, 170524.64, 170524.64), ci=c(403226.7, 403226.7, 403226.7, 403226.7, 403226.7, 403226.7))
# headTail(dummy)

# plot the data

		black.italic.text<- element_text(family="Arial", face="bold.italic", color="black", size=28)
		black.bold.text<- element_text(family="Arial", face="bold", color="black", size=24)
		italic.text<- element_text(family="Arial", face="italic",color="black", size=24)
		plain.text<- element_text(family="Arial", face="plain",color="black", size=18)
#pdf(file="cells_trt_time.pdf")
	p<-ggplot(cellsSM, aes(x=day, y=median, colour=trt)) +
		geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
		geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
		geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
		#geom_ribbon(aes(fill=trt,ymin=min, ymax=max), alpha=0.05,position=position_dodge(width=0.6)) +
		scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
		#xlim(2,14) +
		theme_bw() +
		theme(legend.position=c(0.08,0.15), legend.title=element_blank())
	p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Symbiont Density (cells cm'^-2*")"))))
	plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
#dev.off()


########### Now load all the Fv/Fm data from all timepoints
## set the working directory for all analyses
setwd("~/Desktop/PUBS/StressMemoryPaper/additionalAnalyses/")

####################################
## load the adjusted Chl data
adj_chl<-read.csv("~/Desktop/PUBS/stressmemorypaper/additionalAnalyses/adj_totalchlSM.csv", header=TRUE, stringsAsFactors=TRUE)
adj_chl$trt=factor(adj_chl$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
## load the adjusted cell data
adj_cells<-read.csv("~/Desktop/PUBS/stressmemorypaper/additionalAnalyses/adj_cellsSM.csv", header=TRUE, stringsAsFactors=TRUE)
adj_cells$trt=factor(adj_cells$trt,levels=c("C","N","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))

## load the fvfm data
dat=read.csv("~/Desktop/PUBS/StressMemoryPaper/additionalAnalyses/allfvfmtp.csv", header=TRUE, stringsAsFactors=TRUE)
## reorder factors to appear as desired
dat$genet=factor(dat$genet,levels=c("C","D","E","F","G","H","I","J"))
dat$trt=factor(dat$trt,levels=c("C","U","LL","LH","HL","HH"), labels=c("C","N","LL","LH","HL","HH"))
#dat$day=factor(dat$day, levels=c("1","2","3","6","10","11")) #day is a numerical value, only do this for timepoint

## examine the complete dataset
names(dat) # what variables are there
str(dat) # look at the data structure
headTail(dat) # look at the first few rows
tail(dat) # look at the last few rows
summary(dat) # look for possible missing values or unbalanced design

# some summary stats
library(plyr)

fvfmStats <- ddply(dat, c("day","trt"), summarise,
                   n = length(fvfm),
                   mean = mean(round(fvfm,3)),
                   median = median(round(fvfm,3)),
                   min = min(round(fvfm,3)),
                   max = max(round(fvfm,3)),
                   iqr = IQR(round(fvfm,3), na.rm=TRUE),
                   sd = sd(round(fvfm,3)),
                   se = sd/ sqrt(n),
                   ci = se * qt(.95/2 + .5, n-1)	# Confidence at the 95% interval
)
fvfmStats
fvfmStats %>% mutate_if(is.numeric, round, 3)
write.csv(fvfmStats,file="summaryFvFm_x_Trt_Timeseries.csv")


## This is to identify the outliers via Mahalanobis distance
library(stats)
library(reshape)
library(psych)

data_wide<-reshape(dat, idvar = "trt_genet", timevar = "day", direction = "wide")
headTail(data_wide)
fvfm<- select(data_wide, fvfm.0, fvfm.1, fvfm.2, fvfm.5, fvfm.10, fvfm.11); dim(fvfm); summary(fvfm)
mahal_fvfm = mahalanobis(fvfm, colMeans(fvfm, na.rm=TRUE), cov(fvfm,use="pairwise.complete.obs"));
mahal_fvfm
## determine a cutoff score
cutoff=qchisq(1-.001, ncol(fvfm))
cutoff
ncol(fvfm) # df
summary(mahal_fvfm < cutoff)  ## NAs
noout_fvfm = subset(fvfm, mahal_fvfm < cutoff) # remove samples with mahal values < the cutoff (> are outlying)
#noout_fvfm ## outliers #NONE  were removed from all timepoints
#THIS DF HAS NO NAs BECAUSE THERE ARE NO OUTLIERS!



####### This plots median +-1sd of cells, total chl and Fv/Fm over time in all trts #######

black.italic.text<- element_text(family="Arial", face="bold.italic", color="black", size=28)
black.bold.text<- element_text(family="Arial", face="bold", color="black", size=24)
italic.text<- element_text(family="Arial", face="italic",color="black", size=24)
plain.text<- element_text(family="Arial", face="plain",color="black", size=18)

# cells
p<-ggplot(cellsSM, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw() +
  theme(legend.position=c(0.08,0.15), legend.title=element_blank())
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Symbiont Density (cells cm'^-2*")"))))
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# total chl
p<-ggplot(totalchlSM, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw() +
  theme(legend.position=c(0.08,0.15), legend.title=element_blank())
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Total Chlorophyll (' ~ mu *'g Chl cm'^-2*") ")))) 
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# adjusted total chl (only day 10 and 11)
p<-ggplot(adj_chl, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.05)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.05)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.05)) +
  scale_x_continuous(breaks=c(9, 10,11,12)) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Total Chlorophyll (' ~ mu *'g Chl cm'^-2*") ")))) 
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# adjusted cells (only day 10 and 11)
p<-ggplot(adj_cells, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.05)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.05)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.05)) +
  scale_x_continuous(breaks=c(9, 10,11,12)) +
  theme_bw() +
  theme(legend.position="none")
p.labs <- p +  labs(x = "Time (days)", y = expression(bold(paste('Symbiont Density (cells cm'^-2*")"))))
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))

# fvfm with all timepoints
p<-ggplot(fvfmStats, aes(x=day, y=median, colour=trt)) +
  geom_point(aes(colour=trt, fill=trt), size=5, position=position_dodge(width=0.7)) +
  geom_line(aes(colour=trt), size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.1, position=position_dodge(width=0.7)) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  theme_bw() +
  theme(legend.position=c(0.08,0.15), legend.title=element_blank())
p.labs <- p +  labs(x = "Time (days)", y = expression(bolditalic(paste('Photochemical Efficiency (F'[V]*'/F'[M]*") " )))  )
plot(p.labs + theme(axis.title=black.bold.text, axis.text=black.bold.text, legend.text=plain.text))
