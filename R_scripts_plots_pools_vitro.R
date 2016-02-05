# Line plots for plerixafor +/- G-CSF experiments
# In vitro data
# Paul Miller, paulhmiller@gmail.com

library(reshape2)
library(ggplot2)
library(grid)

source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")
setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/in_vitro_pools')

dat <- read.csv("plerixafor_pools_vitro.csv",header=T,check.names=F) #,row.names=1)
dat <- dat[c(1:24)]  # Take just the columns you need 
assays <- names(dat[5:24])  # Take assay names
dat <- melt(dat,id.vars=c("Pool","Sample","Timepoint","Group"),variable.name="assay",value.name="value",na.rm=TRUE)  # Convert to tall
dat$value <- log10(dat$value)  # Log10 transform

#Rename group factors
levels(dat$Group)[levels(dat$Group)=="A"]  <-  "P"
levels(dat$Group)[levels(dat$Group)=="B"]  <-  "G+P"
levels(dat$Group)[levels(dat$Group)=="A+B"]  <-  "Baseline"

#Re-order factor levels
dat$Timepoint <- factor(dat$Timepoint,levels=c("BL","Day -1","Day 0","Day 1"))  

# Pool baseline data so that both groups use the same pooled set of baseline data.
# Marrow
BLP <- dat[dat$Sample=="Marrow" & dat$Timepoint=="BL",]
BLP$Group <- "P"
BLG <- dat[dat$Sample=="Marrow" & dat$Timepoint=="BL",]
BLG$Group <- "G+P"
notBL <- dat[dat$Sample=="Marrow" & dat$Timepoint!="BL",]
BM <- rbind(BLP,BLG,notBL)
# Blood
BLP <- dat[dat$Sample=="Blood" & dat$Timepoint=="BL",]
BLP$Group <- "P"
BLG <- dat[dat$Sample=="Blood" & dat$Timepoint=="BL",]
BLG$Group <- "G+P"
notBL <- dat[dat$Sample=="Blood" & dat$Timepoint!="BL",]
PB <- rbind(BLP,BLG,notBL)

# Re-name time points
levels(PB$Timepoint)[levels(PB$Timepoint)=="BL"] <- "D0"
levels(PB$Timepoint)[levels(PB$Timepoint)=="Day -1"] <- "Post-G"
levels(PB$Timepoint)[levels(PB$Timepoint)=="Day 0"] <- "P+4 h"
levels(PB$Timepoint)[levels(PB$Timepoint)=="Day 1"] <- "P+24 h"


PB$Group<-factor(PB$Group,levels=c("P","G+P")) 

PB<-summarySE(PB, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) 
# Above line is unnecessary because there is only one measurement (no error bar); 
# but it seems to be needed to have the P group shapes on top of the G+P shapes. 


### Mobilised Peripheral Blood sourced cells ### 

CD34ul <- PB[PB$assay=="34 per uL of original",]
CD34RAul <- PB[PB$assay=="34p45RAn per uL of original",]
CD3438ul <- PB[PB$assay=="34p38n per uL of original",]
CD34RA38ul <- PB[PB$assay=="34p38n45RAn per uL of original",]
viable <- PB[PB$assay=="viable per uL of original",]

tCFC <- PB[PB$assay=="CFC",]
w3LTCTNC <- PB[PB$assay=="TNC.wk3",]
w3LTC34 <- PB[PB$assay=="CD34.wk3",]
w3LTCIC <- PB[PB$assay=="LTC-IC.wk3",]
w6LTCTNC <- PB[PB$assay=="TNC.wk6",]
w6LTC34 <- PB[PB$assay=="CD34.wk6",]
w6LTCIC <- PB[PB$assay=="LTC-IC.wk6",]

plot <- ggplot(data=CD34ul,aes(x=factor(Timepoint), y=as.numeric(value), group=Group, colour=Group))+
	ylab(expression(paste("CD34 x 10"^"3","/mL"))) +
	coord_cartesian(ylim=c(-0.39794,1.60206)) +
	geom_point(aes(shape=Group),size=3) +
	scale_colour_manual(values=p2) +
	geom_line()+
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8),labels = c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8))) + 
	annotation_logticks(sides = "l", size=0.3) +
	scale_x_discrete("")+ #, labels=c("??"))  +
	guides(colour=F, shape=F)+
	themePM1() 
plot

plot %+% CD34ul +
	ylab(expression(paste("CD34 x 10"^"3","/mL"))) +
	coord_cartesian(ylim=c(-0.39794,1.60206)) +
	geom_point(data=subset(CD34ul, Timepoint == "Post-G") ,size=3, colour = "#4169E1", shape = 17) +
	themePM1() 
ggsave(filename="pool PB CD34.pdf",width=5.75,height=5.75, units="cm")

plot %+% tCFC +
	ylab(expression(paste("CFC x 10"^"3","/mL"))) +
	coord_cartesian(ylim=c(-0.39794,1.60206)) +
	geom_point(data=subset(tCFC, Timepoint == "Post-G") ,size=3, colour = "#4169E1", shape = 17) +
  themePM1() + theme(plot.margin = unit(c(0, 0.2, 0, 0.1+0.12), "cm"))
ggsave(filename="pool PB tCFC.pdf",width=5.75,height=5.75, units="cm")

plot %+% w3LTCIC +
	ylab(expression(paste("LTC-CFC x 10"^"3","/mL"))) +
	coord_cartesian(ylim=c(-2,2)) +
	geom_point(data=subset(w3LTCIC, Timepoint == "Post-G") ,size=3, colour = "#4169E1", shape = 17) +
  themePM1() + theme(plot.margin = unit(c(0, 0.2, 0, 0.1-0.095), "cm"))
ggsave(filename="pool PB w3LTCIC.pdf",width=5.75,height=5.75, units="cm")

plot %+% w6LTCIC +
	ylab(expression(paste("LTC-CFC x 10"^"3","/mL")))+
	geom_point(data=subset(w6LTCIC, Timepoint == "Post-G") ,size=3, colour = "#4169E1", shape = 17) +
	coord_cartesian(ylim=c(-2,2)) +
	  themePM1()
ggsave(filename="pool PB w6LTCIC.pdf",width=5.75,height=5.75, units="cm")



###MOBILIZED BONE MARROW SOURCED CELLS###

BMi <- BM[BM$Sample=="Marrow",]

# Change scale of CFC from per 10e5
BMi[BMi$assay=="CFC","value"] <- BMi[BMi$assay=="CFC","value"]-2

BM <- summarySE(BMi, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) 
# summarySE provides the standard deviation, standard error of the mean, and 
# a (default 95%) confidence interval. #measurevar is the x-axis to only plot 
# means, change plots to plot XXXs only.

# Omit day-1 (due to poor recovery data)
BM <- BM[BM$Timepoint != "Day -1",]

# Re-name time points
levels(BM$Timepoint)[levels(BM$Timepoint)=="BL"] <- "D0"
levels(BM$Timepoint)[levels(BM$Timepoint)=="Day -1"] <- "Post-G"
levels(BM$Timepoint)[levels(BM$Timepoint)=="Day 0"] <- "P+4 h"
levels(BM$Timepoint)[levels(BM$Timepoint)=="Day 1"] <- "P+24 h"

CD34pcs <- BM[BM$assay=="34 percent of prefreeze",]
tCFC <- BM[BM$assay=="CFC",]
w3LTCTNC <- BM[BM$assay=="TNC.wk3",]
w3LTC34 <- BM[BM$assay=="CD34.wk3",]
w3LTCIC <- BM[BM$assay=="LTC-IC.wk3",]
w6LTCTNC <- BM[BM$assay=="TNC.wk6",]
w6LTC34 <- BM[BM$assay=="CD34.wk6",]
w6LTCIC <- BM[BM$assay=="LTC-IC.wk6",]


plot<-ggplot(data=CD34pcs,aes(x=factor(Timepoint), y=as.numeric(value), group=Group, colour=Group))+
	geom_point(aes(shape=Group),size=3) +
	geom_errorbar(data=CD34pcs,aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +
	scale_colour_manual(values=p2) +
	coord_cartesian(ylim=c(-2,1)) +
	geom_line()+
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8),labels = c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8))) + 
	annotation_logticks(sides = "l", size=0.3) +
	scale_x_discrete("")+ 
	guides(colour=F,lty=F,shape=F)+
	ggtitle("") +
	ylab("% CD34") +
	themePM1()  
plot

plot %+% CD34pcs +
	coord_cartesian(ylim=c(-2,1)) +
	ylab("% CD34") +
	themePM1()  
ggsave(filename="pool BM CD34.pdf",width=5.75,height=5.75, units="cm")

plot %+% tCFC +
	coord_cartesian(ylim=c(-1,1)) +
	ylab(expression(paste("CFC per 10"^"3","cells"))) +
  themePM1() + theme(plot.margin = unit(c(0, 0.2, 0, 0.1-0.07), "cm"))
ggsave(filename="pool BM CFC.pdf",width=5.75,height=5.75, units="cm")

plot %+% w3LTCIC +
	coord_cartesian(ylim=c(-1,2)) +
	ylab(expression(paste("LTC-CFC per 10"^"5","cells")))  +
  themePM1() + theme(plot.margin = unit(c(0, 0.2, 0, 0.1-0.155), "cm"))
ggsave(filename="pool BM w3LTCIC.pdf",width=5.75,height=5.75, units="cm")

plot %+% w6LTCIC +
	coord_cartesian(ylim=c(-1,2)) +
	ylab(expression(paste("LTC-CFC per 10"^"5","cells"))) +
  themePM1() + theme(plot.margin = unit(c(0, 0.2, 0, 0.1-0.095), "cm"))
ggsave(filename="pool BM w6LTCIC.pdf",width=5.75,height=5.75, units="cm")
