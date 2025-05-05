# Taylor M. Hatcher, https://github.com/taylorhatcher/WARP2024/blob/main/24hrTPCConstantAnalysis.R
# 24 Hr Constant TPC Analysis for P.rapae
# This analysis compares historical constant feeding relative growth rates at 6hr and 24hr at a constant temperature

# Load required libraries 
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(reshape2)
library(viridis)
library(tidyverse)
library(lubridate)
library(nlme)
library(lme4)
library(sjPlot)
library(car)

colm<- viridis_pal(option = "mako")(8)
cols<- colm[c(2,4,7)]
cols2<- colm[c(3,6)]

#toggle between desktop (y) and laptop (n)
desktop<- "n"

# Load data
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/WARP/projects/TPCconstant/out/")
if(desktop=="n") setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My Drive/Buckley/Work/WARP/projects/TPCconstant/out/")

tpc<- read.csv("PastPresentFilteredConstantTpc2024.csv")

#-----
#instar label
in.lab<- c("4th instar", "5th instar")
tpc$in.lab <- in.lab[match(tpc$instar, c(4,5))]

#set up time classes: 6 and 24 hr studies
tpc$time.class<- NA
#expand windows
tpc$time.class[which(tpc$time>5.5 & tpc$time<7.5)]<-6 #or to 8.5
tpc$time.class[which(tpc$time>21.5 & tpc$time<26)]<-24
#--------
#check time dependence of feeding
tpc1<- tpc[which(tpc$time.per=="past"),]
tpc1$UniID <- factor(tpc1$UniID)

#restrict to potential time limits
mod= lm(rgrlog ~ Mo + poly(temp)*time, data= tpc1[which(tpc1$time.class==6),]) 
mod= lm(rgrlog ~ Mo + poly(temp)*time, data= tpc1[which(tpc1$time.class==24),]) 
anova(mod)

plot_model(mod, type = "pred", terms = c("time", "temp"), show.data=TRUE)

#lme: issues accounting for individual
mod.lmer <- lme(rgrlog ~ Mo + poly(temp)*time, random=~1|mom, data = na.omit(tpc1[which(tpc1$time>5 & tpc1$time<10),]))
mod.lmer <- lme(rgrlog ~ Mo + poly(temp)*time, random=~1|mom, data = na.omit(tpc1[which(tpc1$time>21 & tpc1$time<26),]))
anova(mod.lmer)

plot_model(mod.lmer, type = "pred", terms = c("time", "temp"), show.data=TRUE)

#plot temp and time
ggplot(tpc1[which(tpc1$time>5 & tpc1$time<10),], aes(x = temp, y=gr, color = time)) + 
  geom_point(alpha=0.4, position = position_jitterdodge()) +geom_smooth()

#Model
mod.lmer <- lme(rgrlog ~ Mo + poly(temp)*time*time.per*time.class*instar, random=~1|UniID, data = na.omit(tpc))
anova(mod.lmer)

#plot feeding rate over time to assess time adjustments
tpc.plot <- tpc[which(tpc$time>0),]
tpc.plot$mass.r<- log10(tpc.plot$fw/tpc.plot$Mo)

ggplot(tpc.plot[tpc.plot$time.class==24,], aes(x = time, y = gr, color = time.per)) + #y = mgain
  geom_point(alpha=0.4, position = position_jitterdodge()) +
  facet_grid(temp ~ instar) +xlim(0,40)

ggplot(tpc.plot, aes(x = time, y = gr, color = time.per)) + #y = mgain
  geom_point(alpha=0.4, position = position_jitterdodge()) +
  facet_grid(temp ~ instar) +xlim(0,40)

#---------------------
#PLOT
#Figure 2
#plot temperature sensitivity
tpc.plot <- tpc[which(!is.na(tpc$time.class)),]

#update labels
in.lab<- c("4th instar", "5th instar")
tpc.plot$in.lab <- in.lab[match(tpc.plot$instar, c(4,5))]
hr.lab<- c("6 hour","24 hour")
tpc.plot$hr.lab <- hr.lab[match(tpc.plot$time.class, c(6,24))]
tpc.plot$hr.lab <- factor(tpc.plot$hr.lab, levels=c("6 hour","24 hour"), ordered=TRUE)
tpc.plot$time.per <- factor(tpc.plot$time.per, levels=c("past","current"), ordered=TRUE)

#rgr or gr
tpc.plot$grow= tpc.plot$rgrlog
#tpc.plot$grow= tpc.plot$gr

rgr.plot <- ggplot(tpc.plot, aes( x = temp, y = grow, color = time.per)) +
  geom_point(alpha=0.4, position = position_jitterdodge()) +
  facet_grid(hr.lab ~ in.lab) 

#---
#plot temp means
tpc.agg <- tpc.plot %>%
  group_by(temp, time.per, time.class, instar, hr.lab, in.lab) %>% 
  dplyr::summarise(
    mean = mean(grow, na.rm = TRUE),
    n= length(grow),
    sd = sd(grow, na.rm = TRUE) )
tpc.agg$se= tpc.agg$sd / sqrt(tpc.agg$n)

#restrict to points 
tpc.agg<- tpc.agg[which(tpc.agg$n>5),]

#plotting temp means with error bars
Fig2_growth.plot= rgr.plot + 
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(data=tpc.agg, aes(x=temp, y = mean, fill=time.per), size=3, col="black", pch=21)+
  geom_line(data=tpc.agg, aes(x=temp, y = mean))+
  theme_bw(base_size=16)+xlab("Temperature (°C)")+ylab("Growth rate (g/g/h)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Time period", fill="Time period")+theme(legend.position="bottom")

#add family lines
#not individuals at every temp
#rgr.plot= rgr.plot +
#geom_line(data=tpc.agg.f, aes(x=temp, y = mean, group=mom), linewidth=1)

#-------
#plot time periods together
Fig3_rgrtime.plot <- ggplot(tpc.agg, aes( x = temp, y = mean, color = time.per, lty=factor(time.class))) +
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(size=2) + geom_line()+
  facet_grid(.~ in.lab) +
  theme_bw(base_size=16)+xlab("Temperature (°C)")+ylab("Growth rate (g/g/h)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Time period", fill="Time period", lty="Time (hr)")+theme(legend.position="bottom")+
  geom_hline(yintercept=0, col="darkgray")

#plot 4th and 5th together
rgr45.plot <- ggplot(tpc.agg, aes( x = temp, y = mean, color = time.per, lty=factor(instar))) +
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(size=2) + geom_line()+
  facet_grid(hr.lab ~ .) +
  theme_bw(base_size=16)+xlab("Temperature (°C)")+ylab("Growth rate (g/g/h)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Time period", fill="Time period", lty="Instar")+
  theme(legend.position="bottom")

#----------------
#distribution plots
tpc$time.per<- factor(tpc$time.per, levels=c("past","current", ordered=TRUE))

#initial weights
Fig4_plot.mass<- ggplot(tpc, aes(x=Mo,color=time.per, group=time.per)) + 
  geom_density(aes(fill=time.per), alpha=0.5, adjust=1.8)+
  ylab("Density") +xlab("Mass (mg)")+
  facet_wrap(.~in.lab, scales="free")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(2,6)])+
theme_classic(base_size=16) +theme(legend.position = c(0.9, 0.8))+
  labs(color="Time period", fill="Time period")

#compare distributions
#compare variance
leveneTest(Mo ~ time.per, tpc[tpc$instar==4,])
#compare means
t.test(Mo ~ time.per, data=tpc[tpc$instar==4,], alternative = "two.sided", var.equal = FALSE)
