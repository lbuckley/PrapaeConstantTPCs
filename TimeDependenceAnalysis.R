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
library(mgcv)

colm<- viridis_pal(option = "mako")(8)
cols<- colm[c(2,4,7)]
cols2<- colm[c(3,6)]

#toggle between desktop (y) and laptop (n)
desktop<- "n"

# Load data
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/WARP/projects/TPCconstant/")
if(desktop=="n") setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My Drive/Buckley/Work/WARP/projects/TPCconstant/")

tpc<- read.csv("out/PastPresentFilteredConstantTpc2024.csv")

#-----
#instar label
in.lab<- c("4th instar", "5th instar")
tpc$in.lab <- in.lab[match(tpc$instar, c(4,5))]

#set up time classes: 6 and 24 hr studies
tpc$time.class<- NA
#expand windows
tpc$time.class[which(tpc$time>5.5 & tpc$time<7.5)]<-6 
tpc$time.class[which(tpc$time>21.5 & tpc$time<26)]<-24

#set up broader time class
tpc$time.class.b<- NA
#expand windows
tpc$time.class.b[which(tpc$time>5 & tpc$time<10)]<-6 
tpc$time.class.b[which(tpc$time>=21.5 & tpc$time<26)]<-24 #or 27

#fix hour label
#hr.lab<- c("6 hour","24 hour")
#tpc.plot$hr.lab <- hr.lab[match(tpc.plot$time.class.b, c(6,24))]
#tpc.plot$hr.lab <- factor(tpc.plot$hr.lab, levels=c("6 hour","24 hour"), ordered=TRUE)

#--------
#check time dependence of feeding
tpc1<- tpc[which(tpc$time.per=="past"),]
tpc1$UniID <- factor(tpc1$UniID)

#restrict to potential time limits
mod= lm(rgrlog ~ Mo + poly(temp,3)*time, data= tpc1[which(tpc1$time.class==6),]) 
mod= lm(rgrlog ~ Mo + poly(temp,3)*time, data= tpc1[which(tpc1$time.class==24),]) 
anova(mod)

plot_model(mod, type = "pred", terms = c("time", "temp"), show.data=TRUE)

#lme: issues accounting for individual
mod.lmer <- lme(rgrlog ~ Mo + poly(temp,3)*time, random=~1|mom/ID, data = na.omit(tpc1[which(tpc1$time>5 & tpc1$time<10),]))
mod.lmer <- lme(rgrlog ~ Mo + poly(temp,3)*time, random=~1|mom/ID, data = na.omit(tpc1[which(tpc1$time>21 & tpc1$time<26),]))
anova(mod.lmer)

plot_model(mod.lmer, type = "pred", terms = c("time", "temp"), show.data=TRUE)

#plot temp and time
ggplot(tpc1[which(tpc1$time>5 & tpc1$time<10),], aes(x = temp, y=gr, color = time)) + 
  geom_point(alpha=0.4, position = position_jitterdodge()) +geom_smooth()

#Model
#mod.lmer <- lme(rgrlog ~ Mo + poly(temp,3)*time*time.per*time.class*instar, random=~1|UniID, data = na.omit(tpc))
#anova(mod.lmer)

#plot feeding rate over time to assess time adjustments
tpc.plot <- tpc[which(tpc$time>0),]
tpc.plot$mass.r<- log10(tpc.plot$fw/tpc.plot$Mo)

ggplot(tpc.plot[tpc.plot$time.class==6,], aes(x = time, y = gr, color = time.per)) + #y = mgain
  geom_point(alpha=0.4, position = position_jitterdodge()) +
  facet_grid(temp ~ instar) +xlim(0,40)

ggplot(tpc.plot, aes(x = time, y = gr, color = time.per)) + #y = mgain
  geom_point(alpha=0.4, position = position_jitterdodge()) +
  facet_grid(temp ~ instar) +xlim(0,40)

#----
# #GAM
# mod.gam <- gam(rgrlog ~ Mo + 
#                  s(temp, bs = "cr", k=5, by = interaction(time, time.class, instar))+
#                  time*time.class*instar, #+s(UniID, bs = "re"),                               # Random intercept
#                data = na.omit(tpc),
#                method = "REML")
# summary(mod.gam)

#-----
#Compare TPCs with two time classes

#rgr or gr
tpc.plot$grow= tpc.plot$gr

#narrower and broader period comparisons
tpc.agg <- tpc.plot %>%
  group_by(temp, time.per, time.class, instar, in.lab) %>% 
  dplyr::summarise(
    mean = mean(grow, na.rm = TRUE),
    n= length(grow),
    sd = sd(grow, na.rm = TRUE),
    mean.mass = mean(fw-Mo, na.rm = TRUE),
    sd.mass = sd(fw-Mo, na.rm = TRUE) )

tpc.agg$se= tpc.agg$sd / sqrt(tpc.agg$n)
tpc.agg$se.mass= tpc.agg$sd / sqrt(tpc.agg$n)

tpc.agg.b <- tpc.plot %>%
  group_by(temp, time.per, time.class.b, instar, in.lab) %>% 
  dplyr::summarise(
    mean = mean(grow, na.rm = TRUE),
    n= length(grow),
    sd = sd(grow, na.rm = TRUE),
    mean.mass = mean(fw-Mo, na.rm = TRUE),
    sd.mass = sd(fw-Mo, na.rm = TRUE) )

tpc.agg.b$se= tpc.agg.b$sd / sqrt(tpc.agg.b$n)
tpc.agg.b$se.mass= tpc.agg.b$sd / sqrt(tpc.agg.b$n)
names(tpc.agg.b)[3] <- "time.class"

tpc.agg$time.type<- "narrow"
tpc.agg.b$time.type<- "broad"
tpc.all<- rbind(tpc.agg, tpc.agg.b)

#tpc.agg$mean-tpc.agg.b$mean

#-----
#plot 

#drop NA time class
tpc.all<- tpc.all[which(!is.na(tpc.all$time.class)),]

#COMPARE NARROW AND BROAD DEFINITION
FigSx_rgrtime.timeplot.comp <- ggplot(tpc.all[which(tpc.all$time.per=="past"),], aes( x = temp, y = mean, lty = time.per, color=factor(time.type))) +
  #geom_errorbar(data=tpc.all, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(size=2.5) + geom_line(linewidth=1.25)+
  facet_grid(time.class ~ in.lab) +
  theme_bw(base_size=16)+xlab("Temperature (°C)")+ylab("Growth rate (g/g/h)")+
  scale_color_viridis_d()+
  labs(color="Time type", fill="Time type", lty="Time (hr)")+theme(legend.position="bottom")+
  geom_hline(yintercept=0, col="darkgray")

#--------------
#model differences

#find points omitted
tpc.plot$added<- "both"
tpc.plot$added[which(is.na(tpc.plot$time.class) & !is.na(tpc.plot$time.class.b))]<- "broad"

tpc.comp.n<- tpc.plot[which(!is.na(tpc.plot$time.class)),c("rgrlog","grow","Mo","temp","time","mom","ID","time.per","instar","time.class.b","in.lab")] 
tpc.comp.b<- tpc.plot[which(!is.na(tpc.plot$time.class.b)),c("rgrlog","grow","Mo","temp","time","mom","ID","time.per","instar","time.class.b","in.lab")] 

tpc.comp.n$time.type<- "narrow"
tpc.comp.b$time.type<- "broad"

tpc.comp<- rbind(tpc.comp.n, tpc.comp.b)

#restrict to potential time limits
mod.lmer4 <- lme(grow ~ poly(temp,3)*time.per*Mo, random=~1|time.type/mom/ID, data = na.omit(tpc.comp[which(tpc.comp$time.class.b==6 &tpc.comp$instar==4),]))

mod.lmer5 <- lme(grow ~ poly(temp,3)*time.per*Mo, random=~1|time.type/mom/ID, data = na.omit(tpc.comp[which(tpc.comp$time.class.b==6 &tpc.comp$instar==5),]))

summary(mod.lmer4)
anova(mod.lmer4)
sigma(mod.lmer4)
VarCorr(mod.lmer4)

#----
#Table S1

#Save anova
tablesS1<- rbind(as.data.frame(anova(mod.lmer4)), as.data.frame(anova(mod.lmer5)) )

colnames(tablesS1)[3:4]<- c("F","p")
tablesS1$sig<-""
tablesS1$sig[tablesS1$p<0.05]<-"*"
tablesS1$sig[tablesS1$p<0.01]<-"**"
tablesS1$sig[tablesS1$p<0.001]<-"***"
tablesS1$F= round(tablesS1$F,1)
tablesS1$p= round(tablesS1$p,4)

#Table S1
write.csv(tablesS1, "./figures/tablesS1_time_window.csv")

#------
#plot data points

rgr.plot.nb <- ggplot(tpc.plot[which(tpc.plot$time.per=="past" & !is.na(tpc.plot$time.class.b) ),], aes( x = temp, y = grow, color = added)) +
  geom_point(alpha=0.4, position = position_jitterdodge()) +
  facet_grid(time.class.b ~ in.lab, scales="free_y") 

#-----------------
#aes(fill=time.type), pch=21, 
#scale_fill_manual(values=colm[c(4,7)])+

#check differences
tpc.s<- tpc[which(!is.na(tpc$time.class.b)),]
tpc.s2<- tpc.s[which(tpc.s$instar==5 & tpc.s$temp==29),]

#data by time 
ggplot(tpc.s[which(tpc.s$time.per=="past" & tpc.s$time.class==6),], aes(x = time, y = gr, color = factor(temp))) + #y = mgain
  geom_point(alpha=0.4, position = position_jitterdodge()) + geom_vline(xintercept = 27)

#check outliers
ggplot(tpc.s[tpc.s$time.class.b==6,], aes(x = time, y = gr, color = time.per)) + #y = mgain
  geom_point(alpha=0.4, position = position_jitterdodge()) 

tpc.agg.b[which(tpc.agg.b$time.per=="past" & tpc.agg.b$time.class.b==24 & tpc.agg.b$instar==4),]

