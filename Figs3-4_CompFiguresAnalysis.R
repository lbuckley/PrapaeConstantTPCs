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

# Load data
tpc<- read.csv("data/PastPresentFilteredConstantTpc2024.csv")

#-----
#instar label
in.lab<- c("4th instar", "5th instar")
tpc$in.lab <- in.lab[match(tpc$instar, c(4,5))]

#set up time classes: 6 and 24 hr studies
tpc$time.class<- NA
#expand windows
tpc$time.class[which(tpc$time==0)]<-0 
tpc$time.class[which(tpc$time>5 & tpc$time<10)]<-6 
tpc$time.class[which(tpc$time>21.5 & tpc$time<26)]<-24

# tpc.ch<- tpc[which(tpc$instar==4 & tpc$temp==35 & tpc$time.per=="past"),]
# ggplot(tpc.ch, aes( x = time, y = rgrlog, color = time.class))+
#   geom_point()+xlim(0,10)

fams<- unique(paste(tpc$mom, tpc$ID  ))

#---
#plot temperature sensitivity
tpc.plot <- tpc
#tpc.plot <- tpc[which(!is.na(tpc$time.class)),]

#update labels
in.lab<- c("4th instar", "5th instar")
tpc.plot$in.lab <- in.lab[match(tpc.plot$instar, c(4,5))]
hr.lab<- c("6 hour","24 hour")
tpc.plot$hr.lab <- hr.lab[match(tpc.plot$time.class, c(6,24))]
tpc.plot$hr.lab <- factor(tpc.plot$hr.lab, levels=c("6 hour","24 hour"), ordered=TRUE)

#rename time period
tpc.plot$time.per <- c("1999","2024")[match(tpc.plot$time.per, c("past","current"))]
tpc.plot$time.per <- factor(tpc.plot$time.per, levels=c("1999","2024"), ordered=TRUE)

#relative growth rate
tpc.plot$rgr=  (tpc.plot$fw/tpc.plot$Mo) / tpc.plot$time
tpc.plot$rgr2=  tpc.plot$mgain/tpc.plot$Mo/tpc.plot$time

#absolute growth rate
tpc.plot$agr= tpc.plot$mgain/tpc.plot$time

plot(tpc.plot$rgr, tpc.plot$agr)
plot(tpc.plot$gr, tpc.plot$agr) #gr is absolute growth rate

#select growth rate
tpc.plot$grow= tpc.plot$rgrlog
#tpc.plot$grow= tpc.plot$agr

#family counts
tpc.agg<- tpc.plot 
tpc.agg2 <- tpc.agg %>%
  group_by(temp, time.per, time.class, instar, in.lab, mom, hr.lab) %>% 
  dplyr::summarise(
    n= length(grow)
   )

tpc.agg2 <- tpc.agg %>%
  group_by(temp, time.per, time.class, instar, in.lab, hr.lab) %>% 
  dplyr::summarise(
    n= length(grow)
  )

#---------------------
#PLOT
#Figure 4

rgr.plot <- ggplot(tpc.plot[which(!is.na(tpc.plot$hr.lab)),], aes( x = temp, y = grow, color = time.per)) +
  geom_point(alpha=0.4, 
             position = position_jitterdodge(jitter.width = 5,
                                             jitter.height = 0,
                                             dodge.width = 0.75)) +
  facet_grid(hr.lab ~ in.lab, scales="free_y") 

#tpc.plot[which(tpc.plot$temp==41 & tpc.plot$time.per=="recent" & tpc.plot$instar==4),]

#---
#plot temp means
tpc.agg<- tpc.plot #[which(!is.na(tpc.plot$time.class)),]
tpc.agg <- tpc.agg %>%
  group_by(temp, time.per, time.class, instar, hr.lab, in.lab) %>% 
  dplyr::summarise(
    mean = mean(grow, na.rm = TRUE),
    n= length(grow),
    sd = sd(grow, na.rm = TRUE),
    mean.mass = mean(fw-Mo, na.rm = TRUE),
    sd.mass = sd(fw-Mo, na.rm = TRUE) )

tpc.agg$se= tpc.agg$sd / sqrt(tpc.agg$n)
tpc.agg$se.mass= tpc.agg$sd.mass / sqrt(tpc.agg$n)

#restrict to points 
tpc.agg<- tpc.agg[which(tpc.agg$n>5),]
tpc.agg.mass<- tpc.agg
tpc.agg<- tpc.agg[which(!is.na(tpc.agg$hr.lab)),]

#plotting temp means with error bars
Fig4A_growth.plot= rgr.plot + 
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(data=tpc.agg, aes(x=temp, y = mean, fill=time.per), size=3, col="black", pch=21)+
  geom_line(data=tpc.agg, aes(x=temp, y = mean))+
  theme_bw(base_size=16)+
  ylab("Relative growth rate (RGR, log10 mg/mg/h)")+ 
  #ylab("Absolute growth rate (AGR, mg/h)")+
  xlab("Temperature (°C)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Year", fill="Year")+theme(legend.position="none")+
  ylim(-0.015, 0.06)+
  geom_hline(yintercept=0, col="darkgray")

#add family lines
#not individuals at every temp
#rgr.plot= rgr.plot +
#geom_line(data=tpc.agg.f, aes(x=temp, y = mean, group=mom), linewidth=1)

#----
#Model
mod.lmer4 <- lme(gr ~ poly(temp,3)*time.per*time.class*Mo, random=~1|mom/ID, data = na.omit(tpc.plot[tpc.plot$instar==4,]))
anova(mod.lmer4)
mod.lmer5 <- lme(gr ~ poly(temp,3)*time.per*time.class*Mo, random=~1|mom/ID, data = na.omit(tpc.plot[tpc.plot$instar==5,]))
anova(mod.lmer5)

sigma(mod.lmer4)
sigma(mod.lmer5)

#Save anova
tables2<- rbind(as.data.frame(anova(mod.lmer4)), as.data.frame(anova(mod.lmer5)) )

colnames(tables2)[3:4]<- c("F","p")
tables2$sig<-""
tables2$sig[tables2$p<0.05]<-"*"
tables2$sig[tables2$p<0.01]<-"**"
tables2$sig[tables2$p<0.001]<-"***"
tables2$F= round(tables2$F,1)
tables2$p= round(tables2$p,4)

#Table 2
write.csv(tables2, "figures/Tables2_comp_growth.csv")

#compare across instars for text
mod.lmer45 <- lme(gr ~ poly(temp,3)*time.per*time.class*instar, random=~1|mom/ID, data = na.omit(tpc.plot))
anova(mod.lmer45)

#-------
#plot time periods together
Fig4B_rgrtime.plot <- ggplot(tpc.agg, aes( x = temp, y = mean, color = time.per, lty=factor(time.class))) +
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(size=2.5) + geom_line(linewidth=1.25)+
  facet_grid(.~ in.lab) +
  theme_bw(base_size=16)+xlab("Temperature (°C)")+
  #ylab("RGR (log10 mg/mg/h)")+
  ylab("AGR (mg/h)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Year", fill="Year", lty="Time (hr)")+theme(legend.position="bottom")+
  geom_hline(yintercept=0, col="darkgray")

#plot 4th and 5th together
rgr45.plot <- ggplot(tpc.agg, aes( x = temp, y = mean, color = time.per, lty=factor(instar))) +
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(size=2) + geom_line()+
  facet_grid(hr.lab ~ .) +
  theme_bw(base_size=16)+xlab("Temperature (°C)")+ylab("Growth rate (mg/mg/h)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Year", fill="Year", lty="Instar")+
  theme(legend.position="bottom")

#-------
#Mass Plot 
tpc.agg.mass$per.temp<- paste(tpc.agg.mass$time.per, tpc.agg.mass$temp, sep="_")
#restrict time
tpc.agg.mass<- tpc.agg.mass[which(!is.na(tpc.agg.mass$time.class)),]
#fix time 0 error bars
tpc.agg.mass$se.mass[which(tpc.agg.mass$time.class==0)]<- 0

FigS4_mass.plot= ggplot(data=tpc.agg.mass, aes(x=time.class, y = mean.mass, color=factor(temp), group=factor(per.temp), lty=time.per)) +
  geom_line(size=1.25)+
  geom_errorbar(data=tpc.agg.mass, aes(x=time.class, y=mean.mass, ymin=mean.mass-se.mass, ymax=mean.mass+se.mass), width=0, col="black")+
  geom_point(size=3)+
  theme_bw(base_size=16)+xlab("Time (h)")+ylab("Mass gain (mg)")+
  scale_color_viridis_d(option="plasma")+
  labs(color="Temperature (°C)", lty="Year")+theme(legend.position="bottom")+
  facet_wrap(.~in.lab, scales="free_y") +xlim(0,24)

#----------------
#Figure 3. distribution plots
tpc$time.per <- c("1999","2024")[match(tpc$time.per, c("past","current"))]
tpc$time.per <- factor(tpc$time.per, levels=c("1999","2024"), ordered=TRUE)

#initial weights
Fig3_plot.mass<- ggplot(tpc, aes(x=Mo,color=time.per, group=time.per)) + 
  geom_density(aes(fill=time.per), alpha=0.5, adjust=1.8)+
  ylab("Density") +xlab("Mass (mg)")+
  facet_wrap(.~in.lab, scales="free")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(2,6)])+
theme_bw(base_size=16) +theme(legend.position = c(0.9, 0.8))+
  labs(color="Year", fill="Year")

#compare distributions
#compare variance
leveneTest(Mo ~ time.per, tpc[tpc$instar==5,])
#compare means
t.test(Mo ~ time.per, data=tpc[tpc$instar==5,], alternative = "two.sided", var.equal = FALSE)
#unequal variance using Welch modification to the degrees of freedom

#Wilcoxon Rank Sum Tests
wilcox.test(Mo ~ time.per, data=tpc[tpc$instar==5,], alternative = "two.sided")

#------------
#write out plots
#save figures 
pdf("figures/Fig3_mass.pdf",height = 6, width = 8)
Fig3_plot.mass
dev.off()

design <- "AA
            AA
            AA
            AA
            BB
             BB"

pdf("figures/Fig4_relative_growthplot.pdf",height = 10, width = 7)
#pdf("figures/FigS3_absolute_growthplot.pdf",height = 10, width = 7)
Fig4A_growth.plot + Fig4B_rgrtime.plot +
  plot_layout(design=design)+plot_annotation(tag_levels = 'A')
dev.off()

#supplementary mass plot
pdf("figures/FigS4_masstime.pdf",height = 6, width = 8)
FigSx_mass.plot
dev.off()

#------------
#Analysis by dates

tpc.plot$doy.group <- NA
#1999: 209-228; 249-251; 298-308; 349-354
#2024: 122-137; 253-270
tpc.plot$doy.group[tpc.plot$Jdate %in% c(122:137)] <- 130
tpc.plot$doy.group[tpc.plot$Jdate %in% c(209:228)] <- 220
tpc.plot$doy.group[tpc.plot$Jdate %in% c(249:251)] <- 250
tpc.plot$doy.group[tpc.plot$Jdate %in% c(253:270)] <- 260
tpc.plot$doy.group[tpc.plot$Jdate %in% c(298:308)] <- 304
tpc.plot$doy.group[tpc.plot$Jdate %in% c(349:354)] <- 350

doy.plot <- ggplot(tpc.plot[which(tpc.plot$time.per==2024),], aes( x = temp, y = grow, color = factor(doy.group))) +
  geom_point(alpha=0.2, 
             position = position_jitterdodge(jitter.width = 7,
                                             jitter.height = 0,
                                             dodge.width = 0.75)) +
  facet_grid(hr.lab ~ in.lab, scales="free_y") 

#model with date
mod.lmer5 <- lme(gr ~ poly(temp,3)*Jdate*time.class*Mo, random=~1|mom/ID, data = na.omit(tpc.plot[which(tpc.plot$instar==5 & tpc.plot$time.per==2024),]))
anova(mod.lmer5)
