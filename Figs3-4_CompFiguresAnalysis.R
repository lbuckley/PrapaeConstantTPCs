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
library(ggeffects)

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

n_families <- tpc.agg %>%
  group_by(time.per, instar, in.lab) %>%
  summarise(n_moms = n_distinct(mom), .groups = "drop")
# 1999: 24 4th and 21 5th moms
# 2024: 40 4th and 37 5th moms

n_ind <- tpc.agg %>%
  group_by(time.per, instar, in.lab) %>%
  summarise(n_ind = n_distinct(UniID), .groups = "drop")
# 1999: 106 4th and 127 5th ind
# 2024: 800 4th and 561 5th ind

#dates
n_doy <- tpc.agg[tpc.agg$time.class==6,] %>%
  group_by(time.per, instar, in.lab, Jdate) %>%
  summarise(doy = unique(Jdate))

n_doy[n_doy$time.per==2024 & n_doy$instar==5,]$doy

# 1999: 4th Aug-Sep on 216, 224, 225, 249; 5th Nov-Dec on 298, 300, 305, 349, 351 
# 2024: 4th on 122 (May 2) 123 135 253 254 255 263 264 265 266; 5th on 123 129 135 136 254 255 256 264 265 266 268 269
#some May, mostly Sept

#range number of individuals per family
range_ind <- tpc.agg %>%
  group_by(time.per, instar, in.lab, temp, mom) %>%
  summarise(n_ind = n_distinct(UniID), 
            .groups = "drop")

range_ind <- range_ind %>%
  group_by(time.per) %>%
  summarise(mean_ind= mean(n_ind),
            med_ind= median(n_ind),
            min_ind= min(n_ind),
            max_ind= max(n_ind),
            .groups = "drop")            

#per treatment
range_ind <- tpc.agg %>%
  group_by(time.per, instar, in.lab, temp) %>%
  summarise(n_ind = n_distinct(UniID), 
            .groups = "drop")

range_ind <- range_ind %>%
  group_by(time.per) %>%
  summarise(mean_ind= mean(n_ind),
            med_ind= median(n_ind),
            min_ind= min(n_ind),
            max_ind= max(n_ind),
            .groups = "drop")   

tpc.agg2 <- tpc.agg %>%
  group_by(temp, time.per, time.class, instar, in.lab, hr.lab) %>% 
  dplyr::summarise(
    n= length(grow)
  )

#sample sizes per family and treatment, perhaps in a supplementary table.

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

# Plot effects to examine interactions
preds1 <- ggpredict(mod.lmer4, terms = c("temp [all]", "Mo","time.per"))
preds2 <- ggpredict(mod.lmer4, terms = c("temp [all]", "time.class", "Mo"))
preds3 <- ggpredict(mod.lmer4, terms = c("temp [all]", "Mo"))
preds4 <- ggpredict(mod.lmer4, terms = c("Mo [all]","time.per"))


preds5.5 <- ggpredict(mod.lmer5, terms = c("temp [all]","time.class", "time.per"))
preds2.5 <- ggpredict(mod.lmer5, terms = c("temp [all]","time.class", "Mo"))
preds4.5 <- ggpredict(mod.lmer5, terms = c("Mo [all]","time.per"))

#4th instar
plot.eff1<- plot(preds1) +
  labs(
    title  = "Interaction of temperature, mass, and year",
    x      = "Temperature (°C)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Initial mass (mg)"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d()

plot.eff2<- plot(preds2) +
  labs(
    title  = "Interaction of temperature, mass, and duration",
    x      = "Temperature (°C)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Duration (h)"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d() #option = "turbo"

plot.eff3<- plot(preds3) +
  labs(
    title  = "Interaction of temperature and mass",
    x      = "Temperature (°C)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Initial mass (mg)"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d()

plot.eff4<- plot(preds4) +
  labs(
    title  = "Interaction of mass and year",
    x      = "Initial mass (mg)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Year"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d()

#5th instar
plot.eff5.5<- plot(preds5.5) +
  labs(
    title  = "Interaction of temperature, duration, and year",
    x      = "Temperature (°C)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Duration (h)"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d()

plot.eff2.5<- plot(preds2.5) +
  labs(
    title  = "Interaction of temperature, mass, and duration",
    x      = "Temperature (°C)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Duration (h)"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d() #option = "turbo"

plot.eff4.5<- plot(preds4.5) +
  labs(
    title  = "Interaction of mass and year",
    x      = "Initial mass (mg)",
    y      = "Predicted relative growth rate (RGR, log10 mg/mg/h)",
    colour = "Year"   
  )+
  theme(axis.title = element_text(size = 14)) +
  scale_color_viridis_d()+scale_fill_viridis_d()

#-----------
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

#plot effects
preds <- ggpredict(mod.lmer45, terms = c("temp [all]", "time.per","time.class"))
preds <- ggpredict(mod.lmer45, terms = c("temp [all]", "time.per","instar"))
plot(preds)

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

#plot_model(rgr45.plot, type = "pred", terms = c("grow_t6hour","Mo"))

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
#Figure 3. mass distribution plots
tpc$time.per <- c("1999","2024")[match(tpc$time.per, c("past","current"))]
tpc$time.per <- factor(tpc$time.per, levels=c("1999","2024"), ordered=TRUE)

#initial weights
Fig3_plot.mass.dist<- ggplot(tpc, aes(x=Mo,color=time.per, group=time.per)) + 
  geom_density(aes(fill=time.per), alpha=0.5, adjust=1.8)+
  ylab("Density") +xlab("Mass (mg)")+
  facet_wrap(.~in.lab, scales="free")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(2,6)])+
theme_bw(base_size=16) +theme(legend.position = c(0.9, 0.8))+
  labs(color="Year", fill="Year")

#---
#Plot as box plots

Fig3A_plot.mass<- ggplot(tpc, aes(x=time.per, y=Mo,color=time.per, group=time.per)) + 
  geom_violin(quantile.linetype="solid")+
#  geom_density(aes(fill=time.per), alpha=0.5, adjust=1.8)+
  ylab("Mass (mg)") +xlab("Year")+
  facet_grid(in.lab~., scales="free", switch = "y")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(2,6)])+
  theme_bw(base_size=16) +theme(legend.position = c(0.7, 0.8))+
  theme(legend.position = "none")

#means and percent change
tpc %>%
  group_by(time.per, instar) %>%
  summarise(mean_value = mean(Mo, na.rm = TRUE),
            sd_value   = sd(Mo, na.rm = TRUE),
            n          = sum(!is.na(Mo)),
            se_value   = sd_value / sqrt(n),
            .groups    = "drop"
            )

round((1-11.9/12.9)*100, 1)
round((1-48.1/49.2)*100, 1)
#--------
#compare distributions
#compare variance
leveneTest(Mo ~ time.per, tpc[tpc$instar==5,])
#compare means
t.test(Mo ~ time.per, data=tpc[tpc$instar==5,], alternative = "two.sided", var.equal = FALSE)
#unequal variance using Welch modification to the degrees of freedom

#Wilcoxon Rank Sum Tests
wilcox.test(Mo ~ time.per, data=tpc[tpc$instar==5,], alternative = "two.sided")

#models accounting for family
#add Jdate
mass.lmer4 <- lme(Mo ~ time.per, random=~1|mom/ID, data=tpc[tpc$instar==4,])
anova(mass.lmer4)

mass.lmer5 <- lme(Mo ~ time.per, random=~1|mom/ID, data=tpc[tpc$instar==5,])
anova(mass.lmer5)

sigma(mod.lmer4)
sigma(mod.lmer5)

#plot_model(mass.lmer5, type = "pred", terms = c("time.per", "Jdate"))

#--------
#Growth at 6 vs 24 hrs
#Check role of mass

tpc_wide<- tpc.plot[!is.na(tpc.plot$hr.lab),]
#subtract one from 24hr measures
tpc_wide$hr.lab<- gsub(" ","",tpc_wide$hr.lab)
tpc_wide[tpc_wide$hr.lab=="24hour","Jdate"]<- tpc_wide[tpc_wide$hr.lab=="24hour","Jdate"]-1


tpc_wide <- tpc_wide %>%
  pivot_wider(
    id_cols = c(UniID, mom, ID, Jdate, Mo, temp, time.per, instar, in.lab),
    names_from  = hr.lab,
    values_from = grow,
    names_prefix = "grow_t"
  )

#just intermediate temperatures
tpc_wide<- tpc_wide[tpc_wide$temp %in% c(17, 23, 29, 35),]

# [tpc_wide$instar==4,]
Fig3B_plot.rgrmass_4th <- ggplot(tpc_wide[tpc_wide$instar==4,], aes( x = grow_t6hour, y = grow_t24hour, fill = Mo, color=time.per)) +
  geom_point(shape=21, alpha=0.4, size=3, stroke = 0.7) +
  #facet_grid(temp ~ ., scales="free_y") +
  scale_fill_viridis_c()+
  scale_color_manual(values=cols2)+
  ylab("24 hour Growth rate") +xlab("6 hour Growth rate (mg/mg/h)")+
  theme_bw(base_size=16) + #theme(legend.position = c(0.7, 0.8))+
  labs(fill="Mass (mg)", color="Year")+ #+xlim(-0.01, 0.07)+ylim(-0.01, 0.07)
  geom_abline(color="darkgray") #+annotate("text", label = "4th instar", x = -0.002, y = 0.035)

Fig3B_plot.rgrmass_5th <- ggplot(tpc_wide[tpc_wide$instar==5,], aes( x = grow_t6hour, y = grow_t24hour, fill = Mo, color=time.per)) +
  geom_point(shape=21, alpha=0.4, size=3, stroke = 0.7) +
  #facet_grid(in.lab ~ ., scales="free_y") +
  scale_fill_viridis_c()+
  scale_color_manual(values=cols2)+
  ylab("24 hour Growth rate") +xlab("6 hour Growth rate (mg/mg/h)")+
  theme_bw(base_size=16) + #theme(legend.position = c(0.7, 0.8))+
  labs(fill="Mass (mg)", color="Year")+ #+xlim(-0.01, 0.07)+ylim(-0.01, 0.07)
  geom_abline(color="darkgray")+ #+annotate("text", label = "5th instar", x = -0.005, y = 0.028)+
  guides(color = "none")

#check role of initial mass
mod.lmer4 <- lme(grow_t24hour ~ grow_t6hour*time.per*Mo, random=~1|mom/ID, data = na.omit(tpc_wide[tpc_wide$instar==4,]))
anova(mod.lmer4)

mod.lmer5 <- lme(grow_t24hour ~ grow_t6hour*time.per*Mo, random=~1|mom/ID, data = na.omit(tpc_wide[tpc_wide$instar==5,]))
anova(mod.lmer5)

plot_model(mod.lmer4, type = "pred", terms = c("grow_t6hour","Mo"))

#correlations
cor(tpc_wide[tpc_wide$instar==4,]$grow_t6hour, tpc_wide[tpc_wide$instar==4,]$grow_t24hour, use="complete.obs")
cor(tpc_wide[tpc_wide$instar==5,]$grow_t6hour, tpc_wide[tpc_wide$instar==5,]$grow_t24hour, use="complete.obs")

#------------
#write out plots
#save figures 

design <- "AABB
            AACC"

pdf("figures/Fig3_mass.pdf",height = 8, width = 10)
Fig3A_plot.mass +Fig3B_plot.rgrmass_4th +Fig3B_plot.rgrmass_5th+
  plot_layout(design=design)+plot_annotation(tag_levels = 'A')
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
#pdf("figures/FigS4_masstime.pdf",height = 6, width = 8)
#FigSx_mass.plot
#dev.off()

#supplementary effect plot
pdf("figures/FigSx_model4th.pdf",height = 10, width = 10)
plot.eff1 +plot.eff2 +plot.eff3 +plot.eff4
dev.off()

pdf("figures/FigSx_model5th.pdf",height = 10, width = 10)
plot.eff2.5 / (plot.eff5.5 +plot.eff4.5)
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

# Generate predictions at representative Jdate values
preds <- ggpredict(mod.lmer5, terms = c("temp [all]", "Jdate [meansd]","time.class"))
plot(preds)

#------------
#Plot as in Kingsolver 2000 PBZ, growth rate (g/g/h) 

tpc.plot$rgr=  (tpc.plot$fw/tpc.plot$Mo) / tpc.plot$time

#plot temp means
tpc.agg2<- tpc.plot #[which(!is.na(tpc.plot$time.class)),]
tpc.agg2 <- tpc.agg2 %>%
  group_by(temp, time.per, time.class, instar, hr.lab, in.lab) %>% 
  dplyr::summarise(
    mean = mean(rgr, na.rm = TRUE),
    n= length(rgr),
    sd = sd(rgr, na.rm = TRUE),
    mean.mass = mean(fw-Mo, na.rm = TRUE),
    sd.mass = sd(fw-Mo, na.rm = TRUE) )

tpc.agg2$se= tpc.agg2$sd / sqrt(tpc.agg2$n)
tpc.agg2$se.mass= tpc.agg2$sd.mass / sqrt(tpc.agg2$n)

#restrict to points 
tpc.agg2<- tpc.agg2[which(tpc.agg2$n>5),]
tpc.agg2<- tpc.agg[which(!is.na(tpc.agg2$hr.lab)),]

#plot 4th and 5th together
rgr45.plot <- ggplot(tpc.agg, aes( x = temp, y = mean, color = time.per, lty=factor(instar))) +
  geom_errorbar(data=tpc.agg, aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  geom_point(size=2) + geom_line()+
  facet_grid(hr.lab ~ .) +
  theme_bw(base_size=16)+xlab("Temperature (°C)")+ylab("Growth rate (mg/mg/h)")+
  scale_color_manual(values=cols2)+scale_fill_manual(values=colm[c(4,7)])+
  labs(color="Year", fill="Year", lty="Instar")+
  theme(legend.position="bottom")

pdf("figures/FigSx_ComparePBZ2000.pdf",height = 6, width = 6)
rgr45.plot
dev.off()

