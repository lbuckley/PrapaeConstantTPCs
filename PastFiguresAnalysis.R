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
library(ggridges)

colm<- viridis_pal(option = "mako")(8)
cols<- colm[c(2,4,7)]
cols2<- colm[c(3,6)]

#toggle between desktop (y) and laptop (n)
desktop<- "n"

# Load data
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/WARP/projects/TPCconstant/out/")
if(desktop=="n") setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My Drive/Buckley/Work/WARP/projects/TPCconstant/out/")

tpc<- read.csv("PastPresentFilteredConstantTpc2024.csv")

#==========================
#INITIAL DATA

#rgr by time
#4th and 5th instar

tpc1<- tpc[which(tpc$time.per=="past"),]
#restrict to less that 55 hours
#tpc1<- tpc1[which(tpc1$time<55),]

#time classes
ggplot(tpc1, aes(x = time)) + geom_density()

#set up time classes
tpc1$time.class<- NA
tpc1$time.class[which(tpc1$time<90)]<- 80
tpc1$time.class[which(tpc1$time<65)]<- 54
tpc1$time.class[which(tpc1$time<40)]<- 24
tpc1$time.class[which(tpc1$time<15)]<- 6
tpc1$time.class<- factor(tpc1$time.class, levels=c(6, 24, 54, 80), ordered=TRUE)

#specify metric
tpc1$grow<- tpc1$rgrlog
#tpc1$grow<- tpc1$gr

#mean across time classes
tpc.agg <- tpc1 %>%
  group_by(temp, time.per, time.class, instar) %>% 
  dplyr::summarise(
    mean = mean(grow, na.rm = TRUE),
    n= length(grow),
    sd = sd(grow, na.rm = TRUE) )
tpc.agg$se= tpc.agg$sd / sqrt(tpc.agg$n)

#restrict to points with 6 measurements
tpc.agg<- tpc.agg[-which(tpc.agg$n<6),]

#instar label
in.lab<- c("4th instar", "5th instar")
tpc.agg$in.lab <- in.lab[match(tpc.agg$instar, c(4,5))]

#plot
tpc.agg$instar= factor(tpc.agg$instar)

Fig1_time.plot<- ggplot(tpc.agg, aes(x = temp, y = mean, color = time.class)) + 
  geom_point() + geom_line()+
  geom_errorbar(aes(x=temp, y=mean, ymin=mean-se, ymax=mean+se), width=0, col="black")+
  facet_grid(. ~ in.lab) +
  theme_bw(base_size=16) +theme(legend.position = "bottom")+
  xlab("Temperature (°C)")+ylab("Growth rate (g/g/hr)")+
  labs(color="Time (hr)")+scale_color_viridis_d()

#model
mod= lm(rgrlog ~ Mo + poly(temp)*time.class*instar, data= tpc1) 
#by instar
mod= lm(rgrlog ~ Mo + poly(temp)*time.class, data= tpc1[tpc1$instar==5,]) 
anova(mod)

plot_model(mod, type = "pred", terms = c("time.class", "temp"), show.data=TRUE)

#lme: issues accounting for individual
mod.lmer <- lme(rgrlog ~ Mo + poly(temp)*time.class*instar, random=~1|UniID, data = na.omit(tpc1))
anova(mod)

#Kingsolver et al. notes
#Because growth during the 5th instar is approximately isometric (rather than allometric or exponential) we modeled mass on an arithmetic (rather than log) scale

#Time 0-24h of the 5th instar for all test temperatures, Period 2 represents Time 24-48h or Time 30-54h depending on the test temperature, and Period 3 

#We used linear mixed-effects models (function lme in R library nlme) to analyze larval mass (m). We considered time (t) and test temperature (T) as fixed effects; time was modeled as a 2nd order (P. rapae)
#Individual was included as a random (intercept) effect in the model.

#==================================================
#TEMP ANALYSIS

#OPERATIVE TEMPS
#toggle between desktop (y) and laptop (n)
desktop<- "n"
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared drives/TrEnCh/Projects/WARP/Projects/PrapaeGardenExpt/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared drives/TrEnCh/Projects/WARP/Projects/PrapaeGardenExpt/")

tdat<- read.csv("./data/PrapaeGardenTemps_WARP.csv")

#drop 2024 air temperature
tdat<- tdat[-which(tdat$T=="Logger3.T4.shadedT"),]

#average across sensors
tdat.mean <- tdat %>%
  group_by(dt, Date, Time, Year) %>%
  summarise(Tmean = mean(value, na.rm = TRUE), n=length(value))

#running average across time periods
#1,2,6,12, 24 hr 
#ggrides

#-----
#plot distributions for overlapping data, #1997: 224-238; 1999: 227:237
Tdist.plot <- ggplot(tdat.mean[which(tdat.mean$dt>223 & tdat.mean$dt<238),], aes(x=Tmean, color=factor(Year), fill=factor(Year), group=factor(Year))) +  geom_density(alpha=0.5)+
  scale_fill_viridis_d() +scale_color_viridis_d()+
  theme_classic(base_size = 18)+
  labs(x = "Temperature (°C)", color = "Year", fill="Year") 

#restrict to daylight 223 238 Aug 11-26; Aug 11: 6,8:30; Aug 26: 6:20-8
#tdat.mean$dec.dt<- tdat.mean$dt - floor(tdat.mean$dt)
#tdat.day<- tdat.mean[which(tdat.mean$dec.dt>0.25 & tdat.mean$dec.dt<0.85),]
#tdat.day<- tdat.day[which(tdat.day$dt>223 & tdat.day$dt<238),]

Tdist.day.plot <- ggplot(tdat.mean, aes(x=Tmean, color=factor(Year), fill=factor(Year), group=factor(Year))) +  geom_density(alpha=0.5)+
  scale_fill_viridis_d() +scale_color_viridis_d()+
  theme_classic(base_size = 18)+ xlim(0,40)+
  labs(title="doy 223-238, daytime" , x = "Temperature (°C)", color = "Year", fill="Year") 
#hourly max, min
tdat.mean.hr <- tdat %>%
  group_by(Date, hr, Year) %>%
  summarise(Tmin = min(value, na.rm = TRUE), Tmax = max(value, na.rm = TRUE), Tmean = mean(value, na.rm = TRUE), n=length(value), .groups = 'drop')

Tdist.hr.plot <- ggplot(tdat.mean.hr, aes(x=Tmax, color=factor(Year), fill=factor(Year), group=factor(Year))) +  geom_density(alpha=0.5)+
  scale_fill_viridis_d() +scale_color_viridis_d()

#------------------------
#ENVI TEMPS

#GHCND date
#find stations: https://ncics.org/portfolio/monitor/ghcn-d-station-data/ 
# SEATTLE SAND POINT WEATHER FORECAST OFFICE, WA US (USW00094290)

if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared drives/TrEnCh/Projects/WARP/Projects/PrapaeGardenExpt/data/GHCNdata/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared drives/TrEnCh/Projects/WARP/Projects/PrapaeGardenExpt/data/GHCNdata/")

t.dat<- read.csv("USW00094290_2025.csv")
t.dat$site="Seattle"

t.dat$tmin= t.dat$TMIN /10 #divide by ten
t.dat$tmax= t.dat$TMAX /10

t.dat$month= round(month(as.POSIXlt(t.dat$DATE)))
t.dat$year= year(as.POSIXlt(t.dat$DATE))
#restrict to growing season
t.dat= t.dat[which(t.dat$month %in% c(4,5,6,7,8,9)),] 

#code season
t.dat$season<- NA
t.dat$season[which(t.dat$month %in% c(4:6))] ="spring"
t.dat$season[which(t.dat$month %in% c(7:9))] ="summer"

#restrict years
t.dat= t.dat[which(t.dat$year %in% c(1994:2024)),]
t.dat1= t.dat[which(t.dat$year %in% c(1990:1999)),]
t.dat1$period="initial"
t.dat2= t.dat[which(t.dat$year %in% c(2015:2024)),]
t.dat2$period="recent"
#combine
t.dat= rbind(t.dat1,t.dat2)

#dtr
dtr=function(T_max, T_min, t=7:18){
  gamma= 0.44 - 0.46* sin(0.9 + pi/12 * t)+ 0.11 * sin(0.9 + 2 * pi/12 * t);   # (2.2) diurnal temperature function
  T = T_max*gamma + T_min - T_min*gamma
  return(T)
}

#hourly?
temps= sapply(t.dat$tmax, FUN="dtr", T_min=t.dat$tmin)
temps= as.data.frame(t(temps))
temps$period= t.dat$period
temps$site= t.dat$site
#temps$season= t.dat$season
temps$month= t.dat$month

#to long format
temps1<- melt(temps, id.vars=c("period","site","month"))

#density plot
month.plot<- ggplot(temps1, aes(x=value, y=month, color=factor(month), fill=factor(month), lty=period))+
  geom_density_ridges(lwd=1.2, alpha=0.5)+
  scale_color_viridis_d(option="mako")+scale_fill_viridis_d(option="mako")+
  #scale_color_manual(values=rev(cols))+scale_fill_manual(values=rev(cols))+
  xlim(0,42)+
  #ylim(5.9, 9.5)+
  xlab("Temperature (°C)") +ylab("Month")+ 
  guides(fill="none", color="none")+ labs(lty="Period")+
  theme_classic(base_size = 18)+theme(legend.position = c(0.9, 0.8))

#------------
#Test increasing incidence of warm temperatures
#make variable whether temperature hot
mod1<- lm(value~month*period, data=temps1)

temps1$o30<- 0
temps1$o30[temps1$temp>=25]<- 1

mod1 <- glm(o30~month+period, family=binomial, data=temps1) 

summary(mod1)
anova(mod1)







