#Past feeding and temp data

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
library(zoo)
library(TrenchR)
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
  theme_bw(base_size=18) +theme(legend.position = "bottom")+
  xlab("Temperature (°C)")+ylab("Growth rate (g/g/hr)")+
  labs(color="Time (hr)")+scale_color_viridis_d()+
  xlim(11,41)

#make time numeric
tpc1$time.n <- as.numeric(tpc1$time.class)

#model
mod= lm(rgrlog ~ Mo + poly(temp)*time.class*instar, data= tpc1) 
#by instar
mod= lm(rgrlog ~ Mo + poly(temp)*time.class, data= tpc1[tpc1$instar==5,]) 
anova(mod)

plot_model(mod, type = "pred", terms = c("time.class", "temp"), show.data=TRUE)


#lme: issues accounting for individual ~1|UniID
mod.lmer <- lme(rgrlog ~ Mo + poly(temp)*time.class*instar, random=~1|mom/ID, data = na.omit(tpc1))
anova(mod.lmer)

#by instar
mod.lmer <- lme(rgrlog ~ Mo + poly(temp)*time.n, random=~1|mom/ID, data = na.omit(tpc1[tpc1$instar==5,]))
summary(mod.lmer)
anova(mod.lmer)
#4th: time.n            -0.00292623 0.000242389
#5th: time.n            -0.00305272 0.000212852

mod= lm(rgrlog ~ Mo + poly(temp)*time.class, data= tpc1[tpc1$instar==5,]) 

mod= lm(rgrlog ~ Mo + poly(temp)*time.class, data= tpc1[tpc1$instar==5,]) 
#----------------

#GAM
mod.gam <- gam(rgrlog ~ Mo + 
                 s(temp, bs = "cr", k=5, by = interaction(time.class, instar))+
               time.class*instar, #+s(UniID, bs = "re"),                               # Random intercept
               data = na.omit(tpc1),
               method = "REML")
summary(mod.gam)

#Kingsolver et al. notes
#Because growth during the 5th instar is approximately isometric (rather than allometric or exponential) we modeled mass on an arithmetic (rather than log) scale

#Time 0-24h of the 5th instar for all test temperatures, Period 2 represents Time 24-48h or Time 30-54h depending on the test temperature, and Period 3 

#We used linear mixed-effects models (function lme in R library nlme) to analyze larval mass (m). We considered time (t) and test temperature (T) as fixed effects; time was modeled as a 2nd order (P. rapae)
#Individual was included as a random (intercept) effect in the model.

#==================================================
#TEMP ANALYSIS

#OPERATIVE TEMPS
#toggle between desktop (y) and laptop (n)
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/Shared drives/TrEnCh/Projects/WARP/Projects/PrapaeGardenExpt/")
if(desktop=="n") setwd("/Users/lbuckley/Google Drive/Shared drives/TrEnCh/Projects/WARP/Projects/PrapaeGardenExpt/")

tdat<- read.csv("./data/PrapaeGardenTemps_WARP.csv")

#drop 2024 air temperature
tdat<- tdat[-which(tdat$T=="Logger3.T4.shadedT"),]

#doy
tdat$doy= floor(tdat$dt)

#average across sensors
tdat.mean <- tdat %>%
  group_by(dt, Date, Time, Year, hr, doy) %>%
  summarise(Tmean = mean(value, na.rm = TRUE), n=length(value))

#hourly max, min
tdat.mean.hr <- tdat.mean %>%
  group_by(Date, hr, Year, doy) %>%
  summarise(Tmin = min(Tmean, na.rm = TRUE), Tmax = max(Tmean, na.rm = TRUE), Tmean = mean(Tmean, na.rm = TRUE), n=length(Tmean), .groups = 'drop')

#check dates
#c(unique(tdat.mean.hr[which(tdat.mean.hr$Year==2023),"doy"]))
#1997: 224-238
#1999: 209-217, 227-237
#2023: 190-202, 212-226
#2024: 173-240

#restrict to overlapping time periods 227-237
tdat<- tdat.mean.hr[which(tdat.mean.hr$doy %in% c(227:237)), ]
#just 199 and 2024
tdat<- tdat[which(tdat$Year %in% c(1999,2024)), ]

#order by time
tdat<- tdat[order(tdat$Year, tdat$doy, tdat$hr), ]

#-----
#running average across time periods

# Calculate running averages for different time windows
tave <- tdat %>%
  # Calculate rolling means for different windows
  mutate(
    # 2-hour rolling average
    t2h = rollmean(Tmean, k = 2, fill = NA, align = "right"),
    # 6-hour rolling average
    t6h = rollmean(Tmean, k = 6, fill = NA, align = "right"),
    # 12-hour rolling average
    t12h = rollmean(Tmean, k = 12, fill = NA, align = "right"),
    # 24-hour rolling average
    t24h = rollmean(Tmean, k = 24, fill = NA, align = "right"),
    # 12-hour rolling average
    t54h = rollmean(Tmean, k = 54, fill = NA, align = "right"),
    # 24-hour rolling average
    t80h = rollmean(Tmean, k = 80, fill = NA, align = "right")
  )

#to long format
tave<- na.omit(melt(tave[,c("hr","Year","doy","t2h","t6h", "t12h","t24h","t54h","t80h")], id.vars=c("hr","Year","doy")))

#define time window
tave$hours<- gsub("t","",tave$variable)
tave$hours<- gsub("h","",tave$hours)
tave$hours<- factor(tave$hours, levels=c(2,6,12,24,54,80), ordered=TRUE)

#plot distributions
hr.plot.op<- ggplot(tave, aes(x=value, y=hours, color=factor(Year), fill=factor(Year) ))+
  geom_density_ridges(lwd=1.2, alpha=0.5)+
  scale_color_manual(values=cols2)+scale_fill_manual(values=cols2)+
  theme_bw(base_size=18) +theme(legend.position = "bottom")+
  xlab("Operative temperature (°C)")+ylab("Time average (hr)")+
  labs(color="Year", fill="Year") + 
  xlim(11,35)

tpc.agg$hours <- tpc.agg$time.class

# #plot distributions with tpcs on top
# ttplot<- ggplot()+
#   #add densities
#   geom_density(data= tave, aes(x=value, color=factor(Year), fill=factor(Year)))+
#   #add lines
#   geom_line(data= tpc.agg[which(tpc.agg$instar==5),], aes(x=temp, y=mean*10, color=factor(Year), group=hours))+
#   #facet
#   facet_wrap(.~hours, scales="free_y")

#analysis
tave$hours.n<- as.numeric(tave$hours)
mod<- lm(value~Year*hours.n, data=tave)
anova(mod)

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

#dtr over time
t.dat$dtr<- t.dat$tmax-t.dat$tmin
t.dtr <- t.dat %>%
  group_by(year) %>% 
  dplyr::summarise(dtrm = mean(dtr, na.rm = TRUE))
  
ggplot(t.dtr, aes(x=year, y=dtrm))+geom_line()+geom_smooth(method="lm") 
    
#-------

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
dtr=function(x, t=1:24){
  T_max= x[1]
  T_min= x[2]
  gamma= 0.44 - 0.46* sin(0.9 + pi/12 * t)+ 0.11 * sin(0.9 + 2 * pi/12 * t);   # (2.2) diurnal temperature function
  T = T_max*gamma + T_min - T_min*gamma
  return(T)
}

#hourly
temps=apply(t.dat[,c("tmax","tmin")], FUN="dtr", MARGIN=1)
temps= as.data.frame(t(temps))
temps$period= t.dat$period
temps$site= t.dat$site
temps$season= t.dat$season
temps$month= t.dat$month
temps$date= t.dat$DATE

#to long format
colnames(temps)[1:24]=1:24
thr<- na.omit(melt(temps, id.vars=c("period","site","month","date","season"), value.name ="Tmean",variable.name = "hr"))

#order by time
thr<- thr[order(thr$date, thr$hr), ]

# Calculate running averages for different time windows
tave <- thr %>%
  # Calculate rolling means for different windows
  mutate(
    # 2-hour rolling average
    t2h = rollmean(Tmean, k = 2, fill = NA, align = "right"),
    # 6-hour rolling average
    t6h = rollmean(Tmean, k = 6, fill = NA, align = "right"),
    # 12-hour rolling average
    t12h = rollmean(Tmean, k = 12, fill = NA, align = "right"),
    # 24-hour rolling average
    t24h = rollmean(Tmean, k = 24, fill = NA, align = "right"),
    # 12-hour rolling average
    t54h = rollmean(Tmean, k = 54, fill = NA, align = "right"),
    # 24-hour rolling average
    t80h = rollmean(Tmean, k = 80, fill = NA, align = "right")
  )

#to long format
tave<- na.omit(melt(tave[,c("period","hr","season","t2h","t6h", "t12h","t24h","t54h","t80h")], id.vars=c("hr","period","season")))

#define time window
tave$hours<- gsub("t","",tave$variable)
tave$hours<- gsub("h","",tave$hours)
tave$hours<- factor(tave$hours, levels=c(2,6,12,24,54,80), ordered=TRUE)

#plot distributions
hr.plot.ws<- ggplot(tave, aes(x=value, y=hours, color=factor(period), fill=factor(period) ))+
  geom_density_ridges(lwd=1.2, alpha=0.5)+
  #facet_wrap(.~season)+
  scale_color_manual(values=cols2)+scale_fill_manual(values=cols2)+
  theme_bw(base_size=18) +theme(legend.position = "bottom", axis.title.y=element_blank())+
  xlab("Environmental temperature (°C)")+ylab("")+ #ylab("Time average (hr)")+
  labs(color="Period", fill="Period") +
  xlim(0,35)

#analysis
tave$hours.n<- as.numeric(tave$hours)
mod<- lm(value~period*hours.n, data=tave)
anova(mod)

#------------
#write out plot

if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/WARP/projects/TPCconstant/")
if(desktop=="n") setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My Drive/Buckley/Work/WARP/projects/TPCconstant/")

design <- "AA
            AA
            BC
            BC
             BC"

#save figure 
pdf("./figures/Fig1.pdf",height = 8, width = 9)
Fig1_time.plot + hr.plot.op +hr.plot.ws +plot_layout(design = design)+plot_annotation(tag_levels = 'A')
dev.off()

#------------
#Test increasing incidence of warm temperatures
#make variable whether temperature hot
mod1<- lm(value~hours*period, data=tave)

tave$o30<- 0
tave$o30[tave$value>=30]<- 1

mod1 <- glm(o30~hours*period, family=binomial, data=tave) 

summary(mod1)
anova(mod1, test="Chisq")







