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
library(TrenchR)

colm<- viridis_pal(option = "mako")(8)
cols<- colm[c(2,4,7)]
cols2<- colm[c(3,6)]

#toggle between desktop (y) and laptop (n)
desktop<- "n"

if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/WARP/projects/TPCconstant/Data/")
if(desktop=="n") setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My Drive/Buckley/Work/WARP/projects/TPCconstant/Data/")

#### read in past raw data 
tpc1 = read.csv("PrapaeW.1999.ConstantTempTPCs.4thinstar.jul2021.xlsx - data.csv")
tpc1$instar=4 # identify instar

#assemble data
# filter out and include only caterpillars who were active
  tpc1 <- tpc1 %>% filter(active == "yes")

tpc2 <- read.csv("PrapaeW.1999.ConstantTempTPCs.5thinstar.jul2021.xlsx - data.csv")
names(tpc2) <- names(tpc1)
tpc2$instar <- 5 # Identify instar

# Combine past datasets using rbind
tpc.p <- rbind(tpc1, tpc2)
tpc.p$time.per = "past"

#time: time since start of experiment (h)
#duration: duration of trial (h)
#M0	mass at start of experiment (time = 0)	mg
#Mi	mass at start of current trial	mg
#iw	initial mass (start of trial)	mg
#fw	final mass (end of trial)	mg
#mgain	mass gain during trial	mg

#----------------------
#### load in recent 2024  Constant TPC data
tpc.c = read.csv("2024PrapaeConstantTPCCombineddata.csv", skip = 1)

#drop dead
tpc.c<- tpc.c[which(!is.na(tpc.c$M0)),]
tpc.c<- tpc.c[-which(tpc.c$fw=="dead"),]
tpc.c<- tpc.c[-which(tpc.c$fw==""),]
## CHECK HUGE FW

#fix format issues
tpc.c$Date[which(tpc.c$Date=="9/11+434:469/2024")] <- "11-Sep-24"
tpc.c$Date[which(tpc.c$Date=="15-May-23")] <- "15-May-24"

tpc.c$fw[which(tpc.c$fw=="14:51")] <- 14.51
tpc.c$fw[which(tpc.c$fw=="14:34")] <- 14.34
tpc.c$fw[which(tpc.c$fw=="11:05")] <- 11.05
tpc.c$fw<- as.numeric(tpc.c$fw)

#Calculate mass gained in each time treatment
tpc.c$mgain= tpc.c$fw - tpc.c$M0
#FIX: CHECK big negative and positive values

# Paste date in t.in and t.out column to make duration calculation easier
tpc.c$t.in <- paste(tpc.c$Date, tpc.c$t.in, sep = " ")
tpc.c$t.out <- paste(tpc.c$Date, tpc.c$t.out, sep = " ")

# Paste mom and individual together to create UniID
tpc.c <- tpc.c %>%
  mutate(UniID = paste(Female, Individual, sep = " "))

# Convert to POSIXct with correct formatting
tpc.c <- tpc.c %>%
  mutate(
    t.in = as.POSIXct(t.in, format = "%d-%b-%y %H:%M", tz = "UTC"),
    t.out = as.POSIXct(t.out, format = "%d-%b-%y %H:%M", tz = "UTC")
  )%>%
  group_by(UniID)%>% 
  mutate(
    first_t_out = first(t.out),
    duration = as.numeric(difftime(t.out, t.in, units = "hours")),
    time = as.numeric(difftime(t.out, first_t_out, units = "hours")))

#FIX
#CHECK time NAs
#CHECK current times that deviate from 6 and 24 hrs
tpc.c[which(tpc.c$duration>30),]

# Make sure that new data follows naming of old data sets
tpc.c$mom= tpc.c$Female
tpc.c$ID= tpc.c$Individual

tpc.c$time.per = "current"
tpc.c$Mo = tpc.c$M0

#add doy of recent
tpc.c$Jdate= gsub("May","5", tpc.c$Date)
tpc.c$Jdate= gsub("Sep","9", tpc.c$Jdate)
tpc.c$Jdate= day_of_year(tpc.c$Jdate, format = "%d-%m-%y")

#----------
#combine historic and current
tpc.ps= tpc.p[,c("UniID","mom","ID","temp", "active","instar","time","duration","Mo","fw","time.per","mgain","Jdate")]
tpc.cs= tpc.c[,c("UniID","mom","ID","temp", "active","instar","time","duration","Mo","fw","time.per","mgain","Jdate")]

# Ensure data types match for both datasets before combining
tpc.ps$mom <- as.character(tpc.ps$mom)

# Combine past and current datasets while preserving time.per
tpc <- rbind(tpc.cs, tpc.ps) 

#align active
tpc$active[which(tpc$active %in% c("yes","y?"))]<- "y"
tpc$active[which(tpc$active %in% c("no","n?"))]<- "n"

#==========================
#analysis
# Filter caterpillars to include the ones who were active
tpc <- tpc %>% filter(active %in% c("y"))

#drop initial zero time estimates
#tpc<- tpc[which(tpc$time>0),]

# calculate relative growth rate using logarithmic scale for 1999 past data set
#Kingsolver and Gomulkiewicz 2003
#“We have quantified short-term (6 hr) growth rates of P. rapae caterpillars at different temperatures in terms of the relative growth rate (RGR), defined as RGR =[ln(mf / mi )]/(tf - ti ), where mi is initial mass at time ti , and mf is final mass at time tf . Figure 1 shows short-term RGR as a function of temperature for 4th and 5th instar P. rapae from Seattle WA (Kingsolver, 2000).”
#paper plot matches log10 not ln
#tpc$rgrlog= log10(tpc$fw/tpc$Mo)/tpc$time 
#OR with conversion to mg: 
tpc$rgrlog= (log10(tpc$fw)-log10(tpc$Mo))/tpc$time 

# calculate relative growth rate using arithmetic scale 
#tpc$rgrarith = (tpc$fw/tpc$Mo) / tpc$time 

#growth rate
tpc$gr= tpc$mgain/tpc$time 

# Save data frame to new Csv
if(desktop=="y") setwd("/Users/laurenbuckley/Google Drive/My Drive/Buckley/Work/WARP/projects/TPCconstant/out/")
if(desktop=="n") setwd("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My Drive/Buckley/Work/WARP/projects/TPCconstant/out/")

write.csv(tpc, "PastPresentFilteredConstantTpc2024.csv")
