# Taylor M. Hatcher
# Analysing data for constant tpc data sets from 1999 and 2024

# Install necessary packages
install.packages("lme4")

# Load libraries
library("lme4")

# Fit the linear mixed-effects model for constant tpc data

tpcvismodel <- read.csv("PastPresentFilteredConstantTpc2024.csv")

tpcvismodel_filtered<- tpcvismodel %>%
  filter(rgrlog > 0 & !is.na(rgrlog))



historicdata <- lmer(
  formula = rgrlog ~ poly(dur, 2) + factor(temp) + (1 | mom),
  data = tpcvismodel_filtered,
  REML = FALSE
)

# Summary of the model 
summary(model)

# Check model 
plot(model)