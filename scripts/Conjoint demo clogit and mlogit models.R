# Load packages
library(survival)
library(mlogit)

# Load 'Travel Mode Choice Data' example from AER package
data("TravelMode", package = "AER")

# Fit model conditional logistic regression model
# and multinomial logit model
travel.mod1 <- clogit(as.numeric(choice)-1 ~ wait + gcost + mode + strata(individual), data = TravelMode)
travel.mod2 <- mlogit(choice ~ wait + gcost, data = TravelMode, choice = "choice", shape = "long", id.var = "individual", alt.var = "mode")

# Gives the same result
# We can use clogit or mlogit to fit 
summary(travel.mod1)
summary(travel.mod2)
