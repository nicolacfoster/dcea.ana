

# IMPUTATION AND ANALYSES MODELS
# Author: Dr Nicola Foster, UCL, LSHTM
# Project: ASCENT DCEA

install.packages("pan")
install.packages("ICC")

#load packages
require(mice)
require(lattice)
require(pan)
require(purrr)
require(missMethods)
require(ggplot2)
require(GGally)
require(tidyverse)
require(VIM)
require(mice)
require(miceadds)
require(dplyr)
require(foreign)
require(survey)
require(broom.mixed)
require(multilevel)
require(haven)
require(labelled)
require(coda)
require(R2OpenBUGS)
require(R2jags)
require(brms)
require(survey)
require(ggthemes)
require(extrafont)
require(remotes)
require(sandwich)
require(ggstatsplot)
require(howManyImputations)


# SETS UP THE ANALYTICAL ENVIRONMENT
###############################################################################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/project_ASCENT/DATA-repository/FINAL")

set.seed(42)
data1 <- read_dta("ASCENT_DCEA_2024_1213.dta")
attach(data1)

# display different missing data patterns. Look at the number of different missing data patterns in the dataframe
# and which patterns occur most frequently in the data
md.pattern(data1)
# develop and evaluate a missingness parameter to explore the relationship between missingness and other variables
histogram(~ DALY_miss | is.na(DALY_miss), data=data1)

# look at the intraclass correlation for incomplete variables SEP_miss, DALY_miss and csttotaptient_miss
ICC1(aov(DALY_miss~facilityID))
#ICC1(aov(as.factor(SEP_miss)~facilityID))
ICC1(aov(csttotaptient_miss~facilityID))


# IMPUTATION WITHOUT CLUSTERING
###############################################################################

ini <- mice(data1, maxit=0)
meth <- ini$meth
meth

meth[c(3, 5, 6, 7)] <- "norm"
meth

pred <- ini$pred
pred

pred[, "DALY_miss"] <- 0
pred[, "csttotaptient_miss"] <- 0
pred

imp1 <- mice(data1, meth=meth, pred=pred, print=FALSE)

# summary of the imputed dataset
summary(complete(imp1))
summary(data1)

# comparing the ICCs of the variables first imputed with those in the incomplete dataset
data.frame(vars = names(data1[c(7,12,13)]),
           observed = c(ICC1(aov(DALY_miss~facilityID, data1)),
                        ICC1(aov(SEP_miss~facilityID, data1)),
                        ICC1(aov(csttotaptient_miss~facilityID, data1))),
           norm = c(ICC1(aov(DALY_miss~facilityID, complete(imp1))),
                    ICC1(aov(SEP_miss~facilityID, complete(imp1))),
                    ICC1(aov(csttotaptient_miss~facilityID, complete(imp1)))))



# IMPUTATION WITH CLUSTERING
###############################################################################
# now do the imputation again but with the cluster variable

pred[, "facilityID"] <- 0

imp2 <- mice(data1, meth=meth, pred=pred, print=FALSE)

data.frame(vars=names(data1[c(7,12,13)]),
           observed=c(ICC1(aov(DALY_miss~facilityID, data1)),
                      ICC1(aov(SEP_miss~facilityID, data1)),
                      ICC1(aov(csttotaptient_miss~facilityID, data1))),
           norm=c(ICC1(aov(DALY_miss~facilityID, complete(imp1))),
                  ICC1(aov(SEP_miss~facilityID, complete(imp1))),
                  ICC1(aov(csttotaptient_miss~facilityID, complete(imp1)))),
           normclass=c(ICC1(aov(DALY_miss~facilityID, complete(imp2))),
                       ICC1(aov(SEP_miss~facilityID, complete(imp2))),
                       ICC1(aov(csttotaptient_miss~facilityID, complete(imp2)))))

# In including the class variable during estimation for the imputation, we are using a fixed effects approach
# approximate to formulating separate regression models for each facility and imputing within facilities from 
# these models

# check convergence of the imputations - inspecting trace lines
plot(imp2, c("DALY_miss",
             "SEP_miss",
             "csttotaptient_miss"))

# Add another 10 iterations to inspect trace lines again
imp3 <- mice.mids(imp2, maxit=10)
plot(imp3, c("DALY_miss", "SEP_miss", "csttotaptient_miss"))

# Seems okay. Now add another 20 iterations to confirm this - hairy caterpillars!!
imp4 <- mice.mids(imp3, maxit=20, print=FALSE)
plot(imp4, c("DALY_miss", "SEP_miss", "csttotaptient_miss"))

# Plot the densities of observed and imputed data
densityplot(imp2, ~DALY_miss)
densityplot(imp2, ~csttotaptient_miss)
densityplot(imp2, ~DALY_miss | .imp)
# imp2 overestimates magnitudes of costs

# View the imputed dataset
complete(imp2, 1) [1:17, ]

# Do a final imputation of the dataset, using predictive mean matching,
# including all variables as predictors

imp5 <- mice(data1)
densityplot(imp5, ~DALY_miss)
densityplot(imp5, ~csttotaptient_miss)
# comparing the pmm(imp5) against norm(imp2), pmm fits the shape of the data better as samples from the observed values

# Now, compare the ICCs from each of the datasets
data.frame(vars=names(data1[c(7,12,13)]),
           observed=c(ICC1(aov(DALY_miss~facilityID, data1)),
                      ICC1(aov(SEP_miss~facilityID, data1)),
                      ICC1(aov(csttotaptient_miss~facilityID, data1))),
           norm=c(ICC1(aov(DALY_miss~facilityID, complete(imp1))),
                  ICC1(aov(SEP_miss~facilityID, complete(imp1))),
                  ICC1(aov(csttotaptient_miss~facilityID, complete(imp1)))),
           normclass=c(ICC1(aov(DALY_miss~facilityID, complete(imp2))),
                       ICC1(aov(SEP_miss~facilityID, complete(imp2))),
                       ICC1(aov(csttotaptient_miss~facilityID, complete(imp2)))),
           pmm=c(ICC1(aov(DALY_miss~facilityID, complete(imp5))),
                 ICC1(aov(SEP_miss~facilityID, complete(imp5))),
                 ICC1(aov(csttotaptient_miss~facilityID, complete(imp5)))),
           orig=c(ICC1(aov(DALY_miss~as.factor(facilityID),data1)),
                  ICC1(aov(SEP_miss~as.factor(facilityID), data1)),
                  ICC1(aov(csttotaptient_miss~as.factor(facilityID), data1))))


# OTHER FLEXIBLE MULTILEVEL APPROACHES
###############################################################################
# 5 DIFFERENT APPROACHES TO CONSIDER
# 2l.norm: Imputes univariate missing data using a two-level normal model with heterogeneous within group variances
# 2l.pan: Imputes univariate missing data using a two-level normal model with homogeneous within group variances
# 2lonly.mean: Imputes the mean of the class within the class
## 2lonly.norm: Imputes univariate missing data at level 2 using Bayesian linear regression analysis
## 2lonly.pmm: Imputes univariate missing data at level 2 using predictive mean matching

# when specifying the predictor matrix, -2 denotes the class variable, a value of 1 indicates a fixed effect and a
# value of 2 indicates a random effect

data1$facilityID <- as.numeric(data1$facilityID)
data1$study_id <- as.factor(data1$study_id)
data1$DALY_miss <- as.integer(data1$DALY_miss)
data1$class <- as.factor(data1$facilityID)

## CHECK THIS!
#### USING 2L.NORM ############################################################
pred["DALY_miss",] <- c(2, -2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
meth <- c("", "", "", "", "", "", "", "", "", "", "", "2l.norm", "", "", "", "", "")
imp6 <- mice(data1, meth=meth, pred=pred, print=FALSE)

densityplot(imp6, ~csttotaptient_miss)

# the plot is similar to those with pmm and class as a fixed effect, 
# plots the density curves
densityplot(data1, ~csttotaptient_miss) #true data
lines(densityplot(imp5, ~csttotaptient_miss), col = "red", lwd=2) #2l.norm
lines(density(complete(imp6)$csttotaptient_miss), col = "green", lwd=2) #PMM

# consider convergence
plot(imp6)
# run more imputations, convergence not apparent
imp6b <- mice.mids(imp6, maxit=10, print=FALSE)
plot(imp6b)

#### USING 2L.PAN #############################################################
# if group variances are homogenous - use 2l.pan to impute
pred["DALY_miss",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 2, 2, 2, 2, 2)
meth <- c("", "", "", "", "", "", "", "", "", "", "", "2l.pan", "", "", "", "", "")
imp7 <- mice(data1, meth=meth, pred=pred, print=FALSE)


# REVIEW DATASET AND IMPUTE EACH VARIABLE BASED ON ITS CHARACTERISTICS
# DALY_MISS [13] <- 2L.NORM
# SEP_MISS [7] <- 2LOGREG
# CSTTOTAPTIENT_MISS [12] <- 2L.NORM

pred["DALY_miss",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2) #2l.norm
pred["SEP_miss",] <- c(0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2) #2lonly.mean
pred["csttotaptient_miss",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2) #2l.norm

meth <- c("", "", "", "", "", "", "2l.norm", "", "", "", "", "2l.norm", "2l.norm", "", "", "", "")
imp8 <- mice(data1, meth=meth, pred=pred, print=FALSE)

# repeat the imputations using pmm for everything
pmmdata <- data1
pmmdata$class <- as.factor(data1$facilityID)
imp9 <- mice(pmmdata, m=5, print=FALSE)

# plot the densities
densityplot(imp9)
plot(imp9)



