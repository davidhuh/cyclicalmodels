################################################################################
## Modeling Cyclical Patterns in Daily College Drinking Data with Many Zeroes ##
##                                                                            ##
##  Authors: Huh, D., Kaysen, D., & Atkins, D.C.                              ##
## Software: R 3.1.3                                                          ##
##  Updated: 09 April 2015
################################################################################

######## Install glmmADMB (if necessary) ########

## Installation instructions at http://glmmadmb.r-forge.r-project.org/

## Install Method #1
#install.packages("glmmADMB", repos="http://r-forge.r-project.org",
#                 type="source")

## Method #2 (if #1 doesn't work)
#install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                 getOption("repos")),
#                 type="source")

## Load glmmADMB package
library(glmmADMB)


########## Load Data ##########################

dash.df <- read.csv("dash.csv")


########## Prepare Data for Analysis ############

dash.df <- within(dash.df, {
  ## ID variable must be a factor for glmmADMB
  id <- factor(id)
  
  ## Divide outcome into binary and truncated count components
  y.bin <- ifelse(drinks > 0, 1, 0)
  y.cnt <- ifelse(drinks > 0, drinks, NA)

  ## Create cyclical variables;
  cycPha <- sin(2*pi*(dayofwk/7))
  cycAmp <- cos(2*pi*(dayofwk/7))
  
  ## Create saturated set of dummy variables
  Tue <- as.numeric(dayofwk==1)
  Wed <- as.numeric(dayofwk==2)
  Thu <- as.numeric(dayofwk==3)
  Fri <- as.numeric(dayofwk==4)
  Sat <- as.numeric(dayofwk==5)
  Sun <- as.numeric(dayofwk==6)
})

## Create separate datasets for the logistic and truncated count regressions
dash.bin.df <- dash.df[!is.na(dash.df$y.bin), ]
dash.cnt.df <- dash.df[!is.na(dash.df$y.cnt), ]


########## Multilevel Hurdle Negative Binomial Model -- No Time Predictors ########

## Logistic submodel w/ Saturated set of dummy variables ###
##  -- Random intercepts
hunb.cyc.bin.0 <- glmmadmb(y.bin ~ Tue + Wed + Thu + Fri + Sat + Sun + (1|id),
                           data=dash.bin.df,
                           family="binomial",
                           corStruct="diag",
                           verbose = TRUE)

### Truncated Negative Binomial (NB2) submodel w/ Cyclical variables ###
##  -- Random intercepts and slopes
hunb.cyc.cnt.0 <- glmmadmb(y.cnt ~ cycAmp + cycPha + (cycAmp + cycPha|id),
                           data=dash.cnt.df,
                           family="truncnbinom",
                           corStruct="diag",
                           verbose=TRUE)

summary(hunb.cyc.bin.0); summary(hunb.cyc.cnt.0)


########## Multilevel Hurdle Negative Binomial Model -- with Time Predictor ######

## Logistic submodel w/ Saturated dummy variables ###
hunb.cyc.bin.1 <- glmmadmb(y.bin ~ dmqsoc*(Tue + Wed + Thu + Fri + Sat + Sun) + (1|id),
                           data=dash.bin.df,
                           family="binomial",
                           corStruct="diag",
                           verbose = TRUE)

## Truncated Negative Binomial (NB2) submodel w/ Cyclical variables ###
hunb.cyc.cnt.1 <- glmmadmb(y.cnt ~ dmqsoc*(cycAmp+cycPha) + (cycAmp+cycPha|id),
                           data=dash.cnt.df,
                           family="truncnbinom",
                           corStruct="diag",
                           verbose=TRUE)

summary(hunb.cyc.bin.1); summary(hunb.cyc.cnt.1)
