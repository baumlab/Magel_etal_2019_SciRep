
# Effects of bleaching-associated mass coral mortality on reef structural complexity across a gradient of local human disturbance

# Authors: Jennifer M.T. Magel [1], John H.R. Burns [2], Ruth D. Gates [2], Julia K. Baum [1,2]
# Institution: [1] Department of Biology, University of Victoria, Victoria, British Columbia, V8P 5C2, Canada
# [2] Hawai‘i Institute of Marine Biology, Kāne‘ohe, Hawai‘i, 96744, USA
# Corresponding Author: Jennifer Magel, Email: jenn.magel@gmail.com


# Script to fit linear mixed-effects models for coral reef structural complexity metrics and perform model selection 


##############################

# Load necessary packages
library(arm)
library(nlme)
library(MuMIn)
library(knitr)

# Set your working directory
# Make sure that this contains the "KI_complexity_data.csv" data file
setwd("C:/Users/...") # If on a PC
setwd("/Users/...") # If on a Mac

# Load the data
ki_SC <- read.csv("KI_complexity_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Change variables to the proper format
# Change Site to a factor
ki_SC$Site <- as.factor(ki_SC$Site)
# Change Year to a factor
ki_SC$Year <- as.factor(ki_SC$Year)
# Change heat_stress to a factor, order levels
ki_SC$heat_stress <- factor(ki_SC$heat_stress, levels = c('Before', 'After'))
# Change hum_dist to a factor, order levels
ki_SC$hum_dist <- factor(ki_SC$hum_dist, levels = c('Very Low', 'Medium', 'Very High'))

# Create new data frame with standardized predictor variables
ki_SCz <- ki_SC
ki_SCz$dbranching <- rescale(ki_SC$dbranching)
ki_SCz$dplating <- rescale(ki_SC$dplating)
ki_SCz$dmassive <- rescale(ki_SC$dmassive)


##############################

## Calculate mean, standard deviation, and standard error for each structural complexity metric

# Rugosity
rug <- aggregate(rugosity ~ heat_stress, ki_SC, mean)
rug.sd <- aggregate(rugosity ~ heat_stress, ki_SC, function(x) sd(x))
rug.se <- aggregate(rugosity ~ heat_stress, ki_SC, function(x) sd(x)/sqrt(length(x)))

# Terrain ruggedness
vrm <- aggregate(terrain_rug ~ heat_stress, ki_SC, mean)
vrm.sd <- aggregate(terrain_rug ~ heat_stress, ki_SC, function(x) sd(x))
vrm.se <- aggregate(terrain_rug ~ heat_stress, ki_SC, function(x) sd(x)/sqrt(length(x)))

# Curvature
curv <- aggregate(abs_curv ~ heat_stress, ki_SC, mean)
curv.sd <- aggregate(abs_curv ~ heat_stress, ki_SC, function(x) sd(x))
curv.se <- aggregate(abs_curv ~ heat_stress, ki_SC, function(x) sd(x)/sqrt(length(x)))


##############################

## Fit models and perform model selection for each of the three strutural complexity metrics (surface rugosity, terrain ruggedness, and absolute curvature)


############################
##### SURFACE RUGOSITY #####
############################

# Fit the global model using the standardized predictors
# Including interaction between heat stress and human disturbance, and a random factor of Site
modelr <- lme(rugosity ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
              data = ki_SCz, method = "ML")

# Check goodness of fit of global model
r.squaredGLMM(modelr)

# Fit additional global models with a varIdent variance structure to see if this improves model fit
# Variance structure for year
modelrx <- lme(rugosity ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
               weights = varIdent(form = ~1 | Year), data = ki_SCz, method = "ML")
# Variance structure for human disturbance
modelrz <- lme(rugosity ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
               weights = varIdent(form = ~1 | hum_dist), data = ki_SCz, method = "ML")

# Check whether variance structure improves model fit
AICc(modelr, modelrx, modelrz)

# Inclusion of a variance structure for year improves model fit for the rugosity models


# Fit models containing all possible combinations of variables
model.setr <- dredge(modelrx,  extra = "r.squaredGLMM") # Add R squared values to the model selection output

# Create model selection table
model_tabler <- model.sel(model.setr)
names(model_tabler) <- c("(Int)", "Branch", "Mass", "Plat", "Heat", "Dist", "Year*Dist", "R2m", "R2c",
                         "df", "LL", "AICc", "delta", "weight")
kable(model_tabler, digits = 3)

# Create the top model set (all models with delta AICc < 4)
model.setr4 <- get.models(model.setr, subset = delta < 4)

# Perform model averaging of the top model set
avg_modelr <- model.avg(model.setr4)

# View model-averaged parameter estimates for each predictor
summary(avg_modelr)

# Produce 95% confidence intervals for the model-averaged parameter estimates, calculated using the zero method
confint(avg_modelr, full = TRUE)


##############################
##### TERRAIN RUGGEDNESS #####
##############################

# Fit the global model using the standardized predictors
# Including interaction between heat stress and human disturbance, and a random factor of Site
modelv <- lme(terrain_rug ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
              data = ki_SCz, method = "ML")

# Check goodness of fit of global model
r.squaredGLMM(modelv)

# Fit additional global models with a varIdent variance structure to see if this improves model fit
# Variance structure for year
modelvx <- lme(terrain_rug ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
               weights = varIdent(form = ~1 | Year), data = ki_SCz, method = "ML")
# Variance structure for human disturbance
modelvz <- lme(terrain_rug ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
               weights = varIdent(form = ~1 | hum_dist), data = ki_SCz, method = "ML")

# Check whether variance structure improves model fit
AICc(modelv, modelvx, modelvz)

# Inclusion of a variance structure does not improve model fit for the terrain ruggedness models


# Fit models containing all possible combinations of variables
model.setv <- dredge(modelv,  extra = "r.squaredGLMM")

# Create model selection table
model_tablev <- model.sel(model.setv)
names(model_tablev) <- c("(Int)", "Branch", "Mass", "Plat", "Heat", "Dist", "Year*Dist", "R2m", "R2c", 
                         "df", "LL", "AICc", "delta", "weight")
kable(model_tablev, digits = 3)

# Create the top model set (all models with delta AICc < 4)
model.setv4 <- get.models(model.setv, subset = delta < 4)

# Perform model averaging of the top model set
avg_modelv <- model.avg(model.setv4)

# View model-averaged parameter estimates for each predictor
summary(avg_modelv)

# Produce 95% confidence intervals for the model-averaged parameter estimates
confint(avg_modelv, full = TRUE)


##############################
##### ABSOLUTE CURVATURE #####
##############################

# Fit the global model using the standardized predictors
# Including interaction between heat stress and human disturbance, and a random factor of Site
modelc <- lme(abs_curv ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
              data = ki_SCz, method = "ML")

# Check goodness of fit of global model
r.squaredGLMM(modelc)

# Fit additional global models with a varIdent variance structure to see if this improves model fit
# Variance structure for year
modelcx <- lme(abs_curv ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
               weights = varIdent(form = ~1 | Year), data = ki_SCz, method = "ML")
# Variance structure for human disturbance
modelcz <- lme(abs_curv ~ heat_stress * hum_dist + dbranching + dplating + dmassive, random = ~1 | Site, 
               weights = varIdent(form = ~1 | hum_dist), data = ki_SCz, method = "ML")

# Check whether variance structure improves model fit
AICc(modelc, modelcx, modelcz)

# Inclusion of a variance structure for heat stress improves model fit for the curvature models


# Fit models containing all possible combinations of variables
model.setc <- dredge(modelcx,  extra = "r.squaredGLMM")

# Create model selection table
model_tablec <- model.sel(model.setc)
names(model_tablec) <- c("(Int)", "Branch", "Mass", "Plat", "Heat", "Dist", "Year*Dist", "R2m", "R2c", 
                         "df", "LL", "AICc", "delta", "weight")
kable(model_tablec, digits = 3)

# Create the top model set (all models with delta AICc < 4)
model.setc4 <- get.models(model.setc, subset = delta < 4)

# Perform model averaging of the top model set
avg_modelc <- model.avg(model.setc4)

# View model-averaged parameter estimates for each predictor
summary(avg_modelc)

# Produce 95% confidence intervals for the model-averaged parameter estimates
confint(avg_modelc, full = TRUE)
