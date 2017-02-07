###############################################################################

### ABVEGROUND BIOMASS BASED ON CHAVE 2014 ALLOMETRY
### PROPAGATES UNCERTAINTY

###############################################################################

###PACKAGES
library('ggplot2')
library('dplyr')
source("r_code/functions.R")
library('R2jags')
library('reshape2')
source("r_code/generalBAYESIANcode.R")


###############################################################################

### READ IN DATASET FROM CHAVE 2014
chave_data <- read.csv("data/Chave_GCB_Direct_Harvest_Data.csv",header = TRUE)
str(chave_data)

### READ IN GEOLOCATIONS E VARIABLE
chave_geo <- read.csv("outputs/chave_sites_E.csv",header = TRUE)
str(chave_geo)
### rename locality=site
colnames(chave_geo)[1] <- "Site"
head(chave_geo)

### ADD E PARAMETER TO DIRECT HAVEST DATA 
chave_data2 <- merge(chave_data,chave_geo[,c("Site","E")],by="Site",all=FALSE)
head(chave_data2)

### PREP DATA FOR CHAVE MODEL
nobs <- nrow(chave_data2)
E <- as.numeric(chave_data2$E)
WSG <- as.numeric(chave_data2$Wood.specific.gravity)
DBH <- as.numeric(chave_data2$DBH.cm.)
AGB <- as.numeric(chave_data2$Dry.total.AGB.kg.)
HEIGHT <- as.numeric(chave_data2$Total.height.m.)
HEIGHT_measured <- as.numeric(chave_data2$Total.height.m.)


###############################################################################
###############################################################################


# AAPLY CHAVE DIAMETER-HEIGHT ALLOMETRY AND AGB ALLOMETRY TO ESTIMATE (1ST) HEIGHT AND (2ND) AGB FOR SURINAME 
# DATASET


#READ IN SURINAME DATA
suriname_data <- read.csv("outputs/kabo_clean_longNOV22.csv")
str(suriname_data)


#ADD E PARAMETER
#su_E_parameter <- read.csv("outputs/Kabo_chave_E_parameter.csv") 
#str(su_E_parameter)
#colnames(su_E_parameter)[2] <- "PLOT"


#MERGE E DATAFILE WITH DATASET BASED ON PLOT ID
#suriname_data2 <- merge(suriname_data[,c("PLOT","census","BOT_NAME","DIAM","gwd_wd")],su_E_parameter[,c("PLOT","E")],by="PLOT",all = FALSE)

#PREP DATA FOR JAGS
SU_DBH <- as.numeric(suriname_data$DIAM)
SU_WSG <- as.numeric(suriname_data$gwd_wd)
SU_E <- as.numeric(suriname_data$E)
SU_nobs <-nrow(suriname_data) 

cat ("model {
     
     #uniformed priors
     Beta0 ~ dnorm(0,1e-06)
     Beta1 ~ dnorm(0,1e-06)
     Beta2 ~ dnorm(0,1e-06)

     Beta_agb ~ dgamma(0.001,0.001) #constrain to positive; negative biomass values not realistic
     Beta_pow ~ dgamma(0.001,0.001)
     
     #parameter for gamma distribution
     shape ~ dunif(0, 100) #shape parameter for gamma dist
     shape2 ~ dunif(0, 100)
     
     #diameter-height allometry: based on Chave dataset
     
     for (i in 1:nobs) {
     
     log(mu[i]) <- Beta0 - E[i] + Beta1 * log(DBH[i]) - Beta2 * (log(DBH[i]))^2
     HEIGHT[i] ~ dgamma(shape,shape/mu[i]) #MEASURED HEIGHT WITH SAMPLING UNCERTAINTY
     
     mu_agb[i] <- Beta_agb * pow((WSG[i] * (DBH[i]^2) * HEIGHT_measured[i]),Beta_pow)
     AGB[i] ~ dgamma(shape2,shape2/mu_agb[i]) #MEASURED AGB
     
     
     
     }
     
     # predict tree height for suriname data and use in chave eq.4
     
     for (i in 1:SU_nobs) {
     
     log(hei_mu_SU[i]) <- Beta0 - SU_E[i] + Beta1 * log(SU_DBH[i]) - Beta2 * (log(SU_DBH[i]))^2
     hei_SU[i] ~ dgamma(shape,shape/hei_mu_SU[i]) #height estimate with sampling uncertainty
     
     
     #estimate biomass based on Chave allometry eq. 4 
     
     mu_agb_SU[i] <- Beta_agb * pow((SU_WSG[i] * (SU_DBH[i]^2) * hei_SU[i]),Beta_pow)
     AGB_SU[i] ~ dgamma(shape2,shape2/mu_agb_SU[i]) #MEASURED AGB
     
     
     }
     
     
     
     }",fill=TRUE,file="outputs/CHAVE_SU_HEIGHT_AGB.txt")

### MODEL PARAMETERS TO SAVE
params=c("Beta0","Beta1","Beta2","Beta_agb","Beta_pow","shape2","shape","AGB","hei_SU","AGB_SU")

### SEND DATA TO JACS
jags.data_uncertainty<-list(nobs=nobs,DBH=DBH,HEIGHT=HEIGHT,HEIGHT_measured=HEIGHT_measured,AGB=AGB,WSG=WSG,E=E,SU_nobs=SU_nobs,SU_DBH=SU_DBH,SU_WSG=SU_WSG,SU_E=SU_E)



### RUN MODEL AND CALL RESULTS
CHAVE_HEIGHT_AGB_SU<-jags(data=jags.data_uncertainty,  inits=NULL,parameters.to.save=params, "outputs/CHAVE_SU_HEIGHT_AGB.txt",  n.chains=3, n.iter=100000, n.burnin=10000, n.thin=200,DIC=TRUE)

save(CHAVE_HEIGHT_AGB_SU,file="outputs/CHAVE_HEIGHT_AGB_SU.RData")

#
#PLOT 1: Rhat (Model Convergence)

hist(Rhat(CHAVE_HEIGHT_AGB_SU),breaks=1000,main="Model Convergence")
text(1.005,0.8,paste("Max Rhat=",round(max(Rhat(CHAVE_HEIGHT_AGB_SU)),2)))
bads<-length(which(Rhat(CHAVE_HEIGHT_AGB_SU)>1.1))/length(Rhat(CHAVE_HEIGHT_AGB_SU))
text(1.005,0.7,paste("% Rhat >1.1 =",round(bads,5)))

### COMPARE PARAMETERS FROM MODEL OUTPUT TO CHAVE ALLOMETRY ESTIMATES
Q(CHAVE_HEIGHT_AGB_SU,"Beta2")

#check chain convergence
plot(as.mcmc(CHAVE_HEIGHT_AGB_SU))
#extract out parameter results from model
#coda summary format
outJG<-as.mcmc(CHAVE_HEIGHT_AGB_SU)
summary(outJG)


###############################################################################
###############################################################################