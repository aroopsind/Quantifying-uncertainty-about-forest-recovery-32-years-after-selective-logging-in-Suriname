######################################################################

#GROWTH,RECRUITMENT AND MORTALITY ANALYSIS FOR ACS
#INFORMED PRIORS

######################################################################
###load libraries
library('ggplot2')
library('dplyr')
source("r_code/functions.R")
library('R2jags')
library('reshape2')
source("r_code/generalBAYESIANcode.R")



###import clean kabo_data; long format
stem_data <- read.csv("outputs/kabo_clean_longNOV22.csv")

##read in census dates
dates <- read.csv("outputs/census_dates_ALL.csv")
##drop unneeded census
dates2 <- dates %>%
      select(c(PLOT,census83_00,census00_12))


###drop 1978,1980, censuses
stem_data2 <- stem_data[which(stem_data$census %in% c(1983,2000,2012)),]

###add census number (ordinal numbering)
stem_data2$num.census <- census.num(stem_data2$census)

###for loop to identify when a new tree enters the >15cm diameter class at each census
stem_data2$RECRUIT <- NA

for (i in 1:nrow(stem_data2)) {
      
      stem_data2$RECRUIT[i] = if (stem_data2$num.census[i] == 1) {
            ("A")
      }
      else {
            previous.censuses = c (1:3)[which(c(1:3) < stem_data2$num.census[i])]
            checkers=unique(stem_data2[stem_data2$num.census %in% previous.censuses,]$rec.no)
            recruit = stem_data2$rec.no[i] %in% checkers
            if (recruit == TRUE) {("A")} else { ("I")}
            
            
            
      }
}

#write.csv(stem_data2,"outputs/stem_data_clean_recruits.csv",row.names = FALSE)
#stem_data2=read.csv("outputs/stem_data_clean_recruits.csv",header = TRUE)


#check code worked
recruits=stem_data2[which(stem_data2$RECRUIT =="I"),]; head(recruits);unique(recruits$RECRUIT)

recruits %>%
      group_by(PLOT,census) %>%
      summarise(freq=n())

#number of trees that survived from previous census (1983,2000) compared to 2012
rec_id_2012<- unique(stem_data2[which(stem_data2$census == 2012),]$rec.no)
rec_other <- unique(stem_data2[which(!stem_data2$census == 2012),]$rec.no)
length(rec_id_2012 %in% rec_other)


###ACS by individual stems based on measured diameter @15 CM
#height allometry
stem_data2$height_base <- exp(0.893 - stem_data2$E + 0.760 * (log(15)) - 0.0340 *(log(15))^2)


#biomass allometry/converted to carbon and Mg
stem_data2$ACS_baseline <- ((0.0673 * (stem_data2$gwd_wd * (15^2)*stem_data2$height_base)^0.976)*0.47)/1000
#write.csv(stem_data2,"outputs/stem_data_acs_table.csv",row.names=FALSE)

#################################################################################
#PROCESSES UNDERLYING NET ACS CHANGE:
#CLARK 2001: APPROACH 2
#STAND INCREMENT = (SUM OF ACS @T2 - SUM OF ACS @T1) +
#                  (SUM OF BIOMASS OF TREE THAT DIED IN THE INTERVAL) -
#                 [BIOMASS OF A MINIMUM SIZE TREE X NUMBER OF NEW TREES]
#################################################################################

#1. ACS STOCKS (1983,2000,2012)
#2. ACS MORTALITY (2000,2012), 
#3. ACS increment of recruits = ACS of same tree @ 15 cm DBH * NUMBER OF RECRUITS

#################################################################################

#1. Initial ACS (t0,t1, & t2)
###stripping away all the trees alive at census t0=1983;t1=2000; and t2=2012; 
stems_t0 <- stem_data2[which(stem_data2$census==1983),] ; stems_t0_id <- unique(stems_t0$rec.no)
stems_t1 <- stem_data2[which(stem_data2$census==2000),] ; stems_t1_id <- unique(stems_t1$rec.no)
stems_t2 <- stem_data2[which(stem_data2$census==2012),] ; stems_t2_id <- unique(stems_t2$rec.no) #use v2, that includes missed trees


###trees recorded in 1983; missed in 2000 and reappeared in 2012; use recruit status =="A" 
missed_t1 <- stems_t2[which(!stems_t2$rec.no %in% stems_t1_id & stems_t2$RECRUIT=="A"),]
stems_t1_v2=rbind(stems_t1,missed_t1)

#################################################################################
#################################################################################
#2. ACS mortality, trees that die between to-t1 & t1-t2 = ACS of dead tree - ACS  of same tree @ 10cm DBH

#residual trees between 1983 to 2000
census_2000_residual <- stems_t1[which(stems_t1$RECRUIT == "A"),]
#trees recorded in 1983 and died by 2000
census_2000_mort <- stems_t0[which(!stems_t0$rec.no %in% census_2000_residual$rec.no),]


#summarise by plot
census_2000_mort %>%
      group_by(PLOT) %>%
      summarise(ACS_mort_2000=sum(ACS))



#residual trees between 2000 to 2012
census_2012_residual <- stems_t2[which(stems_t2$RECRUIT == "A"),]

#trees recorded in 2000 and died by 2012
census_2012_mort <- stems_t1_v2[which(!stems_t1_v2$rec.no %in% census_2012_residual$rec.no),]


#summarise by plot
census_2012_mort %>%
      group_by(PLOT) %>%
      summarise(ACS_mort_12=sum(ACS))

#################################################################################

#3. ACS increment of recruits = ACS of same tree @ 15 cm DBH * NUMBER OF RECRUITS


recruits=stem_data2[which(stem_data2$RECRUIT =="I"),]; head(recruits);unique(recruits$RECRUIT)


ACS_recruit <-  recruits %>%
      group_by(PLOT,census) %>%
      summarise(ACS_recruit=sum(ACS_baseline))




#################################################################################
#1. residual standing carbon stocks


#summarise ACS by plot for 1983
stem_1983_plot <- as.data.frame(stems_t0 %>%
                                      group_by(PLOT) %>%
                                      summarise(ACS1=sum(ACS)))

#summarise ACS by plot for 2000
stem_2000_plot <- as.data.frame(stems_t1_v2 %>%
                                      group_by(PLOT) %>%
                                      summarise(ACS2=sum(ACS)))

#summarise ACS by plot for 2012
stem_2012_plot <- as.data.frame (stems_t2 %>%
                                       group_by(PLOT) %>%
                                       summarise(ACS3=sum(ACS)))

### combine residual carbon stocks
residual_ACS <- Reduce(function(...) merge(..., all=TRUE), list(stem_1983_plot, stem_2000_plot,stem_2012_plot))

#################################################################################

#2. mortality between censuses

# dead trees between 1983 - 2000
mort1 <- as.data.frame(census_2000_mort %>%
                             group_by(PLOT) %>%
                             summarise(mort_2000=sum(ACS)))

# dead trees between 2000 - 2013
mort2 <- as.data.frame(census_2012_mort %>%
                             group_by(PLOT) %>%
                             summarise(mort_2012=sum(ACS)))

### combine mortality
mortality_ACS <- Reduce(function(...) merge(..., all=TRUE), list(mort1,mort2))

#################################################################################

#3. recruit between censuses

# recruited ACS 1983 - 2000
ACS_recruit1 <- as.data.frame(ACS_recruit[which(ACS_recruit$census==2000),])
names(ACS_recruit1)[2:3] <- c("census2000","recruit_2000")

# recruited ACS 2000- 2012
ACS_recruit2 <- as.data.frame(ACS_recruit[which(ACS_recruit$census==2012),])
names(ACS_recruit2)[2:3] <- c("census2012","recruit_2012")


### combine mortality
recruited_ACS <- Reduce(function(...) merge(..., all=TRUE), list(ACS_recruit1,ACS_recruit2))


#################################################################################
#STAND INCREMENT = (SUM OF ACS @T2 - SUM OF ACS @T1) +
#                  (SUM OF BIOMASS OF TREE THAT DIED IN THE INTERVAL) -
#                 [BIOMASS OF A MINIMUM SIZE TREE X NUMBER OF NEW TREES]


npp_increment <- Reduce(function(...) merge(..., all=TRUE), list(residual_ACS,mortality_ACS,recruited_ACS,dates2))

npp_increment$npp_1 <- ((npp_increment$ACS2 - npp_increment$ACS1) + (npp_increment$mort_2000) - (npp_increment$recruit_2000))/npp_increment$census83_00

npp_increment$npp_2 <- ((npp_increment$ACS3 - npp_increment$ACS2) + (npp_increment$mort_2012) - (npp_increment$recruit_2012))/npp_increment$census00_12

#################################################################################

### modify for JAGS
npp_increment2 <- npp_increment %>%
      select(c(PLOT, npp_1, npp_2))

#convert to long format
npp_long <- melt(npp_increment2,
                 id.vars=c("PLOT"),
                 measure.vars=c("npp_1","npp_2"),
                 variable.name="npp",
                 value.name="carbon_increment"
)

#add logging treatment based on plot
npp_long$m3 <- log_intensity_m3(plot_id = npp_long$PLOT)


############################################################################


#NPP CORRECTION USING INFORMED PRIORS MODEL

#data for jags model
nobs=as.numeric(nrow(npp_long))
npp_measured=as.numeric((npp_long$carbon_increment))
plot_id <- as.factor(npp_long$PLOT)
log_m3 <- as.numeric(npp_long$m3)

##########################################################################################

#priors from Johnson supplemental 2016
johnson_wp <- data.frame(site_id=c("ELD-01","ELD-02","ELD-03","ELD-04","FMH-01","IWO-21","IWO-22","NOU-01","NOU-02","NOU-03","NOU-04","NOU-05","NOU-06","NOU-07","NOU-08","NOU-09","NOU-10","NOU-11","NOU-12","NOU-13","NOU-14","NOU-15","NOU-16","NOU-17","NOU-18","NOU-19","NOU-20","NOU-21","NOU-22","PAR-20","PAR-21","PAR-22","PAR-23","PAR-24","PAR-26","PAR-27","PAR-28","PAR-29","RIO-01","RIO-02","SCR-05"), C_ha_yr=c(5.25,3.46,4.30,3.79,3.91,2.51,2.35,5.30,4.35,4.78,4.84,3.70,4.13,3.47,3.85,3.79,3.39,4.17,3.60,3.18,2.79,2.57,3.04,4.06,3.29,3.24,4.23,3.13,3.08,2.22,2.16,3.32,3.11,4.20,2.48,2.29,2.50,2.95,3.92,3.99,3.16),
                         W_l=c(3.45,0.22,8.48,2.58,0.42,1.20,0.77,2.38,1.95,3.94,2.53,4.77,2.56,2.35,2.57,3.70,3.36,1.63,2.26,1.80,4.44,2.35,3.86,1.03,0.81,3.05,1.16,0.78,3.99,16.69,3.45,3.90,0.99,7.43,4.25,3.42,1.13,0.72,4.55,2.87,1.87))
                             
wp_mu <- mean(johnson_wp$C_ha_yr)
wp_sd <- sd(johnson_wp$C_ha_yr)  
sigShape = gammaShRaFromMeanSD( mean=sd(johnson_wp$C_ha_yr) , sd=sd(johnson_wp$C_ha_yr))$shape 
sigRate = gammaShRaFromMeanSD( mean=sd(johnson_wp$C_ha_yr) , sd=sd(johnson_wp$C_ha_yr))$rate 
sd(johnson_wp$C_ha_yr)/sqrt(length(johnson_wp$C_ha_yr))
##########################################################################################

#JAGS MODEL WITH INFORMED PRIORS
cat("model { 
    
    #uninformed priors for npp
    npp_mu ~ dnorm(wp_mu, 1/wp_sd^2) #informed priors
    tau <- pow(sigma, -2) #precision of residuals
    sigma ~ dgamma(sigShape, sigRate) #informed priors
    
    sd_plot ~ dunif(0,10) #plot level variance
    tau_plot <- 1/(sd_plot*sd_plot)
    
    #priors for treatment
    logtrt~dnorm(0, .0001)
    
    #plot level random effects
    for (i in 1:12) # for each plot
    {
    plot_Eff[i] ~ dnorm (0,tau_plot) # random plot effects
    }
    
    
    for (i in 1:nobs) # 
    {
    
    mu[i] <- npp_mu + plot_Eff[plot_id[i]] + logtrt * log_m3[i]
    npp_measured[i] ~ dnorm (mu[i],tau) # random plot effects
    
    #posterior predictive distribution
    pred_npp[i]~dnorm(mu[i],tau)   
    }
    
    
    
    }", fill=TRUE, file="outputs/npp_model_InformedPriors.txt")

#params to save
params=c("npp_mu","tau","plot_Eff","logtrt","pred_npp","sigma")

#jags data
jags.data<-list(nobs=nobs,npp_measured=npp_measured,plot_id=plot_id,log_m3=log_m3,wp_mu=wp_mu,wp_sd=wp_sd,sigRate=sigRate,sigShape=sigShape)

draws(n.iter =100000 ,n.burnin =20000,n.thin = 200 ,n.chains = 3)

#Run the model and call the results
npp_model<-jags(data=jags.data,  inits=NULL,parameters.to.save=params, "outputs/npp_model_InformedPriors.txt",  n.chains=3, n.iter=300000, n.burnin=20000, n.thin=2000,DIC=TRUE)

#save(npp_model,file="outputs/npp_model_informedPriors.RData")
#load("outputs/npp_model_informedPriors8July.RData")

#check chain convergence
plot(as.mcmc(npp_model))
#extract out parameter results from model
par(mfrow=c(1,1))
pars<-gipar(npp_model)

# to get DIC or specify DIC=TRUE in jags() or do the following
dic.samples(npp_model$model, n.iter=1000, type="pD")

#coda summary format
outJG<-as.mcmc.list(npp_model$BUGSoutput)
summary(outJG)


#######################################################################
#######################################################################

#PLOT 1: Rhat (Model Convergence)

hist(Rhat(npp_model),breaks=1000,main="Model Convergence")
text(1.008,0.8,paste("Max Rhat=",round(max(Rhat(npp_model)),2)))
bads<-length(which(Rhat(npp_model)>1.1))/length(Rhat(npp_model))
text(1.013,0.7,paste("% Rhat >1.1 =",round(bads,5)))

#######################################################################
#######################################################################

### CARBON INCREMENT PREDICTIONS FOR TABLE 2

npp_pred <- gipar(npp_model)$pred_npp

### add column names with plot id
colnames(npp_pred) <- plot_id

### convert to long format
npp_pred_long <- melt(npp_pred)
colnames(npp_pred_long) <- c("draws","plot","npp_acs")
### add logging category
npp_pred_long$log.intensity <- log_intensity_cat(npp_pred_long$plot)

### summarise by logging treatments
table2_npp <- npp_pred_long %>%
      group_by(log.intensity) %>%
      summarise(`0.025%`=quantile(npp_acs, probs=0.025),
                 `50%`=quantile(npp_acs, probs=0.5),
                 `0.975%`=quantile(npp_acs, probs=0.975),
                 avg=mean(npp_acs))

write.csv(table2_npp,"outputs/table2_npp.csv",row.names = FALSE)

npp_pred_logg <- npp_pred_long[which(!npp_pred_long$log.intensity == "control"),]
mean(npp_pred_logg$npp_acs)

#######################################################################
#######################################################################

### RESIDUAL GROWTH; RECRUITMENT; MORTALITY FOR TABLE 2
head(stem_data2)

### RESIDUAL GROWTH

##convert to wide format; residual growth between censuses
stem_data_wide <- dcast(stem_data2,rec.no + PLOT ~ census, value.var="ACS")

##add census dates
stem_data_wide2 <- merge(stem_data_wide,dates2,by="PLOT",all=TRUE)

##increments between censuses
stem_data_wide2$increment1 <- (stem_data_wide2$`2000` - stem_data_wide2$`1983`)/stem_data_wide2$census83_00
stem_data_wide2$increment2 <- (stem_data_wide2$`2012` - stem_data_wide2$`2000`)/stem_data_wide2$census00_12

##convert to long format and summarise by plot
stem_increment <- melt(stem_data_wide2,id.vars="PLOT",measure.vars=c("increment1","increment2"),variable.name="census_interval",value.name="ACS")
head(stem_increment)

###summarise by plot
stem_increment_plot <- stem_increment %>%
      group_by(PLOT,census_interval) %>%
      summarise(ACS_increment=sum(ACS,na.rm=TRUE))

###add logging category
stem_increment_plot$log_cat <- log_intensity_cat(stem_increment_plot$PLOT)

stem_increment_plot %>%
      group_by(log_cat) %>%
      summarise(mean_increment=mean(ACS_increment),se_increment=sd(ACS_increment)/sqrt(length(ACS_increment)))

#######################################################################

### RECRUIMENT

head(stem_data2)

stem_recruit <- stem_data2[which(stem_data2$RECRUIT=="I"),]
stem_recruit$ingrowth <- stem_recruit$ACS - stem_recruit$ACS_baseline


###summarise by plot
stem_ingrowth_plot <- stem_recruit%>%
      group_by(PLOT,census) %>%
      summarise(ACS_ingrowth=sum(ingrowth,na.rm=TRUE))



### convert to wide format
stem_recruit_wide <- dcast(stem_ingrowth_plot ,PLOT ~ census, value.var="ACS_ingrowth")

## add census interval
stem_recruit_wide2 <- merge(stem_recruit_wide, dates2,by="PLOT",all=TRUE)

stem_recruit_wide2$ingrowth1 <- stem_recruit_wide2$`2000`/stem_recruit_wide2$census83_00
stem_recruit_wide2$ingrowth2 <- stem_recruit_wide2$`2012`/stem_recruit_wide2$census00_12

##convert to long format and summarise by plot
stem_recruit_long <- melt(stem_recruit_wide2,id.vars="PLOT",measure.vars=c("ingrowth1","ingrowth2"),variable.name="census_interval",value.name="ACS")
head(stem_recruit_long)

###add logging category
stem_recruit_long$log_cat <- log_intensity_cat(stem_recruit_long$PLOT)

stem_recruit_long %>%
      group_by(log_cat) %>%
      summarise(mean_increment=mean(ACS),se_increment=sd(ACS)/sqrt(length(ACS)))


#######################################################################

### MORTALITY

head(mortality_ACS)

## add census interval
mortality_ACS2 <- merge(mortality_ACS,dates2,by="PLOT",all=TRUE)
mortality_ACS2$mort1 <- mortality_ACS2$mort_2000/mortality_ACS2$census83_00
mortality_ACS2$mort2 <- mortality_ACS2$mort_2012/mortality_ACS2$census00_12

### convert to long format
mortality_long <- melt(mortality_ACS2, id.vars="PLOT",measure.vars=c("mort1","mort2"),variable.name="mort_interval",value.name="mort_ACS")
###add logging category
mortality_long$log_cat <- log_intensity_cat(mortality_long$PLOT)

mortality_long %>%
      group_by(log_cat) %>%
      summarise(mean_increment=mean(mort_ACS),se_increment=sd(mort_ACS)/sqrt(length(mort_ACS)))

mortality_logg <- mortality_long[which(mortality_long$log_cat !="control"),]
mortality_control <- mortality_long[which(mortality_long$log_cat =="control"),]
mean(mortality_control$mort_ACS)
##################################################################

#### NPP FROM CLARK 2 FOR CONTROL PLOTS
npp_long$log_cat <- log_intensity_cat(npp_long$PLOT)
npp_long %>%
      group_by(log_cat) %>%
      summarise(avg=mean(carbon_increment))

npp_control <- npp_long[which(npp_long$log_cat=="control"),]
mean(npp_control$carbon_increment)
