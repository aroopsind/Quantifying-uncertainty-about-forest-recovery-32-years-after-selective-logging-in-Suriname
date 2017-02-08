##############################################################################################################
#COMMERCIAL TIMBER STOCK RECOVERY MODELS

###load packages and functions
library('ggplot2')
library('dplyr')
source("r_code/functions.R")
library('R2jags')
library('reshape2')
source("r_code/generalBAYESIANcode.R")


##############################################################################################################
##############################################################################################################

### Read in data and process ACS as weighted to pre-census values

comm_volume <- read.csv("outputs/commercial_HARV_VOL23NOV_35cm.csv")

#empty dataframe to store weighted values
comm_weight <- matrix(NA,12,6)

#weight recovery by pre-logging values; 
for (i in 1:12){
      for (j in 1:6) {
            comm_weight[i,j] <- comm_volume[i,j + 1]/comm_volume[i, 2]
      }
}


comm_weight2 <- as.data.frame(cbind(comm_weight,comm_volume$PLOT))#add plot id
colnames(comm_weight2)[1:7] <- c("1978","1980","1981","1983","2000","2012","PLOT") # add colnames

#recast in long format
comm_data_long <- melt(comm_weight2,
                       id.vars="PLOT",
                       measure.vars = c("1978","1980","1981","1983","2000","2012"),
                       variable.name = "census",
                       value.name = "comm_recovery")

### KEEP ONLY LOGGED PLOTS/DROP PRE-LOGGED CENSUS
comm_data_logPLOTS <- comm_data_long[which(!comm_data_long$PLOT %in% c(41,42,43)),]
comm_data_logPLOTS2 <- comm_data_logPLOTS[which(!comm_data_logPLOTS$census %in% c(1978)),]


##############################################################################################################
##############################################################################################################

### CONTROL PLOTS
comm_vol_control <- comm_volume [which(comm_volume$PLOT %in% c(41,42,43)),][c(1:3),c(1,5,6,7)]

#empty matrix to store values
control_weighted <- matrix(NA,3,3)

for (i in 1:3) {
      for (j in 1:3) {
            control_weighted[i,j] <- comm_vol_control[i,j+1]/comm_vol_control[i, 2] #weighted by 1983 plot levels
      }
}


comm_control2 <- as.data.frame(cbind(control_weighted,comm_vol_control$PLOT)) #add plot

colnames(comm_control2 )[1:4] <- c("1983","2000","2012","PLOT")  # add colnames


#recast in long format
comm_control_long <- melt(comm_control2,
                          id.vars="PLOT",
                          measure.vars = c("1983","2000","2012"),
                          variable.name = "census",
                          value.name = "comm_recovery")

##############################################################################################################
##############################################################################################################
### PROCESS DATA FOR JAGS

#combine control and logged plots
comm_data_recovery <- rbind(comm_data_logPLOTS2 ,comm_control_long)

### add time since logged
comm_data_recovery$tsl <- time.logged(census.yr = comm_data_recovery$census, plot.id = comm_data_recovery$PLOT)


#add logging intensity
comm_data_recovery$log_intensity <- log_intensity_m3(comm_data_recovery$PLOT)


plot_id <- as.factor(comm_data_recovery$PLOT) #unique plot id
time.since.logged <- as.numeric(comm_data_recovery$tsl) #time since logged (years)
timb_recovery <- as.numeric(comm_data_recovery$comm_recovery) #comm. volume (m3 ha)
logged_m3 <- as.numeric(comm_data_recovery$log_intensity) #basal area logged (m2 ha)
Nobs <- nrow(comm_data_recovery) #numb. observations

##############################################################################################################
##############################################################################################################

### Model timber recovery: response weighted recovery by pre-logging census data
#covariates: (1) logging intensity - fixed effect (continuous variable; volume (m3 ha))
#            (2) time since logged - fixed effect (continuous variable - year)
#            (3) plot level random effects
#            (4) interaction tsl * logging intensity
#response variable weighted by commercial volume in 1978

cat("model { #commercial timber stock recovery
    
    #Uninformative priors
    sd_plot ~ dunif(0,10) #variance of plot to plot variation
    intercept ~ dnorm(0,1e-06) #intercept term (baseline;unlogged plots)
    
    tau ~ dgamma (0.001,0.001) #precision of residuals
    tau_plot <- 1/(sd_plot*sd_plot)
    
    logtrt ~ dnorm(0,1e-06) #uninformative prior for logging effect
    tsl ~ dnorm(0,1e-06) #uninformative prior for time since logged
    beta3 ~ dnorm(0,1e-06) #uninformative prior for interaction term

    ### Mariginal and conditional R2
    Vfixed <- (sd(pred_timb_recovery))^2 
    Vresidual<-1/tau # get the variance of residuals
    Vrandom<-1/tau_plot  # get the variance of random plot effect
    marginalR2<-Vfixed/(Vfixed + Vrandom + Vresidual) # calculate marginalR2
    conditionalR2<-(Vrandom + Vfixed)/(Vfixed + Vrandom + Vresidual) # calculate conditional R2
    
    for (i in 1:12) # for each plot
    {
    plot_Eff[i] ~ dnorm (0,tau_plot) # random plot effects
    }
    
    
    for (i in 1:Nobs) 
    { #timb. stocks recovery for each plot at each census interval post logging weighted by pre-logging volumes
    
    #likelihood
    mu[i] <- intercept + logtrt * logged_m3[i] + tsl * time.since.logged[i] + plot_Eff[plot_id[i]] + beta3 * logged_m3[i]*time.since.logged[i]
    
    #observed volume recovery
    timb_recovery[i] ~ dnorm(mu[i],tau)
    
    #posterior predictive distribution
    pred_timb_recovery[i]~dnorm(mu[i],tau)      
    residuals[i] <- (timb_recovery[i] - pred_timb_recovery[i])
    
    }
    }", fill=TRUE, file="outputs/kabo_model_timber.txt")



#params to save
params=c("intercept","logtrt","tsl","plot_Eff","tau","tau_plot","pred_timb_recovery","beta3","marginalR2","conditionalR2","residuals")

#jags data
jags.data<-list(logged_m3=logged_m3,time.since.logged=time.since.logged,plot_id=plot_id,timb_recovery=timb_recovery,Nobs=Nobs)


#Run the model and call the results
kabo_model_timber<-jags(data=jags.data,  inits=NULL,parameters.to.save=params, "outputs/kabo_model_timber.txt",  n.chains=3, n.iter=500000, n.burnin=25000, n.thin=2500,DIC=TRUE)
save(kabo_model_timber,file="outputs/kabo_model_timber.RData")

#check chain convergence
plot(as.mcmc(kabo_model_timber))
#extract out parameter results from model
par(mfrow=c(1,1))
pars<-gipar(kabo_model_timber)

# to get DIC or specify DIC=TRUE in jags() or do the following
dic.samples(kabo_model_timber$model, n.iter=1000, type="pD")

#coda summary format
outJG<-as.mcmc.list(kabo_model_timber$BUGSoutput)
summary(outJG)


#######################################################################
#######################################################################

#PLOT 1: Rhat (Model Convergence)

hist(Rhat(kabo_model_timber),breaks=1000,main="",xlab="",ylab="")
text(1.015,1.5,paste("Max Rhat=",round(max(Rhat(kabo_model_timber)),2)))
bads<-length(which(Rhat(kabo_model_timber)>1.1))/length(Rhat(kabo_model_timber))
text(1.015,1.3,paste("% Rhat >1.1 =",round(bads,5)))

#######################################################################


### MODEL PREDICTIONS

#time since logged
tsl.seq=seq(from=0.1,to=32,length.out = 100)
nsim=nrow(pars$intercept)


### CONTROL PREDICTIONS
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))

for (i in 1:100) { 
      temp <- rnorm(nsim,pars$intercept + pars$logtrt * 0 + (pars$beta3 * 0 * tsl.seq[i]) + pars$tsl * tsl.seq[i],sd=1/sqrt(pars$tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      
}


predict_control <- yquantile
predict_control_mean <- yhat
predict_contr <- as.data.frame(cbind(predict_control ,predict_control_mean,tsl.seq))
colnames(predict_contr)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")



### LOW INTENSITY
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco_15 <- (matrix(NA,100,1))

for (i in 1:100) { 
      temp <- rnorm(nsim,pars$intercept + pars$logtrt * 15.0 + (pars$beta3 * 15.0 * tsl.seq[i]) +pars$tsl * tsl.seq[i],sd=1/sqrt(pars$tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco_15[i,] <- quantile(temp,c(0.18))
}


predict_15m3 <- yquantile
predict_15m3_mean <- yhat
predict_15m3 <- as.data.frame(cbind(predict_15m3 ,predict_15m3_mean,tsl.seq))
colnames(predict_15m3)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")


### MEDIUM INTENSITY
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco_23 <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 15,23, & 46 volume (m3 ha)
      temp <- rnorm(nsim,pars$intercept + pars$logtrt * 23.0 + (pars$beta3 * 23.0 * tsl.seq[i]) +pars$tsl * tsl.seq[i],sd=1/sqrt(pars$tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco_23[i,] <- quantile(temp,c(0.19))
}

predict_23m3 <- yquantile
predict_23m3_mean <- yhat
predict_23m3 <- as.data.frame(cbind(predict_23m3 ,predict_23m3_mean,tsl.seq))
colnames(predict_23m3)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")

### HIGH INTENSITY
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco_46 <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 15,23, & 46 volume (m3 ha)
      temp <- rnorm(nsim,pars$intercept + pars$logtrt * 46.0 + (pars$beta3 * 46.0 * tsl.seq[i]) +pars$tsl * tsl.seq[i],sd=1/sqrt(pars$tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco_46[i,] <- quantile(temp,c(0.30))
}

predict_46m3 <- yquantile
predict_46m3_mean <- yhat
predict_46m3 <- as.data.frame(cbind(predict_46m3 ,predict_46m3_mean,tsl.seq))
colnames(predict_46m3)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")



### PLOT OUT PREDICTIONS

png(file = "outputs/REVISIONS2/timber_high.png", bg = "transparent", type = c("cairo"), width=2000, height=2000, res=300)

###Plot out predictions
ggplot() + 
      #geom_line(data = predict_contr, aes(x=tsl, y=mean),stat="identity",colour="darkgreen",size=1) +
      #geom_ribbon(data=predict_contr,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="darkgreen") +
      #geom_line(data = predict_15m3, aes(x=tsl, y=mean),stat="identity",colour="royalblue2",size=1) +
      #geom_ribbon(data=predict_15m3,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="royalblue2") +
      #geom_line(data = predict_23m3, aes(x=tsl, y=mean),stat="identity",colour="gold3",size=1) +
      #geom_ribbon(data=predict_23m3,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="gold3") +
      geom_line(data = predict_46m3, aes(x=tsl, y=mean),stat="identity",colour="#e66101",size=1) +
      geom_ribbon(data=predict_46m3,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="#e66101") +
      geom_hline(aes(yintercept=1),size=2,linetype=2,colour="black") +
      scale_y_continuous(limits=c(0.1,1.8),breaks=c(0.2,0.6,1.0,1.4,1.8)) +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      theme(axis.text.x=element_text(size=14)) +
      theme(axis.text.y=element_text(size=14)) +
      theme(legend.text = element_text(size = 14)) +
      ylab(expression("")) +
      xlab(expression("")) +
      #theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_blank())
#theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
#theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
#xlab (expression(bold("Time since logged (years)"))) 


dev.off()

##############################################################################
##############################################################################
### COEFFICIENT PLOT

coef_points_model <- list(LOG_INTENSITY = pars$logtrt, TIME_SINCE_LOGGED =pars$tsl,INTERACTION=pars$beta3)
pplotCOEF(coef_points_model)

##############################################################################
##############################################################################
### SURINAME POLICY IMPLEMENTATION
tsl.seq=seq(from=0.1,to=70,length.out = 100)
nsim=nrow(pars$intercept)


yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))


for (i in 1:100) { 
      temp <- rnorm(nsim,pars$intercept + pars$logtrt * 25.0 + (pars$beta3 * 25.0 * tsl.seq[i]) + pars$tsl * tsl.seq[i],sd=1/sqrt(pars$tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco[i,] <- quantile(temp,c(0.33))
}

policy_predictions <- yquantile
policy_pred_mean <- yhat
predict_policy <- as.data.frame(cbind(policy_predictions ,policy_pred_mean,tsl.seq))
colnames(predict_policy)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")


###Plot out policy predictions
ggplot() + 
      geom_line(data = predict_policy, aes(x=tsl, y=mean),stat="identity",colour="darkgreen",size=1) +
      geom_ribbon(data=predict_policy,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="darkgreen") +
      geom_hline(aes(yintercept=1),size=2,linetype=2,colour="black") +
      #scale_y_continuous(limits=c(0.3,1.3),breaks=c(0.3,0.7,1.0,1.3)) +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      theme(axis.text.x=element_text(size=14)) +
      theme(axis.text.y=element_text(size=14)) +
      theme(legend.text = element_text(size = 14)) 

##############################################################################

### RESIDUALs VERSES FITTED VALUES

plot(1,pch="",xlim=c(0,2.5),ylim=c(-1.5,1.5),ylab="",xlab="")

for(i in 1:570) {
      points(pars$residuals[i,]~comm_data_recovery$comm_recovery,col="blue")
      
}
abline(h=0,col = "lightgray", lty=3,lwd=5)

##############################################################################

### MARGINAL AND CONDITIONAL R SQUARE


hist(pars$marginalR2,ylab="",xlab="",main="")
mean(pars$marginalR2); quantile(pars$marginalR2, probs = c(0.025,0.975)) #  variance explained by fixed effects
hist(pars$conditionalR2 - pars$marginalR2,ylab="",xlab="",main="") #variance explained by random effects
mean(pars$conditionalR2 - pars$marginalR2); quantile(pars$conditionalR2 - pars$marginalR2, probs = c(0.025,0.975)) 

hist(pars$conditionalR2,ylab="",xlab="",main="") #variance explained by both fixed and random effects
mean(pars$conditionalR2); quantile(pars$conditionalR2, probs = c(0.025,0.975)) 
