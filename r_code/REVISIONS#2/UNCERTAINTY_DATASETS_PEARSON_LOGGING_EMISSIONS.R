### PROCESS DATASETS FOR ALLOMETRIC UNCERTAINTY

###############################################################################

###PACKAGES
library('ggplot2')
library('dplyr')
source("r_code/functions.R")
library('R2jags')
library('reshape2')
source("r_code/generalBAYESIANcode.R")
library('data.table')

###############################################################################


###READ IN MODEL OUTPUTS
load("outputs/CHAVE_HEIGHT_AGB_SU.RData")
k_draws <- gipar(CHAVE_HEIGHT_AGB_SU)

###EXTRACT OUT BIOMASS ESIMATES
agb_pred <- as.data.frame(t(k_draws[["AGB_SU"]])) #1350 draws
dim(agb_pred)

### CONVERT TO CARBON AND Mg
carbon <- matrix(NA,nrow=nrow(agb_pred),ncol=ncol(agb_pred))

for (i in 1:nrow(carbon)) {
      for (j in 1:ncol(carbon)) {
            carbon[i,j] <- (agb_pred[i,j]/1000) * 0.47
      }
}

### save datafile
### write.csv(carbon,"outputs/carbon_pred.csv",row.names = FALSE)



### LOAD THE SURINAME DATA
suriname_data <- read.csv("outputs/kabo_clean_long.csv")
str(suriname_data)


### MERGE AGB PREDICT WITH SURINAME DATA
suriname_pred <- cbind(suriname_data[,c("PLOT","census")],carbon)
dim(suriname_pred)


### SUMMARISE BY PLOT AND CENSUS
carbon_sum <- suriname_pred %>%
            group_by(PLOT,census) %>%
            summarise_each(funs(sum))

###############################################################################
### WEIGHT RECOVERY BY PRE-LOGGED VALUES BASED ON PEARSON (2014) LOGGING 
### EMISSIONS FACTOR ASSESSED FOR GUYANA [2.33 MG C M3]

### Drop 1978 and 1980 census [commercial species censused only]
carbon_plots <- as.data.frame(carbon_sum[which(!carbon_sum$census %in% c(1978,1980)),])

### ESTIMATE PRE-LOGGED CARBON STOCKS BASED ON LOGGING INTENSITY
carbon_1981<- as.data.frame(carbon_plots[which(carbon_plots$census %in% c(1981)),])

###1. add logging intensity
carbon_1981$log_m3 <- log_intensity_m3(carbon_1981$PLOT)

###2. create empty dataframe to store initial values
carbon_initial <- matrix(NA,nrow=nrow(carbon_1981),ncol=ncol(agb_pred))

###3. estimate initial/pre-logged carbon based on emissions factor reported in Pearson 2014
for (i in 1:nrow(carbon_initial)) {
      for (j in 1:ncol(carbon_initial)) {
            carbon_initial[i,j] <- (carbon_1981[i,1353] * 1.52) + carbon_1981[i,j+2]
      }
}


###4. add plot id to initial carbon estimate
carbon_intial2 <- cbind(carbon_1981$PLOT,as.data.frame(carbon_initial))
names(carbon_intial2)[1] <- "PLOT"


###5. logged plots only
carbon_plots_logged <- carbon_plots[which(!carbon_plots$PLOT %in% c(41,42,43)),]
###6. drop census id columns
carbon_plots_logged2 <- carbon_plots_logged %>%
      select(-census)

###7. convert to data tables
setDT(carbon_plots_logged2)
setDT(carbon_intial2)

###
names(carbon_intial2) <- colnames(carbon_plots_logged2) #ensure colnames match

##8. create vector to match plot id
xcols <-names(carbon_plots_logged2)[-1]
icols <- paste0("i.",xcols)

#9. join tables and weight by initial pre-logged carbon estimate by plot
rec.rate.logged <- carbon_plots_logged2[carbon_intial2,on="PLOT",Map('/',mget(xcols),mget(icols)),by=.EACHI]


#10. add back census column
rec.rate.logged2 <- cbind(carbon_plots_logged[,c("census")],rec.rate.logged)
names(rec.rate.logged2)[1] <- "census"

###############################################################################
### PROCESS CONTROL PLOTS SEPARATELY: FIRST CENSUS = 1983
carbon_plots_control <- carbon_plots[which(carbon_plots$PLOT %in% c(41,42,43)),]
### drop census id
carbon_plots_control2 <- carbon_plots_control %>%
      select(-census)

### initial carbon
carbon_control_initial <- carbon_plots_control[which(carbon_plots_control$census %in% c(1983)),]
## drop census id
carbon_control_initial2 <- carbon_control_initial %>%
      select(-census)

### convert to data tables
setDT(carbon_plots_control2)
setDT(carbon_control_initial2)


## create vector to match plot id
xcols <-names(carbon_plots_control2)[-1]
icols <- paste0("i.",xcols)

### join tables and weight by initial pre-logged carbon estimate by plot
rec.rate.control <- carbon_plots_control2[carbon_control_initial2,on="PLOT",Map('/',mget(xcols),mget(icols)),by=.EACHI]

#10. add back census column
rec.rate.control2 <- cbind(carbon_plots_control[,c("census")],rec.rate.control)
names(rec.rate.control2)[1] <- "census"

##############################################################################
## combined logged and control weighted datatables

rec.rate.weighted <- as.data.frame(rbind(rec.rate.control2,rec.rate.logged2))

### add time since logged
rec.rate.weighted$tsl <- time.logged(census.yr = rec.rate.weighted$census, plot.id = rec.rate.weighted$PLOT)

#add logging intensity; volume of timber harvest
rec.rate.weighted$log_intensity <- log_intensity_m3(rec.rate.weighted$PLOT)

##############################################################################
##############################################################################

### JAGS MODEL CODE

### Model: response weighted ACS recovery 
#covariates: (1) logging intensity - fixed effect (continuous variable; volume_m3 (ha)
#            (2) time since logged - fixed effect (continuous variable - year)
#            (3) plot level random effects
#            (4) interaction term between time since logged and logging intensity

cat("model { #ACS recovery rate
    
    #Uninformative priors
    sd_plot ~ dunif(0,10) #variance of plot to plot variation
    intercept ~ dnorm(0,1e-06) #intercept term (baseline;control plots)
    
    tau ~ dgamma (0.001,0.001) #precision of residuals
    tau_plot <- 1/(sd_plot*sd_plot)
    
    logtrt ~ dnorm(0,1e-06) #uninformative prior for logging effect
    tsl ~ dnorm(0,1e-06) #uninformative prior for time since logged
    beta3 ~ dnorm(0,1e-06) #uninformative prior for interaction term


    ### Mariginal and conditional R2
    Vfixed <- (sd(pred_ACS_recovery))^2 
    Vresidual<-1/tau # get the variance of residuals
    Vrandom<-1/tau_plot  # get the variance of random plot effect
    marginalR2<-Vfixed/(Vfixed + Vrandom + Vresidual) # calculate marginalR2
    conditionalR2<-(Vrandom + Vfixed)/(Vfixed + Vrandom + Vresidual) # calculate conditional R2
    
    for (i in 1:12) # for each plot
    {
    plot_Eff[i] ~ dnorm (0,tau_plot) # random plot effects
    }
    
    
    for (i in 1:Nobs) 
    { #ACS recovery for each plot at each census interval post logging weighted by 1983 ACS values from control plots
    
    #likelihood
    mu[i] <- intercept + logtrt * logged_m3[i] + tsl * time.since.logged[i] + plot_Eff[plot_id[i]] + beta3 * logged_m3[i]*time.since.logged[i]
    
    #observed volume recovery
    ACS_recovery[i] ~ dnorm(mu[i],tau)
    
    #posterior predictive distribution
    pred_ACS_recovery[i]~dnorm(mu[i],tau)
    # Residuals
    residuals[i] <-  ACS_recovery[i] - pred_ACS_recovery[i]


        }
    }", fill=TRUE, file="outputs/kabo_model1_ACS.txt")




##############################################################################

###JAGS data
plot_id <- as.factor(rec.rate.weighted$PLOT) # plot id
time.since.logged <- as.numeric(rec.rate.weighted$tsl) #time since logging
ACS_recovery <- rec.rate.weighted %>%
      select(-c(PLOT,census,tsl,log_intensity))# weighted carbon recovery rates
logged_m3 <- as.numeric(rec.rate.weighted$log_intensity) #volume_m3 logged
Nobs <- nrow(rec.rate.weighted) # no. observations


# PARAMS TO SAVE
params=c("intercept","logtrt","tsl","plot_Eff","tau","tau_plot","pred_ACS_recovery","beta3","marginalR2","conditionalR2","residuals")
# NUMBER OF POSTERIOR DRAWS
post.draws <- draws(n.iter =300000 ,n.burnin =20000,n.thin = 2000 ,n.chains = 3) ###=1350


### OBJECTS TO STORE MODEL RUNS PARAMETER ESTIMATES
recov_intercept<-vector("list",length=post.draws) # length in number of posterior draws
recov_logtrt<-vector("list",length=post.draws)
recov_beta3 <- vector("list",length=post.draws)
recov_tau <- vector("list",length=post.draws)
recov_pred_ACS <- vector("list",length=post.draws)
recov_tsl<-vector("list",length=post.draws)
marginalR2 <- vector("list",length=post.draws)
conditionalR2<-vector("list",length=post.draws)
residual <- vector("list",length=post.draws)

### OBJECT TO STORE MODEL DIC VALUES
#dic_value <- vector("list",length = post.draws)#posterior draws) ??



### FOR LOOP FOR TO RUN JAGS MODEL FOR EACH UNCERTAINTY DATASETS

for(i in 1:ncol(ACS_recovery)) { # N DRAWS OF THE UNCERTAINTY
      
      #jags data
      jags.data<-list(logged_m3=logged_m3,time.since.logged=time.since.logged,plot_id=plot_id,ACS_recovery=ACS_recovery[,i],Nobs=Nobs)
      
      
      #Run the model and call the results
      kabo_model1_ACS<-jags(data=jags.data,  inits=NULL,parameters.to.save=params, "outputs/kabo_model1_ACS.txt", n.chains=3,n.iter=200000,n.burnin=20000,n.thin=2000)
      
      #ONLY STORE MODELS WITH GOOD CONVERGENCE
      if(any(Rhatty(kabo_model1_ACS)>1.1)) {next}
      
      else{
      recov_pred_ACS[[i]]<-gipar(kabo_model1_ACS)$pred_ACS_recovery
      recov_intercept[[i]]<-gipar(kabo_model1_ACS)$intercept
      recov_logtrt[[i]]<-gipar(kabo_model1_ACS)$logtrt
      recov_tsl[[i]]<-gipar(kabo_model1_ACS)$tsl
      recov_beta3[[i]]<-gipar(kabo_model1_ACS)$beta3
      recov_tau[[i]]<-gipar(kabo_model1_ACS)$tau
      marginalR2[[i]]<-gipar(kabo_model1_ACS)$marginalR2
      conditionalR2[[i]]<-gipar(kabo_model1_ACS)$conditionalR2
      residual[[i]] <- gipar(kabo_model1_ACS)$residuals
      
      
            
      }
      
}


pred_ACS <- do.call("rbind",recov_pred_ACS); write.csv(pred_ACS,"outputs/pred_ACS.csv",row.names=FALSE)
intercept <- do.call("rbind",recov_intercept); write.csv(pred_ACS,"outputs/intercept.csv",row.names=FALSE)
log_effect <- do.call("rbind",recov_logtrt); write.csv(pred_ACS,"outputs/log_effect.csv",row.names=FALSE)
log_tsl <- do.call("rbind",recov_tsl); write.csv(pred_ACS,"outputs/log_tsl.csv",row.names=FALSE)
interact_X <- do.call("rbind",recov_beta3); write.csv(pred_ACS,"outputs/interact_X.csv",row.names=FALSE)
model_tau <- do.call("rbind",recov_tau); write.csv(pred_ACS,"outputs/model_tau.csv",row.names=FALSE)
marginalR2_output <- do.call("rbind",marginalR2); write.csv(marginalR2_output,"outputs/marginalR2_output.csv",row.names=FALSE)
conditionalR2_output<- do.call("rbind",conditionalR2); write.csv(conditionalR2_output,"outputs/conditionalR2_output.csv",row.names=FALSE)
residual_output<- do.call("rbind",residual); write.csv(conditionalR2_output,"outputs/residual_output.csv",row.names=FALSE)

### READ IN MODEL OUTPUTS
#pred_ACS <- read.csv("outputs/pred_ACS.csv")
#intercept <- read.csv("outputs/intercept.csv")
#log_effect <- read.csv("outputs/log_effect.csv")
#log_tsl <- read.csv("outputs/log_tsl.csv")
#interact_X <- read.csv("outputs/interact_X.csv")
#model_tau <- read.csv("outputs/model_tau.csv")
#residual_output <- read.csv("outputs/residual_output.csv")
#marginalR2_output <- read.csv("outputs/marginalR2_output.csv")
#conditionalR2_output <- read.csv("outputs/conditionalR2_output.csv")


post_quant <- matrix(NA,nrow=Nobs,ncol=3) # empty dataframe to store 0.025 - 97.75

for (i in 1:Nobs){
      post_quant[i,1] <- quantile(pred_ACS[,i],0.025)
      post_quant[i,2] <- mean(pred_ACS [,i],0.5)
      post_quant[i,3] <- quantile(pred_ACS [,i],0.975)
}
colnames(post_quant) <- c("lower","model_mean","upper")


#time since logged
tsl.seq=seq(from=0.1,to=70,length.out = 100)
nsim=nrow(log_effect)

### CONTROL PREDICTIONS
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 1,2, & 3 volume_m3 (m2 ha)
      temp <- rnorm(nsim,intercept + log_effect * 0.0 + (interact_X * 0.0 * tsl.seq[i]) +log_tsl * tsl.seq[i],sd=1/sqrt(model_tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco[i,] <- quantile(temp,c(0.46))
}

control_predictions <- yquantile
control_pred_mean <- yhat
predict_contr <- as.data.frame(cbind(control_predictions ,control_pred_mean,tsl.seq))
colnames(predict_contr)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")


### LOW INTENSITY

yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 1,2, & 3 volume_m3 (m2 ha)
      temp <- rnorm(nsim,intercept + log_effect * 15.0 + (interact_X * 15.0 * tsl.seq[i]) +log_tsl * tsl.seq[i],sd=1/sqrt(model_tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco[i,] <- quantile(temp,c(0.52))
}

low_predictions <- yquantile
low_pred_mean <- yhat
predict_low <- as.data.frame(cbind(low_predictions ,low_pred_mean,tsl.seq))
colnames(predict_low)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")

### MEDIUM INTENSITY

yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 1,2, & 3 volume_m3 (m2 ha)
      temp <- rnorm(nsim,intercept + log_effect * 23.0 + (interact_X * 23.0 * tsl.seq[i]) +log_tsl * tsl.seq[i],sd=1/sqrt(model_tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco[i,] <- quantile(temp,c(0.59))
}

med_predictions <- yquantile
med_pred_mean <- yhat
predict_med <- as.data.frame(cbind(med_predictions ,med_pred_mean,tsl.seq))
colnames(predict_med)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")

### HIGH INTENSITY
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 1,2, & 3 volume_m3 (m2 ha)
      temp <- rnorm(nsim,intercept + log_effect * 46.0 + (interact_X * 46.0 * tsl.seq[i]) +log_tsl * tsl.seq[i],sd=1/sqrt(model_tau)) 
      yquantile[i,] <- quantile(temp,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp)
      prob_no_reco[i,] <- quantile(temp,c(0.75))
}

high_predictions <- yquantile
high_pred_mean <- yhat
predict_high <- as.data.frame(cbind(high_predictions ,high_pred_mean,tsl.seq))
colnames(predict_high)[1:5] <- c("q0.025","q0.5","q0.975","mean","tsl")



#png(file = "outputs/REVISIONS2/high_uncertainty.png", bg = "transparent", type = c("cairo"), width=2000, height=2000, res=300)

###Plot out predictions
ggplot() + 
      #geom_line(data = predict_contr, aes(x=tsl, y=mean),stat="identity",colour="darkgreen",size=1) +
      #geom_ribbon(data=predict_contr,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="darkgreen") +
      #geom_line(data = predict_low, aes(x=tsl, y=mean),stat="identity",colour="royalblue2",size=1) +
      #geom_ribbon(data=predict_low,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="royalblue2") +
      #geom_line(data = predict_med, aes(x=tsl, y=mean),stat="identity",colour="gold3",size=1) +
      #geom_ribbon(data=predict_med,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="gold3") +
      geom_line(data = predict_high, aes(x=tsl, y=mean),stat="identity",colour="#e66101",size=1) +
      geom_ribbon(data=predict_high,aes(x=tsl,ymin=q0.025,ymax=q0.975),alpha=0.2,fill="#e66101") +
      geom_hline(aes(yintercept=1),size=2,linetype=2,colour="black") +
      scale_y_continuous(limits=c(0.3,1.3),breaks=c(0.3,0.7,1.0,1.3)) +
      theme_bw() +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      theme(axis.text.x=element_text(size=14)) +
      theme(axis.text.y=element_text(size=14)) +
      #theme(axis.title.y = element_text (size=14)) +
      #theme(axis.title.x = element_text (size=14)) +
      #theme(strip.text.x = element_text(size = 14)) +
      theme(legend.text = element_text(size = 14)) +
      ylab(expression("")) +
      xlab(expression("")) + 
      theme(axis.text.x = element_blank()) +
      theme(axis.text.y = element_blank())
#xlab (expression(bold("Time since logged (years)"))) 
#theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) 

#dev.off()

#time since logged
tsl.seq=seq(from=0.1,to=70,length.out = 100)
nsim=nrow(log_effect)

### CONTROL PREDICTIONS
yquantile <- matrix(NA,100,3)
yhat <- (matrix(NA,100,1))
prob_no_reco <- (matrix(NA,100,1))

for (i in 1:100) { #fixed logging intensity 1,2, & 3 volume_m3 (m2 ha)
      temp <- rnorm(nsim,intercept + log_effect * 25.0 + (interact_X * 25.0 * tsl.seq[i]) +log_tsl * tsl.seq[i],sd=1/sqrt(model_tau)) 
      temp2 <- temp * 184.95
      yquantile[i,] <- quantile(temp2,c(0.025,0.5,0.975))
      yhat[i,] <- mean(temp2)
      prob_no_reco[i,] <- quantile(temp2,c(0.46))
}



##############################################################################
##############################################################################
### COEFFICIENT PLOT

coef_points_model <- list(LOG_INTENSITY = log_effect, TIME_SINCE_LOGGED =log_tsl,INTERACTION=interact_X)
pplotCOEF(coef_points_model)

##############################################################################

### RESIDUALs VERSES FITTED VALUES

### observed values
rec.rate.weighted_fitted <- rec.rate.weighted %>%
      select(-c(census,PLOT,tsl,log_intensity))
rec.rate.weighted_colmean <- apply(rec.rate.weighted_fitted,1,mean)

### truncate residual outputs; too many draws
residual_output2 <- residual_output[1:570,]

plot(1,pch="",xlim=c(0,1.5),ylim=c(-1.0,1.0),ylab="",xlab="")

for(i in 1:570) {
      points(residual_output2[i,]~rec.rate.weighted_colmean,col="blue")
      
}
abline(h=0,col = "lightgray", lty=3,lwd=5)

##############################################################################
##############################################################################

### MARGINAL AND CONDITIONAL R SQUARE


hist(marginalR2_output,ylab="",xlab="",main="")
mean(marginalR2_output); quantile(marginalR2_output, probs = c(0.025,0.975)) #  variance explained by fixed effects
hist(conditionalR2_output - marginalR2_output,ylab="",xlab="",main="") #variance explained by random effects
mean(conditionalR2_output - marginalR2_output); quantile(conditionalR2_output - marginalR2_output, probs = c(0.025,0.975)) 

hist(conditionalR2_output,ylab="",xlab="",main="") #variance explained by both fixed and random effects
mean(conditionalR2_output); quantile(conditionalR2_output, probs = c(0.025,0.975)) 

##############################################################################
##############################################################################

### ESTIMATE CARBON LOSSES IMMEDIATELY AFTER LOGGING

ACS_loss <- rec.rate.logged2[which(rec.rate.logged2$census == 1981),]
str(ACS_loss)

### add logging intensity
ACS_loss$log_cat <- log_intensity_cat(ACS_loss$PLOT)

### remove plot and census columns
ACS_loss2 <- ACS_loss %>%
      select(-c(PLOT,census))

### summarise by logging intensity
ACS_loss_logging <- ACS_loss2 %>%
      group_by(log_cat) %>%
      summarise_each(funs(mean))

apply(ACS_loss_logging[,2:1351],1,mean)
