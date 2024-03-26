# This code recreates the figures and analyses from See et al. 2024. Decadal increases in carbon uptake offset by respiratory losses across northern permafrost ecosystems

library(lmerTest)
library(MuMIn)
library(ggplot2)
library(sjPlot)
library(nlme)
library(brms)
library(bayestestR)
library(cv)


##################################################################################################################
#                    FIGURE TWO                        
# #########################################################################################################################
# ##BASED ON BELSHE ANALYSIS CODE http://bbolker.github.io/mixedmodels-misc/ecostats_chap.html

Ann<-read.csv("AnnualFluxes.csv")
Ann$mat<-as.numeric(Ann$mat)
NEEall<-subset(Ann,is.na(nee)==F)
NEEall$meas_year<-as.numeric(NEEall$meas_year)

##################################
# ANNUAL PERMAFROST NEE TIMESERIES
BelNEE1<-lme(nee~meas_year,
            data=subset(NEEall,permafrost=="Yes"), method = "REML",
            random = ~1 +meas_year|study_id,
            correlation=corCAR1(form=~meas_year|study_id),
            control=list(maxIter=100000, niterEM=100000)
) 
summary(BelNEE1)
tab_model(BelNEE1)
intervals(BelNEE1,which="fixed")
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
#If "cv" function begins throwing errors, restarting R, masking all other packages besides "lme" and "cv," and rerun.
cv(BelNEE1,clusterVariables="study_id",data=subset(NEEall,permafrost=="Yes"))### Leave-one-out cross validation

#  Bayesian refit of annual permafrsot NEE model
summary(BelNEE1)$sigma
priors <- c(prior(normal(-3371, 3090), class = Intercept), # intercept prior based on model summary
            prior(normal(1.67, 1.54), class = b), # slope prior based on model summary
            prior(normal(48.7, 10), class = sigma) # population variance based on model
)
BelNEE1_bayes<-brm(nee~meas_year + (1 + meas_year|study_id)+ar(p=1),
                    cores = 8, chains = 8, iter=20000, warmup = 12000, 
                    data=subset(NEEall,permafrost=="Yes"), prior=priors,
                    control = list(adapt_delta=0.9999))
summary(BelNEE1_bayes)
p_map(BelNEE1_bayes)
plot(BelNEE1_bayes)

#QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelNEE1,type=c("p","smooth")),
             plot(BelNEE1,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelNEE1,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelNEE1,abline=c(0,1)                    ))

####################################
#ANNUAL NONPERMAFROST NEE TIMESERIES
BelNEE2<-lme(nee~meas_year,
            data=subset(NEEall, permafrost=="No"), method = "REML",
            random = ~1 +meas_year|study_id,
            correlation=corCAR1(form=~meas_year|study_id),
            control=list(maxIter=100000, niterEM=100000)
) 
summary(BelNEE2)
tab_model(BelNEE2)
intervals(BelNEE2,which="fixed")
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelNEE2,clusterVariables="study_id",data=subset(NEEall,permafrost=="No"))

#  Bayesian refit of annual nonpermafrost NEE model
summary(BelNEE2)$sigma
priors <- c(prior(normal(9608, 4928), class = Intercept), # intercept prior based on model
            prior(normal(-4.8, 2.46), class = b), # slope prior based on model
            prior(normal(57.9, 10), class = sigma) # population variance based on model
)
BelNEE2_bayes<-brm(nee~meas_year + (1 + meas_year|study_id)+ar(p=1),
                    cores = 8, chains = 8, iter=20000, warmup = 12000, 
                    data=subset(NEEall,permafrost=="No"), prior=priors,
                    control = list(adapt_delta=0.9999))
summary(BelNEE2_bayes)
p_map(BelNEE2_bayes)
plot(BelNEE2_bayes)

#QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelNEE2,type=c("p","smooth")),
             plot(BelNEE2,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelNEE2,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelNEE2,abline=c(0,1)                    ))

################
# FIGURE 2b NEE
nonperm<-expand.grid(meas_year=seq(1995,2022,by=0.1),permafrost="No")
nonperm$pred <- predict(BelNEE2, nonperm, level = 0)
Designmat <- model.matrix(formula(BelNEE2)[-2], nonperm)
predvar <- diag(Designmat %*% vcov(BelNEE2) %*% t(Designmat))
nonperm$SE <- sqrt(predvar)
nonperm$SE2 <- sqrt(predvar+BelNEE2$sigma^2)

perm<-expand.grid(meas_year=seq(1995,2022,by=0.1),permafrost="Yes")
perm$pred <- predict(BelNEE1, perm, level = 0)
Designmat <- model.matrix(formula(BelNEE1)[-2], perm)
predvar <- diag(Designmat %*% vcov(BelNEE1) %*% t(Designmat))
perm$SE <- sqrt(predvar)
perm$SE2 <- sqrt(predvar+BelNEE1$sigma^2)

newdat<-as.data.frame(rbind(nonperm,perm))

# Graph Theme
longtheme<-theme(
  plot.title=element_text(size=10,hjust=.5, vjust=0),
    panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_rect(color="black", fill="white"),
  panel.border=element_blank(),
  legend.position = "none",#c(.9,.9),
  legend.background = element_rect(colour=NA), 
  legend.key=element_rect(fill="white"),
  axis.text=element_text(size=9, color='black'),
  #axis.title.x=element_blank(),
  axis.title.y=element_text(size=9, color='black'))


pd <- position_dodge(width=0.4)
Fig2b <- ggplot(newdat,aes(x=meas_year,y=pred,colour=permafrost))+
  geom_point(position=pd,size = .8) +scale_color_manual(values = c("No"="orange",
                                                                  "Yes"="turquoise3") )+
  geom_linerange(aes(ymin=pred-1.96*SE,ymax=pred+1.96*SE,size = 0.2), position=pd,alpha=0.05)+
  #geom_point(data=NEEall,(aes(meas_year,nee,colour=permafrost,size = 0.15)))+
  xlab('Year')+ylab(bquote(Annual~CO[2]~Flux~"(g"~C~m^-2~")"))+
  xlim(1994,2022)+ylim(-210,110)+ggtitle("Annual NEE")+guides(size = FALSE)+
  longtheme
Fig2b
ggsave("PaperFigs/2b.pdf", Fig2b,width = 88, height = 88, dpi=300, units="mm")

#######################################
# SUMMERTIME PERMAFROST NEE TIMESERIES

GSall<-read.csv("SummerFluxes.csv")#Load summertime fluxes
NEEall_GS<-subset(GSall,is.na(nee)==F)

BelNEE1_GS<-lme(nee~meas_year,
             data=subset(NEEall_GS,permafrost=="Yes"), method = "REML",
             random = ~1 +meas_year|study_id,
             correlation=corCAR1(form=~meas_year|study_id),
             control=list(maxIter=10000, niterEM=10000)
) 
summary(BelNEE1_GS)
tab_model(BelNEE1_GS)
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
### Leave-one-out cross validation increased max iterations to maxIter=100000 and niterEM=100000 to ensure convergence 
cv(BelNEE1_GS,clusterVariables="study_id",data=subset(NEEall_GS,permafrost=="Yes"))

#  Bayesian refit of summertime permafrost NEE model
summary(BelNEE1_GS)$sigma
priors <- c(prior(normal(6026, 1652), class = Intercept), # intercept prior based on model
            prior(normal(-3, 0.82), class = b), # slope prior based on model
            prior(normal(43.7, 10), class = sigma) # population variance based on model
)

BelNEE1_GS_bayes<-brm(nee~meas_year + (1 + meas_year|study_id)+ar(p=1),
                   cores = 8, chains = 8, 
                   iter=20000, warmup = 12000, 
                   data=subset(NEEall_GS,permafrost=="Yes"), prior=priors,
                   control = list(adapt_delta=0.9999))
summary(BelNEE1_GS_bayes)
p_map(BelNEE1_GS_bayes)
plot(BelNEE1_GS_bayes)

#QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelNEE1_GS,type=c("p","smooth")),
             plot(BelNEE1_GS,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelNEE1_GS,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelNEE1_GS,abline=c(0,1)                    ))

##########################################
# SUMMERTIME NON-PERMAFROST NEE TIMESERIES
BelNEE2_GS<-lme(nee~meas_year,
             data=subset(NEEall_GS,permafrost=="No"), method = "REML",
             random = ~1 +meas_year|study_id,
             correlation=corCAR1(form=~meas_year|study_id),
             control=list(maxIter=10000, niterEM=10000)
)
summary(BelNEE2_GS)
tab_model(BelNEE2_GS)
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelNEE2_GS,clusterVariables="study_id",data=subset(NEEall_GS,permafrost=="No"))

#Bayesian refit of summertime non-permafrost NEE model
summary(BelNEE2_GS)$sigma
priors <- c(prior(normal(5141, 1840), class = Intercept), # intercept prior based on model
            prior(normal(-2.6, 0.92), class = b), # slope prior based on model
            prior(normal(43.6, 10), class = sigma) # population variance based on model
)

BelNEE2_GS_bayes<-brm(nee~meas_year + (1 + meas_year|study_id)+ar(p=1),
                   cores = 8, chains = 8, 
                   iter=20000, warmup = 12000, 
                   data=subset(NEEall_GS,permafrost=="No"), prior=priors,
                   control = list(adapt_delta=0.999))
summary(BelNEE2_GS_bayes)
p_map(BelNEE2_GS_bayes)
plot(BelNEE2_GS_bayes)

#QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelNEE2_GS,type=c("p","smooth")),
             plot(BelNEE2_GS,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelNEE2_GS,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelNEE2_GS,abline=c(0,1)                    ))


################
# FIGURE 2a NEE
nonperm_GS<-expand.grid(meas_year=seq(1992,2022,by=0.1),permafrost="No")
nonperm_GS$pred <- predict(BelNEE2_GS, nonperm_GS, level = 0)
Designmat <- model.matrix(formula(BelNEE2_GS)[-2], nonperm_GS)
predvar <- diag(Designmat %*% vcov(BelNEE2_GS) %*% t(Designmat)) 
nonperm_GS$SE <- sqrt(predvar) 
nonperm_GS$SE2 <- sqrt(predvar+BelNEE2_GS$sigma^2)

perm_GS<-expand.grid(meas_year=seq(1992,2022,by=0.1),permafrost="Yes")
perm_GS$pred <- predict(BelNEE1_GS, perm_GS, level = 0)
Designmat <- model.matrix(formula(BelNEE1_GS)[-2], perm_GS)
predvar <- diag(Designmat %*% vcov(BelNEE1_GS) %*% t(Designmat)) 
perm_GS$SE <- sqrt(predvar) 
perm_GS$SE2 <- sqrt(predvar+BelNEE1_GS$sigma^2)

newdat_GS<-as.data.frame(rbind(nonperm_GS,perm_GS))

library(ggplot2)
pd <- position_dodge(width=0.4)
Fig2a <- ggplot(newdat_GS,aes(x=meas_year,y=pred,colour=permafrost))+
  geom_point(position=pd,size = .8) +scale_color_manual(values = c("No"="orange",
                                                                   "Yes"="turquoise3") )+
  geom_linerange(aes(ymin=pred-1.96*SE,ymax=pred+1.96*SE,size = 0.2), position=pd,alpha=0.05)+
  #geom_point(data=NEEall,(aes(meas_year,nee,colour=permafrost,size = 0.15)))+
  xlab('Year')+ylab(bquote(Summer~CO[2]~Flux~"(g"~C~m^-2~")"))+
  xlim(1994,2022)+ylim(-210,110)+ggtitle("Summer NEE")+guides(size = FALSE)+
  longtheme
Fig2a
ggsave("PaperFigs/2a.pdf", Fig2a,width = 88, height = 88, dpi=300, units="mm")

######################
###  COMBINED FIGURE 2
fig2<-ggarrange(Fig2a, Fig2b,ncol = 2)
fig2
ggsave("PaperFigs/Figure2.jpeg", fig2,width = 180, height = 88, dpi=300, units="mm")




##################################################################################################################
#                    FIGURE THREE                       
# #########################################################################################################################

####################################
##  ANNUAL PERMAFROST GPP TIMESERIES
GPPall<-subset(NEEall,is.na(gpp)==F)

BelGPP1<-lme(gpp~meas_year,
         data=subset(GPPall,permafrost=="Yes"), method = "REML",
         random = ~1 +meas_year|study_id,
         correlation=corCAR1(form=~meas_year|study_id),
         control=list(maxIter=100000, niterEM=100000)
) 
summary(BelGPP1)
tab_model(BelGPP1)
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelGPP1,clusterVariables="study_id",data=subset(GPPall,permafrost=="Yes"))### Leave-one-out cross validation
  
##Bayesian refit of annual permafrost GPP model
summary(BelGPP1)$sigma
priors <- c(prior(normal(12131, 5817), class = Intercept), # intercept prior based on model
            prior(normal(-6.27, 2.9), class = b), # slope prior based on model
            prior(normal(77.37, 10), class = sigma) # population variance based on model
    )
BelGPP1_bayes<-brm(gpp~meas_year + (1 + meas_year|study_id)+ar(p=1),
                        cores = 8, chains = 8, iter=20000, warmup = 12000, 
                        data=subset(GPPall,permafrost=="Yes"), prior=priors,
                        control = list(adapt_delta=0.9999))
summary(BelGPP1_bayes)
p_map(BelGPP1_bayes)
plot(BelGPP1_bayes)
    
##QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
    grid.arrange(plot(BelGPP1,type=c("p","smooth")),
                 plot(BelGPP1,sqrt(abs(resid(.)))~fitted(.),
                      # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                      type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
                 ## "sqrt(abs(resid(x)))"),
                 plot(BelGPP1,resid(.,type="pearson")~meas_year,
                      type=c("p","smooth")),
                 qqnorm(BelGPP1,abline=c(0,1)                    ))

#######################################
# ANNUAL NON_PERMAFROST GPP TIMESERIES
BelGPP2<-lme(gpp~meas_year,
                 data=subset(GPPall,permafrost=="No"), method = "REML",
                 random = ~1 +meas_year|study_id,
                 correlation=corCAR1(form=~meas_year|study_id),
                 control=list(maxIter=100000, niterEM=10000)
    )
summary(BelGPP2)
tab_model(BelGPP2)
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelGPP2,clusterVariables="study_id",data=subset(GPPall,permafrost=="No"))### Leave-one-out cross validation

##Bayesian refit of annual non-permafrost GPP model
summary(BelGPP2)$sigma
  priors <- c(prior(normal(12711, 11487), class = Intercept), # intercept prior based on model
              prior(normal(-6.7, 5.7), class = b), # slope prior based on model
              prior(normal(155, 10), class = sigma) # population variance based on model
  )
BelGPP2_bayes<-brm(gpp~meas_year + (1 + meas_year|study_id)+ar(p=1),
                   cores = 8, chains = 8, iter=20000, warmup = 12000, 
                   data=subset(GPPall,permafrost=="No"), prior=priors,
                   control = list(adapt_delta=0.9999))
summary(BelGPP2_bayes)
p_map(BelGPP2_bayes)
plot(BelGPP2_bayes)    
    
##QQ Plots of Frequentist Model Residuals   
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
    grid.arrange(plot(BelGPP2,type=c("p","smooth")),
                 plot(BelGPP2,sqrt(abs(resid(.)))~fitted(.),
                      # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                      type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
                 ## "sqrt(abs(resid(x)))"),
                 plot(BelGPP2,resid(.,type="pearson")~meas_year,
                      type=c("p","smooth")),
                 qqnorm(BelGPP2,abline=c(0,1)                    ))
###########
# FIGURE 3c

nonperm_GPP_Annual<-expand.grid(meas_year=seq(1995,2022,by=0.1),permafrost="No")
nonperm_GPP_Annual$pred <- predict(BelGPP2, nonperm_GPP_Annual, level = 0)
Designmat <- model.matrix(formula(BelGPP2)[-2], nonperm_GPP_Annual)
predvar <- diag(Designmat %*% vcov(BelGPP2) %*% t(Designmat))
nonperm_GPP_Annual$SE <- sqrt(predvar)
nonperm_GPP_Annual$SE2 <- sqrt(predvar+BelGPP2$sigma^2)

perm_GPP_Annual<-expand.grid(meas_year=seq(1995,2022,by=0.1),permafrost="Yes")
perm_GPP_Annual$pred <- predict(BelGPP1, perm_GPP_Annual, level = 0)
Designmat <- model.matrix(formula(BelGPP1)[-2], perm_GPP_Annual)
predvar <- diag(Designmat %*% vcov(BelGPP1) %*% t(Designmat))
perm_GPP_Annual$SE <- sqrt(predvar)
perm_GPP_Annual$SE2 <- sqrt(predvar+BelGPP1$sigma^2)

newdat_GPP_Annual<-as.data.frame(rbind(nonperm_GPP_Annual,perm_GPP_Annual))
    
pd <- position_dodge(width=0.4)
Fig3c <- ggplot(newdat_GPP_Annual,aes(x=meas_year,y=pred,colour=permafrost))+
    geom_point(position=pd,size = .8) +scale_color_manual(values = c("No"="orange",
                                                                       "Yes"="turquoise3") )+
    geom_linerange(aes(ymin=pred-1.96*SE,ymax=pred+1.96*SE,size = 0.2), position=pd,alpha=0.08)+
    #geom_point(data=NEEall,(aes(meas_year,nee,colour=permafrost,size = 0.15)))+
    xlab('Year')+ylab(bquote(Annual~CO[2]~Flux~"(g"~C~m^-2~")"))+
    xlim(1994,2022)+ylim(-1200,-200)+ggtitle("Annual GPP")+guides(size = FALSE)+
    longtheme
Fig3c
ggsave("PaperFigs/3c.pdf", Fig3c,width = 88, height = 88, dpi=300, units="mm")
 

########################################
##  SUMMERTIME PERMAFROST GPP TIMESERIES    
GPPall_GS<-subset(GSall,is.na(gpp)==F)
    
BelGPP1_GS<-lme(gpp~meas_year,
             data=subset(GPPall_GS,permafrost=="Yes"), method = "REML",
             random = ~1 +meas_year|study_id,
             correlation=corCAR1(form=~meas_year|study_id),
             control=list(maxIter=100000, niterEM=100000)
    ) 
summary(BelGPP1_GS)
tab_model(BelGPP1_GS)
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelGPP1_GS,clusterVariables="study_id",data=subset(GPPall_GS,permafrost=="Yes"))

##Bayesian refit of  summertime permafrost GPP model
summary(BelGPP1_GS)$sigma
priors <- c(prior(normal(13382, 4187), class = Intercept), # intercept prior based on model
          prior(normal(-6.8, 2), class = b), # slope prior based on model
          prior(normal(85.3, 10), class = sigma) # population variance based on model
    )
BelGPP1_GS_bayes<-brm(gpp~meas_year + (1 + meas_year|study_id)+ar(p=1),
                   cores = 8, chains = 8, iter=20000, warmup = 12000, 
                   data=subset(GPPall_GS,permafrost=="Yes"), prior=priors,
                   control = list(adapt_delta=0.9999))
summary(BelGPP1_GS_bayes)
p_map(BelGPP1_GS_bayes)
plot(BelGPP1_GS_bayes)
 
##QQ Plots of Frequentist Model Residuals 
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelGPP1_GS,type=c("p","smooth")),
             plot(BelGPP1_GS,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelGPP1_GS,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelGPP1_GS,abline=c(0,1)                    ))
    
############################################
##  SUMMERTIME NON-PERMAFROST GPP TIMESERIES    
BelGPP2_GS<-lme(gpp~meas_year,
           data=subset(GPPall_GS,permafrost=="No"), method = "REML",
           random = ~1 +meas_year|study_id,
           correlation=corCAR1(form=~meas_year|study_id),
           control=list(maxIter=100000, niterEM=100000)
    ) 
summary(BelGPP2_GS)
tab_model(BelGPP2_GS)
intervals(BelGPP2_GS,which="fixed")
cv(BelGPP2_GS,clusterVariables="study_id",data=subset(GPPall_GS,permafrost=="No"))### Failed to converge even at max itereations= 100,000
# a<-subset(GPPall_GS,permafrost=="No") #calculate MSE for model in same manner as cv package (fixed effects only)
# a$pred<- 2874.605+(a$meas_year*-1.643)
# sqrt(mse(a$gpp,a$pred))

##Bayesian refit of  summertime non-permafrost GPP model
summary(BelGPP2_GS)$sigma
priors <- c(prior(normal(2875, 6717), class = Intercept), # intercept prior based on model
            prior(normal(-1.6, 3.3), class = b), # slope prior based on model
            prior(normal(105, 10), class = sigma) # population variance based on model
  )
BelGPP2_GS_bayes<-brm(gpp~meas_year + (1 + meas_year|study_id)+ar(p=1),
           cores = 8, chains = 8, iter=20000, warmup = 12000, 
           data=subset(GPPall_GS,permafrost=="No"), prior=priors,
           control = list(adapt_delta=0.9999))
summary(BelGPP2_GS_bayes)
p_map(BelGPP2_GS_bayes)
plot(BelGPP2_GS_bayes)    

#QQ Plots of Frequentist Model Residuals 
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelGPP2_GS,type=c("p","smooth")),
         plot(BelGPP2_GS,sqrt(abs(resid(.)))~fitted(.),
        # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
        type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
       ## "sqrt(abs(resid(x)))"),
         plot(BelGPP2_GS,resid(.,type="pearson")~meas_year,
        type=c("p","smooth")),
         qqnorm(BelGPP2_GS,abline=c(0,1)                    ))
    
###########
# FIGURE 3a
nonperm_GPP_GS<-expand.grid(meas_year=seq(1992,2022,by=0.1),permafrost="No")
nonperm_GPP_GS$pred <- predict(BelGPP2_GS, nonperm_GPP_GS, level = 0)
Designmat <- model.matrix(formula(BelGPP2_GS)[-2], nonperm_GPP_GS)
predvar <- diag(Designmat %*% vcov(BelGPP2_GS) %*% t(Designmat)) 
nonperm_GPP_GS$SE <- sqrt(predvar) 
    #nonperm$SE2 <- sqrt(predvar+BelGPP2$sigma^2)
    
perm_GPP_GS<-expand.grid(meas_year=seq(1992,2022,by=0.1),permafrost="Yes")
perm_GPP_GS$pred <- predict(BelGPP1_GS, perm_GPP_GS, level = 0)
Designmat <- model.matrix(formula(BelGPP1_GS)[-2], perm_GPP_GS)
predvar <- diag(Designmat %*% vcov(BelGPP1_GS) %*% t(Designmat)) 
perm_GPP_GS$SE <- sqrt(predvar) 
#perm_GPP_GS$SE2 <- sqrt(predvar+BelGPP1$sigma^2)

newdat_GPP_GS<-as.data.frame(rbind(nonperm_GPP_GS,perm_GPP_GS))
    
pd <- position_dodge(width=0.4)
Fig3a <- ggplot(newdat_GPP_GS,aes(x=meas_year,y=pred,colour=permafrost))+
    geom_point(position=pd,size = .8) +scale_color_manual(values = c("No"="orange",
                   "Yes"="turquoise3") )+
    geom_linerange(aes(ymin=pred-1.96*SE,ymax=pred+1.96*SE,size = 0.2), position=pd,alpha=0.08)+
      #geom_point(data=NEEall,(aes(meas_year,nee,colour=permafrost,size = 0.15)))+
    xlab('Year')+ylab(bquote(Summer~CO[2]~Flux~"(g"~C~m^-2~")"))+
    xlim(1994,2022)+ylim(-620,20)+ggtitle("Summer GPP")+guides(size = FALSE)+
    longtheme
Fig3a
ggsave("PaperFigs/3a.pdf", Fig3c,width = 88, height = 88, dpi=300, units="mm")

    
############################################
##  ANNUAL PERMAFROST RECO TIMESERIES 
RecoAll<-subset(NEEall,is.na(reco)==F)

BelReco1<-lme(reco~meas_year,
             data=subset(RecoAll,permafrost=="Yes"), method = "REML",
             random = ~1 +meas_year|study_id,
             correlation=corCAR1(form=~meas_year|study_id),
             control=list(maxIter=100000, niterEM=100000)) #exact model structure without weights, because already aggregated at site level in ABC/Santa Barb.
summary(BelReco1)
tab_model(BelReco1)
intervals(BelReco1,which="fixed")
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelReco1,clusterVariables="study_id",data=subset(RecoAll,permafrost=="Yes"))##Failure to converge
 # a<-subset(RecoAll,permafrost=="Yes") #calculate MSE for model in same manner as cv package (fixed effects only)
 # a$pred<- -11758.338+(a$meas_year*6.079)
 # sqrt(mse(a$reco,a$pred))

##Bayesian refit of  annual permafrost Reco model
summary(BelReco1)$sigma
priors <- c(prior(normal(-11758, 6457), class = Intercept), # intercept prior based on model
            prior(normal(6.1, 3.2), class = b), # slope prior based on model
            prior(normal(87, 10), class = sigma) # population variance based on model
)

BelReco1_bayes<-brm(reco~meas_year + (1 + meas_year|study_id)+ar(p=1),
                    cores = 8, chains = 8, iter=20000, warmup = 12000, 
                    data=subset(RecoAll,permafrost=="Yes"), prior=priors,
                    control = list(max_treedepth = 15,adapt_delta=0.9999))
summary(BelReco1_bayes)
p_map(BelReco1_bayes)
plot(BelReco1_bayes)  

# QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelReco1,type=c("p","smooth")),
             plot(BelReco1,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelReco1,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelReco1,abline=c(0,1)                    ))

#######################################
# ANNUAL NON-PERMAFROST RECO TIMESERIES 
BelReco2<-lme(reco~meas_year,
              data=subset(RecoAll,permafrost=="No"), method = "REML",
              random = ~1 +meas_year|study_id,
              correlation=corCAR1(form=~meas_year|study_id),
              control=list(maxIter=1000000, niterEM=1000000)) 
summary(BelReco2)
tab_model(BelReco2)
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
### one or more LOOCV models fail to converge
cv(BelReco2,clusterVariables="study_id",data=subset(RecoAll,permafrost=="No"))
# a<-subset(RecoAll,permafrost=="No") 
# a$pred<- -7629.582+(a$meas_year*4.156)
# sqrt(mse(a$reco,a$pred))

##Bayesian refit of  annual non-permafrost Reco model
summary(BelReco2)$sigma
priors <- c(prior(normal(-7630, 9267), class = Intercept), # intercept prior based on model
            prior(normal(4.2, 4.6), class = b), # slope prior based on model
            prior(normal(147, 10), class = sigma) # population variance based on model
)
BelReco2_bayes<-brm(reco~meas_year + (1 + meas_year|study_id)+ar(p=1),
                   cores = 8, chains = 8, iter=20000, warmup = 12000, 
                   data=subset(RecoAll,permafrost=="No"), prior=priors,
                   control = list(max_treedepth = 15,adapt_delta=0.9999))
summary(BelReco2_bayes)
p_map(BelReco2_bayes)
plot(BelReco2_bayes)  

# QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelReco2,type=c("p","smooth")),
             plot(BelReco2,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelReco2,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelReco2,abline=c(0,1)                    ))

############
#FIGURE 3d 
nonperm_Reco_Annual<-expand.grid(meas_year=seq(1995,2022,by=0.1),permafrost="No")
nonperm_Reco_Annual$pred <- predict(BelReco2, nonperm_Reco_Annual, level = 0)
Designmat <- model.matrix(formula(BelReco2)[-2], nonperm_Reco_Annual)
predvar <- diag(Designmat %*% vcov(BelReco2) %*% t(Designmat))
nonperm_Reco_Annual$SE <- sqrt(predvar)
nonperm_Reco_Annual$SE2 <- sqrt(predvar+BelReco2$sigma^2)

perm_Reco_Annual<-expand.grid(meas_year=seq(1995,2022,by=0.1),permafrost="Yes")
perm_Reco_Annual$pred <- predict(BelReco1, perm_Reco_Annual, level = 0)
Designmat <- model.matrix(formula(BelReco1)[-2], perm_Reco_Annual)
predvar <- diag(Designmat %*% vcov(BelReco1) %*% t(Designmat))
perm_Reco_Annual$SE <- sqrt(predvar)
perm_Reco_Annual$SE2 <- sqrt(predvar+BelGPP1$sigma^2)

newdat_Reco_Annual<-as.data.frame(rbind(nonperm_Reco_Annual,perm_Reco_Annual))

Fig3d <- ggplot(newdat_Reco_Annual,aes(x=meas_year,y=pred,colour=permafrost))+
  geom_point(position=pd,size = .8) +scale_color_manual(values = c("No"="orange",
                                                                   "Yes"="turquoise3") )+
  geom_linerange(aes(ymin=pred-1.96*SE,ymax=pred+1.96*SE,size = 0.2), position=pd,alpha=0.08)+
  #geom_point(data=NEEall,(aes(meas_year,nee,colour=permafrost,size = 0.15)))+
  xlab('Year')+ylab(bquote(Annual~CO[2]~Flux~"(g"~C~m^-2~")"))+
  xlim(1994,2022)+ylim(200,1200)+ggtitle(bquote(Annual~R[eco]))+guides(size = FALSE)+
  longtheme
Fig3d
ggsave("PaperFigs/3d.pdf", Fig3d,width = 88, height = 88, dpi=300, units="mm")


############################################
##  SUMMERTIME PERMAFROST RECO TIMESERIES 
RecoAll_GS<-subset(GSall,is.na(reco)==F)
BelReco1_GS<-lme(reco~meas_year,
              data=subset(RecoAll_GS,permafrost=="Yes"), method = "REML",
              random = ~1 +meas_year|study_id,
              correlation=corCAR1(form=~meas_year|study_id),
              control=list(maxIter=100000, niterEM=100000)) 
summary(BelReco1_GS)  
tab_model(BelReco1_GS)
intervals(BelReco1_GS,which="fixed")
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
cv(BelReco1_GS,clusterVariables="study_id",data=subset(RecoAll_GS,permafrost=="Yes"))# One or more models failed to converge at 100,00 iterations
a<-subset(RecoAll_GS,permafrost=="Yes") 
a$pred<- -6001.526+(a$meas_year*3.076)
sqrt(mse(a$reco,a$pred))

##Bayesian refit of  summertime permafrost Reco model
summary(BelReco1_GS)$sigma
priors <- c(prior(normal(-6001, 2193), class = Intercept), # intercept prior based on model
            prior(normal(3.1, 1.1), class = b), # slope prior based on model
            prior(normal(48.5, 10), class = sigma) # population variance based on model
)
BelReco1_GS_bayes<-brm(reco~meas_year + (1 + meas_year|study_id)+ar(p=1),
                    cores = 8, chains = 8, iter=20000, warmup = 12000, 
                    data=subset(GPPall_GS,permafrost=="Yes"), prior=priors,
                    control = list(adapt_delta=0.9999))
summary(BelReco1_GS_bayes)
p_map(BelReco1_GS_bayes)
plot(BelReco1_GS_bayes)    

# QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelReco1_GS,type=c("p","smooth")),
             plot(BelReco1_GS,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelReco1_GS,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelReco1_GS,abline=c(0,1)                    ))

############################################
##  SUMMERTIME NON_PERMAFROST RECO TIMESERIES 
BelReco2_GS<-lme(reco~meas_year,
              data=subset(RecoAll_GS,permafrost=="No"), method = "REML",
              random = ~1 +meas_year|study_id,
              correlation=corCAR1(form=~meas_year|study_id),
              control=list(maxIter=100000, niterEM=100000)) 
summary(BelReco2_GS)    
tab_model(BelReco2_GS)
intervals(BelReco2_GS,which="fixed")
#retrieve MSE for full model and for LOOCV model (calculate square root for RMSE in table S5)
###Convergence failure at 100,000 iterations 
cv(BelReco2_GS,clusterVariables="study_id",data=subset(RecoAll_GS,permafrost=="No"))
# a<-subset(RecoAll_GS,permafrost=="No") 
# a$pred<- -1222.7028+(a$meas_year*0.7772)
# sqrt(mse(a$reco,a$pred))

##Bayesian refit of  summertime non-permafrost Reco model
summary(BelReco2_GS)$sigma
priors <- c(prior(normal(-1223, 4788), class = Intercept), # intercept prior based on model
            prior(normal(0.8, 2.4), class = b), # slope prior based on model
            prior(normal(72.5, 10), class = sigma) # population variance based on model
)

BelReco2_GS_bayes<-brm(reco~meas_year + (1 + meas_year|study_id)+ar(p=1),
                    cores = 8, chains = 8, iter=20000, warmup = 12000, 
                    data=subset(GPPall_GS,permafrost=="No"), prior=priors,
                    control = list(adapt_delta=0.9999))
summary(BelReco2_GS_bayes)
p_map(BelReco2_GS_bayes)
plot(BelReco2_GS_bayes)    

# QQ Plots of Frequentist Model Residuals
colvec <- c("#ff1111","#007eff") ## second colour matches lattice default
grid.arrange(plot(BelReco2_GS,type=c("p","smooth")),
             plot(BelReco2_GS,sqrt(abs(resid(.)))~fitted(.),
                  # col=ifelse(mc2$Site=="Toolik, AK",colvec[1],colvec[2]),
                  type=c("p","smooth"),ylab=expression(sqrt(abs(resid)))),
             ## "sqrt(abs(resid(x)))"),
             plot(BelReco2_GS,resid(.,type="pearson")~meas_year,
                  type=c("p","smooth")),
             qqnorm(BelReco2_GS,abline=c(0,1)                    ))

############
# FIGURE 3b
nonperm_Reco_GS<-expand.grid(meas_year=seq(1992,2022,by=0.1),permafrost="No")
nonperm_Reco_GS$pred <- predict(BelReco2_GS, nonperm_Reco_GS, level = 0)
Designmat <- model.matrix(formula(BelReco2_GS)[-2], nonperm_Reco_GS)
predvar <- diag(Designmat %*% vcov(BelReco2_GS) %*% t(Designmat)) 
nonperm_Reco_GS$SE <- sqrt(predvar) 
nonperm_Reco_GS$SE2 <- sqrt(predvar+BelReco2_GS$sigma^2)

perm_Reco_GS<-expand.grid(meas_year=seq(1992,2022,by=0.1),permafrost="Yes")
perm_Reco_GS$pred <- predict(BelReco1_GS, perm_Reco_GS, level = 0)
Designmat <- model.matrix(formula(BelReco1_GS)[-2], perm_Reco_GS)
predvar <- diag(Designmat %*% vcov(BelReco1_GS) %*% t(Designmat)) 
perm_Reco_GS$SE <- sqrt(predvar) 
perm_Reco_GS$SE2 <- sqrt(predvar+BelReco1_GS$sigma^2)

newdat_Reco_GS<-as.data.frame(rbind(nonperm_Reco_GS,perm_Reco_GS))

library(ggplot2)
pd <- position_dodge(width=0.4)
Fig3b <- ggplot(newdat_Reco_GS,aes(x=meas_year,y=pred,colour=permafrost))+
  geom_point(position=pd,size = .8) +scale_color_manual(values = c("No"="orange",
                                                                   "Yes"="turquoise3") )+
  geom_linerange(aes(ymin=pred-1.96*SE,ymax=pred+1.96*SE,size = 0.2), position=pd,alpha=0.08)+
  #geom_point(data=NEEall,(aes(meas_year,nee,colour=permafrost,size = 0.15)))+
  xlab('Year')+ylab(bquote(Summer~CO[2]~Flux~"(g"~C~m^-2~")"))+
  xlim(1994,2022)+ylim(20,620)+ggtitle(bquote(Summer~R[eco]))+guides(size = FALSE)+
  longtheme
Fig3b
ggsave("PaperFigs/3b.pdf", Fig3b,width = 88, height = 88, dpi=300, units="mm")

#####################
#COMBINED FIGURE 3
fig3<-ggarrange(Fig3a, Fig3b,Fig3c, Fig3d,ncol = 2,nrow=2)
fig3
ggsave("PaperFigs/Figure3.jpeg", fig3,width = 180, height = 180, dpi=300, units="mm")


# #########################################################################################################################
#                     FIGURE FOUR  
# #########################################################################################################################

library(lmerTest)
library(ggplot2)
library(nlme)
ABCraw<-read.csv("MonthlyFluxes.csv", header=T) #Read in monthly values


##################################
# Calculate GPP time series models by month

ABCraw$meas_year<-as.numeric(ABCraw$meas_year)

Borslopes<-data.frame(matrix(ncol = 7,nrow=12))
names(Borslopes)<-c("nee","neeSE","neeP","month","Perm","N","monLet")
Arcslopes<-data.frame(matrix(ncol = 7,nrow=12))
names(Arcslopes)<-c("nee","neeSE","neeP","month","Perm","N","monLet")
for (i in 1:12){
  mnthly=subset(ABCraw,interval_month==i)
  mnthly<-subset(mnthly,is.na(nee)==F)
  nee.fullmodel <-lme(nee~meas_year,
                      data=subset(mnthly,permafrost=="No"), method = "REML",na.action = na.exclude,
                      random = ~1 |study_id,
                      #correlation=corCAR1(form=~meas_year|study_id),
                      control=list(maxIter=100000, niterEM=100000)
  )
  Borslopes[i,1:6]<-c(coef(summary(nee.fullmodel))["meas_year","Value"],coef(summary(nee.fullmodel))["meas_year","Std.Error"],coef(summary(nee.fullmodel))["meas_year","p-value"],i,'No Permafrost',nrow(subset(mnthly,permafrost=="No")))
  nee.fullmodel_P  <-lme(nee~meas_year,
                         data=subset(mnthly,permafrost=="Yes"), method = "REML",
                         random = ~1 |study_id,
                         #random = ~1 +meas_year|study_id,
                         #correlation=corCAR1(form=~meas_year|study_id),
                         control=list(maxIter=100000, niterEM=100000)
  )
  Arcslopes[i,1:6]<-c(coef(summary(nee.fullmodel_P))["meas_year","Value"],coef(summary(nee.fullmodel_P))["meas_year","Std.Error"],coef(summary(nee.fullmodel_P))["meas_year","p-value"],i,'Permafrost',nrow(subset(mnthly,permafrost=="Yes")))
}
Borslopes$monLet<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
Arcslopes$monLet<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
slopes<-rbind(Borslopes,Arcslopes)
slopes$nee<-as.numeric(slopes$nee)
slopes$neeSE<-as.numeric(slopes$neeSE)
slopes$neeP<-as.numeric(slopes$neeP)
slopes$month<-as.numeric(slopes$month)
slopes$monLet <- reorder(slopes$monLet, slopes$month)

#write.csv(slopes,"fig2_nee_int.csv")
#slopes<-read.csv("fig2_nee_int.csv")
slopes$monLet <- reorder(slopes$monLet, slopes$month)

pNEE <- ggplot(slopes, aes(x=monLet, y=nee, fill=Perm)) + scale_fill_manual(values = c("No Permafrost"="orange",
                                                                                       "Permafrost"="turquoise3") )+
  geom_bar(stat="identity", position=position_dodge(0.9),colour="black") +
  geom_errorbar(aes(ymin=nee-neeSE, ymax=nee+neeSE), width=.2,position=position_dodge(0.9)
  )
Fig4a<-pNEE +theme_bw()+ xlab("")+ylab("")+# ggtitle("Change in NEE")+ ylab(bquote(CO[2]~Flux~"(g"~C~m^-2~yr^-1~")") ) +
  theme(legend.background = element_rect(fill="transparent"), 
        legend.key = element_rect(colour="transparent"), 
        legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        plot.title = element_text(hjust = 0.5)
  )
Fig4a
# ggsave("PaperFigs/4a.pdf", Fig4a)
 #ggsave("Final Figs/R/4a.jpg", Fig4a,width = 88, height = 60, dpi=300, units="mm")

##################################
# GPP time series models by month
Borslopes_GPP<-data.frame(matrix(ncol = 7,nrow=12))
names(Borslopes_GPP)<-c("gpp","gppSE","gppP","month","Perm","N","monLet")
Arcslopes_GPP<-data.frame(matrix(ncol = 7,nrow=12))
names(Arcslopes_GPP)<-c("gpp","gppSE","gppP","month","Perm","N","monLet")
for (i in 1:12){
  mnthly=subset(ABCraw,interval_month==i)
  mnthly<-subset(mnthly,is.na(gpp)==F)
  gpp.fullmodel <-lme(gpp~meas_year,
                      data=subset(mnthly,permafrost=="No"), method = "REML",
                      random = ~1 |study_id,
                      #correlation=corCAR1(form=~meas_year|study_id),
                      control=list(maxIter=100000, niterEM=100000)
  )
  Borslopes_GPP[i,1:6]<-c(coef(summary(gpp.fullmodel))["meas_year","Value"],coef(summary(gpp.fullmodel))["meas_year","Std.Error"],coef(summary(gpp.fullmodel))["meas_year","p-value"],i,'No Permafrost',nrow(subset(mnthly,permafrost=="No")))
  gpp.fullmodel_P  <-lme(gpp~meas_year,
                         data=subset(mnthly,permafrost=="Yes"), method = "REML",
                         random = ~1 |study_id,
                         #correlation=corCAR1(form=~meas_year|study_id),
                         control=list(maxIter=1000000, niterEM=1000000)
  )
  Arcslopes_GPP[i,1:6]<-c(coef(summary(gpp.fullmodel_P))["meas_year","Value"],coef(summary(gpp.fullmodel_P))["meas_year","Std.Error"],coef(summary(gpp.fullmodel_P))["meas_year","p-value"],i,'Permafrost',nrow(subset(mnthly,permafrost=="Yes")))
}
Borslopes_GPP$monLet<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
Arcslopes_GPP$monLet<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
slopes_GPP<-rbind(Borslopes_GPP,Arcslopes_GPP)
slopes_GPP$gpp<-as.numeric(slopes_GPP$gpp)
slopes_GPP$gppSE<-as.numeric(slopes_GPP$gppSE)
slopes_GPP$gppP<-as.numeric(slopes_GPP$gppP)
slopes_GPP$month<-as.numeric(slopes_GPP$month)
slopes_GPP$monLet <- reorder(slopes_GPP$monLet, slopes_GPP$month)
slopes_GPP$Perm <- ifelse(slopes_GPP$Perm=="No Permafrost","Non-permafrost", slopes_GPP$Perm)
write.csv(slopes_GPP,"fig2_gpp_int.csv")
# slopes_GPP<-read.csv("fig2_gpp_int.csv")

slopes_GPP$monLet <- reorder(slopes_GPP$monLet, slopes_GPP$month)
pGPP <- ggplot(slopes_GPP, aes(x=monLet, y=gpp, fill=Perm)) + scale_fill_manual(values = c("Non-permafrost"="orange",
                                                                                           "Permafrost"="turquoise3") )+
  geom_bar(stat="identity", position=position_dodge(0.9),colour="black") +
  geom_errorbar(aes(ymin=gpp-gppSE, ymax=gpp+gppSE), width=.2,position=position_dodge(0.9)
  )+ ylim(-3.25,1)
Fig4b<-pGPP +theme_bw() + xlab("")+ylab("")+#+ylab(bquote(CO[2]~Flux~"(g"~C~m^-2~yr^-1~")")) + xlab("")+ ggtitle("Change in GPP")+
  theme(legend.background = element_rect(fill="transparent"), 
        legend.key = element_rect(colour="transparent"), 
        legend.title = element_blank(), 
        legend.position = c(0.25, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        plot.title = element_text(hjust = 0.5)
  )
Fig4b
# ggsave("PaperFigs/4b.pdf", Fig4b)
# ggsave("Final Figs/R/4b.jpg", Fig4b,width = 88, height = 60, dpi=300, units="mm")

##################################
# Reco time series models by month
Borslopes_Reco<-data.frame(matrix(ncol = 7,nrow=12))
names(Borslopes_Reco)<-c("reco","recoSE","recoP","month","Perm","N","monLet")
Arcslopes_Reco<-data.frame(matrix(ncol = 7,nrow=12))
names(Arcslopes_Reco)<-c("reco","recoSE","recoP","month","Perm","N","monLet")
for (i in 1:12){
  mnthly=subset(ABCraw,interval_month==i)
  mnthly<-subset(mnthly,is.na(reco)==F)
  reco.fullmodel <-lme(reco~meas_year,
                       data=subset(mnthly,permafrost=="No"), method = "REML",
                       random = ~1|study_id,
                       #correlation=corCAR1(form=~meas_year|study_id),
                       control=list(maxIter=10000000, niterEM=10000,msTol=1e-6)
  )
  Borslopes_Reco[i,1:6]<-c(coef(summary(reco.fullmodel))["meas_year","Value"],coef(summary(reco.fullmodel))["meas_year","Std.Error"],coef(summary(reco.fullmodel))["meas_year","p-value"],i,'No Permafrost',nrow(subset(mnthly,permafrost=="No")))
  reco.fullmodel_P  <-lme(reco~meas_year,
                          data=subset(mnthly,permafrost=="Yes"), method = "REML",
                          random = ~1 |study_id,
                          #correlation=corCAR1(form=~meas_year|study_id),
                          control=list(maxIter=10000000, niterEM=10000,msTol=1e-6)
  )
  Arcslopes_Reco[i,1:6]<-c(coef(summary(reco.fullmodel_P))["meas_year","Value"],coef(summary(reco.fullmodel_P))["meas_year","Std.Error"],coef(summary(reco.fullmodel_P))["meas_year","p-value"],i,'Permafrost',nrow(subset(mnthly,permafrost=="Yes")))
}
Borslopes_Reco$monLet<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
Arcslopes_Reco$monLet<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
slopes_Reco<-rbind(Borslopes_Reco,Arcslopes_Reco)
slopes_Reco$reco<-as.numeric(slopes_Reco$reco)
slopes_Reco$recoSE<-as.numeric(slopes_Reco$recoSE)
slopes_Reco$recoP<-as.numeric(slopes_Reco$recoP)
slopes_Reco$month<-as.numeric(slopes_Reco$month)
slopes_Reco$monLet <- reorder(slopes_Reco$monLet, slopes_Reco$month)
write.csv(slopes_Reco,"fig2_reco_int.csv")
#slopes_Reco<-read.csv("fig2_reco_int.csv")

#slopes_Reco$monLet <- reorder(slopes_Reco$monLet, slopes_Reco$month)
p <- ggplot(slopes_Reco, aes(x=monLet, y=reco, fill=Perm)) + scale_fill_manual(values = c("No Permafrost"="orange",
                                                                                          "Permafrost"="turquoise3") )+
  geom_bar(stat="identity", position=position_dodge(0.9),colour="black") +
  geom_errorbar(aes(ymin=reco-recoSE, ymax=reco+recoSE), width=.2,position=position_dodge(0.9)
  ) + ylim(-1,3.25)
Fig4c<-p +theme_bw()+xlab("")+ylab("")+#ylab(bquote(CO[2]~Flux~"(g"~C~m^-2~yr^-1~")")) + xlab("")+ggtitle(bquote(Change~"in"~R[eco]))+
  theme(legend.background = element_rect(fill="transparent"), 
        legend.key = element_rect(colour="transparent"), 
        legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        plot.title = element_text(hjust = 0.5)
  )
Fig4c
# ggsave("PaperFigs/4c.pdf", Fig4c)
 #ggsave("Final Figs/R/4c.jpg", Fig4c,width = 88, height = 60, dpi=300, units="mm")

####################
# COMBINED FIGURE 4
fig4<-ggarrange(Fig4a, Fig4b,Fig4c,ncol = 1)
fig4
#ggsave("PaperFigs/Figure4.jpeg", fig4,width = 88, height = 190, dpi=300, units="mm")


####################################################################################################################
#        FIGURE 5
####################################################################################################################
library(ggplot2)
library(sjPlot)
library(gridExtra) 
library(nlme)
library(boot)
library (emmeans)
c<-read.csv("TempSlopes.csv")

###########################
#NEE-TEMP by Factor
cc<-subset(c,is.na(neeT_slp)==F) #Get Rid of NA values for slopes of temperature vs NEE. 

#PERMAFROST COMPARISON
permNEE1<-summary(lm(neeT_slp~permafrost-1, weights = 1/neeT_SE, data = cc))
permNEE1
table(cc$permafrost)#Get sample sizes 
permNEE<-emmeans(lm(neeT_slp~permafrost-1, weights = 1/neeT_SE, data = cc), pairwise ~ permafrost)
#BIOME NEE COMPARISON
biomeNEE1<-summary(lm(neeT_slp~biome-1, weights = 1/neeT_SE, data = cc))
biomeNEE1
tab_model(biomeNEE1)
table(cc$biome)#Get sample sizes 
biomeNEE<-emmeans(lm(neeT_slp~biome-1, weights = 1/neeT_SE, data = cc),pairwise~biome)
biomeNEE
#Soil CN COMPARISON
NlimNEE1<-summary(lm((neeT_slp)~Nlim5-1, weights = 1/neeT_SE, data = cc))
NlimNEE1
tab_model(NlimNEE1)
table(cc$Nlim5)#Get sample sizes 
NlimNEE<-emmeans(lm(neeT_slp~Nlim5-1, weights = 1/neeT_SE, data = cc),pairwise~Nlim5)
NlimNEE

fxNEE<-data.frame(matrix(ncol = 6,nrow=6))
names(fxNEE)<-c("var","factor","estimate","lower.CL","upper.CL","pvalue")
fxNEE[1,]<-c("unfrozen","perm",summary(permNEE$emmeans)[1,2],summary(permNEE$emmeans)[1,5],summary(permNEE$emmeans)[1,6],summary(permNEE$contrasts)[1,6])
fxNEE[2,]<-c("permafrost","perm",summary(permNEE$emmeans)[2,2],summary(permNEE$emmeans)[2,5],summary(permNEE$emmeans)[2,6],summary(permNEE$contrasts)[1,6])
fxNEE[3,]<-c("boreal","biome",summary(biomeNEE$emmeans)[1,2],summary(biomeNEE$emmeans)[1,5],summary(biomeNEE$emmeans)[1,6],summary(biomeNEE$contrasts)[1,6])
fxNEE[4,]<-c("tundra","biome",summary(biomeNEE$emmeans)[2,2],summary(biomeNEE$emmeans)[2,5],summary(biomeNEE$emmeans)[2,6],summary(biomeNEE$contrasts)[1,6])
fxNEE[5,]<-c("C limited","nuts",summary(NlimNEE$emmeans)[1,2],summary(NlimNEE$emmeans)[1,5],summary(NlimNEE$emmeans)[1,6],summary(NlimNEE$contrasts)[1,6])
fxNEE[6,]<-c("N limited","nuts",summary(NlimNEE$emmeans)[2,2],summary(NlimNEE$emmeans)[2,5],summary(NlimNEE$emmeans)[2,6],summary(NlimNEE$contrasts)[1,6])
fxNEE$var<-factor(fxNEE$var,levels=c("unfrozen","permafrost","tundra","boreal","C limited","N limited","dry","wet"))
fxNEE$estimate<-as.numeric(fxNEE$estimate)
fxNEE$lower.CL<-as.numeric(fxNEE$lower.CL)
fxNEE$upper.CL<-as.numeric(fxNEE$upper.CL)

fxNEE$var<-reorder(fxNEE$var,c(6,5,4,3,2,1))#MAKE NEE FIGURE
NEEfx<-ggplot(fxNEE, aes(x=estimate, y=var)) +geom_point(size=4,col="chartreuse3")+
  geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL), width=.1,size = 1.3,col="chartreuse3") +
  geom_vline(xintercept=0, linetype='dotted', col = 'black',size=0.75)+labs(x ="", y = "")+scale_x_continuous(limits=c(-0.4,0.4),breaks=seq(-0.4,0.4,0.2))
NEErr<-NEEfx+theme_bw() +theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"),
                               legend.title = element_blank(),
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               #axis.text.y=element_text(colour="black"),
                               axis.text.y=element_blank(),
                               plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                               panel.background = element_rect(fill='transparent'), #transparent panel bg
                               plot.background = element_rect(fill='transparent', color=NA)
)
NEErr
#ggsave("PaperFigs/5a_long.pdf", plot=g0,height = 5,width = 3,units = "in")
#ggsave("PaperFigs/5a_long.jpg", plot=g0,height = 4.5,width = 3,units = "in")

#####################
#GPP-TEMP BY FACTOR
dd<-subset(c,is.na(gppT_slp)==F) #Get Rid of NA values for GPP inthis portion of analysis

#PERMAFROST COMPARISON
permGPP1<-summary(lm(gppT_slp~permafrost-1, weights = 1/gppT_SE, data = dd))
permGPP1
table(dd$permafrost)#Get sample sizes 
permGPP<-emmeans(lm(gppT_slp~permafrost-1, weights = 1/gppT_SE, data = dd), pairwise ~ permafrost)
permGPP
#BIOME GPP COMPARISON
biomeGPP1<-summary(lm(gppT_slp~biome-1, weights = 1/gppT_SE, data = dd))
biomeGPP1
table(dd$biome)#Get sample sizes 
biomeGPP<-emmeans(lm(gppT_slp~biome-1, weights = 1/gppT_SE, data = dd),pairwise~biome)
biomeGPP
#Soil CN COMPARISON
NlimGPP1<-summary(lm((gppT_slp)~Nlim5-1, weights = 1/gppT_SE, data = dd))
NlimGPP1
table(dd$Nlim5)#Get sample sizes 
NlimGPP<-emmeans(lm(gppT_slp~Nlim5-1, weights = 1/gppT_SE, data = dd),pairwise~Nlim5)
NlimGPP

fxGPP<-data.frame(matrix(ncol = 6,nrow=6))
names(fxGPP)<-c("var","factor","estimate","lower.CL","upper.CL","pvalue")
fxGPP[1,]<-c("unfrozen","perm",summary(permGPP$emmeans)[1,2],summary(permGPP$emmeans)[1,5],summary(permGPP$emmeans)[1,6],summary(permGPP$contrasts)[1,6])
fxGPP[2,]<-c("permafrost","perm",summary(permGPP$emmeans)[2,2],summary(permGPP$emmeans)[2,5],summary(permGPP$emmeans)[2,6],summary(permGPP$contrasts)[1,6])
fxGPP[3,]<-c("boreal","biome",summary(biomeGPP$emmeans)[1,2],summary(biomeGPP$emmeans)[1,5],summary(biomeGPP$emmeans)[1,6],summary(biomeGPP$contrasts)[1,6])
fxGPP[4,]<-c("tundra","biome",summary(biomeGPP$emmeans)[2,2],summary(biomeGPP$emmeans)[2,5],summary(biomeGPP$emmeans)[2,6],summary(biomeGPP$contrasts)[1,6])
fxGPP[5,]<-c("C limited","nuts",summary(NlimGPP$emmeans)[1,2],summary(NlimGPP$emmeans)[1,5],summary(NlimGPP$emmeans)[1,6],summary(NlimGPP$contrasts)[1,6])
fxGPP[6,]<-c("N limited","nuts",summary(NlimGPP$emmeans)[2,2],summary(NlimGPP$emmeans)[2,5],summary(NlimGPP$emmeans)[2,6],summary(NlimGPP$contrasts)[1,6])
fxGPP$var<-factor(fxGPP$var,levels=c("unfrozen","permafrost","tundra","boreal","C limited","N limited","dry","wet"))
fxGPP$estimate<-as.numeric(fxGPP$estimate)
fxGPP$lower.CL<-as.numeric(fxGPP$lower.CL)
fxGPP$upper.CL<-as.numeric(fxGPP$upper.CL)

#MAKE GPP FIGURE
fxGPP$var<-reorder(fxGPP$var,c(6,5,4,3,2,1))

GPPfx<-ggplot(fxGPP, aes(x=estimate, y=var)) +geom_point(size=4,col="darkgoldenrod1")+
  geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL), width=.1,size = 1.3,col="darkgoldenrod1") +
  geom_vline(xintercept=0, linetype='dotted', col = 'black',size=0.75)+labs(x ="", y = "")+scale_x_continuous(limits=c(-0.6,0.6),breaks=seq(-0.6,0.6,0.3))
GPPrr<-GPPfx+theme_bw() +theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"),
                               legend.title = element_blank(),
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               #axis.text.y=element_text(colour="black"),
                               axis.text.y=element_blank(),
                               plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                               panel.background = element_rect(fill='transparent'), #transparent panel bg
                               plot.background = element_rect(fill='transparent', color=NA)
)
GPPrr
#ggsave("PaperFigs/5b_long.jpg",plot=g0,height = 4.5,width = 3,units = "in")

#####################
#Reco-TEMP BY FACTOR
ee<-subset(c,is.na(recoT_slp)==F) #Get Rid of Reco NA values for this portion of analysis

#PERMAFROST COMPARISON
permReco1<-summary(lm(recoT_slp~permafrost-1, weights = 1/recoT_SE, data = ee))
permReco1
table(ee$permafrost)#Get sample sizes 
permReco<-emmeans(lm(recoT_slp~permafrost-1, weights = 1/recoT_SE, data = ee), pairwise ~ permafrost)
permReco
#BIOME Reco COMPARISON
biomeReco1<-summary(lm(recoT_slp~biome-1, weights = 1/recoT_SE, data = ee))
biomeReco1
table(ee$biome)#Get sample sizes 
biomeReco<-emmeans(lm(recoT_slp~biome-1, weights = 1/recoT_SE, data = ee),pairwise~biome)
biomeReco
#Soil CN COMPARISON
NlimReco1<-summary(lm((recoT_slp)~Nlim5-1, weights = 1/recoT_SE, data = ee))
NlimReco1
table(ee$Nlim5)#Get sample sizes 
NlimReco<-emmeans(lm(recoT_slp~Nlim5-1, weights = 1/recoT_SE, data = ee),pairwise~Nlim5)
NlimReco

fxReco<-data.frame(matrix(ncol = 6,nrow=6))
names(fxReco)<-c("var","factor","estimate","lower.CL","upper.CL","pvalue")
fxReco[1,]<-c("unfrozen","perm",summary(permReco$emmeans)[1,2],summary(permReco$emmeans)[1,5],summary(permReco$emmeans)[1,6],summary(permReco$contrasts)[1,6])
fxReco[2,]<-c("permafrost","perm",summary(permReco$emmeans)[2,2],summary(permReco$emmeans)[2,5],summary(permReco$emmeans)[2,6],summary(permReco$contrasts)[1,6])
fxReco[3,]<-c("boreal","biome",summary(biomeReco$emmeans)[1,2],summary(biomeReco$emmeans)[1,5],summary(biomeReco$emmeans)[1,6],summary(biomeReco$contrasts)[1,6])
fxReco[4,]<-c("tundra","biome",summary(biomeReco$emmeans)[2,2],summary(biomeReco$emmeans)[2,5],summary(biomeReco$emmeans)[2,6],summary(biomeReco$contrasts)[1,6])
fxReco[5,]<-c("C limited","nuts",summary(NlimReco$emmeans)[1,2],summary(NlimReco$emmeans)[1,5],summary(NlimReco$emmeans)[1,6],summary(NlimReco$contrasts)[1,6])
fxReco[6,]<-c("N limited","nuts",summary(NlimReco$emmeans)[2,2],summary(NlimReco$emmeans)[2,5],summary(NlimReco$emmeans)[2,6],summary(NlimReco$contrasts)[1,6])
fxReco$var<-factor(fxReco$var,levels=c("unfrozen","permafrost","tundra","boreal","C limited","N limited","dry","wet"))
fxReco$estimate<-as.numeric(fxReco$estimate)
fxReco$lower.CL<-as.numeric(fxReco$lower.CL)
fxReco$upper.CL<-as.numeric(fxReco$upper.CL)
#MAKE Reco FIGURE
fxReco$var<-reorder(fxReco$var,c(6,5,4,3,2,1))
Recofx<-ggplot(fxReco, aes(x=estimate, y=var)) +geom_point(size=4,col="aquamarine2")+
  geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL), width=.1,size = 1.3,col="aquamarine2") +
  geom_vline(xintercept=0, linetype='dotted', col = 'black',size=0.75)+labs(x ="", y = "")+scale_x_continuous(limits=c(-0.6,0.6),breaks=seq(-0.6,0.6,0.3))
Recorr<-Recofx+theme_bw() +theme(legend.background = element_rect(fill="transparent"),
                                 legend.key = element_rect(colour="transparent"),
                                 legend.title = element_blank(),
                                 legend.position = "none",
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 axis.text = element_text(size = 9),
                                 axis.text.x=element_text(colour="black"),
                                 #axis.text.y=element_text(colour="black"),
                                 axis.text.y=element_blank(),
                                 plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                                 panel.background = element_rect(fill='transparent'), #transparent panel bg
                                 plot.background = element_rect(fill='transparent', color=NA) #transparent plot bg
)
Recorr
#ggsave("PaperFigs/5c_long.jpg", plot=g0,height = 4.5,width = 3,units = "in")

fig5<-ggarrange(NEErr, NULL,GPPrr,NULL, Recorr,nrow = 1, widths = c(1, -0.05,1,-0.05, 1))
fig5
ggsave("PaperFigs/Figure5.jpeg", fig5,width = 150, height = 90, dpi=300, units="mm")


####################################################################################################################
#        FIGURE 6
####################################################################################################################
c$H2Obal<-c$aetmean-c$pptmean # Calculate difference between summertime transpiration and precipitation 

#################################
#  WATER USE VS TEMP-FLUX SLOPES

BAL_NEE<-lm(neeT_slp~H2Obal, weights = 1/neeT_SE,data=c)
summary(BAL_NEE)
tab_model(BAL_NEE)

BAL_reco<-lm(recoT_slp~H2Obal, weights = 1/neeT_SE,data=c)
summary(BAL_reco)
tab_model(BAL_reco)

BAL_gpp<-lm(gppT_slp~H2Obal, weights = 1/neeT_SE,data=c)
summary(BAL_gpp)
tab_model(BAL_gpp)

## Water balance FIGURES
fig6a<-ggplot(c, aes(x=H2Obal, y=neeT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/neeT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                            "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  geom_smooth(aes(weight = 1/neeT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Summer water use index (mm)", y = "Temperature effect on NEE \n (standardized slope)")
fig6a

fig6a<-fig6a +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
                               plot.margin = unit(c(1,2,1,1), "mm")
)
fig6a
ggsave("Final Figs/R/6a_NEE.jpg", fig6a, width = 80, height = 60, dpi=300, units="mm")

#FIGURE 6B
fig6b<-ggplot(c, aes(x=H2Obal, y=gppT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/gppT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                            "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  geom_smooth(aes(weight = 1/gppT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Summer water use index (mm)", y = "Temperature effect on GPP \n (standardized slope)") 
fig6b

fig6b<-fig6b +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
                               plot.margin = unit(c(1,2,1,1), "mm")
)
fig6b
ggsave("Final Figs/R/6b_GPP.jpg", fig6b, width = 80, height = 60, dpi=300, units="mm")

#FIGURE 6C
fig6c<-ggplot(c, aes(x=H2Obal, y=recoT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/recoT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                             "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  geom_smooth(aes(weight = 1/recoT_SE), method = lm, size = 1,colour="black")+
  labs(x ="",y="")
  #labs(x ="Summer water use index (mm)", y = expression(atop("Temperature effect on R"[eco], paste("(standardized slope)")))) 
fig6c
fig6c<-fig6c +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
                               plot.margin = unit(c(1,2,1,1), "mm")
)
fig6c
#ggsave("PaperFigs/h2obal_Reco.jpg", fig6c)
ggsave("Final Figs/R/6c_GPP.jpg", fig6c, width = 80, height = 60, dpi=300, units="mm")

#FIGURE 6
fig6<-ggarrange(fig6a,fig6b,fig6c,ncol = 1)
fig6
#ggsave("PaperFigs/Figure6.jpeg", fig6,width = 88, height = 190, dpi=300, units="mm")

####################################################################################################################
#        EXTENDED DATA FIGURES
####################################################################################################################


#####################################################################################################################
# Extended Data Fig 1
## Models using soil C:N ratio (top 5 cm) to predict standardized slopes of and summer temperature vs fluxes. 
CN_GPP<-lm(gppT_slp~CN_5, weights = 1/neeT_SE,data=c)
summary(CN_GPP)
tab_model(CN_GPP)

CN_R<-lm(recoT_slp~CN_5, weights = 1/neeT_SE,data=c)
summary(CN_R)
tab_model(CN_R)

fig8a<-ggplot(c, aes(x=CN_5, y=gppT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/recoT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                             "Yes"="turquoise3") )+ 
  geom_smooth(aes(group=NA),colour="black",method = lm)+
  geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  labs(x ="Soil C:N ratio", y = "Temperature effect on GPP") 
fig8a<-fig8a +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                             legend.key = element_rect(colour="transparent"), 
                             legend.title = element_blank(), 
                             legend.position = "none",
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.text = element_text(size = 9),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"),
)
fig8a
ggsave("Final Figs/R/8a.jpg", fig8a, width = 80, height = 80, dpi=300, units="mm")

fig8b<-ggplot(c, aes(x=CN_5, y=recoT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/recoT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                             "Yes"="turquoise3") )+ 
  geom_smooth(aes(group=NA),colour="black",method = lm)+
  geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  labs(x ="Soil C:N ratio", y = "Temperature effect on Reco") 
fig8b<-fig8b +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                             legend.key = element_rect(colour="transparent"), 
                             legend.title = element_blank(), 
                             legend.position = "none",
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.text = element_text(size = 9),
                             axis.text.x=element_text(colour="black"),
                             axis.text.y=element_text(colour="black"),
)
fig8b
ggsave("Final Figs/R/8b.jpg", fig8b, width = 80, height = 80, dpi=300, units="mm")

#FIGURE 8
fig8<-ggarrange(fig8a, fig8b,ncol = 2)
fig8


###############################################################################################################################
# Extended Data Fig 2
## Models using 30 year mean summer AET and precipitation to predict standardized slopes of and summer temperature vs fluxes. 
AET_NEE<-lm(neeT_slp~aetmean, weights = 1/neeT_SE,data=c)
summary(AET_NEE)
tab_model(AET_NEE)

AET_reco<-lm(recoT_slp~aetmean, weights = 1/neeT_SE,data=c)
summary(AET_reco)
tab_model(AET_reco)

AET_gpp<-lm(gppT_slp~aetmean, weights = 1/neeT_SE,data=c)
summary(AET_gpp)
tab_model(AET_gpp)

MAP_NEE<-lm(neeT_slp~pptmean, weights = 1/neeT_SE,data=c)
summary(MAP_NEE)
tab_model(MAP_NEE)

MAP_reco<-lm(recoT_slp~pptmean, weights = 1/neeT_SE,data=c)
summary(MAP_reco)
tab_model(MAP_reco)

MAP_gpp<-lm(gppT_slp~pptmean, weights = 1/neeT_SE,data=c)
summary(MAP_gpp)
tab_model(MAP_gpp)

####FIGURE S2
fig9a<-ggplot(c, aes(x=aetmean, y=neeT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/neeT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                            "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  geom_smooth(aes(weight = 1/neeT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "")
  #labs(x ="Mean summer AET (mm)", y = "Temperature effect on NEE") 
fig9a<-fig9a +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
)
fig9a
ggsave("Final Figs/R/9a.jpg", fig9a, width = 60, height = 60, dpi=300, units="mm")

fig9b<-ggplot(c, aes(x=pptmean, y=neeT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/neeT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                            "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  geom_smooth(aes(weight = 1/neeT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Mean summer precipitation (mm)", y = "Temperature effect on NEE") 
fig9b<-fig9b +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
)
fig9b
ggsave("Final Figs/R/9b.jpg", fig9b, width = 60, height = 60, dpi=300, units="mm")


fig9c<-ggplot(c, aes(x=aetmean, y=gppT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/gppT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                            "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  #geom_smooth(aes(weight = 1/gppT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Mean summer AET (mm)", y = "Temperature effect on GPP") 
fig9c<-fig9c +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
)
fig9c
ggsave("Final Figs/R/9c.jpg", fig9c, width = 60, height = 60, dpi=300, units="mm")

fig9d<-ggplot(c, aes(x=pptmean, y=gppT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/gppT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                            "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  #geom_smooth(aes(weight = 1/neeT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Mean summer precipitation (mm)", y = "Temperature effect on GPP") 
fig9d<-fig9d +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
)
fig9d
ggsave("Final Figs/R/9d.jpg", fig9d, width = 60, height = 60, dpi=300, units="mm")

fig9e<-ggplot(c, aes(x=aetmean, y=recoT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/recoT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                             "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  #geom_smooth(aes(weight = 1/recoT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Mean summer AET (mm)", y = "Temperature effect on Reco") 
fig9e<-fig9e +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
)
fig9e
ggsave("Final Figs/R/9e.jpg", fig9e, width = 60, height = 60, dpi=300, units="mm")

fig9f<-ggplot(c, aes(x=pptmean, y=recoT_slp)) + 
  geom_point(aes(fill = permafrost, size=1/recoT_SE), colour="black", shape=21)+scale_fill_manual(values = c("No"="orange",
                                                                                                             "Yes"="turquoise3") )+
  scale_size(range = c(2,7)) +geom_hline(yintercept=0, linetype='dotted', col = 'black',size=0.75)+
  #geom_smooth(aes(weight = 1/neeT_SE), method = lm, size = 1,colour="black")+
  labs(x ="", y = "") 
  #labs(x ="Mean summer precipitation (mm)", y = "Temperature effect on Reco") 
fig9f<-fig9f +theme_bw()+theme(legend.background = element_rect(fill="transparent"),
                               legend.key = element_rect(colour="transparent"), 
                               legend.title = element_blank(), 
                               legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.text = element_text(size = 9),
                               axis.text.x=element_text(colour="black"),
                               axis.text.y=element_text(colour="black"),
)
fig9f
ggsave("Final Figs/R/9f.jpg", fig9f, width = 60, height = 60, dpi=300, units="mm")

### EXTENDED DATA FIGURE 2
fig9<-ggarrange(fig7a, fig7b,fig7c, fig7d, fig7e, fig7f, ncol = 3,nrow=2)
fig9
#ggsave("PaperFigs/Figure9.jpeg", fig7,width = 250, height = 180, dpi=300, units="mm")



