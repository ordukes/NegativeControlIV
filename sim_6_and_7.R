#Simulation code for Experiments 6-7 (Table 2) in supplementary material.

install.packages("ivreg")
install.packages("ivtools")
install.packages("geex")
install.packages("gmm")
install.packages("drgee")
install.packages("remotes")
install.packages("grattanCharts")

library(ivreg)
library(ivtools)
library(geex)
library(gmm)
library(drgee)
library(remotes)
library(grattanCharts)

rm(list=ls())

#setwd(" ") 
source("nciv_fun.R")

#Run simulations 6 and 7.

nem.sim <- function(n,nsim,seed,experiment){
  
  res<-matrix(NA,nsim,16)
  
  for(i in 1:nsim){
    
    set.seed(i*seed)
    U<-rbinom(n,1,0.5)
    C1<-rbinom(n,1,0.5)
    C2<-rnorm(n)
    S<-rbinom(n,1,plogis(-0.5+C1+0.6*C2+0.5*C1*C2))
    Z<-rbinom(n,1,plogis(0.25*C1-0.25*C2+0.5*C1*C2))

    if (experiment==6){ 
      A<-rbinom(n,1,S*plogis(1-1.5*Z-0.75*C1-0.3*C2-0.5*C1*C2+1*U)) 
      Y<-rnorm(n,1+U+A*S+0.5*Z*S*(1-0.4*C1-0.4*C2+0.5*C1*C2)+S*(0.5*C1+0.5*C2+0.5*C1*C2)+(0.5*C1+0.5*C2-0.5*C1*C2))
    }
    if (experiment==7){ 
      A<-rbinom(n,1,S*plogis(1-0.25*Z-0.75*C1-0.3*C2-0.5*C1*C2+1*U)) 
      Y<-rnorm(n,1+U+A*S+Z*(1-0.4*C1-0.4*C2+0.5*C1*C2)+S*(0.5*C1+0.5*C2+0.5*C1*C2)+(0.5*C1+0.5*C2-0.5*C1*C2))
    }
    
    simdata=data.frame(Y,A,Z,C1,C2,S)
    
    iv_f=Z~C1+C2+I(C1*C2);pop_f=S~C1+C2+I(C1*C2);odds_f=~1;tport_f=~C1+C2+I(C1*C2);om_f=Y~C1+C2+I(C1*C2)
    se_f=~C1+C2+I(C1*C2);ps_f=A~Z+C1+C2+I(C1*C2)
    
    subsim<-simdata[S==1,]
    dr_fit<-drgee(oformula = formula(Y~Z+C1+C2+I(C1*C2)),
                  eformula = formula(A~Z+C1+C2+I(C1*C2)),iaformula=formula(~1),
                  olink = "identity", elink = "logit",
                  data = subsim, estimation.method = "dr")
    res[i,1]<-dr_fit$coefficients
    res[i,2]<-dr_fit$vcov
    
    nciv_tsls_fit<-nem_tsls(simdata,tport_formula=tport_f,outcome_formula=om_f,seff_formula=se_f,teff_formula=~1)
    res[i,3]<-coef(nciv_tsls_fit)[length(coef(nciv_tsls_fit))]
    res[i,4]<-vcov(nciv_tsls_fit)[nrow(vcov(nciv_tsls_fit)),ncol(vcov(nciv_tsls_fit))]
    
    nciv_z_fit<-nem_g_Z(simdata,iv_formula=iv_f,odds_formula = odds_f,tport_formula=tport_f,teff_formula=~1)
    res[i,5]<-coef(nciv_z_fit)[length(coef(nciv_z_fit))]
    res[i,6]<-vcov(nciv_z_fit)[nrow(vcov(nciv_z_fit)),ncol(vcov(nciv_z_fit))]
    #
    nciv_s_fit<-nem_g_S(simdata,pop_formula=pop_f,odds_formula = odds_f,seff_formula=se_f,teff_formula=~1)
    res[i,7]<-coef(nciv_s_fit)[length(coef(nciv_s_fit))]
    res[i,8]<-vcov(nciv_s_fit)[nrow(vcov(nciv_s_fit)),ncol(vcov(nciv_s_fit))]
    #
    nciv_ipw_fit<-nem_ipw(simdata,iv_formula = iv_f,odds_formula = odds_f,pop_formula=pop_f,teff_formula=~1)
    res[i,9]<-coef(nciv_ipw_fit)[length(coef(nciv_ipw_fit))]
    res[i,10]<-vcov(nciv_ipw_fit)[nrow(vcov(nciv_ipw_fit)),ncol(vcov(nciv_ipw_fit))]
    
    nciv_mr_fit<-nem_mr(simdata,iv_formula=iv_f,pop_formula=pop_f,tport_formula=tport_f,
                        outcome_formula=om_f,teff_formula=~1,odds_formula=~1,seff_formula=se_f,
                        propensity_formula=ps_f,speff=FALSE)
    res[i,11]<-coef(nciv_mr_fit)[length(coef(nciv_mr_fit))]
    res[i,12]<-vcov(nciv_mr_fit)[nrow(vcov(nciv_mr_fit)),ncol(vcov(nciv_mr_fit))]
    
    nciv_eff_fit<-nem_mr(simdata,iv_formula=iv_f,pop_formula=pop_f,tport_formula=tport_f,
                         outcome_formula=om_f,teff_formula=~1,odds_formula=~1,seff_formula=se_f,
                         propensity_formula=ps_f,speff=TRUE)
    res[i,13]<-coef(nciv_eff_fit)[length(coef(nciv_mr_fit))]
    res[i,14]<-vcov(nciv_eff_fit)[nrow(vcov(nciv_mr_fit)),ncol(vcov(nciv_mr_fit))]
    
    fitA.ZC<-glm(formula=A~Z+C1+C2+I(C1*C2), data=subsim)
    fitY.AC<-glm(formula=Y~A+C1+C2+I(C1*C2), data=subsim)
    iv_fit<- ivglm(estmethod="ts", fitX.LZ=fitA.ZC, fitY.LX=fitY.AC, data=subsim,
                   ctrl=TRUE)
    res[i,15]<-iv_fit$est[2]
    res[i,16]<-iv_fit$vcov[2,2]
    
    cat(i,round(c(mean(res[1:i,1]),mean(res[1:i,3]),
                  mean(res[1:i,5]),mean(res[1:i,7]),
                  mean(res[1:i,9]),mean(res[1:i,11]),
                  mean(res[1:i,13]),mean(res[1:i,15])),digits=3),"\n")
  }
  return(res)
}

res6<-nem.sim(n=5000,nsim=2000,seed=320,experiment=6)
res7<-nem.sim(n=5000,nsim=2000,seed=44,experiment=7) 

truth=1
mean(res6[,1]-truth)
mean(res7[,1]-truth)
mean(res6[,1]-1.96*sqrt(res6[,2])<truth & res6[,1]+1.96*sqrt(res6[,2])>truth)
mean(res7[,1]-1.96*sqrt(res7[,2])<truth & res7[,1]+1.96*sqrt(res7[,2])>truth)

mean(res6[,15]-truth)
mean(res7[,15]-truth)
mean(res6[,15]-1.96*sqrt(res6[,2])<truth & res6[,15]+1.96*sqrt(res6[,2])>truth)
mean(res7[,15]-1.96*sqrt(res7[,2])<truth & res7[,15]+1.96*sqrt(res7[,2])>truth)

#Create Table 2.

seq<-c(3,5,7,9,11,13)
n_est<-6;n_exp<-2;n_stat<-3
result_mat<-matrix(NA,n_exp,n_est*n_stat)
for(z in seq){
  q<-3*((z-1)/2)-2
  
  result_mat[1,q]<-mean(res6[,z]-truth)
  result_mat[1,q+1]<-sqrt(var(res6[,z]))
  result_mat[1,q+2]<-mean(res6[,z]-1.96*sqrt(res6[,z+1])<truth & res6[,z]+1.96*sqrt(res6[,z+1])>truth)
  
  result_mat[2,q]<-mean(res7[,z]-truth)
  result_mat[2,q+1]<-sqrt(var(res7[,z]))
  result_mat[2,q+2]<-mean(res7[,z]-1.96*sqrt(res7[,z+1])<truth & res7[,z]+1.96*sqrt(res7[,z+1])>truth)
  }
result_mat_1<-round(result_mat,digits=2)[,-13:-15] #Exclude inefficient MR
categ<-rep(c("Bias","ESE","Cov"),times=n_est)
experiment<-c(1:5)
result_mat_N<-as.data.frame(cbind(experiment,result_mat_1))
colnames(result_mat_N)<-c("Exp","$\\hat{\\psi}_{TSLS}$__Bias","$\\hat{\\psi}_{TSLS}$__SD","$\\hat{\\psi}_{TSLS}$__Cov",
                          "$\\hat{\\psi}_{g-Z}$__Bias","$\\hat{\\psi}_{g-Z}$__SD","$\\hat{\\psi}_{g-Z}$__Cov",
                          "$\\hat{\\psi}_{g-S}$__Bias","$\\hat{\\psi}_{g-S}$__SD","$\\hat{\\psi}_{g-S}$__Cov",
                          "$\\hat{\\psi}_{g-IPW}$__Bias","$\\hat{\\psi}_{g-IPW}$__SD","$\\hat{\\psi}_{g-IPW}$__Cov",
                          "$\\hat{\\psi}_{MR-eff}$__Bias","$\\hat{\\psi}_{MR-eff}$__SD","$\\hat{\\psi}_{MR-eff}$__Cov")
print_2heading_xtable(.data=result_mat_N,booktabs=FALSE, separator = "__",xtable.dots = list(label = "Kang2"),table.placement="htbp",sanitize.text.function=identity)

