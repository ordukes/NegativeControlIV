#Functions for implementing estimators and obtaining standard errors.

nem_ipw_fun <- function(theta,data){
  Y=cbind(data$Y);A=cbind(data$A);S=cbind(data$S);Z=cbind(data$Z)
  cov_iv=data.matrix(cbind(data$cov_iv));cov_pop=data.matrix(cbind(data$cov_pop))
  cov_odds=data.matrix(cbind(data$cov_odds));cov_te=data.matrix(cbind(data$cov_te))
  cov_iv_S=data.matrix(cbind(cov_odds*as.vector(S),cov_iv))
  cov_pop_Z=data.matrix(cbind(cov_odds*as.vector(Z),cov_pop))
  
  iv_pos<- 1:ncol(cov_iv)
  init_odds_pos_1<-(max(iv_pos) + 1):(max(iv_pos) + (ncol(cov_odds)))
  pop_pos <- (max(init_odds_pos_1) + 1):(max(init_odds_pos_1) + (ncol(cov_pop)))
  init_odds_pos_2<-(max(pop_pos) + 1):(max(pop_pos) + (ncol(cov_odds)))
  odds_pos<-(max(init_odds_pos_2) + 1):(max(init_odds_pos_2) + (ncol(cov_odds)))
  te_pos <- (max(odds_pos) + 1):(max(odds_pos) + (ncol(cov_te)))
  
  iv_pred <- plogis(cov_iv%*% theta[iv_pos]+cov_odds%*%theta[init_odds_pos_1]*S)
  pop_pred <- plogis(cov_pop%*%theta[pop_pos]+cov_odds%*%theta[init_odds_pos_2]*Z)
  pop_pred_0 <- plogis(cov_pop%*%theta[pop_pos]);iv_pred_0 <- plogis(cov_iv %*% theta[iv_pos])
  iv_pred_1 <- plogis(cov_odds%*%theta[odds_pos]+cov_iv%*%theta[iv_pos])
  delta<-(1+((1-pop_pred_0)*iv_pred_0)/(pop_pred_0*iv_pred_1))^(-1)
  C_l<-exp(cov_odds%*%theta[odds_pos])*iv_pred_0*pop_pred_0 +
    (1-iv_pred_0)*pop_pred_0+iv_pred_0*(1-pop_pred_0)+(1-iv_pred_0)*(1-pop_pred_0)
  phi<-(Z*S)/((exp(cov_odds%*%theta[odds_pos])*iv_pred_0*pop_pred_0)/C_l)+
    ((1-Z)*(1-S))/(((1-iv_pred_0)*(1-pop_pred_0))/C_l)-
    ((1-Z)*S)/(((1-iv_pred_0)*pop_pred_0)/C_l)-(Z*(1-S))/((iv_pred_0*(1-pop_pred_0))/C_l)

  m1<-cov_iv_S*as.vector((Z-iv_pred))
  m2<-cov_pop_Z*as.vector((S-pop_pred))
  m3<-cov_odds*as.vector((S-delta)*(Z-plogis(cov_odds%*%theta[odds_pos]*S+cov_iv%*%theta[iv_pos])))
  m4<-cov_te*as.vector(phi*(Y-cov_te%*%theta[te_pos]*A*S))
  return(cbind(m1,m2,m3,m4))
}

nem_ipw <- function(data,iv_formula,pop_formula,odds_formula,teff_formula){
  S=cbind(data$S);Z=cbind(data$Z)
  cov_iv<-model.matrix(iv_formula,data=data);cov_pop<-model.matrix(pop_formula,data=data)
  cov_odds<-model.matrix(odds_formula,data=data);cov_te<-model.matrix(teff_formula,data=data)
  nparms <- ncol(cov_iv) + ncol(cov_pop) + 3*ncol(cov_odds) + ncol(cov_te) 
  iv_mod<-glm(Z~-1+cov_iv+as.matrix(cov_odds*S),family=binomial,data=data)
  pop_mod<-glm(S~-1+cov_pop+as.matrix(cov_odds*Z),family=binomial,data=data)
  odds_fit<-drgee(oformula = iv_formula,eformula = pop_formula,iaformula = odds_formula,
                  olink = "logit", elink = "logit",data = data, estimation.method = "dr")
  init <- c(coef(iv_mod),coef(pop_mod),coef(odds_fit),rep(0,ncol(cov_te)))
  gmm_data<-list(Y=as.numeric(data$Y),A=as.numeric(data$A),S=as.numeric(data$S),Z=as.numeric(data$Z),cov_iv=cov_iv,cov_pop=cov_pop,cov_odds=cov_odds,cov_te=cov_te)
  gmm(nem_ipw_fun, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-ncol(cov_te))),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb")
}

nem_g_fun_Z <- function(theta,data){
  Y=cbind(data$Y);A=cbind(data$A);S=cbind(data$S);Z=cbind(data$Z)
  cov_iv=data.matrix(cbind(data$cov_iv));cov_tp=data.matrix(cbind(data$cov_tp))
  cov_odds=data.matrix(cbind(data$cov_odds));cov_te=data.matrix(cbind(data$cov_te))
  cov_iv_S=data.matrix(cbind(cov_odds*as.vector(S),cov_iv))
  
  iv_pos<- 1:ncol(cov_iv)
  odds_pos<- (max(iv_pos) + 1):(max(iv_pos) + (ncol(cov_odds)))
  tp_pos<-(max(odds_pos) + 1):(max(odds_pos) + (ncol(cov_tp)))
  te_pos <- (max(tp_pos) + 1):(max(tp_pos) + (ncol(cov_te)))
  
  iv_pred <- plogis(cov_iv%*% theta[iv_pos]+cov_odds%*%theta[odds_pos]*S)
  tp_pred <- (cov_tp %*% theta[tp_pos])*Z
  m1<-cov_iv_S*as.vector(Z-iv_pred)
  m2<-cov_tp*as.vector((1-S)*(Z-iv_pred)*(Y-tp_pred))
  m3<-cov_te*as.vector(S*(Z-iv_pred)*(Y-tp_pred-cov_te%*%theta[te_pos]*A))
  return(cbind(m1,m2,m3))
}

nem_g_Z <- function(data,iv_formula,odds_formula,tport_formula,teff_formula){
  S=as.numeric(data$S);Z=as.numeric(data$Z);Y=as.numeric(data$Y)
  cov_iv<-model.matrix(iv_formula,data=data); cov_odds<-model.matrix(odds_formula,data=data);
  cov_tp<-model.matrix(tport_formula,data=data);cov_te<-model.matrix(teff_formula,data=data)
  nparms <- ncol(cov_iv) + ncol(cov_odds) + ncol(cov_tp) + ncol(cov_te)
  
  iv_mod<-glm(Z~-1+cov_iv+data.matrix(cov_odds*S),family=binomial,data=data)
  iv_pred<-predict(iv_mod,type="response")
  g_est_Z<-solve(t(cov_tp*(Z-iv_pred)*(1-S))%*%(cov_tp*Z))%*%(t(cov_tp*(Z-iv_pred)*(1-S))%*%Y)
  init <- c(coef(iv_mod),g_est_Z,rep(0,ncol(cov_te)))

  gmm_data<-list(Y=as.numeric(data$Y),A=as.numeric(data$A),S=as.numeric(data$S),Z=as.numeric(data$Z),cov_iv=cov_iv,cov_odds=cov_odds,cov_tp=cov_tp,cov_te=cov_te)
  gmm(nem_g_fun_Z, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-ncol(cov_te))),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb")
}

nem_g_fun_S <- function(theta,data){
  Y=cbind(data$Y);A=cbind(data$A);S=cbind(data$S);Z=cbind(data$Z)
  cov_pop=data.matrix(cbind(data$cov_pop));cov_se=data.matrix(cbind(data$cov_se))
  cov_odds=data.matrix(cbind(data$cov_odds));cov_te=data.matrix(cbind(data$cov_te))
  cov_pop_Z=data.matrix(cbind(cov_odds*as.vector(Z),cov_pop))

  pop_pos<- 1:ncol(cov_pop)
  odds_pos<- (max(pop_pos) + 1):(max(pop_pos) + (ncol(cov_odds)))
  se_pos<-(max(odds_pos) + 1):(max(odds_pos) + (ncol(cov_se)))
  te_pos <- (max(se_pos) + 1):(max(se_pos) + (ncol(cov_te)))
  
  pop_pred <- plogis(cov_pop%*% theta[pop_pos]+cov_odds%*%theta[odds_pos]*Z)
  se_pred <- (cov_se %*% theta[se_pos])*S
  m1<-cov_pop_Z*as.vector(S-pop_pred)
  m2<-cov_se*as.vector((S-pop_pred)*(Y-se_pred-cov_te%*%theta[te_pos]*A*S))
  m3<-cov_te*as.vector(Z*(S-pop_pred)*(Y-se_pred-cov_te%*%theta[te_pos]*A*S))
  return(cbind(m1,m2,m3))
}

nem_g_S <- function(data,pop_formula,odds_formula,seff_formula,teff_formula){
  S=as.numeric(data$S);Z=as.numeric(data$Z);Y=as.numeric(data$Y);A=as.numeric(data$A)
  cov_pop<-model.matrix(pop_formula,data=data); cov_odds<-model.matrix(odds_formula,data=data);
  cov_se<-model.matrix(seff_formula,data=data); cov_te<-model.matrix(teff_formula,data=data)
  nparms <- ncol(cov_pop) + ncol(cov_odds) + ncol(cov_se) + ncol(cov_te)
  
  pop_mod<-glm(S~-1+cov_pop+as.matrix(cov_odds*Z),family=binomial,data=data)
  pop_pred<-predict(pop_mod,type="response")
  cov_se_te<-data.matrix(cbind(cov_se,cov_te*Z)); cov_se_te_d<-data.matrix(cbind(cov_se*S,cov_te*A*S))
  g_est_S<-solve(t(cov_se_te*(S-pop_pred))%*%(cov_se_te_d))%*%(t(cov_se_te*(S-pop_pred))%*%Y)
  init <- c(coef(pop_mod),g_est_S)
  
  gmm_data<-list(Y=as.numeric(data$Y),A=as.numeric(data$A),S=as.numeric(data$S),Z=as.numeric(data$Z),cov_pop=cov_pop,cov_odds=cov_odds,cov_se=cov_se,cov_te=cov_te)
  gmm(nem_g_fun_S, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-ncol(cbind(cov_se,cov_te)))),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb")
}
 
nem_tsls_fun <- function(theta,data){
  Y=cbind(data$Y);A=cbind(data$A);S=cbind(data$S);Z=cbind(data$Z)
  cov_tp=data.matrix(cbind(data$cov_tp));cov_om=data.matrix(cbind(data$cov_om))
  cov_se=data.matrix(cbind(data$cov_se));cov_te=data.matrix(cbind(data$cov_te))

  om_pos<-1:ncol(cov_om)
  tp_pos <- (max(om_pos) + 1):(max(om_pos) + (ncol(cov_tp)))
  se_pos<-(max(tp_pos) + 1):(max(tp_pos) + (ncol(cov_se)))
  te_pos <- (max(se_pos) + 1):(max(se_pos) + (ncol(cov_te)))
  
  om <- cov_om%*%theta[om_pos];tp_pred <-(cov_tp %*% theta[tp_pos])*Z
  se_pred <- (cov_se %*% theta[se_pos])*S
  m1<-cov_om*as.vector((Y-tp_pred-om)*(1-S))
  m2<-data.matrix(cov_tp*as.vector(Z))*as.vector((Y-tp_pred-om)*(1-S))
  m3<-cov_se*as.vector(Y-tp_pred-om-se_pred-cov_te%*%theta[te_pos]*A*S)
  m4<-data.matrix(cov_te*as.vector(Z*S))*as.vector(Y-tp_pred-om-se_pred-cov_te%*%theta[te_pos]*A*S)
  return(cbind(m1,m2,m3,m4))
}

nem_tsls <- function(data,tport_formula,seff_formula,outcome_formula,teff_formula){
  cov_tp<-model.matrix(tport_formula,data=data);cov_om<-model.matrix(outcome_formula,data=data)
  cov_se<-model.matrix(seff_formula,data=data);cov_om<-model.matrix(outcome_formula,data=data)
  cov_te<-model.matrix(teff_formula,data=data);Z=as.numeric(data$Z)
  nparms <- ncol(cov_tp) + ncol(cov_om) + ncol(cov_se) + ncol(cov_te)
  
  lmod_Y0<-lm(Y~-1+cov_om+cbind(cov_tp*as.vector(Z)),data=data,subset=S==0)
  gamma<-(cov_tp%*%coef(lmod_Y0)[(ncol(cov_om)+1):(ncol(cov_om)+ncol(cov_tp))])*Z
  om<-cov_om%*%coef(lmod_Y0)[1:ncol(cov_om)]
  lmod_Y1<-ivreg(I(Y-gamma)~-1+cbind(cov_se*S)+cbind(cov_te*A)|-1+cbind(cov_se*S)+cbind(cov_te*Z),data=data,offset=om)
  init <- c(coef(lmod_Y0),coef(lmod_Y1))
  
  gmm_data<-list(Y=as.numeric(data$Y),A=as.numeric(data$A),S=as.numeric(data$S),Z=as.numeric(data$Z),cov_tp=cov_tp,cov_om=cov_om,cov_se=cov_se,cov_te=cov_te)
  gmm(nem_tsls_fun, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-ncol(cov_te))),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb")
}

nem_mr_fun <- function(theta,data){
  Y=cbind(data$Y);A=cbind(data$A);S=cbind(data$S);Z=cbind(data$Z);speff=data$speff
  cov_iv=data.matrix(cbind(data$cov_iv));cov_pop=data.matrix(cbind(data$cov_pop))
  cov_tp=data.matrix(cbind(data$cov_tp));cov_om=data.matrix(cbind(data$cov_om))
  cov_te=data.matrix(cbind(data$cov_te));cov_ps=data.matrix(cbind(data$cov_ps));
  cov_ps_0=data.matrix(cbind(data$cov_ps_0));cov_ps_1=data.matrix(cbind(data$cov_ps_1))
  cov_odds=data.matrix(cbind(data$cov_odds));cov_se=data.matrix(cbind(data$cov_se))
  
  ps_pos <- 1:ncol(cov_ps)
  iv_pos<- (max(ps_pos) + 1):(max(ps_pos) + (ncol(cov_iv)))
  init_odds_pos_1<-(max(iv_pos) + 1):(max(iv_pos) + (ncol(cov_odds)))
  pop_pos <- (max(init_odds_pos_1) + 1):(max(init_odds_pos_1) + (ncol(cov_pop)))
  init_odds_pos_2<-(max(pop_pos) + 1):(max(pop_pos) + (ncol(cov_odds)))
  odds_pos<-(max(init_odds_pos_2) + 1):(max(init_odds_pos_2) + (ncol(cov_odds)))
  om_pos<- (max(odds_pos) + 1):(max(odds_pos) + (ncol(cov_om)))
  init_tp_pos <- (max(om_pos) + 1):(max(om_pos) + (ncol(cov_tp)))
  tp_pos <- (max(init_tp_pos) + 1):(max(init_tp_pos) + (ncol(cov_tp)))
  se_pos <- (max(tp_pos) + 1):(max(tp_pos) + (ncol(cov_se)))
  init_te_pos <- (max(se_pos ) + 1):(max(se_pos) + (ncol(cov_te)))
  te_pos <- (max(init_te_pos) + 1):(max(init_te_pos) + (ncol(cov_te)))
  
  ps_pred <- plogis(cov_ps%*%theta[ps_pos]) 
  iv_pred_i <- plogis(cov_iv%*% theta[iv_pos]+cov_odds%*%theta[init_odds_pos_1]*S)
  pop_pred_i <- plogis(cov_pop%*%theta[pop_pos]+cov_odds%*%theta[init_odds_pos_2]*Z)
  pop_pred_0 <- plogis(cov_pop%*%theta[pop_pos]);iv_pred_0 <- plogis(cov_iv %*% theta[iv_pos])
  iv_pred_1 <- plogis(cov_odds%*%theta[odds_pos]+cov_iv%*%theta[iv_pos])
  delta<-(1+((1-pop_pred_0)*iv_pred_0)/(pop_pred_0*iv_pred_1))^(-1)
  C_l<-exp(cov_odds%*%theta[odds_pos])*iv_pred_0*pop_pred_0 +
    (1-iv_pred_0)*pop_pred_0+iv_pred_0*(1-pop_pred_0)+(1-iv_pred_0)*(1-pop_pred_0)
  phi<-(Z*S)/((exp(cov_odds%*%theta[odds_pos])*iv_pred_0*pop_pred_0)/C_l)+
    ((1-Z)*(1-S))/(((1-iv_pred_0)*(1-pop_pred_0))/C_l)-
    ((1-Z)*S)/(((1-iv_pred_0)*pop_pred_0)/C_l)-(Z*(1-S))/((iv_pred_0*(1-pop_pred_0))/C_l)
  
  C_Z<-exp(cov_odds%*%theta[odds_pos]*S)*iv_pred_0+(1-iv_pred_0)
  iv_pred<-(exp(cov_odds%*%theta[odds_pos]*S)*iv_pred_0)/C_Z
  C_S<-exp(cov_odds%*%theta[odds_pos]*Z)*pop_pred_0+(1-pop_pred_0)
  pop_pred<-(exp(cov_odds%*%theta[odds_pos]*Z)*pop_pred_0)/C_S
  om <- cov_om%*%theta[om_pos]
  tp_pred_i <- (cov_tp %*% theta[init_tp_pos])*Z
  tp_pred <- (cov_tp %*% theta[tp_pos])*Z
  se_pred <- (cov_se %*% theta[se_pos])*S
  
  if(speff==FALSE){ weight<-rep(1,dim(cov_te)[1]) } else {
    weight_num<-plogis(cov_ps_1%*%theta[ps_pos])-plogis(cov_ps_0%*%theta[ps_pos])
    weight_denom<-1/((exp(cov_odds%*%theta[odds_pos])*iv_pred_0*pop_pred_0)/C_l)+
      1/(((1-iv_pred_0)*(1-pop_pred_0))/C_l)+1/(((1-iv_pred_0)*pop_pred_0)/C_l)+
      1/((iv_pred_0*(1-pop_pred_0))/C_l)
    weight<-weight_num/ weight_denom
  }
  
  m1<-cov_ps*as.vector((A-ps_pred)*S)
  m2<-cov_iv*as.vector((Z-iv_pred_i))
  m3<-data.matrix(cov_odds*as.vector(S))*as.vector((Z-iv_pred_i))
  m4<-cov_pop*as.vector((S-pop_pred_i))
  m5<-data.matrix(cov_odds*as.vector(Z))*as.vector((S-pop_pred_i))
  m6<-cov_odds*as.vector((S-delta)*(Z-plogis(cov_odds%*%theta[odds_pos]*S+cov_iv%*%theta[iv_pos])))
  m7<-cov_om*as.vector((Y-om-tp_pred_i)*(1-S))
  m8<-data.matrix(cov_tp*as.vector(Z))*as.vector((Y-om-tp_pred_i)*(1-S))
  m9<-cov_tp*as.vector((1-S)*(Z-iv_pred)*(Y-tp_pred-om))
  m10<-cov_se*as.vector((S-pop_pred)*(Y-tp_pred-om-se_pred-cov_te%*%theta[init_te_pos]*A*S))
  m11<-cov_te*as.vector(Z*(S-pop_pred)*(Y-tp_pred-om-se_pred-cov_te%*%theta[init_te_pos]*A*S))
  m12<-cov_te*as.vector(weight*phi*(Y-tp_pred-om-se_pred-cov_te%*%theta[te_pos]*A*S))

  return(cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12))
}

nem_mr <- function(data,iv_formula,pop_formula,tport_formula,outcome_formula,teff_formula,propensity_formula,odds_formula,seff_formula,speff){
  cov_iv<-model.matrix(iv_formula,data=data);cov_pop<-model.matrix(pop_formula,data=data)
  cov_tp<-model.matrix(tport_formula,data=data);cov_om<-model.matrix(outcome_formula,data=data)
  cov_se<-model.matrix(seff_formula,data=data);cov_te<-model.matrix(teff_formula,data=data);
  cov_odds<-model.matrix(odds_formula,data=data)
  cov_ps<-model.matrix(propensity_formula,data=data);Z=as.numeric(data$Z);S=as.numeric(data$S);Y=as.numeric(data$Y);A=as.numeric(data$A)
  data_Z_0<-data; data_Z_0[, colnames(data_Z_0)== "Z"]=0;data_Z_1<-data; data_Z_1[, colnames(data_Z_1)== "Z"]=1
  cov_ps_0<-model.matrix(propensity_formula,data=data_Z_0);cov_ps_1<-model.matrix(propensity_formula,data=data_Z_1)
  
  ps_mod<-glm(propensity_formula,family=binomial,data=data,subset=S==1)
  ps<-rep(NA,length(A));ps[S==1]<-predict(ps_mod,type="response")
  iv_mod<-glm(Z~-1+cov_iv+as.matrix(cov_odds*S),family=binomial,data=data)
  pop_mod<-glm(S~-1+cov_pop+as.matrix(cov_odds*Z),family=binomial,data=data)
  odds_fit<-drgee(oformula = iv_formula,eformula = pop_formula,iaformula = odds_formula,
                  olink = "logit", elink = "logit",data = data, estimation.method = "dr")
  iv_pred_0<-plogis(cov_iv%*%coef(iv_mod)[1:ncol(cov_iv)])
  C_Z<-exp(cov_odds%*%coef(odds_fit)*S)*iv_pred_0+(1-iv_pred_0)
  iv_pred<-as.vector((exp(cov_odds%*%coef(odds_fit)*S)*iv_pred_0)/C_Z)
  pop_pred_0<-plogis(cov_pop%*%coef(pop_mod)[1:ncol(cov_pop)])
  C_S<-exp(cov_odds%*%coef(odds_fit)*Z)*pop_pred_0+(1-pop_pred_0)
  pop_pred<-as.vector((exp(cov_odds%*%coef(odds_fit)*Z)*pop_pred_0)/C_S)

  lmod_Y0<-lm(Y~-1+cov_om+cbind(cov_tp*as.vector(Z)),data=data,subset=S==0)
  om<-cov_om%*%coef(lmod_Y0)[1:ncol(cov_om)]; gamma<-(cov_tp%*%coef(lmod_Y0)[(ncol(cov_om)+1):(ncol(cov_om)+ncol(cov_tp))])*Z
  g_est_Z<-solve(t(cov_tp*(Z-iv_pred)*(1-S))%*%(cov_tp*Z))%*%(t(cov_tp*(Z-iv_pred)*(1-S))%*%(Y-om))
  cov_se_te<-data.matrix(cbind(cov_se,cov_te*Z)); cov_se_te_d<-data.matrix(cbind(cov_se*S,cov_te*A*S))
  g_est_S_A<-solve(t(cov_se_te*(S-pop_pred))%*%(cov_se_te_d))%*%(t(cov_se_te*(S-pop_pred))%*%(Y-gamma-om))
  
  #NOTE: g-est used as initial value.
  init <- c(coef(ps_mod),coef(iv_mod),coef(pop_mod),coef(odds_fit),coef(lmod_Y0),g_est_Z,g_est_S_A,g_est_S_A[ncol(cov_se_te)-ncol(cov_te)+1])
  nparms <- length(init)

  gmm_data<-list(Y=as.numeric(data$Y),A=as.numeric(data$A),S=as.numeric(data$S),Z=as.numeric(data$Z),
                 cov_iv=cov_iv,cov_pop=cov_pop,cov_tp=cov_tp,cov_om=cov_om,cov_te=cov_te,cov_odds=cov_odds,
                 cov_ps=cov_ps,cov_se=cov_se,cov_ps_0=cov_ps_0,cov_ps_1=cov_ps_1,speff=speff)
  gmm(nem_mr_fun, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-ncol(cov_te))),eqConstFullVcov=TRUE,
       crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=10000)
}

