library("ivmodel")
library("expm")
library("tictoc")
library("xtable")
library("glmnet")
library("sisVIVE")
setwd("C:/Users/Hippo/TSE/Instruments Selection/sim_2022")
source("functions_v5.R")




##### Large sample size

n<-400
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1/sqrt(n),1/sqrt(n),1/n,1/n)
true_z<-c(0,1,0,1,0,1)
alpha_true<-c(0,0,0)

kc<-0.7
B<-20
coverage<-0.95

test_z<-mat_z(kz)

diag1<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag1[2,])
max(diag1[2,])
min(diag1[4,])
max(diag1[4,])
mean(diag1[3,],digits=3)


ns<-2000
set1_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=F,K=20,normalize=F,package=F)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=F,K=20,normalize=F)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,1,0,0,0,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set1_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set1_diag_2sls)



n<-4000
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1/sqrt(n),1/sqrt(n),1/n,1/n)
true_z<-c(0,1,0,1,0,1)
alpha_true<-c(0,0,0)

kc<-0.7
B<-40
coverage<-0.95

test_z<-mat_z(kz)
diag2<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag2[2,])
max(diag2[2,])
min(diag2[4,])
max(diag2[4,])
mean(diag2[3,],digits=3)


ns<-2000
set2_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=F,K=20,normalize=F,package=F)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=F,K=20,normalize=F)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,1,0,0,0,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set2_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set2_diag_2sls)



##### Strong setting

n<-400
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1,1,1,1)
true_z<-c(0,1,0,1,0,1)
alpha_true<-c(1,1,1)

kc<-0.7
B<-20
coverage<-0.95

test_z<-mat_z(kz)
diag3<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag3[2,])
max(diag3[2,])
min(diag3[4,])
max(diag3[4,])
mean(diag3[3,],digits=3)


ns<-2000
set3_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=F,K=20,normalize=F,package=F)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=F,K=20,normalize=F)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,0,1,0,1,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set3_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set3_diag_2sls)



n<-4000
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1,1,1,1)
true_z<-c(0,1,0,1,0,1)
alpha_true<-c(1,1,1)

kc<-0.7
B<-40
coverage<-0.95

test_z<-mat_z(kz)
diag4<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag4[2,])
max(diag4[2,])
min(diag4[4,])
max(diag4[4,])
mean(diag4[3,],digits=3)


ns<-2000
set4_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=F,K=20,normalize=F,package=F)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=F,K=20,normalize=F)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,0,1,0,1,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set4_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set4_diag_2sls)




##### Strong favorable setting


n<-400
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1,1,1,1)
true_z<-c(0,0,0,1,0,1)
alpha_true<-c(1,3)

kc<-0.7
B<-20
coverage<-0.95


test_z<-mat_z(kz)
diag5<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag5[2,])
max(diag5[2,])
min(diag5[4,])
max(diag5[4,])
mean(diag5[3,],digits=3)


ns<-2000
set5_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=1:100,K=20,normalize=T,package=T)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=1:100,K=20,normalize=T)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,1,1,0,1,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set5_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set5_diag_2sls)


n<-4000
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1,1,1,1)
true_z<-c(0,0,0,1,0,1)
alpha_true<-c(1,3)

kc<-0.7
B<-40
coverage<-0.95


test_z<-mat_z(kz)
diag6<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag6[2,])
max(diag6[2,])
min(diag6[4,])
max(diag6[4,])
mean(diag6[3,],digits=3)


ns<-2000
set6_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=1:100,K=20,normalize=T,package=T)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=1:100,K=20,normalize=T)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,1,1,0,1,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set6_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set6_diag_2sls)



##### General setting

n<-400
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1/sqrt(n),1/sqrt(n),1/n,1/n)
true_z<-c(0,1,0,1,0,1)
alpha_true<-c(1,1,1)

kc<-0.7
B<-20
coverage<-0.95

test_z<-mat_z(kz)
diag7<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag7[2,])
max(diag7[2,])
min(diag7[4,])
max(diag7[4,])
mean(diag7[3,],digits=3)


ns<-2000
set7_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=F,K=20,normalize=F,package=F)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=F,K=20,normalize=F)
  
  risk_3<-risk_est(y,x,z,kc,B,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,0,0,0,0,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set7_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set7_diag_2sls)




n<-4000
kz<-6
sig_u<-1
sig_v<-1
rho<-0.5
bet<-2
omega<-matrix(c(sig_u,rho,rho,sig_v),2,2)
sq_omega<-expm::sqrtm(omega)
sig_z<-diag(rep(0.9,kz))+matrix(0.1,kz,kz)
sq_sig_z<-expm::sqrtm(sig_z)
c<-0.5
pi<-c*c(1,1,1/sqrt(n),1/sqrt(n),1/n,1/n)
true_z<-c(0,1,0,1,0,1)
alpha_true<-c(1,1,1)

kc<-0.7
B<-40
coverage<-0.95

test_z<-mat_z(kz)
diag8<-diag_theor(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true)
min(diag8[2,])
max(diag8[2,])
min(diag8[4,])
max(diag8[4,])
mean(diag8[3,],digits=3)


ns<-2000
set8_diag_2sls<-matrix(0,8*9,ns)


for (i in 1:ns){
  set.seed(i)
  
  SB_cv<-as.numeric(sample_boot(n,kc,B,type="cv"))
  uv<-matrix(rnorm(n*2),n,2)%*%sq_omega
  u<-uv[,1]
  v<-uv[,2]
  z<-matrix(rnorm(n*kz),n,kz)%*%sq_sig_z
  x<-z%*%pi+v
  y<-x*bet+z[,true_z==T]%*%alpha_true+u
  
  riskdn<-risk_dn(y,x,z,test_z,type_est="2sls",type_sig="mallows",inc=F,oracle=T)
  min_dn<-test_z[,which.min(riskdn)]
  
  riskan<-risk_ad(y,x,z,test_z,type_est="2sls")
  min_an<-test_z[,which.min(riskan)]
  
  lasso_iv<-risk_lasso(y=y,x=x,z=z,lambda=F,K=20,normalize=F,package=F)
  adalasso_iv<-risk_adalasso(y,x,z,lambda=F,K=20,normalize=F)
  
  risk_3<-risk_est(y,x,z,kc,SB_cv,"2sls",test_z,inc=T,cv=1)
  min_3<-test_z[,apply(risk_3,2, function(e) which.min(e))]
  
  oracle<-c(1,0,0,0,0,0)
  
  diag_risks<-cbind(t(diag_est(y,x,z,bet,coverage,
                               cbind(min_dn,min_an,lasso_iv,adalasso_iv,min_3,
                                     oracle),
                               type_est="2sls",type_ci="AR",
                               inc=T)))
  set8_diag_2sls[,i]<-as.numeric(diag_risks)
  
  print(i)
  
}

diag_estf(set8_diag_2sls)



