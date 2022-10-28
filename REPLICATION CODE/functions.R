
##### Estimators

### IV estimator

est<-function(y,x,z, type="2sls"){
  z<-as.matrix(z)
  n<-dim(z)[1]
  s<-dim(z)[2]
  zx<-t(z)%*%x
  zy<-t(z)%*%y
  zz<-t(z)%*%z
  if (type=="2sls"){
    if (s>1){
      num<-as.numeric(t(zx)%*%solve(zz)%*%zy)
      den<-as.numeric(t(zx)%*%solve(zz)%*%zx)
      b_est<-num/den
    } else if (s==1){
      b_est<-zy/zx
    }
    
  } else if (type=="jive"){
    if (s>1){
      num1<-as.numeric(t(zx)%*%solve(zz)%*%zy)
      den1<-as.numeric(t(zx)%*%solve(zz)%*%zx)
      d<-apply(z,1,function(e) t(e)%*%solve(zz)%*%e)
      num2<-as.numeric(t(x*d)%*%y)
      den2<-as.numeric(t(x*d)%*%x)
      
      b_est<-(num1-num2)/(den1-den2)
    } else if (s==1){
      d<-apply(z,1,function(e) e^2/zz)
      b_est<-(zx*zy/zz-t(x*d)%*%y)/as.numeric(zx^2/zz-t(x*d)%*%x)
    }
  } else if (type=="b2sls"){
    if (s>1){
      num<-as.numeric(t(zx)%*%solve(zz)%*%zy)-as.numeric(t(x)%*%y)*(1-(s-2)/n)*(s-2)/n
      den<-as.numeric(t(zx)%*%solve(zz)%*%zx)-as.numeric(t(x)%*%x)*(1-(s-2)/n)*(s-2)/n
      b_est<-num/den
    } else if (s==1){
      num<-as.numeric(t(zx)%*%solve(zz)%*%zy)-as.numeric(t(x)%*%y)*(1-(s-2)/n)*(s-2)/n
      den<-as.numeric(t(zx)%*%solve(zz)%*%zx)-as.numeric(t(x)%*%x)*(1-(s-2)/n)*(s-2)/n
      b_est<-num/den
    }
  }
  return(b_est)
}

### Estimator of asymptotic variance of IV estimator

b_sig<-function(y,x,z,type="2sls"){
  z<-as.matrix(z)
  n<-dim(z)[1]
  s<-dim(z)[2]
  zx<-t(z)%*%x
  zy<-t(z)%*%y
  zz<-t(z)%*%z
  if (type=="2sls"){
    if (s>1){
      num<-as.numeric(t(zx)%*%solve(zz)%*%zy)
      den<-as.numeric(t(zx)%*%solve(zz)%*%zx)
      b_est<-num/den
      err<-sum((y-x*as.numeric(b_est))^2)/(n*den)
    } else if (s==1){
      b_est<-zy/zx
      err<-sum((y-x*as.numeric(b_est))^2)/(n*(zx^2/zz))
    } 
  } else if (type=="jive"){
    if (s>1){
      num1<-as.numeric(t(zx)%*%solve(zz)%*%zy)
      den1<-as.numeric(t(zx)%*%solve(zz)%*%zx)
      d<-apply(z,1,function(e) t(e)%*%solve(zz)%*%e)
      num2<-as.numeric(t(x*d)%*%y)
      den2<-as.numeric(t(x*d)%*%x)
      
      b_est<-(num1-num2)/(den1-den2)
      err<-sum((y-x*b_est)^2)/(n*abs(den1-den2))
    } else if (s==1){
      d<-sapply(z,function(e) e^2/zz)
      b_est<-(zx*zy/zz-as.numeric(t(x*d)%*%y))/as.numeric(zx^2/zz-as.numeric(t(x*d)%*%x))
      err<-sum((y-x*as.numeric(b_est))^2)/(n*abs(zx^2/zz-as.numeric(t(x*d)%*%x)))
    }
  } else if(type=="b2sls"){
    if (s>1){
      num<-as.numeric(t(zx)%*%solve(zz)%*%zy)-as.numeric(t(x)%*%y)*(1-(s-2)/n)*(s-2)/n
      den<-as.numeric(t(zx)%*%solve(zz)%*%zx)-as.numeric(t(x)%*%x)*(1-(s-2)/n)*(s-2)/n
      b_est<-num/den
      err<-sum((y-x*as.numeric(b_est))^2)/(n*abs(den))
    } else if (s==1){
      num<-as.numeric(t(zx)%*%solve(zz)%*%zy)-as.numeric(t(x)%*%y)*(1-(s-2)/n)*(s-2)/n
      den<-as.numeric(t(zx)%*%solve(zz)%*%zx)-as.numeric(t(x)%*%x)*(1-(s-2)/n)*(s-2)/n
      b_est<-num/den
      err<-sum((y-x*as.numeric(b_est))^2)/(n*abs(den))
    }
  }
  
  return(err)
  
}


##### Confidence intervals

### Confidence interval based on Gaussian asymptotics

ci_an<-function(b_est,sig,coverage){
  conf<-c(b_est+qnorm((1-coverage)/2)*sqrt(sig),b_est+qnorm(1-(1-coverage)/2)*sqrt(sig))
  return(conf)
}

### Weak identification robust confidence interval

ci_robust<-function(y,x,z,coverage,type_ci="AR"){
  
  iv_test<-ivmodel::ivmodel(Y=y,D=x,Z=z,intercept=F)
  
  if (type_ci=="AR"){
    ar_ci<-ivmodel::AR.test(iv_test,beta0=0,1-coverage)$ci
    if (dim(ar_ci)[1]>1){
      cif<-c(-Inf,Inf)
    } else if (dim(ar_ci)[1]==1 & is.na(ar_ci[1])==1){
      cif<-NA
    } else if (dim(ar_ci)[1]==1 & is.na(ar_ci[1])==0){
      cif<-ar_ci
    }
  } else if (type_ci=="CLR"){
    
    clr_ci<-ivmodel::CLR(iv_test,beta0=0,1-coverage)$ci
    if (dim(clr_ci)[1]>1){
      cif<-c(-Inf,Inf)
    } else if (dim(clr_ci)[1]==1 & is.na(clr_ci[1])==1){
      cif<-NA
    } else if (dim(clr_ci)[1]==1 & is.na(clr_ci[1])==0){
      cif<-clr_ci
    }
  }
  
  return(cif)
}


##### Building the risks

### From the paper

## R_EXO

R1<-function(y,x,z,b_est){
  z<-as.matrix(z)
  n<-dim(z)[1]
  s<-dim(z)[2]
  sigma_z<-t(z)%*%z/n
  if (s>1){
    risk<-as.numeric(t(y-x*b_est)%*%z%*%solve(sigma_z)%*%t(z)%*%(y-x*b_est)/n^2)/var(x)^2/(s+1)
  } else if (s==1){
    risk<-as.numeric(t(y-x*b_est)%*%z*t(z)%*%(y-x*b_est)/n^2/sigma_z)/var(x)^2
  }
  
  return(risk)
}

## R_PMSE

R2<-function(y,x,z,b_est){
  z<-as.matrix(z)
  n<-dim(z)[1]
  s<-dim(z)[2]
  zz<-t(z)%*%z
  zx<-t(z)%*%x
  if (s>1){
    pi_z<-solve(zz)%*%zx
  } else if (s==1){
    pi_z<-zx/zz
  } 
  risk<-sum((y-z%*%pi_z%*%b_est)^2)/n
  return(risk)
}

## R_WMSE

R3<-function(y,x,z,b_est){
  z<-as.matrix(z)
  n<-dim(z)[1]
  s<-dim(z)[2]
  sigma_z<-t(z)%*%z/n
  if (s>1){
    d<-apply(z,1,function(e) t(e)%*%solve(sigma_z)%*%e)
  } else if (s==1){
    d<-apply(z,1,function(e) e^2/sigma_z)
  }
  
  err<-sapply(1:n, function(e) (y[e]-x[e]*b_est)^2)
  risk<-sum(d*err)/n/var(x)
  return(risk)
}

## R_MSE

R4<-function(y,x,z,b_est){
  n<-length(y)
  risk<-sum((y-x*b_est)^2)/n/var(x)
  return(risk)
}

### Final risk function

Risk<-function(y,x,z,b_est,type="1"){
  z<-as.matrix(z)
  if (type=="1"){
    r<-R1(y,x,z,b_est)
  } else if (type=="2"){
    r<-R2(y,x,z,b_est)
  } else if (type=="3"){
    r<-R3(y,x,z,b_est)
  } else if (type=="4"){
    r<-R4(y,x,z,b_est)
  } else if (type=="J"){
    r<-R1(y,x,z,b_est)
  }
  return(r)
}


##### Resampling

### Bootstrap sample

resample_bt<-function(n,kc){
  ind<-sample(1:n,n,replace=TRUE)
  return(ind)
}

### Cross-validated sample

resample_cv<-function(n,kc){
  ind_pre<-sample(1:n,n*kc,replace=FALSE)
  ind<-c(ind_pre,(1:n)[-ind_pre])
  return(ind)
}

### Disjoint out-of-bag bootstrap sample

resample_oob<-function(n,kc){
  ind_pre<-sample(1:n,n*kc,replace=FALSE)
  ind1<-sample(ind_pre,n*kc,replace=TRUE)
  ind2<-sample((1:n)[-ind_pre],replace=TRUE)
  ind<-c(ind1,ind2)
  return(ind)
}

### Out-of-bag bootstrap sample

resample_oob2<-function(n,kc){
  ind_pre<-sample(1:n,n*kc,replace=TRUE)
  ind<-c(ind_pre,(1:n)[-ind_pre])
  return(ind)
}

### Final resample function

resample<-function(n,kc,type="bt"){
  if (type=="bt"){
    ind<-sample(1:n,n,replace=TRUE)
  } else if (type=="cv"){
    ind_pre<-sample(1:n,n*kc,replace=FALSE)
    ind<-c(ind_pre,(1:n)[-ind_pre])
  } else if (type=="oob"){
    ind_pre<-sample(1:n,n*kc,replace=FALSE)
    ind1<-sample(ind_pre,n*kc,replace=TRUE)
    ind2<-sample((1:n)[-ind_pre],replace=TRUE)
    ind<-c(ind1,ind2)
  } else if (type=="oob2"){
    ind_pre<-sample(1:n,n*kc,replace=TRUE)
    ind<-c(ind_pre,(1:n)[-ind_pre])
  }
  return(ind)
}

### Resampling B times with proportion kc for training

sample_boot<-function(n,kc,B,type="bt"){
  if (type=="cv2"){
    ind<-matrix(1:n,n/B,B)
  } else {
    ind<-sapply(1:B, function(e) resample(n,kc,type))
  }
  return(ind)
}

##### Risk Estimators

### Apparent risk Estimator

risk_app<-function(y,x,z,type_est="2sls",type_risk="1"){
  b_est<-as.numeric(est(y,x,z,type_est))
  r<-Risk(y,x,z,b_est,type_risk)
  return(r)
}

### Bootstrap risk estimator

risk_bt<-function(y,x,z,SB,type_est="2sls",type_risk="1"){
  n<-length(y)
  nB<-length(SB)
  SB2<-matrix(SB,n,nB/n)
  B<-dim(SB2)[2]
  r_vec<-apply(SB2,2, function(e) risk_app(y[e],x[e],z[e,],type_est,type_risk))
  r<-sum(r_vec)/B
  return(r)
}

### Risk estimator using cross-validation or out-of-bag

risk_out<-function(y,x,z,kc,SB,type_est="2sls",type_risk="1"){
  n<-length(y)
  nB<-length(SB)
  SBf<-matrix(SB,n,nB/n)
  z<-as.matrix(z)
  B<-dim(SBf)[2]
  SB1<-SBf[(1:(n*kc)),]
  SB2<-SBf[((n*kc+1):n),]
  b_est<-apply(SB1,2, function(e) est(y[e],x[e],z[e,],type_est))
  r_vec<-sapply(1:B, function(e) Risk(y[SB2[,e]],x[SB2[,e]],z[SB2[,e],],b_est[e],type_risk))
  r<-sum(r_vec)/B
  return(r)
}

risk_out2<-function(y,x,z,kc,B,type_est="2sls",type_risk="1"){
  n<-length(y)
  SBf<-matrix(resample(n,kc,"cv"),n/B,B)
  SB1<-SBf[(1:(n*kc/B)),]
  SB2<-SBf[((n*kc/B+1):(n/B)),]
  b_est<-apply(SB1,2, function(e) est(y[e],x[e],z[e,],type_est))
  r_vec<-sapply(1:B, function(e) Risk(y[SB2[,e]],x[SB2[,e]],z[SB2[,e],],b_est[e],type_risk))
  r<-sum(r_vec)/B
  return(r)
}


##### Risk estimators in practice

### Risk estimators R_EXO, R_PMSE and R_MSE at the same time

risk_est_bit<-function(y,x,z,kc,B,SB,type_est,test_z,inc=T,cv=1){
  
  r<-rep(0,3)
  z<-as.matrix(z)
  kz<-dim(z)[2]
  vx<-var(x)
  
  if (inc==T & length(test_z)>sum(test_z)){
    zs<-as.matrix(z[,test_z==T])
    zbars<-z[,test_z==F]
    res<-lm(cbind(y,x,zs)~zbars+0)$residuals
    ys<-res[,1]
    xs<-res[,2]
    zs<-as.matrix(res[,3:(sum(test_z)+2)])
  } else {
    ys<-y
    xs<-x
    zs<-as.matrix(z[,test_z==T])
  }
  
  if (cv==1){
    r[1]<-risk_out(ys,xs,zs,kc,SB,type_est,type_risk="1")
    r[2]<-risk_out(ys,xs,zs,kc,SB,type_est,type_risk="2")
    r[3]<-risk_out(ys,xs,zs,kc,SB,type_est,type_risk="4")
  } else if (cv==0){
    r[1]<-risk_app(ys,xs,zs,type_est,type_risk="1")
    r[2]<-risk_app(ys,xs,zs,type_est,type_risk="2")
    r[3]<-risk_app(ys,xs,zs,type_est,type_risk="4")
  } else if (cv==2){
    r[1]<-risk_bt(ys,xs,zs,SB,type_est,type_risk="1")
    r[2]<-risk_bt(ys,xs,zs,SB,type_est,type_risk="2")
    r[3]<-risk_bt(ys,xs,zs,SB,type_est,type_risk="4")
  } else if (cv==3){
    r[1]<-risk_out2(ys,xs,zs,kc,B,type_est,type_risk="1")
    r[2]<-risk_out2(ys,xs,zs,kc,B,type_est,type_risk="2")
    r[3]<-risk_out2(ys,xs,zs,kc,B,type_est,type_risk="4")
  }
  
  return(r)
}

risk_est<-function(y,x,z,kc,B,SB,type_est,test_z,inc=T,cv=1){
  
  test_z<-as.matrix(test_z)
  n_test<-dim(test_z)[2]
  
  if (n_test>1){
    r<-t(sapply(1:n_test, function(e) risk_est_bit(y,x,z,kc,B,SB,type_est,test_z[,e],inc,cv)))
  } else if (n_test==1){
    r<-risk_est_bit(y,x,z,kc,B,SB,type_est,test_z,inc,cv)
  }
  
  return(r)
  
}

### Donald Newey (2001)

crit_mse_dn<-function(y,x,z,sig_v,pre_H,type_sig="cv"){
  
  z<-as.matrix(z)
  v<-lm(x~z+0)$residuals/pre_H
  
  if (type_sig=="cv"){
    d<-apply(z,1, function(e) t(e)%*%solve(t(z)%*%z)%*%e)
    R<-mean((v/(1-d))^2)
  } else if (type_sig=="mallows"){
    n<-length(y)
    s<-dim(as.matrix(z))[2]
    R<-mean(v^2)+sig_v*2*s/n
  }
  return(R)
}

risk_dn_bit<-function(y,x,z,test_z,type_sig="cv",type_est="2sls",inc=T,oracle=T){
  
  n<-length(y)
  z<-as.matrix(z)
  
  if (oracle==T){
    pre_est<-est(y,x,z[,1],type_est)
    pre_H<-mean(lm(x~z[,1]+0)$fitted.values^2)/n
  } else if (oracle==F){
    pre_est<-est(y,x,z[,test_z==T],type_est)
    pre_H<-mean(lm(x~z[,test_z==T]+0)$fitted.values^2)/n
  }
  
  pre_u<-y-x*as.numeric(pre_est)
  pre_v<-lm(x~z+0)$residuals/pre_H
  
  sig_u<-t(pre_u)%*%pre_u/n
  sig_v<-t(pre_v)%*%pre_v/n
  rho<-t(pre_u)%*%pre_v/n
  
  zs<-z[,test_z==T]
  
  if (inc==T & length(test_z)>sum(test_z)){
    zbars<-z[,test_z==F]
    res<-lm(cbind(y,x,zs)~zbars+0)$residuals
    y<-res[,1]
    x<-res[,2]
    zs<-res[,3:(sum(test_z)+2)]
  }
  
  if (type_est=="jive" || type_est=="b2sls"){
    
    RK<-crit_mse_dn(y,x,z[,test_z==1],sig_v,pre_H,type_sig)
    SK<-sig_u*(RK+rho*sum(test_z)/(n*sig_u))
    
  } else if (type_est=="2sls"){
    
    RK<-crit_mse_dn(y,x,z[,test_z==1],sig_v,pre_H,type_sig)
    
    SK1<-rho*sum(test_z)^2/n
    SK2<-sig_u*(RK-sig_v*sum(test_z)/n)
    SK<-SK1+SK2
    
  }
  
  return(SK)
}

risk_dn<-function(y,x,z,test_z,type_est="2sls",type_sig="cv",inc=T,oracle=T){
  
  test_z<-as.matrix(test_z)
  n_test<-dim(test_z)[2]
  
  if (n_test>1){
    r<-t(sapply(1:n_test, function(e) risk_dn_bit(y,x,z,test_z[,e],type_sig,type_est,inc,oracle)))
  } else if (n_test==1){
    r<-risk_dn_bit(y,x,z,test_z,type_sig,type_est,inc,oracle)
  }
  return(r)
  
}


### Andrews (1999)

risk_ad_bit<-function(y,x,z,test_z,type_est="2sls"){
  
  n<-length(y)
  z<-as.matrix(z)
  
  zs<-as.matrix(z[,test_z==T])
  
  est_b<-as.numeric(est(y,x,zs,type=type_est))
  uh<-(y-x*est_b)
  uhz<-t(uh)%*%zs
  if (sum(test_z)>1){
    zu2<-t(sapply(1:length(y), function(e) zs[e,]*uh[e]^2))
  } else {
    zu2<-sapply(1:length(y), function(e) zs[e,]*uh[e]^2)
  }  
  zzu2<-t(zu2)%*%zs
  J<-uhz%*%solve(zzu2)%*%t(uhz)
  Jf<-J-log(n)*(sum(test_z)-1)
  
  return(Jf)
}


risk_ad<-function(y,x,z,test_z,type_est="2sls"){
  test_z<-as.matrix(test_z)
  n_test<-dim(test_z)[2]
  
  if (n_test>1){
    r<-t(sapply(1:n_test, function(e) risk_ad_bit(y,x,z,test_z[,e],type_est)))
  } else if (n_test==1){
    r<-risk_ad_bit(y,x,z,test_z,type_est)
  }
  return(r)
}



### Kang (2016) Post-Lasso

risk_lasso<-function(y,x,z,lambda=F,K=10,normalize=F,package=T){
  if (package==T){
    fit<-cv.sisVIVE(y,x,z,lambdaSeq=lambda,K,intercept=F,normalize=normalize)
    LASSO_iv<-as.numeric(fit$alpha==0)
  } else if (package==F){
    n<-length(y)
    dh<-lm(x~z+0)$fitted.values
    zy<-lm(y~z+0)$fitted.values
    yf<-lm(zy~dh+0)$residuals
    zf<-lm(z~dh+0)$residuals
    lambda_cv<-cv.glmnet(x=zf*n,y=yf*n,nfolds=K,alpha=1,intercept=F,normalize=normalize)$lambda.min
    fit<-glmnet(x=zf*n,y=yf*n,lambda=lambda_cv,alpha=1,intercept=F,normalize=normalize)
    alpha_f<-fit$beta
    LASSO_iv<-as.numeric(alpha_f==0)
  }
  if (sum(LASSO_iv)==0){
    LASSO_iv<-rep(1,kz)
  }
  return(LASSO_iv)
}

### Windmeijer (2019) Post-Adaptive-Lasso

risk_adalasso<-function(y,x,z,lambda=F,K=10,normalize=F){
  
  z<-as.matrix(z)
  kz<-dim(z)[2]
  n<-dim(z)[1]
  
  par_yz<-solve(t(z)%*%z)%*%t(z)%*%y
  par_xz<-solve(t(z)%*%z)%*%t(z)%*%x
  med_est<-median(sapply(1:kz,function(e) par_yz[e]/par_xz[e]))
  alpha_pre<-par_yz-par_xz*med_est
  alpha_pre[alpha_pre==0]<-10^{-6}
  xz<-lm(x~z+0)$fitted.values
  zd<-lm(z~xz+0)$residuals
  zf<-sapply(1:kz, function(e) zd[,e]*abs(alpha_pre[e]))
  
  if (lambda==F){
    lambda_cv<-cv.glmnet(x=zf*n,y=y*n,nfolds=K,alpha=1,intercept=F,normalize=normalize)$lambda.min
  } else {
    lambda_cv<-cv.glmnet(x=zf*n,y=y*n,lambda=lambda,alpha=1,intercept=F,normalize=normalize)$lambda.min
  }
  fit<-glmnet(x=zf*n,y=y*n,lambda=lambda_cv,alpha=1,intercept=F,normalize=normalize)
  alpha_r<-fit$beta
  alpha_f<-sapply(1:kz, function(e) alpha_r[e,]/abs(alpha_pre[e]))/n
  ADALASSO_iv<-as.numeric(alpha_f==0)
  if (sum(ADALASSO_iv)==0){
    ADALASSO_iv<-rep(1,kz)
  }
  return(ADALASSO_iv)
}

##### All possible combinations of IVs

mat_z<-function(kz){
  pos_val<-c(0,1)
  arr_pos_val<-lapply(numeric(kz), function(e) pos_val)
  mat_pos<-t(as.matrix(expand.grid(arr_pos_val)))[,-1]
  return(mat_pos)
}

##### Diagnostics empirical

diag_est_bit<-function(y,x,z,bet,coverage,test_z,type_est="2sls",type_ci="AR",inc=T){
  
  if (inc==T & length(test_z)>sum(test_z)){
    zbars<-z[,test_z==F]
    zs<-z[,test_z==T]
    res<-lm(cbind(y,x,zs)~zbars+0)$residuals
    y<-res[,1]
    x<-res[,2]
    zs<-as.matrix(res[,3:(sum(test_z)+2)])
    b_est<-as.numeric(est(y,x,zs,type_est))
    sig_est<-b_sig(y,x,zs,type_est)
    abs_bias<-abs(b_est-bet)
    sq_bias<-(b_est-bet)^2
    cian<-ci_an(b_est,sig_est,coverage)
    cov_an<-(cian[1]<=bet & bet<=cian[2])
    len_an<-abs(cian[2]-cian[1])
    ciri<-ci_robust(y,x,zs,coverage,type_ci)
    if (is.na(ciri)==0){
      len_ri<-abs(ciri[2]-ciri[1])
      if (len_ri==Inf){
        cov_ri<-1
        good_ri<-0
      } else {
        cov_ri<-c(ciri[1]<=bet & bet<=ciri[2])
        good_ri<-1
      }
    } else {
      len_ri<-NA
      cov_ri<-NA
      good_ri<-0
    }
  } else {
    b_est<-as.numeric(est(y,x,z[,test_z==T],type_est))
    sig_est<-b_sig(y,x,z[,test_z==T],type_est)
    abs_bias<-abs(b_est-bet)
    sq_bias<-(b_est-bet)^2
    cian<-ci_an(b_est,sig_est,coverage)
    cov_an<-(cian[1]<=bet & bet<=cian[2])
    len_an<-abs(cian[2]-cian[1])
    ciri<-ci_robust(y,x,z[,test_z==T],coverage,type_ci)
    if (is.na(ciri)==0){
      len_ri<-abs(ciri[2]-ciri[1])
      if (len_ri==Inf){
        cov_ri<-1
        good_ri<-0
      } else {
        cov_ri<-c(ciri[1]<=bet & bet<=ciri[2])
        good_ri<-1
      }
    } else {
      len_ri<-NA
      cov_ri<-NA
      good_ri<-0
    }
  }
  
  return(c(b_est,sqrt(sig_est),abs_bias,sq_bias,cov_an,len_an,cov_ri,len_ri,good_ri,sum(test_z)))
  
}

diag_est<-function(y,x,z,bet,coverage,test_z,type_est="2sls",type_ci="AR",inc=T){
  
  de<-apply(test_z,2, function(e) diag_est_bit(y,x,z,bet,coverage,e,type_est,type_ci,inc))
  
  return(de)
  
}

diag_estf<-function(set_diag,ri=F){
  
  intqr<-apply(set_diag[1:8,],1,function(e) IQR(e))
  med_abs_bias<-apply(set_diag[17:24,],1, function(e) median(e))
  med_sq_bias<-apply(set_diag[25:32,],1, function(e) median(e))
  cov_an<-apply(set_diag[33:40,],1, function(e) mean(e))
  med_len_an<-apply(set_diag[41:48,],1, function(e) median(e))
  cov_r<-apply(set_diag[49:56,],1, function(e) mean(e))
  med_len_r<-apply(set_diag[57:64,],1, function(e) median(e[is.na(e)==0 & e!=Inf]))
  good_r<-apply(set_diag[65:72,],1, function(e) mean(e))
  avg_iv<-apply(set_diag[73:80,],1, function(e) mean(e))
  
  
  mat_diag<-cbind(intqr,med_abs_bias,med_sq_bias,cov_an,med_len_an,cov_r,med_len_r,good_r,avg_iv)
  rownames(mat_diag)<-c("DN","AN","Post-Lasso","Post-adalasso","R_EXO","R_PMSE","R_MSE","Oracle")
  colnames(mat_diag)<-c("Int Qr R","Med Abs Bias",
                        "Med Sq Bias","Cover An","Med Len An","Cover R","Med Len R","% finite CI","Avg Num IVs")
  
  
  return(xtable(mat_diag))
  
}

### Notation for estimates with p-values

pval_notation<-function(x,pval){
  if (pval<=0.001) {
    val<-paste(as.character(x),"$^{****}$")
  } else if (pval>0.001 & pval<=0.01){
    val<-paste(as.character(x),"$^{***}$")
  } else if (pval>0.01 & pval<=0.05){
    val<-paste(as.character(x),"$^{**}$")
  } else if (pval>0.05 & pval<=0.1){
    val<-paste(as.character(x),"$^{*}$")
  } else {
    val<-as.character(x)
  }
  return(val)
}

### First stage F-statistic

Fstats<-function(y,x,bet,clus=F,clus_x){
  x<-as.matrix(x)
  if (dim(x)[2]>1){
  uh<-y-x%*%lm(y~x+0)$coefficients
  if (clus==F){
    xuh<-sapply(1:length(x), function(e) x[e]%*%uh[e]^2)
    nume<-t(x)%*%xuh
    xx<-t(x)%*%x
    Fstat<-solve(xx)%*%xx/denom
  } else if (clus==T){
    nume<-matrix(rowSums(sapply(1:max(clus_x),
                     function(e) t(x[clus_x==e,1:7])%*%uh[clus_x==e]%*%t(t(x[clus_x==e,1:7])%*%uh[clus_x==e]))),7,7)
    denom<-sum(x^2)^2
    Fstat<-nume/denom
  }} else {
    uh<-y-x*lm(y~x+0)$coefficients
    if (clus==F){
      xuh<-sapply(1:length(x), function(e) x[e]%*%uh[e]^2)
      nume<-t(x)%*%xuh
      denom<-sum(x^2)^2
      Fstat<-nume/denom
    } else if (clus==T){
      nume<-sum(sapply(1:max(clus_x),
                                  function(e) (t(x[clus_x==e,1:7])%*%uh[clus_x==e])^2))
      denom<-sum(x^2)^2
      Fstat<-nume/denom
    }
  }
  return(Fstat)
}

### Sargan Hansen (1982) J statistic

Jstats<-function(y,x,z,bet,clus=F,clus_x){
  z<-as.matrix(z)
  uh<-y-x*as.numeric(bet)
  if (clus==F){
    zu<-t(uh)%*%z
    if (dim(z)[2]>1){
      zu2<-t(sapply(1:length(y), function(e) z[e,]*uh[e]^2))
      zzu2<-t(zu2)%*%z
      J_iv<-zu%*%solve(zzu2)%*%t(zu)
      p_val_j<-1-pchisq(J_iv,df=dim(z)[2]-1)
    } else {
      J_iv<-0
      p_val_j<-1
    }
  } else if (clus==T){
    
  }
  return(c(J_iv,p_val_j))
}

### Final function for diagnotics

diag_appli<-function(y,x,z,test_z,type_est="2sls",se="het",ivex=2){
  
  if (type_est=="ols"){
    
  } else {
    if (sum(test_z)<length(test_z)){
    
      ys<-lm(y~z[,test_z==F]+0)$residuals
      xs<-lm(x~z[,test_z==F]+0)$residuals
      zs<-as.matrix(lm(z[,test_z==T]~z[,test_z==F]+0)$residuals)
      
    } else {
      ys<-y
      xs<-x
      zs<-z
    }
    
   estim<-format(est(ys,xs,zs,type=type_est),4,nsmall=4)
   sd_est<-format(sqrt(b_sig(ys,xs,zs,type=type_est)),4,nsmall=4)
   p_val<-format(2*(1-pnorm(abs(as.numeric(estim)/as.numeric(sd_est)))),4,nsmall=4)
   iv_ind<-(1:dim(z)[2])[test_z==T]
   iv_ind<-sapply(iv_ind, function(e) e*(e<ivex)+(e+1)*(e>=ivex))
   iv_rem<-paste("{",paste(iv_ind,collapse=";"),"}",sep="")
  
   F_iv<-as.numeric(summary(lm(xs~zs+0))$fstatistic)
   F_iv<-format(round(c(F_iv[1],pf(F_iv[1],F_iv[2],F_iv[3],lower.tail=F)),digits=4),4,nsmall=4)
   J_iv<-format(round(Jstats(ys,xs,zs,estim),digits=4),4,nsmall=4)
   
   out<-c(pval_notation(estim,p_val),sd_est,iv_rem,pval_notation(F_iv[1],F_iv[2]),pval_notation(J_iv[1],J_iv[2]))
   names(out)<-c("estimator","standard deviation","judge used as iv","F statistic","J statistic")
  }
  
  return(out)
}

### Diagnostics for OLS

diag_est_ols<-function(y,x,bet,coverage){
  
  n<-length(y)
  b_ols<-as.numeric(t(x)%*%y/t(x)%*%x)
  var_ols<-sum((y-x*b_ols)^2)/n/as.numeric((t(x)%*%x))
  cian<-c(b_ols+qnorm((1-coverage)/2)*sqrt(var_ols),b_ols+qnorm(1-(1-coverage)/2)*sqrt(var_ols))
  abs_bias<-abs(b_ols-bet)
  sq_bias<-(b_ols-bet)^2
  len_an<-abs(cian[2]-cian[1])
  cov_an<-(cian[1]<=bet & bet<=cian[2])
  
  return(c(b_ols,sqrt(var_ols),abs_bias,sq_bias,cov_an,len_an))
}

##### Theoretical diagnostics for simulations

### Function to pick set S

A_z<-function(test_z){
  nz<-sum(test_z)
  nkz<-length(test_z)
  A<-matrix(0,nkz,1)
  for (i in 1:nkz){
    b<-rep(0,nkz)
    if (test_z[i]==1){
      b[i]<-1
      A<-cbind(A,b)
    }
  }
  A<-A[,-1]
  return(A)
}

### Concentration parameter bar(S) not included as controls

mu2_ninc<-function(n,test_z,pi,sig_z,sig_v){
  A_S<-as.matrix(A_z(test_z))
  A_barS<-as.matrix(A_z(1-test_z))
  sig_S<-t(A_S)%*%sig_z%*%A_S
  sig_barS<-t(A_barS)%*%sig_z%*%A_barS
  sig_SbarS<-t(A_S)%*%sig_z%*%A_barS
  pi_S<-t(A_S)%*%pi
  pi_barS<-t(A_barS)%*%pi
  pi_Sf<-pi_S+solve(sig_S)%*%sig_SbarS%*%pi_barS
  sig_vS<-sig_v+t(pi_barS)%*%sig_barS%*%pi_barS-t(pi_barS)%*%t(sig_SbarS)%*%solve(sig_S)%*%sig_SbarS%*%pi_barS
  mu2<-n*t(pi_Sf)%*%sig_S%*%pi_Sf/sig_vS/sum(test_z)
  return(mu2)
}

### Concentration parameter bar(S) included as controls

mu2_inc<-function(n,test_z,pi,sig_z,sig_v){
  if (length(test_z)==sum(test_z)){
    mu2<-mu2_ninc(n,test_z,pi,sig_z,sig_v)
  } else {
    A_S<-as.matrix(A_z(test_z))
    A_barS<-as.matrix(A_z(1-test_z))
    sig_S<-t(A_S)%*%sig_z%*%A_S
    sig_barS<-t(A_barS)%*%sig_z%*%A_barS
    sig_SbarS<-t(A_S)%*%sig_z%*%A_barS
    pi_S<-t(A_S)%*%pi
    sig_Sf<-sig_S-sig_SbarS%*%solve(sig_barS)%*%t(sig_SbarS)
    mu2<-n*t(pi_S)%*%sig_Sf%*%pi_S/sig_v/sum(test_z)
  }
  return(mu2)
}

### Bias of OLS

bias_ols<-function(rho,pi,true_z,alpha_true,sig_z,sig_v){
  A_Strue<-as.matrix(A_z(true_z))
  sig_zStrue<-sig_z%*%A_Strue
  bias<-(rho+t(pi)%*%sig_zStrue%*%alpha_true)/(sig_v+t(pi)%*%sig_z%*%pi)
  return(bias)
}

### Bias of 2SLS, bar(S) included as controls

bias_2sls<-function(pi,true_z,test_z,alpha_true,sig_z){
  A_S<-as.matrix(A_z(test_z))
  A_barS<-as.matrix(A_z(1-test_z))
  A_Strue<-as.matrix(A_z(true_z))
  sig_S<-t(A_S)%*%sig_z%*%A_S
  sig_barS<-t(A_barS)%*%sig_z%*%A_barS
  sig_SbarS<-t(A_S)%*%sig_z%*%A_barS
  sig_SStrue<-t(A_S)%*%sig_z%*%A_Strue
  sig_barSStrue<-t(A_barS)%*%sig_z%*%A_Strue
  pi_S<-t(A_S)%*%pi
  if (length(test_z)==sum(test_z)){
    sig_Sf<-sig_S
    sig_SStruef<-sig_SStrue
    bias<-t(pi_S)%*%sig_SStruef%*%alpha_true/(t(pi_S)%*%sig_Sf%*%pi_S)
  } else {
    sig_Sf<-sig_S-sig_SbarS%*%solve(sig_barS)%*%t(sig_SbarS)
    sig_SStruef<-sig_SStrue-sig_SbarS%*%solve(sig_barS)%*%sig_barSStrue
    bias<-t(pi_S)%*%sig_SStruef%*%alpha_true/(t(pi_S)%*%sig_Sf%*%pi_S)
  }
  return(bias)
}

### Final functions for simulations theoretical diagnostics

diag_theor_bit<-function(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true){
  if (type_diag=="ninc"){
    dth<-c(mu2_ninc(n,test_z,pi,sig_z,sig_v),bias_ols(rho,pi,true_z,alpha_true,sig_z,sig_v),
           bias_2sls(pi,true_z,test_z,alpha_true,sig_z),sum(test_z))
  } else if (type_diag=="inc"){
    dth<-c(mu2_ninc(n,test_z,pi,sig_z,sig_v),mu2_inc(n,test_z,pi,sig_z,sig_v),
           bias_ols(rho,pi,true_z,alpha_true,sig_z,sig_v),
           bias_2sls(pi,true_z,test_z,alpha_true,sig_z),sum(test_z))
  }
  return(dth)
}

diag_theor<-function(n,test_z,pi,sig_z,sig_v,rho,type_diag="inc",true_z,alpha_true){
  test_z<-as.matrix(test_z)
  n_test<-dim(test_z)[2]
  if (n_test>1){
    dth<-sapply(1:n_test, function(e) diag_theor_bit(n,test_z[,e],pi,sig_z,sig_v,rho,type_diag,true_z,alpha_true))
  } else if (n_test==1){
    dth<-diag_theor_bit(n,test_z,pi,sig_z,sig_v,rho,type_diag,true_z,alpha_true)
  }
  return(dth)
}


