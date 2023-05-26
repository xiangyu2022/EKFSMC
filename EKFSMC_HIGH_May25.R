resampling.res = function(mm,ww){
  mj = floor(mm*ww)
  index = rep.int(1,mj[1])
  for(j in 2:mm) index=c(index,rep.int(j,mj[j]))
  ww2 = mm*ww-mj
  if(sum(ww2)>0){
    ww2 = ww2/sum(ww2)
    mres = mm-sum(mj)
    indexrex = sample.int(mm,size=mres,replace=T,prob=ww2)
    index = c(index,indexrex)
  }
  return(index)
}
resampling.str = function(mm,ww){
  u <- runif(1)
  cusumw <- cumsum(ww)
  mj <- floor(mm*cusumw-u)+1
  mj[mj<0] <- 0
  index <- rep(1,mj[1])
  for(j in 2:mm) index <- c(index,rep.int(j,mj[j]-mj[j-1]))
  return(index)
}

data_gen=function(pop, I0, S0=pop-I0, N, t, dist_mat, d, lambda, reco_rate, 
                  trans_rate_inside, trans_rate_outside, sigma_I, sigma_y){
  x_S = x_I = x_R = y = matrix(NA, nrow = N, ncol = t+1)
  x_I[,1] = I0
  x_S[,1] = S0
  x_R[,1] = pop - I0 - S0
  Kd = dnorm(1, mean = dist_mat, sd = d)
  diag(Kd) = trans_rate_inside /  trans_rate_outside
  for(i in 1:t){
    ww = rnorm(N, 0, sigma_I)
    x_I[,i+1] = (1-reco_rate)*x_I[,i] + as.vector(matrix(diag(trans_rate_outside*x_S[,i]),ncol=N)%*% Kd %*% x_I[,i]) + ww
    x_S[,i+1] = x_S[,i] - as.vector(matrix(diag(trans_rate_outside*x_S[,i]),ncol=N)%*% Kd %*% x_I[,i]) - ww
    x_R[,i+1] = pop-x_I[,i+1] - x_S[,i+1]
    y[,i+1] = lambda*x_I[,i+1]+rnorm(N, 0, sigma_y)
    y[runif(10)>= 0.7,i+1] = NA
  }
  
  return(list(y=y, x_I=x_I, x_S=x_S, x_R=x_R, 
              reco_rate=reco_rate, trans_rate_inside=trans_rate_inside, 
              trans_rate_outside = trans_rate_outside,
              sigma_I=sigma_I, sigma_y=sigma_y, t=t, N=N, dist_mat, pop=pop))
}

ekf_smc=function(dat, m, resampling, lambda, dist_mat, d, SCORE=TRUE, reco_rate,
                 trans_rate_inside, trans_rate_outside, sigma_I, sigma_y){
  
  # Initialization
  t = dat$t
  y = dat$y
  N = dat$N
  S = Imat = lgW = lgW_lag = array(NA, c(N, m, t+1))
  alp = array(NA, c(N, m, t+1, 2))
  Imat[ , ,1] = matrix(rep(dat$x_I[,1],m), nrow=N, ncol=m, byrow=FALSE)
  S[ , ,1] = matrix(rep(dat$x_S[,1],m), nrow=N, ncol=m, byrow=FALSE)
  lgW[ , ,1] = 0
  lgW_lag[ , ,1] = 0
  alp[ , ,1 ,] = 0
  loglike = rep(0,N)

  mu_I = mu_cond_I = sigma = sigma_cond = array(NA, c(N, m, t+1))
  sigma[ , ,1] = matrix(rep(sigma_I, m), nrow=N, ncol=m, byrow=FALSE)
  mu_I[ , ,1] = matrix(rep(dat$x_I[,1], m), nrow=N, ncol=m, byrow=FALSE)
  Kd = Kd2 = dnorm(1, mean = dist_mat, sd = d)
  diag(Kd) = trans_rate_inside /  trans_rate_outside
  diag(Kd2) = 0

  KdImat = matrix(NA,nrow=N,ncol=m)
  for (j in 1:t) {
    # mu_cond_I[ , ,j+1] = (1 - reco_rate + trans_rate * S[ , ,j] )*mu_I[ , ,j]
    # sigma_cond[ , ,j+1] = sigma[ , ,j] * (1 - reco_rate + trans_rate *S[, ,j] - trans_rate * mu_I[, ,j] )^2 + sigma_I^2
    # mu_cond_I[ , ,j+1] = (1 - reco_rate + trans_rate * apply(Kd,1,sum) * S[ , ,j])*mu_I[ , ,j]
    for (b in 1:m){
      KdImat[,b] = apply(Kd2*Imat[, b,j],2,sum)
    }
    mu_cond_I[ , ,j+1] = (1 - reco_rate + trans_rate_inside * S[ , ,j] )*mu_I[ , ,j] + trans_rate_outside * S[ , ,j] * KdImat
    sigma_cond[ , ,j+1] = sigma[ , ,j] * (1 - reco_rate + trans_rate_inside *S[ , ,j] - trans_rate_inside *mu_I[ , ,j]  - trans_rate_outside * KdImat)^2 + sigma_I^2
    indx = is.na(y[,j+1])
    K_t = sigma_cond[ , ,j+1]/(sigma_cond[ , ,j+1] + sigma_y^2) 
    K_t[indx,] = y[indx,j+1] = 0
    mu_I[ , ,j+1] = mu_cond_I[ , ,j+1] + K_t * (y[ ,j+1] - mu_cond_I[ , ,j+1])
    sigma[ , ,j+1] = sigma_cond[ , ,j+1] * (1 - K_t)

    Imat[ , ,j+1] = mu_I[ , ,j+1] + matrix(rnorm(N*m, 0, sigma_I), nrow=N, ncol=m)
    S[ , ,j+1] = S[ , ,j] + (1-reco_rate)*Imat[ , ,j] - Imat[ , ,j+1] 
    lgW[ , ,j+1] = dnorm(Imat[ , ,j+1], 
                         mean=(1 - reco_rate + trans_rate_inside * S[ , ,j] )*Imat[ , ,j] + trans_rate_outside * S[ , ,j] * KdImat,
                         sd=sigma_I, log=T) + 
                   dnorm(y[,j+1], mean=Imat[ , ,j+1], sd=sigma_y, log=T) -
                   dnorm(Imat[ , ,j+1], mean=mu_I[ , ,j+1], sd=sigma_I, log=T) + lgW[ , ,j]
    # Scaling the weights
    maxlgW = apply(lgW[ , ,j+1], FUN=max, MARGIN=1)
    lgW[ , ,j+1] = lgW[ , ,j+1] - maxlgW
    lgW_lag[ , ,j+1] = lgW[ , ,j+1]
    loglike = loglike + maxlgW
    
    # Gradients
    if(SCORE!=0){
      if(SCORE==2){ # 2016
        for (l in 1:N){
          wtm1<-exp(lgW[l, ,j])
          grad_f<-cbind(S[l, ,j]*Imat[l, ,j]*(Imat[l, ,j+1] - (1 - reco_rate[l] + trans_rate_inside[l] * S[l, ,j])* Imat[l, ,j]), - Imat[l, ,j]*(Imat[l, ,j+1] - (1 - reco_rate[l] + trans_rate_inside[l] * S[l, ,j])* Imat[l, ,j]))/sigma_I^2
          Stm1<-c(sum(wtm1*alp[l, ,j,1]), sum(wtm1*alp[l, ,j,2]))/sum(wtm1)
          alp[l, ,j+1, ]<-lambda*alp[l, ,j, ]+(1-lambda)*matrix(Stm1, ncol=2)[rep(1,m),]+grad_f
        }
      } else if(SCORE==1){ # 2011
        for (l in 1:N){
          grad_f<-cbind(S[l, ,j]*Imat[l, ,j]*(Imat[l, ,j+1] - (1 - reco_rate[l] + trans_rate_inside[l] * S[l, ,j])* Imat[l, ,j]), - Imat[l, ,j]*(Imat[l, ,j+1] - (1 - reco_rate[l] + trans_rate_inside[l] * S[l, ,j])* Imat[l, ,j]))/sigma_I^2
          alp[l, ,j+1,]<-alp[l, ,j,]+grad_f
        }
        }
    }
    
    # Resampling
    if(resampling[j]==1){
      ww=exp(lgW[ , ,j+1])
      sumww=apply(ww,MARGIN=1,FUN=sum)
      ww_prob=ww/sumww
      for (l in 1:N){
        r_index=resampling.res(mm=m, ww=ww_prob[l,])
        Imat[l, ,j+1]=Imat[l,r_index,j+1];S[l, ,j+1]=S[l,r_index,j+1]
        lgW[l, ,j+1]=0
        loglike[l] =loglike[l]+log(mean(ww[l,]))}
      } 
    else if(resampling[j]==2){
      ww=exp(lgW[ , ,j+1])
      sumww=apply(ww,MARGIN=1,FUN=sum)
      ww_prob=ww/sumww
      for (l in 1:N){
      r_index=resampling.str(mm=m, ww=ww_prob[l,])
      Imat[l, ,j+1]=Imat[l,r_index,j+1];S[l, ,j+1]=S[l,r_index,j+1]
      lgW[l, ,j+1]=0
      loglike[l] =loglike[l]+log(mean(ww[l,]))
      }
    }
    
  }
  
  w_T=exp(lgW[ , ,t+1])
  if(resampling[t]==0){
    loglike=loglike+log(apply(w_T,1,mean))
  }
  
  if(SCORE!=0){
    score=matrix(NA,nrow=l,ncol=2)
    score_cond = array(NA, c(l,2,t))
    for (l in 1:N){
    score[l,]=c(sum(w_T[l,]*alp[l,,t+1,1]), sum(w_T[l,]*alp[l,,t+1,2]))/sum(w_T[l,])
    for (k in 1:t){
      score_cond[l, ,k]=c(sum(exp(lgW[l, ,k+1])*alp[l,,k+1,1]), sum(exp(lgW[l, ,k+1])*alp[l, ,k+1,2]))/sum(exp(lgW[l, ,k+1])) - c(sum(exp(lgW[l, ,k])*alp[l, ,k,1]), sum(exp(lgW[l, ,k])*alp[l, ,k,2]))/sum(exp(lgW[l, ,k]))
    }
  } 
  } else{
    score=NULL
    score_cond=NULL
  }
  
  return(list(S=S, Imat=Imat, mu_t_mat=mu_t_mat, mu_cond_mat=mu_cond_mat, lgW=lgW, lgW_lag=lgW_lag, loglike=loglike, score=score, score_cond=score_cond, lambda=lambda))
}

# Data generation
library(purrr)
set.seed(1)
pop = 50000
I0 = 100 + rnorm(10,1:10*10,10)
S0 = pop - I0
N = 10
t = 40
lambda = rep(1,10)
reco_rate = rnorm(N, mean = 1-(0.5)^(1/8), sd = 1e-2)
trans_rate_inside = rnorm(N, mean = reco_rate/2e4, sd = 1e-6)
trans_rate_outside = trans_rate_inside / 2
sigma_I = 10+rnorm(10,1:10,1)
sigma_y = 10+rnorm(10,1:10,1)
euclidean = function(s1, s2) sqrt(sum((s1 - s2)^2))
coordinates = matrix(rdunif(20,a=1,b=10),ncol=2)
dist_mat = matrix(NA, nrow=N, ncol=N)
for (i in 1:N){
  for (j in 1:N){
    dist_mat[i,j] = euclidean(coordinates[i,],coordinates[j,])
  }
}
d = 2
dat = data_gen(pop, I0, S0, N, t, dist_mat, d, lambda, reco_rate, trans_rate_inside, 
               trans_rate_outside, sigma_I, sigma_y)
m = 1000

#EKF-SMC
set.seed(1)
trans_rate=dat$trans_rate*1; reco_rate=dat$reco_rate*1
resu_plpo=ekf_smc(dat=dat, m=m, lambda=0.9, resampling=rep(rep(c(0,1,0,2),c(1,1,1,1)), 10),
                    SCORE=1, seed=NULL, trans_rate=trans_rate, reco_rate=reco_rate, 
                    sigma_I=dat$sigma_I, sigma_y=dat$sigma_y)


par(mar=c(0.5,0.5,0.5,0.5))
par(mfrow=c(3,3))
l=5 # This specifies which region you want to generate plots
for (tobs in 3:40){
  if (tobs%%8==2){
    plot(1,ylim=c(0,0.9),xaxt='n',yaxt='n',cex.main=2.5,ylab="",xlab="")
    legend(x="center", bty = "n",legend=c("proposals", "y_t", "I_t"),  
           col=1:3, pch=c(1, pch_y_I, pch_y_I), text.width=NULL, 
           y.intersp=0.95,cex=2)
  }
  
  w_obs=exp(resu_plpo$lgW_lag[l, ,tobs+1]) 
  w_obs_nl=w_obs/sum(w_obs)
  
  xmin=min(resu_plpo$Imat[l, ,tobs], dat$y[l, tobs], dat$x_I[l, tobs])
  xmax=max(resu_plpo$Imat[l, ,tobs], dat$y[l, tobs], dat$x_I[l, tobs])
  ymin=min(resu_plpo$Imat[l, ,tobs+1], dat$y[l, tobs+1], dat$x_I[l, tobs+1])
  ymax=max(resu_plpo$Imat[l, ,tobs+1], dat$y[l, tobs+1], dat$x_I[l, tobs+1]) 
  crtn=10
  
  pch_y_I=16
  plot(resu_plpo$Imat[l, ,c(tobs,tobs+1)], cex=w_obs_nl*25, pch=1, 
       xlim=c(xmin-crtn, xmax+crtn), ylim=c(ymin-crtn, ymax+crtn), 
       yaxt="n", xaxt="n"
       #main=paste0("y_t, I_t, and estimated I_t (proposals): t=", tobs-1), 
       #sub=paste0("Proposals sized by weights at t+1; weights sum to 1; m=", m)
       )
  points(x=dat$y[l, tobs], y=dat$y[l, tobs+1], col=2, cex=1.2, pch=pch_y_I); 
  points(x=dat$x_I[l, tobs], y=dat$x_I[l, tobs+1], col=3, cex=1.2, pch=pch_y_I); 
  legend(x="bottomright", legend=c(paste("t=",tobs-1,sep="")),cex=1.2)
}


