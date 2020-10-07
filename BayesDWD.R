BayesDWD = function(X,y,burn_in=100,sample_iters = 1000, lambda2=1,est.lam2=FALSE,lam.prior='logUnif', Int = TRUE){
  p = dim(X)[1]
  n = dim(X)[2]
  b.mat = matrix(nrow=sample_iters,ncol=p)
#  prev.b=rep(0,p)
  prob.mat = matrix(nrow=sample_iters,ncol=n) 
  na.ind = is.na(y)
  l.miss=sum(na.ind)
  lambda2.vec=rep(lambda2,sample_iters) 
  if(est.lam2==TRUE){
    lambda2.seq = c(1/128,1/64, 1/32,1/16,1/8,1/4,1/2,1,2,4,8,16,32,64,128)
  #  Norm.prior = c()
    Norm.prior.log = c()
    for(i in 1:length(lambda2.seq)){
      print(i)
      lambda2=lambda2.seq[i]
      zsims = matrix(rnorm(p*10000),nrow=10000,ncol=p)
      zsims.trans = sqrt(1/(n*lambda2))*zsims 
      res = apply(X=t(zsims.trans),MARGIN=2,FUN=prior.unnorm,Xm=X)
#      Norm.prior[i] = mean(res)*(2*pi/(n*lambda2))^(p/2) ##should it be divided by? -- no
   #   res.log = apply(X=t(zsims.trans),MARGIN=2,FUN=prior.unnorm.log,Xm=X)
      Norm.prior.log[i] = log(mean(res))+(p/2)*(log(2*pi)-log(n*lambda2))
    }
 #   Norm.logLambda = approxfun(log(lambda2.seq),log(Norm.prior)) 
    Norm.logLambda = approxfun(log(lambda2.seq),Norm.prior.log) 
    lambda2=1
  }
  if(l.miss>0){
    y[na.ind] = rep(-1,l.miss) + 2*rbinom(l.miss,size=1,prob=0.5)
  }
  b0=0
  b0.vec=b0
  #initialize
  init.res = sdwd(x=t(X[,!na.ind]),y[!na.ind],lambda=0,lambda2=lambda2,standardize=FALSE)
  prev.b = as.vector(init.res$beta)
 # b0=as.numeric(init.res$b0)
  b0=0
  pb = txtProgressBar(min = 0, max = burn_in+sample_iters, initial = 0, style=3)
  for(t in 1:(burn_in+sample_iters)){
    setTxtProgressBar(pb,t)
    var=0.5*var(prev.b)
    U = runif(p)
    b.stars = rnorm(p,prev.b,sqrt(var))
    for(j in 1:p){
      b.star = prev.b
      b.star[j] = b.stars[j]
      r=exp(post.func.l(b0,b.star,y,X,lambda2,n)-post.func.l(b0,prev.b,y,X,lambda2,n))
      if(r >= U[j]) prev.b[j] = b.star[j]
    }
    if(Int==TRUE){
      b0.star = rnorm(1,b0,sqrt(0.25))
      U0 = runif(1)
      r=exp(post.func.l(b0.star,prev.b,y,X,lambda2,n)-post.func.l(b0,prev.b,y,X,lambda2,n))
      if(r >= U0) b0 = b0.star
    }
    if(l.miss>0) 
      y[na.ind] = rep(-1,l.miss) + 2*rbinom(l.miss,size=1,prob=prob.u(prev.b%*%X[,na.ind]))
    if(est.lam2){
      logLambda2.star = log(lambda2)+0.25*rnorm(1)
      lambda2.star = exp(logLambda2.star)
      if(lambda2.star>1/128 && lambda2.star<128){
        if(lam.prior=='logUnif'){
          r=exp(-0.5*n*lambda2.star*sum(prev.b^2)-Norm.logLambda(logLambda2.star)+0.5*n*lambda2*sum(prev.b^2)+Norm.logLambda(log(lambda2)))
        }
        if(lam.prior=='Unif'){
          r=exp(-0.5*n*lambda2.star*sum(prev.b^2)-Norm.logLambda(logLambda2.star)+0.5*n*lambda2*sum(prev.b^2)+Norm.logLambda(log(lambda2))+logLambda2.star-log(lambda2))
        }
        U.lam = runif(1)
        if(r >= U.lam) lambda2= lambda2.star
      }
    }
    if(t>burn_in){
      b.mat[t-burn_in,] = prev.b
      prob.mat[t-burn_in,] = prob.u(b0+prev.b%*%X)
      lambda2.vec[t-burn_in] = lambda2
      b0.vec[t-burn_in]=b0}
    }
prob.est=colMeans(prob.mat)
b.est=colMeans(b.mat)
b0.est=mean(b0.vec)
Results = list(prob.mat=prob.mat,b.mat=b.mat,prob.est=prob.est,b.est=b.est,b0.vec=b0.vec,b0.est=b0.est,lambda2.vec=lambda2.vec)
return(Results)
}

#post.func = function(b,y,X){
#  return(exp(-0.5*n*sum(b^2))*prod(exp(V(y*b%*%X))))
#} 

post.func.l = function(b0,b,y,X,lambda2,n){
  return(-0.5*n*lambda2*sum(b^2) - sum(V(y*(b0+b%*%X))))
} 

#V.unvec = function(u){
#  if(u <= .5) return(1-u)
#  if(u>0.5) return(1/(4*u))
#}

#V = Vectorize(V.unvec)
V = function(u){
  vu = 1-u
  vu[u>0.5]=1/(4*u[u>0.5])
  return(vu)
}

prob.u = function(u){
  return(exp(-V(u))/(exp(-V(u))+exp(-V(-u))))
}

prior.unnorm  = function(b,Xm){
  return(prod(0.5*exp(-V(b%*%Xm))+0.5*exp(-V(-b%*%Xm))))
} 

prior.unnorm.log  = function(b,Xm){
  return(sum(log(0.5*exp(-V(b%*%Xm))+0.5*exp(-V(-b%*%Xm)))))
} 




