#class labels y (length n)
y = c(rep(-1,50),rep(1,50))

#simulation data X (p X n)
n=100
p=50
sig=0.25
mu0 = sig*rnorm(p)
mu1 = sig*rnorm(p)
y = c(rep(-1,50),rep(1,50))
X = matrix(nrow=p,ncol=n)
for(i in 1:n){
  if(y[i]==-1) X[,i] = mu0+rnorm(p)
  if(y[i]==1) X[,i] = mu1+rnorm(p)
}

#center X
for(i in 1:p) X[i,] = X[i,]-mean(X[i,])

#run BayesDWD for lambda2=1
res <- BayesDWD(X,y,lambda2=1)
#weights:
res$b.est
#scores:
res$b0.est+res$b.est%*%X
#probabilities of class membership:
res$prob.est
#each row of res$b.mat gives simulated weights from posterior distribution,
#can use these to create credible intervals, etc.


