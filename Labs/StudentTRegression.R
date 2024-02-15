# Gibbs sampling for regression with student t errors
# Author: Mattias Villani, for the course Advanced Bayesian Learning

#install.packages("mvtnorm")
library(mvtnorm)

# Random number generator for scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

# Gibbs sampler for student-t regression
GibbsTReg <- function(y, X, nu, tau, nIter){
  n = dim(X)[1]
  q = dim(X)[2]
  beta_ = solve(crossprod(X),crossprod(X,y))
  epsilon = y - X%*%beta_
  sigma = sqrt(sum(epsilon^2)/(n-q))
  postDraws = matrix(0,nIter,q+1)
  for (i in seq(1,nIter)){
    
    # Update U
    U = rScaledInvChi2(n, nu + 1, (nu*sigma^2 + epsilon^2)/(nu+1))

    # Update beta
    A = diag(as.vector(1/(sqrt(U))))
    yt = A%*%y
    Xt = A%*%X
    betaHat = solve(crossprod(Xt),crossprod(Xt,yt)) # Least squares estimate
    SigmaPost = solve(crossprod(Xt) + tau^(-2)*diag(q))
    muPost = SigmaPost%*%(crossprod(Xt)%*%betaHat)
    beta_ = t(rmvnorm(1, muPost, SigmaPost))
    postDraws[i,1:q] = beta_
    epsilon = y - X%*%beta_
    
    # Update sigma
    sigma = sqrt(rgamma(1,n*nu/2, rate = (nu/2)*sum(1/U) ) )
    postDraws[i,q+1] = sigma
    
  }
  return(postDraws)
}

# Simulating from the predictive distribution of the student-t regression
PredTReg <- function(XPred, postDraws, nu){
  nPreds = dim(XPred)[1]
  q = dim(XPred)[2]
  nSim = dim(postDraws)[1]
  yPreds = matrix(0, nSim, nPreds)
  for (i in seq(1,nSim)){
    beta_ = postDraws[i,1:q]
    sigma = postDraws[i,q+1]
    yPreds[i,] = XPred%*%beta_ + sigma*rt(nPreds, nu)
  }
  return(yPreds)
}

# Reading data from the course repo
data = read.csv("https://github.com/mattiasvillani/AdvBayesLearnCourse/raw/master/Labs/regression.csv", sep ='\t', header = FALSE)
names(data) <- c("x","y")
X = cbind(1,data[,1])
y = data[,2]
n = length(y)
plot(X[,2],y)

# Simulating from the posterior of beta and sigma
tau = 1
nu = 2
nIter = 10000
postDraws = GibbsTReg(y, X, nu, tau, nIter)
colMeans(postDraws)
apply(postDraws, 2, sd)
par(mfrow = c(2,2))
hist(postDraws[,1], 50, main = "beta0")
hist(postDraws[,2], 50, main = "beta1")
hist(postDraws[,3], 50, main = "sigma")
hist(postDraws[,3]^2*(nu/(nu-2)), 50, main = "Var(y)")

# Normal approx of posterior by optimization
logLike <- function(theta, y, X, nu){
  n = dim(X)[1]
  p = dim(X)[2]
  beta = theta[1:p]
  sigma = exp(theta[p+1])
  epsilon = (y-X%*%beta)/sigma
  return(sum(log(1/sigma) + dt(epsilon, nu, log = TRUE)))
}

logPrior <- function(theta, tau){
  p = length(theta) - 1
  beta = theta[1:p]
  sigma = exp(theta[p+1])
  logPriorPDF = dmvnorm(beta, mean = rep(0, p), sigma = tau^2*diag(p), log = TRUE) + 
    log(1/sigma)
  return(logPriorPDF)
}

logPost <- function(theta, y, X, nu, tau){
  return(logLike(theta, y, X, nu) + logPrior(theta, tau) )
}

beta_ml = solve(crossprod(X))%*%crossprod(X,y)
sigma_ml = sqrt(crossprod(y-X%*%beta_ml)/n)
theta_ml = c(beta_ml,log(sigma_ml))
logPost(theta_ml, y, X, nu, tau)
logLike(theta_ml, y, X, nu)
logPrior(theta_ml, tau)

initVal = theta_ml
OptimResults <- optim(initVal, logPost, gr=NULL, y, X, nu, tau, method=c("BFGS"), 
  control=list(fnscale=-1), hessian=TRUE)
postMode = OptimResults$par
approxPostStd <- sqrt(diag(-solve(OptimResults$hessian)))
apply(postDraws, 2, sd)

logPost(theta_ml, y, X, nu, tau)
logPost(postMode, y, X, nu, tau)
postMode


# Making predictions and plotting predictive distribution
xGrid = seq(-1,1,by =.1)
XPred = cbind(1,xGrid)
yPreds = PredTReg(XPred, postDraws, nu)
predBands = matrix(0,dim(XPred)[1],5)
for (i in 1:length(xGrid)){
  predBands[i,] = quantile(yPreds[,i], c(0.005,0.025,0.5,0.975,0.995))
}
par(mfrow = c(1,1))
plot(X[,2], y, xlab = "x", ylab = "y")
lines(xGrid,predBands[,3], col = "red")
lines(xGrid,predBands[,1], col = "green")
lines(xGrid,predBands[,2], col = "blue")
lines(xGrid,predBands[,4], col = "blue")
lines(xGrid,predBands[,5], col = "green")

