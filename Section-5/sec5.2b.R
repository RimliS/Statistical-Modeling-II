rest <- read.csv("/Users/rimlisengupta/Dropbox/UT Austin/Spring 2018/SDS383D-StatMod2/Section5/Rest.csv",
                 sep = ",", header=TRUE)


#From a simulated dataset
nit <- 5000
n <- 1000

mu1 <- rep(0,nit)
mu2 <- rep(0,nit)
lambda1 <- rep(0,nit)
lambda2 <- rep(0,nit)
d <- matrix(0,n,nit)

y <- scale(rest$Profit)
pi <- 0.5

mu1[1] <- 2
mu2[1] <- 1
lambda1[1] <- 0.5
lambda2[1] <- 1

n1 <- NULL
n2 <- NULL
g <- matrix(0,n,2)
ind0 <- NULL
ind1 <- NULL
sumx1 <- 0
sumx2 <- 0

count <- 0
dnew <- rep(0,n)
for (m in 2:nit) {
  pvect <- pi*dnorm(y,mu1[m-1],sqrt(1/lambda1[m-1]))/
    (pi*dnorm(y,mu1[m-1],sqrt(1/lambda1[m-1]))+(1-pi)*dnorm(y,mu2[m-1],sqrt(1/lambda2[m-1])))  
  
  dnew <- rbinom(n,1,pvect)
  
  n1 <- sum(dnew==0)
  n2 <- sum(dnew==1)
  g <- cbind(y,dnew)
  ind0 <- which(g[,2]==0)
  ind1 <- which(g[,2]==1)
  sumx1 <- sum(g[ind0,1])
  sumx2 <- sum(g[ind1,1])
  
  mu1.new <- rnorm(1,(lambda1[m-1]*sumx1)/(0.01+n1*lambda1[m-1]),sqrt(1/(0.01+n1*lambda1[m-1])))
  mu2.new <- rnorm(1,(lambda2[m-1]*sumx2)/(0.01+n2*lambda2[m-1]),sqrt(1/(0.01+n2*lambda2[m-1])))
  lambda1.new <- rgamma(1,n/2+2,(1+0.5*(sum(y^2)-2*(mu1.new*sumx1+mu2.new*sumx2)
                                       +n1*mu1.new^2+n2*mu2.new^2)))
  lambda2.new <- rgamma(1,n/2+5,(1+0.5*(sum(y^2)-2*(mu1.new*sumx1+mu2.new*sumx2)
                                        +n1*mu1.new^2+n2*mu2.new^2)))
  
  d[,m] <- dnew
  mu1[m] <- mu1.new
  mu2[m] <- mu2.new
  lambda1[m] <- lambda1.new
  lambda2[m] <- lambda2.new
  
  count <- count+1
}


s <- seq(-4,4,len=1000)
density.est <- matrix(0,1000,nit)
fnew <- rep(0,1000)
density.avg <- rep(0,1000)

count=0
for (m in 1:nit) {
  fnew <- pi*dnorm(s,mu1[m],sqrt(1/lambda1[m]))+(1-pi)*dnorm(s,mu2[m],sqrt(1/lambda2[m]))
  density.est[,m] <- fnew
  count=count+1
}
loc <- 1:1000
density.est <- density.est[,-loc]

f <- apply(density.est,1,mean)




#From the given dataset
mu1 <- rep(0,nit)
mu2 <- rep(0,nit)
lambda1 <- rep(0,nit)
lambda2 <- rep(0,nit)
d <- rest$DinnerService

lambda1[1] <- 0.5
lambda2[1] <- 1

y <- scale(rest$Profit)
pi <- 0.5

mu1[1] <- 2
mu2[1] <- 1

n1 <- NULL
n2 <- NULL
g <- matrix(0,n,2)
ind0 <- NULL
ind1 <- NULL
sumx1 <- 0
sumx2 <- 0

count <- 0

for (m in 2:nit) {
  pvect <- pi*dnorm(y,mu1[m-1],sqrt(1/lambda1[m-1]))/
    (pi*dnorm(y,mu1[m-1],sqrt(1/lambda1[m-1]))+(1-pi)*dnorm(y,mu2[m-1],sqrt(1/lambda2[m-1])))  
  
  
  n1 <- sum(d==0)
  n2 <- sum(d==1)
  g <- cbind(y,d)
  ind0 <- which(g[,2]==0)
  ind1 <- which(g[,2]==1)
  sumx1 <- sum(g[ind0,1])
  sumx2 <- sum(g[ind1,1])
  
  mu1.new <- rnorm(1,(lambda1[m-1]*sumx1)/(0.01+n1*lambda1[m-1]),sqrt(1/(0.01+n1*lambda1[m-1])))
  mu2.new <- rnorm(1,(lambda2[m-1]*sumx2)/(0.01+n2*lambda2[m-1]),sqrt(1/(0.01+n2*lambda2[m-1])))
  lambda1[m] <- rgamma(1,n/2+2,(1+0.5*(sum(y^2)-2*(mu1.new*sumx1+mu2.new*sumx2)
                                        +n1*mu1.new^2+n2*mu2.new^2)))
  lambda2[m] <- rgamma(1,n/2+5,(1+0.5*(sum(y^2)-2*(mu1.new*sumx1+mu2.new*sumx2)
                                        +n1*mu1.new^2+n2*mu2.new^2)))
  
  mu1[m] <- mu1.new
  mu2[m] <- mu2.new
  #lambda[m] <- lambda.new
  
  count <- count+1
}


density.est1 <- matrix(0,1000,nit)
fnew1 <- rep(0,1000)
density.avg1 <- rep(0,1000)

count=0
for (m in 1:nit) {
  fnew1 <- pi*dnorm(s,mu1[m],sqrt(1/lambda1[m]))+(1-pi)*dnorm(s,mu2[m],sqrt(1/lambda2[m]))
  density.est1[,m] <- fnew1
  count=count+1
}
loc <- 1:1000
density.est1 <- density.est1[,-loc]

f1 <- apply(density.est1,1,mean)


plot(s,f,main="Density Estimate of the Mixture Model",ty="l")
lines(s,f1,ty="l",col="red")
