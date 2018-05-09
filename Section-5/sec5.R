rest <- read.csv("/Users/rimlisengupta/Dropbox/UT Austin/Spring 2018/SDS383D-StatMod2/Section5/Rest.csv",
                    sep = ",", header=TRUE)


smp_size <- floor(0.5 * nrow(rest))
set.seed(123)

train_ind <- sample(seq_len(nrow(rest)), size = smp_size)

train <- rest[train_ind, ]
test <- rest[-train_ind, ]

train <- na.omit(train[,-1])

lmfit <- lm(Profit ~ factor(DinnerService)+SeatingCapacity,data=train)
y.hat <- predict(lmfit, test)
mse <- function(a,b) mean((a-b)^2)
rmse <- sqrt(mse(train$Profit,y.hat))
#6514.548

#Coefficients:
#  (Intercept)  factor(DinnerService)1         SeatingCapacity  
#52563                   20486                    4711


# length of chain
nit <- 5000
n <- 500
p <- 3

train$SeatingCapacity <- scale(train[,3])
X <- cbind(beta0=rep(1,n),train[,2:3])
X <- as.matrix(X)
K <- diag(p)
L <- diag(n)
mu <- c(0,0,0)

D <- t(X)%*%L%*%X
S <- solve(D+K)
y <- train$Profit
XX.inv <- solve(crossprod(X,X))
M <- X%*%XX.inv%*%t(X)
mu.beta <- as.vector(crossprod(S,(t(X) %*% M %*% y+K%*%mu)))

a <- 2
b <- 1
a.star <- a+(n+p)/2


beta <- matrix(0,p,nit)
w <- rep(0,nit)

library(MASS)
# start chain
count = 0
beta.new <- rep(0,p)
w.new <- rep(0,nit)
w[1] <- 1
w.new[1] <- w[1]

for(m in 2:nit) {
  beta.new <- mvrnorm(1,mu.beta, S*(1/w[m-1]))
  w.new <- rgamma(1,a.star,b+0.5*crossprod((y-X%*%beta.new),(y-X%*%beta.new))
                  +0.5*crossprod((beta.new-mu),(beta.new-mu)))
  
  beta[,m] <- beta.new
  w[m] <- w.new
  count = count + 1
}

col <- 1:1000

beta.hat <- c(mean(beta[1,-col]),mean(beta[2,-col]),mean(beta[3,-col]))
names(beta.hat) <- c("Intercept","Service","Seating")

beta.hat
#Intercept   Service   Seating 
#52427.763 20536.092  4697.139

#yhat <- beta.hat[1]+beta.hat[2]*train$DinnerService+beta.hat[3]*train$SeatingCapacity

e <- as.numeric((crossprod(y-X%*%beta.hat,y-X%*%beta.hat))/(n-p))
#sqrt(e)=6535.175

res <- y-X%*%beta.hat
hist(res,main="Histogram of residuals",xlab="residuals",freq=FALSE,breaks=20)
hist(res,main="Histogram of residuals",xlab="residuals",freq=FALSE,breaks=20)
curve(dnorm(x, mean(res), sd(res)), col="darkblue", lwd=2, add=TRUE)
hist(train$Profit,main="Histogram of profit",xlab="profit",freq=FALSE,breaks=30)
