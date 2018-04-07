library(readr)
library(rstan)
tea <- tea_discipline_oss
tea <- tea[tea$ACTIONS > 0, ]

tea.data <- list(N = nrow(tea),intercept = rep(1, nrow(tea)),X = as.numeric(tea$GRADE), y = tea$ACTIONS)

setwd("/Users/rimlisengupta/Dropbox/UT/Spring-2018/SDS383D-StatMod2/Section3")
fit <- stan(file ='poisson.stan', data = tea.data, chains = 4, iter = 1000)

names(fit) <- c("Intercept", "Grade", "lp")
print(fit)

#4 chains, Number of iterations=1000; warmup=500; thin=1; 
traceplot(fit)
plot(fit)

