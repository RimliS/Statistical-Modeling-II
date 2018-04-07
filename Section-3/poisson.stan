data {
  int<lower=0> N;
  int intercept[N];
  int X[N];
  int<lower=0> y[N];
}

parameters {
  real beta[2];
}

model {
  beta[1] ~ normal(0, 1);
  beta[2] ~ normal(0, 1);
  
  
  for (i in 1:N) y[i] ~ poisson(exp(beta[1] + beta[2] * X[i]));  
}
