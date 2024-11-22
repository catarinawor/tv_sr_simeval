data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
}
transformed data{
  real logbeta_pr;
  real logbeta_pr_sig;

  logbeta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
  logbeta_pr=log(1/pSmax_mean)-0.5*logbeta_pr_sig*logbeta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction

}
parameters {
  real log_a;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this
    
  //variance components  
  real<lower = 0> sigma;
    
}
transformed parameters{
  real b;
  b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter - wide prior
   
  //variance terms
  sigma ~ normal(0,1.0); //half normal on variance (lower limit of zero)
      
  R_S ~ normal(log_a - S*b, sigma);
}
generated quantities{
  real log_lik_oos;
  log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b, sigma);
}