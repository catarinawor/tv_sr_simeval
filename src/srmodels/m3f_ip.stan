data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  array[N] int ii; //index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
}
transformed data{
real logbeta_pr;
real logbeta_pr_sig;

logbeta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
logbeta_pr=log(1/pSmax_mean)-0.5*logbeta_pr_sig*logbeta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction

}
parameters{
  real<lower = 0>  log_a0;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a; //random walk of log_a
  }
  
}  
model{
  //priors
  log_a0 ~ normal(1.5,2.5); //initial productivity - wide prior
  log_b ~ normal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter - wide prior
  a_dev ~ std_normal(); //standardized (z-scales) deviances

  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_a ~ normal(0,1); //half normal on variance (lower limit of zero)
   
 
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma); 
  
}
 generated quantities{
     vector[N] log_lik;
     real S_max;
     vector[L] U_msy;
     vector[L] S_msy;
     
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[ii[n]] - S[n]*b, sigma);
   
    S_max = 1/b;
    U_msy = 1-lambert_w0(exp(1-log_a));
    S_msy = (1-lambert_w0(exp(1-log_a)))/b;
    }
