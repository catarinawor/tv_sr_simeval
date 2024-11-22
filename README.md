# timevar_sr_simeval

This is the code to generate the results presented "Recommendations for estimating and detecting time-varying spawner-recruit dynamics in fish populations". Below are instructions on how to reproduce the results in the paper. 


## Simulation

The main simulation file is: 

- Sim.R

## Estimation

The Estimation routines were run on a slurm cluster. Note that the satn files were run on R 4.1.2, so all the required packages must be comaptible with that version. All the stan and .cpp files are in the src folder. 
 The main R files are:
 - cluster_run_TMB.R
 - cluster_run_stan.R

## Plotting 

R code to generate plots in the main paper:
 - main_plots_bias_prec_base.R
 - main_plots_bias_prec_er.R 
 - main_plots_confmat
 - main_plots_model_selec_lines.R

R code to generate plot in the Supplementary materials
  - sm_plots_compare_mle_mcmc.R
  - sm_plots_compare_mle_mcmc_er.R
  - sm_plots_compare_mle_mcmc_sens_a.R
  - sm_plots_compare_mle_mcmc_smax.R
  - sm_plots_compare_mle_mcmc_sigmed.R
  - sm_plots_compare_mle_mcmc_siglow.R
  - sm_plots_confmat_base.R
  - sm_plots_confmat_er.R
  - sm_plots_confmat_sensa.R
  - sm_plots_confmat_siglow.R
  - sm_plots_confmat_sigmed.R


The simulation-estimation outcomes were not included in this repository due to file size limitations. But they are available upon request. Email: catarina.wor@dfo-mpo.gc.ca or file an issue in this repository. 


