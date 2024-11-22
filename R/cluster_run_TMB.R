#=================================================
# Code used to run simulations-estimation on the 
# GCHPC
#=================================================
#sbase case scenarios

library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/generic/SimPars.csv")

#run example  
tst1<-tmb_func(path=".",
  a=3,
  u=696)

pars<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb <- slurm_apply(tmb_func, pars, jobname = 'TMBrun',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



res <- get_slurm_out(sjobtmb, outtype = 'table', wait = TRUE)


#AFTER JOB IS DONE IMPORT  the results

saveRDS(res[res$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "resbase1.rds")
saveRDS(res[res$scenario%in%simPars$scenario[floor(nrow(simPars)/2+1):nrow(simPars)],], file = "resbase2.rds")
saveRDS(res, file = "resbase.rds")



#============================================================================
#sensitivity a scenarios
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sensitivity/SimPars.csv")

pars_a<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_a <- slurm_apply(tmb_func, pars_a, jobname = 'TMBrun_a',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))



#AFTER JOB IS DONE IMPORT  the results
res_a <- get_slurm_out(sjobtmb_a, outtype = 'table', wait = TRUE)



saveRDS(res_a[res_a$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_a1.rds")
saveRDS(res_a[res_a$scenario%in%simPars$scenario[(floor(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_a2.rds")
saveRDS(res_a, file = "res_a.rds")



#============================================================================
#smax scenarios 
library(rslurm)
library(samEst)
source("R/tmb_func.R")

simPars <- read.csv("data/Smax_sensitivity/SimPars.csv")

pars_smax<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_smax <- slurm_apply(tmb_func, pars_smax, jobname = 'TMBrun_smax',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))


res_smax <- get_slurm_out(sjobtmb_smax, outtype = 'table', wait = TRUE)


saveRDS(res_smax[res_smax$scenario%in%simPars$scenario[seq_len(nrow(simPars)/2)],], file = "res_smax1.rds")
saveRDS(res_smax[res_smax$scenario%in%simPars$scenario[floor(nrow(simPars)/2+1):nrow(simPars)],], file = "res_smax2.rds")
saveRDS(res_smax, file = "res_smax.rds")



#============================================================================

#============================================================================
#sensitivity sigma scenarios low
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

pars_siglow<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_siglow <- slurm_apply(tmb_func, pars_siglow, jobname = 'TMBrun_siglow',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))


res_siglow <- get_slurm_out(sjobtmb_siglow, outtype = 'table', wait = TRUE)



saveRDS(res_siglow[res_siglow$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_siglow1.rds")
saveRDS(res_siglow[res_siglow$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_siglow2.rds")
saveRDS(res_siglow, file = "res_siglow.rds")



#============================================================================
#sensitivity sigma scenarios med
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

pars_sigmed<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_sigmed <- slurm_apply(tmb_func, pars_sigmed, jobname = 'TMBrun_sigmed',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_sigmed <- get_slurm_out(sjobtmb_sigmed, outtype = 'table', wait = TRUE)


saveRDS(res_sigmed[res_sigmed$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_sigmed1.rds")
saveRDS(res_sigmed[res_sigmed$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_sigmed2.rds")
saveRDS(res_sigmed, file = "res_sigmed.rds")



#============================================================================
#ER scenarios
library(rslurm)
library(samEst)
source("R/tmb_func.R")


simPars <- read.csv("data/genericER/SimPars_ER.csv")

pars_ER<-data.frame(path="..",
  a=rep(seq_len(nrow(simPars)),each=1000),
  u=1:1000)


sjobtmb_er <- slurm_apply(tmb_func, pars_ER, jobname = 'TMBrun_ER',
                    nodes = 300, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst"),
                    rscript_path = "/home/Documents/pfmln/results/cluster-tvsimest",
                    libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.1",
                    global_objects=c("simPars"))




res_er <- get_slurm_out(sjobtmb_er, outtype = 'table', wait = TRUE)



saveRDS(res_er[res_er$scenario%in%simPars$scenario[seq_len(round(nrow(simPars)/2))],], file = "res_er1.rds")
saveRDS(res_er[res_er$scenario%in%simPars$scenario[(round(nrow(simPars)/2)+1):nrow(simPars)],], file = "res_er2.rds")
saveRDS(res_er, file = "res_er.rds")
