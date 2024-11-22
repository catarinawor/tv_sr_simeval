#=============================================================
#Run samSim simulations 
#using the harrison chinook data as example
#samSim updates
#Catarina Wor
#=============================================================

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)
#

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
source("R/utils.R")

## Load relevant input data
# Simulation run parameters describing different scenarios


#Run and save simulated data
for(u in 1:6){
  if(u==1){
    print("base")
    simPars <- read.csv("data/generic/SimPars.csv")
    cuPar <- read.csv("data/generic/CUPars.csv")
  }else if(u==2){
    print("base ER")
    simPars <- read.csv("data/genericER/SimPars_ER.csv")
    cuPar <- read.csv("data/genericER/CUPars.csv")
  }else if(u==3){
    print("sensitivity a")
    simPars <- read.csv("data/sensitivity/SimPars.csv")
    cuPar <- read.csv("data/sensitivity/CUPars.csv")   
  }else if(u==4){
    print("sensitivity Smax")
    simPars <- read.csv("data/Smax_sensitivity/SimPars.csv")
    cuPar <- read.csv("data/Smax_sensitivity/CUPars.csv")
  }else if(u==5){
    print("sigma low")
    simPars <- read.csv("data/sigmalow_sensitivity/SimPars.csv")  
    cuPar <- read.csv("data/sigmalow_sensitivity/CUPars_lowsigma.csv")
  }else if(u==6){
    print("sigma med")
    simPars <- read.csv("data/sigmamed_sensitivity/SimPars.csv")  
    cuPar <- read.csv("data/sigmamed_sensitivity/CUPars_medsigma.csv")
  }

  scenNames <- unique(simPars$scenario)
  for(a in seq_len(nrow(simPars))){ 
    print(a)
    genericRecoverySim(simPar=simPars[a,], 
                        cuPar=cuPar, 
                        catchDat=NULL, 
                        srDat=NULL,
                        variableCU=FALSE, 
                        ricPars=NULL, 
                        larkPars=NULL, 
                        cuCustomCorrMat= NULL,
                        outDir="outs", 
                        nTrials=1000, 
                        makeSubDirs=TRUE, 
                        random=FALSE, 
                        uniqueProd=TRUE,
                        uniqueSurv=FALSE)
  
  }

}




