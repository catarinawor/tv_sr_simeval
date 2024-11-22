#===============================================
#Ricker curve scenario plots
# april 2023
#===============================================
library(ggplot2)
library(dplyr)
library(samEst)
library(cowplot)
library(ggpubr)


mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)


#----------------------------------------
#base case                              |
#----------------------------------------
simPar <- read.csv("data/generic/SimPars.csv")
 

simData<-list()
actualSR<-list()
alldat<-list()

for(a in seq_len(nrow(simPar))){

  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                          simPar$nameOM[a],"/",
                          simPar$scenario[a],"/",
                          paste(simPar$nameOM[a],"_", 
                          simPar$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout


  dat<-simData[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar$scenario[a]
  alldat[[a]]<-dat

  S <- seq(0,600000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar$scenario[a])

}

SRdf<-do.call(rbind,actualSR)
datdf<-do.call(rbind,alldat)


SRdf$scenario_f <-factor(SRdf$scenario, levels=c("stationary",  
                                                  "autocorr",
                                                  "sigmaShift",
                                                  "decLinearProd", 
                                                  "regimeProd", 
                                                  "sineProd",  
                                                  "shiftProd",
                                                  "regimeCap", 
                                                  "decLinearCap", 
                                                  "shiftCap", 
                                                  "regimeProdCap", 
                                                  "decLinearProdshiftCap"))   
                                                                  

SRdf$scencode<-dplyr::case_match(SRdf$scenario, 
      "stationary"~"Base1",
      "autocorr"~"Base2",
      "sigmaShift"~"Base3", 
      "decLinearProd"~"Base4",
      "sineProd"~"Base5",
      "regimeProd"~"Base6",
      "shiftProd"~"Base7",
      "decLinearCap"~"Base8",
      "regimeCap"~"Base9",
      "shiftCap"~"Base10", 
      "regimeProdCap"~"Base11",
      "decLinearProdshiftCap"~"Base12"
      )   

SRdf$scencode <-factor(SRdf$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))

datdf$scenario_f <-factor(datdf$scenario, levels=c("stationary","autocorr","sigmaShift",
              "decLinearProd", "regimeProd", "sineProd","shiftProd",
               "regimeCap", "decLinearCap", "shiftCap",
               "regimeProdCap",  "decLinearProdshiftCap"))

datdf$scencode<-dplyr::case_match(datdf$scenario, 
      "stationary"~"Base1",
      "autocorr"~"Base2",
      "sigmaShift"~"Base3", 
      "decLinearProd"~"Base4",
      "sineProd"~"Base5",
      "regimeProd"~"Base6",
      "shiftProd"~"Base7",
      "decLinearCap"~"Base8",
      "regimeCap"~"Base9",
      "shiftCap"~"Base10", 
      "regimeProdCap"~"Base11",
      "decLinearProdshiftCap"~"Base12"
      )   

datdf$scencode <-factor(datdf$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))

SRexample<-  ggplot(SRdf) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scencode)
SRexample    
 
ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/srcurve_basecase.png",
      plot = SRexample, 
      width = 12, height = 6
    )


#select scenarios for main paper
datdfselec<-datdf[datdf$scenario%in%c("autocorr","decLinearProd","sineProd","shiftProd",
    "decLinearCap","shiftCap","regimeProdCap","decLinearProdshiftCap"),]



datdfselec$genericscenario<-dplyr::case_match(datdfselec$scenario, 
     "autocorr"~"stationary",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~fluctuation",
      "shiftProd"~"shift~down+up",
    "decLinearCap"~"linear~decline",
    "shiftCap"~"shift~down",
    "regimeProdCap"~ "shift~-~both",
    "decLinearProdshiftCap"~"mixed~trend"
      )   

datdfselec$paramvary<-dplyr::case_match(datdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
    "decLinearCap"~"S[max]",
    "shiftCap"~"S[max]",
    "regimeProdCap"~ "both",
    "decLinearProdshiftCap"~"both"
      ) 

datdfselec$paramvary<-factor(datdfselec$paramvary, 
    levels=c("stationary","log(alpha)","S[max]","both"))


SRdfselec<-SRdf[SRdf$scenario%in%c("autocorr","decLinearProd","sineProd","shiftProd",
    "decLinearCap","shiftCap","regimeProdCap","decLinearProdshiftCap"),]

SRdfselec$genericscenario<-dplyr::case_match(SRdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~fluctuation",
      "shiftProd"~"shift~down+up",
    "decLinearCap"~"linear~decline",
    "shiftCap"~"shift~down",
    "regimeProdCap"~ "shift~-~both",
    "decLinearProdshiftCap"~"mixed~trend"
      )   

SRdfselec$paramvary<-dplyr::case_match(SRdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
    "decLinearCap"~"S[max]",
    "shiftCap"~"S[max]",
    "regimeProdCap"~ "both",
    "decLinearProdshiftCap"~"both"
      ) 

SRdfselec$paramvary<-factor(SRdfselec$paramvary, 
    levels=c("stationary","log(alpha)","S[max]","both"))


SRexampleselec<-  ggplot(SRdfselec) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    coord_cartesian(ylim = c(-0, 600000))+
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdfselec,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(paramvary~genericscenario, ncol=8,  labeller=label_parsed)+ 
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
SRexampleselec

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/srcurve_selecmain.png",
      plot = SRexampleselec, 
      width = 8, height = 4
    )

paramdf<-reshape2::melt(datdf,id.vars=c("iteration","year","CU","spawners","recruits","obsSpawners","obsRecruits",
                                "ER", "obsER", "targetER",   
                                "sMSY", "sGen", "uMSY", "scenario", "scenario_f","scencode"))

paramdf$simulated<-dplyr::recode(paramdf$scenario, 
      "stationary"="simple",
      "autocorr"="simple",
      "sigmaShift"="simple", 
      "decLinearProd"="rw",
      "sineProd"="rw",
      "regimeProd"="hmm",
      "decLinearCap"="rw",
      "regimeCap"="hmm",
      "shiftCap"="hmm", 
      "shiftProd"="hmm",
      "regimeProdCap"="hmm",
      "decLinearProdshiftCap"="rw"
      )   


paramdf$paramch<-dplyr::recode(paramdf$scenario, 
      "stationary"="simple",
      "autocorr"="simple",
      "sigmaShift"="simple", 
      "decLinearProd"="tva",
      "sineProd"="tva",
      "regimeProd"="tva",
      "shiftProd"="tva",
      "decLinearCap"="tvb",
      "regimeCap"="tvb",
      "shiftCap"="tvb", 
      "regimeProdCap"="both",
      "decLinearProdshiftCap"="both"
      )   


paramdf<-paramdf[paramdf$variable!="beta",]
 
paramdf$scencode<-dplyr::case_match(paramdf$scenario, 
      "stationary"~"Base1",
      "autocorr"~"Base2",
      "sigmaShift"~"Base3", 
      "decLinearProd"~"Base4",
      "sineProd"~"Base5",
      "regimeProd"~"Base6",
      "shiftProd"~"Base7",
      "decLinearCap"~"Base8",
      "regimeCap"~"Base9",
      "shiftCap"~"Base10", 
      "regimeProdCap"~"Base11",
      "decLinearProdshiftCap"~"Base12"
      )   



paramdf$scenario <-factor(paramdf$scenario, levels=c("stationary",
      "autocorr",
      "sigmaShift", 
      "decLinearProd",
      "sineProd",
      "regimeProd",
      "shiftProd",
      "decLinearCap",
      "regimeCap",
      "shiftCap", 
      "regimeProdCap",
      "decLinearProdshiftCap"))

paramdf$scencode <-factor(paramdf$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))

paramtraj<-ggplot(paramdf) +
    geom_line(aes(x=year,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scencode, scales="free")
paramtraj

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/params_basecase.png",
      plot = paramtraj, 
      width = 12, height = 6
    )

paramtraj2<-ggplot(paramdf) +
    geom_line(aes(x=year,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scenario_f, scales="free")
paramtraj2

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/params_basecase_scnname.png",
      plot = paramtraj2, 
      width = 19, height = 6
    )

#selected scenarios


paramdfselec <- paramdf[paramdf$scenario%in%c("autocorr","decLinearProd","sineProd","shiftProd",
    "decLinearCap","shiftCap","regimeProdCap","decLinearProdshiftCap")&
     paramdf$variable!="sigma",]

paramdfselec$variable<-as.character(paramdfselec$variable)
paramdfselec$variable[paramdfselec$variable=="capacity"]<-paste("S[max]")
paramdfselec$variable[paramdfselec$variable=="alpha"]<-paste("log(alpha)")

unique(paramdfselec$variable)

paramdfselec$genericscenario<-dplyr::case_match(paramdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~fluctuation",
      "shiftProd"~"shift~down+up",
    "decLinearCap"~"linear~decline",
    "shiftCap"~"shift~down",
    "regimeProdCap"~ "shift~-~both",
    "decLinearProdshiftCap"~"mixed~trend"
      )   

paramdfselec$paramvary<-dplyr::case_match(paramdfselec$scenario, 
      "autocorr"~".",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
    "decLinearCap"~"S[max]",
    "shiftCap"~"S[max]",
    "regimeProdCap"~ "both",
    "decLinearProdshiftCap"~"both"
      ) 

paramdfselec$paramvary<-factor(paramdfselec$paramvary, 
    levels=c(".","log(alpha)","S[max]","both"))

paramtrajselec<-ggplot(paramdfselec) +
    geom_line(aes(x=year-54,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85,option="A") +
    labs(col = "year") +
    ylab("")+
    theme(legend.position="none")+
    facet_grid(variable~paramvary+genericscenario, scales="free", labeller=label_parsed)
paramtrajselec

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/params_selecmain.png",
      plot = SRexampleselec, 
      width = 8, height = 4
    )


multi.page.scenario <- ggarrange( paramtrajselec,SRexampleselec,
                        nrow = 2, ncol = 1,
                        legend="none",
                        heights=c(1,1),
                        align="v")
multi.page.scenario




#use cowplot



logatrajselec<-ggplot(paramdfselec[paramdfselec$variable=="log(alpha)",]) +
    geom_line(aes(x=year-54,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85,option="A") +
    xlab( "year") +
    ylab(expression(log(alpha)))+
    theme(legend.position="none",axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
    facet_grid(~paramvary+genericscenario, scales="free", labeller=label_parsed)
logatrajselec

yrcol<-viridis::viridis(n = 5)

smaxtrajselec<-ggplot(paramdfselec[paramdfselec$variable=="S[max]",]) +
    geom_line(aes(x=year-54,y=value/1000,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85,option="A") +
    xlab( "year") +
    ylab(expression(paste(S[max]~"(1000s)")))+
    theme(legend.position="none",  strip.text.x = element_blank())+
    facet_grid(~paramvary+genericscenario, scales="free")
smaxtrajselec




SRexampleselec<-  ggplot(SRdfselec) +
    geom_line(aes(x=spawners/1000,y=recruits/1000, col=as.factor(year)),linewidth=2) +
    mytheme + 
    coord_cartesian(ylim = c(0, 570))+
    scale_colour_viridis_d(end=.85) +
    xlab("spawners") +
    ylab("recruits") +
    geom_point(data=datdfselec,aes(x=spawners/1000,y=recruits/1000,col=as.factor(year)),alpha=.5) +
    facet_wrap(paramvary~genericscenario, ncol=8)+ 
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1), strip.text.x = element_blank())
SRexampleselec


plot2 <-plot_grid(logatrajselec, smaxtrajselec, SRexampleselec, ncol=1)
plot2

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/multi.page.scenario_basecase.png",
      plot = plot2, 
      width = 12, height = 6
    )

#----------------------------------------
#equilibrium parameter plots                            |
#----------------------------------------


names(datdf)
eqdf<-datdf

eqdf$sMSY<-smsyCalc(eqdf$alpha,eqdf$beta)
eqdf$uMSY<-umsyCalc(eqdf$alpha)
eqdf$sGen<-unlist(mapply(sGenCalc,a=eqdf$alpha,
          Smsy=eqdf$sMSY, 
          b=eqdf$beta))


pareqdf<-reshape2::melt(eqdf,id.vars=c("iteration","year","CU","spawners","recruits","obsSpawners","obsRecruits",
                                "ER", "obsER", "targetER", "scenario", "scenario_f"))



pareqdf<-pareqdf[pareqdf$variable!="sigma"&,]


 ggplot(pareqdf) +
    geom_line(aes(x=year,y=value), linewidth=2)+
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scenario_f, scales="free")
