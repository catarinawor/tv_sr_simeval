#routines to visualize smulation estimation Results
#Catarina Wor
#November 2024
#============================================

library(ggplot2)
library(gridExtra)
library(dplyr)
source("R/utils.R")

mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)
#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

restmb<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow.rds")
resstan<-readRDS(file = "outs/simest/sigmalow_sensitivity/resstan_siglow.rds")

res<-rbind(restmb,resstan)
res$parameter[res$parameter=="Smax"]<-"smax"
res$parameter[res$parameter=="alpha"]<-"logalpha"


resparam<-res[res$parameter%in%c("logalpha","smax","smsy","sgen","umsy","sigma"),]
resparam$convergence[resparam$parameter=="logalpha"&resparam$mode>40]<-1
resparam$convergence[resparam$parameter=="smax"&resparam$mode>1e8]<-1


convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})
convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
convstatMCMC<-convstat[convstat$x==0&convstat$method=="MCMC",]


allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])

convsum<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){length(unique(x))})

conv_iter<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){(unique(x))})

convsnc<-as.numeric(rownames(convsum))

resl<-list()
for(i in seq_along(convsnc)){

    sel<-conv_iter[convsnc[i],]
    resl[[i]]<-resparam %>% filter(model==sel$model&
                            scenario==sel$scenario&
                            iteration%in%sel$x[[1]])
}

resparam<-as.data.frame(data.table::rbindlist(resl))

df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))


df$col<-factor(df$variable,levels=c("median","mode", "sim"))
unique(df$scenario)
df$scenario<-case_match(
df$scenario,
 "sigmalow_stationary"~"stationary",
 "sigmalow_decLinearProd"~"decLinearProd",        
 "sigmalow_regimeProd"~ "regimeProd",            
 "sigmalow_sineProd"~ "sineProd",             
 "sigmalow_regimeCap"~ "regimeCap",             
 "sigmalow_decLinearCap"~ "decLinearCap",        
 "sigmalow_regimeProdCap"~ "regimeProdCap",
 "sigmalow_shiftCap"~"shiftCap",
 "sigmalow_decLinearProdshiftCap"~"decLinearProdshiftCap")

df$scenario<-factor(df$scenario,levels=c("stationary",
                                        "decLinearProd",
                                        "sineProd", 
                                        "regimeProd",                     
                                        "decLinearCap",
                                        "regimeCap",
                                        "shiftCap",                      
                                        "regimeProdCap",
                                        "decLinearProdshiftCap"  ))

df$scentype<-dplyr::case_match(df$scenario, 
      "stationary"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "regimeProd"~"log(alpha)",
      "decLinearCap"~"S[max]",
      "regimeCap"~"S[max]",
      "shiftCap"~"S[max]", 
      "regimeProdCap"~"both",
      "decLinearProdshiftCap"~"both"
      )   
df$scentype<-factor(df$scentype, levels=c("stationary","log(alpha)",
    "S[max]","both" )  ) 



df$scendesc<-dplyr::case_match(df$scenario, 
      "stationary"~".",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~trend",
      "regimeProd"~"shift~up+down",
      "decLinearCap"~"linear~decline",
      "regimeCap"~"shift~up+down",
      "shiftCap"~"shift~down", 
      "regimeProdCap"~"regime~both",
      "decLinearProdshiftCap"~"trend+shift"
      )   


df$scendesc<-factor(df$scendesc, levels=c(".","linear~decline",
    "sine~trend","shift~up+down","shift~down",
    "regime~both","trend+shift") ) 
     


unique(df$scenario)
df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"  ))



#df_alpha_sim<- df[df$parameter=="logalpha"&df$variable=="sim",]

#df_alpha_est<- df[df$parameter=="logalpha"&df$variable=="mode",]

head(df)


summarydf<-aggregate(df$value,by=list(scenario=df$scenario, 
    parameter=df$parameter,
    variable=df$variable,
    method=df$method, 
    model=df$model,
    by=df$by,
    scentype=df$scentype,
    scendesc=df$scendesc ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975),na.rm=TRUE)})
summarydf<-do.call(data.frame, summarydf)

summarydf_alpha<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="mode",]

summarydf_alpha_sim<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="sim",]


plotalpha<-ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by-54,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim,aes(x=by-54,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.4,2.3))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
ggtitle(expression(paste(sigma,"=0.1"))) +
facet_grid(scentype+scendesc~model, scales="free_y",labeller =  label_parsed)#,
plotalpha
ggsave("figures/compareMCMC_MLE_siglow_alpha.png",
    plot=plotalpha, width = 15,height = 14)

#=======================================================
#b estimates

summarydf_smax<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="mode",]

summarydf_smax_sim<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="sim",]

plotsmax<-ggplot() + 
geom_pointrange(data=summarydf_smax,aes(x=by-54,y= x.50./1000,ymin = x.2.5./1000, 
    ymax = x.97.5./1000, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim,aes(x=by-54,y= x.50./1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max]~"(1000s)")) +
xlab("year") +
ggtitle(expression(paste(sigma,"=0.1"))) +
coord_cartesian(ylim = c(60000,400000)/1000)+ 
facet_grid(scentype+scendesc~model, scales="free_y",labeller =  label_parsed)
plotsmax

ggsave("figures/compareMCMC_MLE_siglow_smax.png",
    plot=plotsmax, width = 15,height = 14)

#=======================================================
#smsy estimates


summarydf_smsy<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="mode",]

summarydf_smsy_sim<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="sim",]



plotsmsy<-ggplot() + 
geom_pointrange(data=summarydf_smsy,aes(x=by-54,y= x.50./1000,ymin = x.2.5./1000, 
    ymax = x.97.5./1000, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim,aes(x=by-54,y= x.50./1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY]~"(1000s)")) +
xlab("year") +
ggtitle(expression(paste(sigma,"=0.1"))) +
coord_cartesian(ylim = c(20000,150000)/1000)+ 
facet_grid(scentype+scendesc~model, scales="free_y",labeller =  label_parsed)
plotsmsy

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sig_low/compareMCMC_MLE_siglow_smsy.png",
    plot=plotsmsy, width = 15,height = 14)
ggsave("figures/compareMCMC_MLE_siglow_smsy.png",
    plot=plotsmsy, width = 15,height = 14)

