#======================================================
#Routines to visualize smulation estimation results
#Catarina Wor
#November 2024
#======================================================


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

#=================================================================================================================
#sensitivity smax
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

restmb<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")

resstansmax1<-readRDS(file = "outs/simest/Smax_sensitivity/resstan_smax1.rds")
resstansmax2<-readRDS(file = "outs/simest/Smax_sensitivity/resstan_smax2.rds")

resstan<-rbind(resstansmax1,resstansmax2)

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"
res$parameter[res$parameter=="alpha"]<-"logalpha"

resparam<-res[res$parameter%in%c("logalpha","smax","sigma","smsy","sgen","umsy"),]

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

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))


summarydf<-aggregate(df$value,by=list(scenario=df$scenario,
    parameter=df$parameter,  
    method=df$method, 
    model=df$model,
    by=df$by, 
    variable=df$variable),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975),na.rm=T)})
summarydf<-do.call(data.frame, summarydf)

summarydf$scenario<-factor(summarydf$scenario,levels=c("regimeSmax025",
                                                       "regimeSmax050",
                                                       "regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax025",
                                                       "trendLinearSmax050",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"  ))


summarydf$typechange<-dplyr::case_match(summarydf$scenario, 
      "regimeSmax025"~"regime",
        "regimeSmax050"~"regime",
        "regimeSmax150"~"regime",     
        "regimeSmax200"~"regime",
        "regimeSmax300"~"regime",
        "trendLinearSmax025"~"trend",
        "trendLinearSmax050"~"trend",
        "trendLinearSmax150"~"trend",
        "trendLinearSmax200"~"trend",
        "trendLinearSmax300"~"trend"
      )   


summarydf$magnitude_ch<-0

summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax025")]<-"S[max]~x~0.25"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax050")]<-"S[max]~x~0.5"    
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax150")]<-"S[max]~x~1.5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax200")]<-"S[max]~x~2.0"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax300")]<-"S[max]~x~3.0"   
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax025")]<-"S[max]~x~0.25"   
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax050")]<-"S[max]~x~0.5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax150")]<-"S[max]~x~1.5" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax200")]<-"S[max]~x~2.0" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax300")]<-"S[max]~x~3.0"   

summarydf$magnitude_ch<-factor(summarydf$magnitude_ch,levels=c("S[max]~x~0.25",
"S[max]~x~0.5",
"S[max]~x~1.5", 
"S[max]~x~2.0", 
"S[max]~x~3.0"))

summarydf$model<-factor(summarydf$model, levels=c("simple", "autocorr", "rwa",  "hmma", "rwb", "hmmb", "rwab", "hmmab"))

summarydf_alpha_sim<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="sim",]

summarydf_alpha<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="mode",]


p_sens_smax_alpha<-ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by-54,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim,aes(x=by-54,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.2,2.0))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(typechange+magnitude_ch~model, scales="free_y", labeller =  label_parsed)
p_sens_smax_alpha

ggsave("figures/compareMCMC_MLE_sens_smax_alpha.png",
    plot=p_sens_smax_alpha, width = 15,height = 18 )





#smax
summarydf_smax_sim<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="sim",]

summarydf_smax<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="mode",]


p_sens_smax_smax<-ggplot() + 
geom_pointrange(data=summarydf_smax,aes(x=by-54,y= x.50./1000,ymin = x.2.5./1000,
 ymax = x.97.5./1000, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim,aes(x=by-54,y= x.50./1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,630000)/1000)+ #570000
mytheme+ 
ylab(expression(S[max]~"(1000s)")) +
xlab("year") +
facet_grid(typechange+magnitude_ch~model, scales="free_y", labeller =  label_parsed)
p_sens_smax_smax
ggsave("figures/compareMCMC_MLE_sens_smax_smax.png",
    plot=p_sens_smax_smax, width = 15,height = 18)


#=======================================================================
#smsy


summarydf_smsy_sim<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="sim",]

summarydf_smsy<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="mode",]


p_sens_smax_smsy<-ggplot() + 
geom_pointrange(data=summarydf_smsy,aes(x=by-54,y= x.50./1000,ymin = x.2.5./1000,
 ymax = x.97.5./1000, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim,aes(x=by-54,y= x.50./1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0,330000)/1000)+ #570000
mytheme+ 
ylab(expression(S[MSY]~"(1000s)")) +
xlab("year") +
facet_grid(typechange+magnitude_ch~model, scales="free_y", labeller =  label_parsed)
p_sens_smax_smsy
ggsave("figures/compareMCMC_MLE_sens_smax_smsy.png",
    plot=p_sens_smax_smsy, width = 15,height = 18)



