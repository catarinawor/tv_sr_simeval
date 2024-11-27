#======================================================
#ER scenarios case full confusion matrices
#Catarina Wor
#November 2022
#======================================================

library(gridExtra)
library(ggplot2)
library(ggpubr)
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
simPar <- read.csv("data/genericER/SimPars_ER.csv")

## Store relevant object names to help run simulation 
res<-readRDS(file = "outs/simest/genericER/res_er.rds")
res<-res[res$convergence==0,]

aic=subset(res,parameter=='AIC'&method=='MLE')
bic=subset(res,parameter=='BIC'&method=='MLE')

unique(aic$scenario)

scn<-factor(unique(aic$scenario), levels=c(
 "lowERLowError",
 "lowERHighError",               
 "highERLowError",
 "highERHighError",              
 "ShiftERLowError",
 "ShiftERHighError",             
"trendERLowError",
"trendERHighError",             
"decLinearProd_ShiftERLowError",
"decLinearProd_highERLowError",
 "incLinearProd_ShiftERLowError",
"incLinearProd_highERLowError", 
"decLinearProd_lowERLowError",
"incLinearProd_lowERLowError"  ) )

EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab")

##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA

cn1<-list()
cn2<-list()

aic_set<-list()
bic_set<-list()
o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
  #head(aic_set[[a]])
  sc1=apply(aic_set[[a]],1,which.min)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
 
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)

}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
"lowERLowError"="stationary",
"lowERHighError"="stationary", 
"highERLowError"="stationary",
"highERHighError"="stationary",
"ShiftERLowError"="stationary",
"ShiftERHighError"="stationary",
 "trendERLowError"="stationary",
 "trendERHighError"="stationary",
"decLinearProd_ShiftERLowError"="dynamic.a", 
"decLinearProd_highERLowError"="dynamic.a", 
 "incLinearProd_ShiftERLowError"="dynamic.a", 
"incLinearProd_highERLowError"="dynamic.a",  
"decLinearProd_lowERLowError"="dynamic.a", 
"incLinearProd_lowERLowError"="dynamic.a") 


conf_matrix$eqem_om <- factor(conf_matrix$eqem_om,
 levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab"))
    
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM




conf_matrix$ERtrend<-case_match(conf_matrix$OM,"decLinearProd_highERLowError"~"highER",
                                   "decLinearProd_ShiftERLowError"~"ShiftER", 
                                   "decLinearProd_lowERLowError"~"lowER",
                                   "incLinearProd_highERLowError"~"highER",  
                                   "incLinearProd_ShiftERLowError"~"ShiftER",
                                   "incLinearProd_lowERLowError"~"lowER",  
                                   "highERHighError"~"highER",               
                                   "highERLowError"~"highER",               
                                   "lowERHighError"~"lowER",                
                                   "lowERLowError"~"lowER",                
                                    "ShiftERHighError"~"ShiftER",              
                                    "ShiftERLowError"~"ShiftER",               
                                    "trendERHighError"~"trendER",              
                                    "trendERLowError"~"trendER")

conf_matrix$ERerror<-case_match(conf_matrix$OM,"decLinearProd_highERLowError"~"LowError",
                                   "decLinearProd_ShiftERLowError"~"LowError", 
                                   "incLinearProd_highERLowError"~"LowError",  
                                   "incLinearProd_ShiftERLowError"~"LowError", 
                                   "decLinearProd_lowERLowError"~"LowError",
                                   "incLinearProd_lowERLowError"~"LowError",
                                   "highERHighError"~"HighError",               
                                   "highERLowError"~"LowError",               
                                   "lowERHighError"~"HighError",                
                                   "lowERLowError"~"LowError",                
                                    "ShiftERHighError"~"HighError",              
                                    "ShiftERLowError"~"LowError",               
                                    "trendERHighError"~"HighError",              
                                    "trendERLowError"~"LowError")

conf_matrix$dynamics<-case_match(conf_matrix$OM,"decLinearProd_highERLowError"~"decLinear",
                                   "decLinearProd_ShiftERLowError"~"decLinear", 
                                   "decLinearProd_lowERLowError"~"decLinear",
                                   "incLinearProd_highERLowError"~"incLinear",  
                                   "incLinearProd_ShiftERLowError"~"incLinear", 
                                   "incLinearProd_lowERLowError"~"incLinear", 
                                   "highERHighError"~"stationary",               
                                   "highERLowError"~"stationary",               
                                   "lowERHighError"~"stationary",                
                                   "lowERLowError"~"stationary",                
                                    "ShiftERHighError"~"stationary",              
                                    "ShiftERLowError"~"stationary",               
                                    "trendERHighError"~"stationary",              
                                    "trendERLowError"~"stationary")


conf_matrix$fullname <- apply( conf_matrix[ , c("ERtrend", "ERerror") ] , 1 , paste , collapse = " " )



conf_matrix$fullname<-factor(conf_matrix$fullname, levels=c(
  "lowER LowError",
  "ShiftER LowError",
  "trendER LowError",
  "highER LowError",
  "lowER HighError",
  "ShiftER HighError",
  "trendER HighError",
  "highER HighError" ) )

conf_matrix_stat<-conf_matrix[conf_matrix$dynamics=="stationary",]

paic1=ggplot(data = conf_matrix_stat, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_stat, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic1

ggsave("figures/AIC_MLE_er.png", 
  plot=paic1, width = 10, height = 8)

conf_matrix_logadyn<-conf_matrix[conf_matrix$dynamics!="stationary",]

conf_matrix_logadyn$OM<-factor(conf_matrix_logadyn$OM,levels=c(
  "decLinearProd_lowERLowError",
  "decLinearProd_ShiftERLowError",
  "decLinearProd_highERLowError", 
  "incLinearProd_lowERLowError", 
  "incLinearProd_ShiftERLowError",
  "incLinearProd_highERLowError" 
))

paic2=ggplot(data = conf_matrix_logadyn, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_logadyn, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic2

ggsave("figures/AIC_MLE_logatrend_er.png", 
  plot=paic2, width = 10, height = 8)



pbic=ggplot(data = conf_matrix_stat, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_stat, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

ggsave("figures/BIC_MLE_er.png",
  plot=pbic, width = 10, height = 8)


pbic2=ggplot(data = conf_matrix_logadyn, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_logadyn, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic2

ggsave("figures/BIC_MLE_logatrend_er.png", 
  plot=pbic2, width = 10, height = 8)





#---------------------------------------------------------------------------------------------
#stan


#===================================================================================
reslfo<-readRDS(file = "outs/simest/genericER/resstanloo_baseER.rds")


conf_matrix_lfo<-expand.grid(EM=EM,OM=scn)
conf_matrix_lfo$LFOmcmc=NA


cn3<-list()
lfomcmc_set<-list()

o=0
for(a in seq_along(scn)){


  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_lfo$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix_lfo$eqem_om <- dplyr::recode(conf_matrix_lfo$OM, 
      "lowERLowError"="stationary",
"lowERHighError"="stationary", 
"highERLowError"="stationary",
"highERHighError"="stationary",
"ShiftERLowError"="stationary",
"ShiftERHighError"="stationary",
 "trendERLowError"="stationary",
 "trendERHighError"="stationary",
"decLinearProd_ShiftERLowError"="dynamic.a", 
"decLinearProd_highERLowError"="dynamic.a", 
 "incLinearProd_ShiftERLowError"="dynamic.a", 
"incLinearProd_highERLowError"="dynamic.a",  
"decLinearProd_lowERLowError"="dynamic.a", 
"incLinearProd_lowERLowError"="dynamic.a"
      )   

conf_matrix_lfo$eqem_om <- factor(conf_matrix_lfo$eqem_om,
 levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab"))
    
conf_matrix_lfo$diag<-conf_matrix_lfo$eqem_om==conf_matrix_lfo$EM



conf_matrix_lfo$ERtrend<-case_match(conf_matrix_lfo$OM,"decLinearProd_highERLowError"~"highER",
                                   "decLinearProd_ShiftERLowError"~"ShiftER", 
                                   "decLinearProd_lowERLowError"~"lowER",
                                   "incLinearProd_highERLowError"~"highER",  
                                   "incLinearProd_ShiftERLowError"~"ShiftER",
                                   "incLinearProd_lowERLowError"~"lowER",  
                                   "highERHighError"~"highER",               
                                   "highERLowError"~"highER",               
                                   "lowERHighError"~"lowER",                
                                   "lowERLowError"~"lowER",                
                                    "ShiftERHighError"~"ShiftER",              
                                    "ShiftERLowError"~"ShiftER",               
                                    "trendERHighError"~"trendER",              
                                    "trendERLowError"~"trendER")

conf_matrix_lfo$ERerror<-case_match(conf_matrix_lfo$OM,"decLinearProd_highERLowError"~"LowError",
                                   "decLinearProd_ShiftERLowError"~"LowError", 
                                   "incLinearProd_highERLowError"~"LowError",  
                                   "incLinearProd_ShiftERLowError"~"LowError", 
                                   "decLinearProd_lowERLowError"~"LowError",
                                   "incLinearProd_lowERLowError"~"LowError",
                                   "highERHighError"~"HighError",               
                                   "highERLowError"~"LowError",               
                                   "lowERHighError"~"HighError",                
                                   "lowERLowError"~"LowError",                
                                    "ShiftERHighError"~"HighError",              
                                    "ShiftERLowError"~"LowError",               
                                    "trendERHighError"~"HighError",              
                                    "trendERLowError"~"LowError")

conf_matrix_lfo$dynamics<-case_match(conf_matrix_lfo$OM,"decLinearProd_highERLowError"~"decLinear",
                                   "decLinearProd_ShiftERLowError"~"decLinear", 
                                   "decLinearProd_lowERLowError"~"decLinear",
                                   "incLinearProd_highERLowError"~"incLinear",  
                                   "incLinearProd_ShiftERLowError"~"incLinear", 
                                   "incLinearProd_lowERLowError"~"incLinear", 
                                   "highERHighError"~"stationary",               
                                   "highERLowError"~"stationary",               
                                   "lowERHighError"~"stationary",                
                                   "lowERLowError"~"stationary",                
                                    "ShiftERHighError"~"stationary",              
                                    "ShiftERLowError"~"stationary",               
                                    "trendERHighError"~"stationary",              
                                    "trendERLowError"~"stationary")


conf_matrix_lfo$fullname <- apply( conf_matrix_lfo[ , c("ERtrend", "ERerror") ],
                     1 , paste , collapse = " " )

conf_matrix_lfo_stat<-conf_matrix_lfo[conf_matrix_lfo$dynamics=="stationary",]

conf_matrix_lfo_stat$fullname<-factor(conf_matrix_lfo_stat$fullname, levels=c(
  "lowER LowError",
  "ShiftER LowError",
  "trendER LowError",
  "highER LowError",
  "lowER HighError", 
  "ShiftER HighError",  
  "trendER HighError",  
 "highER HighError"    
))



pmclfo=ggplot(data =  conf_matrix_lfo_stat, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO")+
  scale_fill_gradient(low="white", high="#009194")  +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix_lfo_stat, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo

ggsave("figures/LFO_MCMC_er.png", 
  plot=pmclfo, width = 10, height = 8)



conf_matrix_lfo_logadyn<-conf_matrix_lfo[conf_matrix_lfo$dynamics!="stationary",]

conf_matrix_lfo_logadyn$OM<-factor(conf_matrix_lfo_logadyn$OM,levels=c(
  "decLinearProd_lowERLowError",
  "decLinearProd_ShiftERLowError",
  "decLinearProd_highERLowError", 
  "incLinearProd_lowERLowError", 
  "incLinearProd_ShiftERLowError",
  "incLinearProd_highERLowError" 
))



pmclfo2=ggplot(data =  conf_matrix_lfo_logadyn, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO")+
  scale_fill_gradient(low="white", high="#009194")  +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix_lfo_logadyn, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo2
ggsave("figures/LFO_MCMC_logatrend_er.png", 
  plot=pmclfo, width = 10, height = 8)



