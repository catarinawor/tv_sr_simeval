#======================================================
#sensitivity sigma=0.3 scenarios case full confusion matrices
#Catarina Wor
#November 2022
#============================================



library(gridExtra)
library(ggplot2)
library(dplyr)
library(cowplot)

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

res<-readRDS(file = "outs/simest/sigmamed_sensitivity/res_sigmed.rds")
res<-res[res$convergence==0,]

res$sigscenario<-"medium"
res$scenario<-case_match(
res$scenario,
 "sigmamed_stationary"~"stationary",
 "sigmamed_decLinearProd"~"decLinearProd",        
 "sigmamed_regimeProd"~ "regimeProd",            
 "sigmamed_sineProd"~ "sineProd",             
 "sigmamed_regimeCap"~ "regimeCap",             
 "sigmamed_decLinearCap"~ "decLinearCap",        
 "sigmamed_regimeProdCap"~ "regimeProdCap",
 "sigmamed_shiftCap"~"shiftCap",
 "sigmamed_decLinearProdshiftCap"~"decLinearProdshiftCap")

aic=subset(res,parameter=='AIC'&method=='MLE')
bic=subset(res,parameter=='BIC'&method=='MLE')


scn<-factor(unique(aic$scenario), levels=c(
   "stationary",
   "decLinearProd",
   "sineProd",
   "regimeProd",
   "decLinearCap",         
   "regimeCap",
   "shiftCap",
   "decLinearProdshiftCap",
   "regimeProdCap"
 ) )


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
  aic_set[[a]]=aic_set[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(aic_set[[a]])
  sc1=apply(aic_set[[a]],1,which.min)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
 
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim

  

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  
  o=max(myseq)

}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )  
conf_matrix$eqem_om<-factor( conf_matrix$eqem_om,levels=c( "stationary",
 "autocorr", "dynamic.a", "regime.a", "dynamic.b", "regime.b", "dynamic.ab", "regime.ab"
  ))      

conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$scentype<-dplyr::case_match(as.character(conf_matrix$OM), 
      "stationary"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "regimeProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
      "decLinearCap"~"Smax",
      "regimeCap"~"Smax",
      "shiftCap"~"Smax", 
      "regimeProdCap"~"both",
      "decLinearProdshiftCap"~"both"
      )   

conf_matrix$scendesc<-dplyr::case_match(as.character(conf_matrix$OM), 
      "stationary"~"",
      "decLinearProd"~"linear decline",
      "sineProd"~"sine trend",
      "regimeProd"~"shift up + down",
      "shiftProd"~"shift down + up",
      "decLinearCap"~"linear decline",
      "regimeCap"~"shift up + down",
      "shiftCap"~"shift down", 
      "regimeProdCap"~"regime both",
      "decLinearProdshiftCap"~"trend & shift"
      )   
 
conf_matrix$fullname <- apply( conf_matrix[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )


conf_matrix$fullname<-factor(conf_matrix$fullname, levels=c(
  " stationary",
  "shift up sigma",
  " autocorrelation",           
  "linear decline log(alpha)",
  "sine trend log(alpha)",  
  "shift up + down log(alpha)", 
  "shift down + up log(alpha)",       
  "linear decline Smax",
  "shift up + down Smax",              
  "shift down Smax",           
  "regime both both",           
 "trend & shift both"         
 ) )

paic_sigmed=ggplot(data =  conf_matrix, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle(expression(paste("AIC", ~sigma, "=0.3")))+
  scale_fill_gradient(low="white", high="#009194")  +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic_sigmed

ggsave("figures/AIC_MLE_sigmed.png",
 plot=paic_sigmed, width = 9,height = 7)



pbic_sigmed=ggplot(data =  conf_matrix, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle(expression(paste("BIC", ~sigma, "=0.3"))) +
  scale_fill_gradient(low="white", high="#009194") +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic_sigmed

ggsave("figures/BIC_MLE_sigmed.png", 
  plot=pbic_sigmed, width = 9,height = 7)






#===============================================================================
 



#---------------------------------------------------------------------------------------------
#stan



#===================================================================================
reslfo<-readRDS(file = "outs/simest/sigmamed_sensitivity/resstanloo_sigmed.rds")


reslfo$scenario<-case_match(
 reslfo$scenario,
 "sigmamed_stationary"~"stationary",
 "sigmamed_decLinearProd"~"decLinearProd",        
 "sigmamed_regimeProd"~ "regimeProd",            
 "sigmamed_sineProd"~ "sineProd",             
 "sigmamed_regimeCap"~ "regimeCap",             
 "sigmamed_decLinearCap"~ "decLinearCap",        
 "sigmamed_regimeProdCap"~ "regimeProdCap",
 "sigmamed_shiftCap"~"shiftCap",             
 "sigmamed_decLinearProdshiftCap"~"decLinearProdshiftCap")

scn<-factor(unique(reslfo$scenario), levels=c(
  "stationary", 
  "decLinearProd",  
  "sineProd",
  "regimeProd",
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  "decLinearProdshiftCap",
  "regimeProdCap"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab")
##Confusion matrices
conf_matrixsigmedlfo<-expand.grid(EM=EM,OM=scn)
conf_matrixsigmedlfo$LFOmcmc=NA

cn3<-list()
lfomcmc_set<-list()

o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(27,8,18,9,24,15,21,12)] #reorder estimation models
  head(lfomcmc_set[[a]])

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrixsigmedlfo$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrixsigmedlfo$eqem_om <- dplyr::recode(conf_matrixsigmedlfo$OM, 
      "stationary"="stationary",
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   


conf_matrixsigmedlfo$eqem_om<-factor(conf_matrixsigmedlfo$eqem_om,levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab"))  
conf_matrixsigmedlfo$diag<-conf_matrixsigmedlfo$eqem_om==conf_matrixsigmedlfo$EM

conf_matrixsigmedlfo$scentype<-dplyr::case_match(as.character(conf_matrixsigmedlfo$OM), 
      "stationary"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "regimeProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
      "decLinearCap"~"Smax",
      "regimeCap"~"Smax",
      "shiftCap"~"Smax", 
      "regimeProdCap"~"both",
      "decLinearProdshiftCap"~"both"
      )   

conf_matrixsigmedlfo$scendesc<-dplyr::case_match(as.character(conf_matrixsigmedlfo$OM), 
      "stationary"~"",
      "decLinearProd"~"linear decline",
      "sineProd"~"sine trend",
      "regimeProd"~"shift up + down",
      "shiftProd"~"shift down + up",
      "decLinearCap"~"linear decline",
      "regimeCap"~"shift up + down",
      "shiftCap"~"shift down", 
      "regimeProdCap"~"regime both",
      "decLinearProdshiftCap"~"trend & shift"
      )   
 
conf_matrixsigmedlfo$fullname <- apply( conf_matrixsigmedlfo[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )


conf_matrixsigmedlfo$fullname<-factor(conf_matrixsigmedlfo$fullname, levels=c(
  " stationary",
  "shift up sigma",
  " autocorrelation",           
  "linear decline log(alpha)",
  "sine trend log(alpha)",  
  "shift up + down log(alpha)", 
  "shift down + up log(alpha)",       
  "linear decline Smax",
  "shift up + down Smax",              
  "shift down Smax",           
  "regime both both",           
 "trend & shift both"         
 ) )


pmclfo=ggplot(data =  conf_matrixsigmedlfo, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle(expression(paste("LFO", ~sigma, "=0.3")))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrixsigmedlfo, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo
ggsave("figures/LFO_MCMC_sigmed.png",
 plot=pmclfo, width = 9,height = 7)


