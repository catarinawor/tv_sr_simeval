#======================================================
#sensitivity log alpha scenarios case full confusion matrices
#Catarina Wor
#November 2024
#============================================


library("ggpubr")
library(gridExtra)
library(ggplot2)
library(cowplot)



mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

#========================================================================================================
#sensitivity a scenario
#read in data

res_a<-readRDS(file = "outs/simest/sensitivity/res_a.rds")

aic_a=subset(res_a, parameter=='AIC'&method=='MLE')
bic_a=subset(res_a, parameter=='BIC'&method=='MLE')


aic_a$mode[aic_a$convergence>0]<-Inf
bic_a$mode[aic_a$convergence>0]<-Inf


scn<-factor(unique(aic_a$scenario), levels=c(
  "trendLinearProd1",
  "trendLinearProd1.3",
  "trendLinearProd2",
  "trendLinearProd5",
  "trendLinearProd7",  
  "trendLinearProd10",
  "regimeProd1",
  "regimeProd1.3",
  "regimeProd2",
  "regimeProd5",
  "regimeProd7",
  "regimeProd10" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab")
##Confusion matrices
conf_matrix_a<-expand.grid(EM=EM,OM=scn)
conf_matrix_a$w_AIC=NA
conf_matrix_a$BIC=NA

cn1<-list()
cn2<-list()

#summarize model selection
aic_set=list()
bic_set=list()

o=0
for(a in seq_along(scn)){

  #head(aic_set[[a]])

   #AIC
  aica<-subset(aic_a,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models


  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_a,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
    head( bic_set[[a]])

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_a$w_AIC[myseq]<-cn1[[a]]
  conf_matrix_a$BIC[myseq]<-cn2[[a]]

  o=max(myseq)

}


conf_matrix_a$eqem_om <- dplyr::recode(conf_matrix_a$OM, 
       "trendLinearProd1"="dynamic.a", 
       "trendLinearProd1.3"="dynamic.a",
       "trendLinearProd2"="dynamic.a", 
       "trendLinearProd5"="dynamic.a",
       "trendLinearProd7"="dynamic.a",
       "trendLinearProd10"="dynamic.a",
       "regimeProd1"="regime.a",
       "regimeProd1.3"="regime.a",
       "regimeProd2"="regime.a",      
       "regimeProd5"="regime.a",      
       "regimeProd7"="regime.a",
       "regimeProd10"="regime.a"
      )   

conf_matrix_a$eqem_om<-factor(conf_matrix_a$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab"         
       ))
conf_matrix_a$diag<-conf_matrix_a$eqem_om==conf_matrix_a$EM


conf_matrix_a$OM2<-dplyr::case_match(conf_matrix_a$OM,
  "trendLinearProd1"~"trend log(a) 1.3 -> 0.03", 
  "trendLinearProd1.3"~"trend log(a) 1.3 -> 0.30",
       "trendLinearProd2"~"trend log(a) 1.3 -> 0.69", 
       "trendLinearProd5"~"trend log(a) 1.3 -> 1.61",
       "trendLinearProd7"~"trend log(a) 1.3 -> 1.95",
       "trendLinearProd10"~"trend log(a) 1.3 -> 2.30",
       "regimeProd1"~"regime log(a) 1.3 -> 0.03", 
       "regimeProd1.3"~"regime log(a) 1.3 -> 0.30",
       "regimeProd2"~"regime log(a) 1.3 -> 0.69",      
       "regimeProd5"~"regime log(a) 1.3 -> 1.61",      
       "regimeProd7"~"regime log(a) 1.3 -> 1.95",
       "regimeProd10"~"regime log(a) 1.3 -> 2.30")

conf_matrix_a$OM2<-factor(conf_matrix_a$OM2, levels=c( "trend log(a) 1.3 -> 0.03", 
 "trend log(a) 1.3 -> 0.30",
  "trend log(a) 1.3 -> 0.69", 
  "trend log(a) 1.3 -> 1.61",
  "trend log(a) 1.3 -> 1.95",
  "trend log(a) 1.3 -> 2.30",
  "regime log(a) 1.3 -> 0.03", 
  "regime log(a) 1.3 -> 0.30",
  "regime log(a) 1.3 -> 0.69",      
  "regime log(a) 1.3 -> 1.61",      
   "regime log(a) 1.3 -> 1.95",
   "regime log(a) 1.3 -> 2.30"))

pa_aic_a=ggplot(data =  conf_matrix_a, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle(expression("AIC Sensitivity"~log(alpha)))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_a, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_aic_a
ggsave("figures/AIC_MLE_sensa.png",
 plot=pa_aic_a, width = 9,height = 7)


pa_bic_a=ggplot(data =  conf_matrix_a, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle(expression("BIC Sensitivity"~log(alpha)))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_a, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_bic_a
ggsave("figures/BIC_MLE_sensa.png", 
  plot=pa_bic_a, width = 9,height = 7)




#LFO
reslfo<-readRDS(file = "outs/simest/sensitivity/resstanloo_a.rds")


##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$LFOmcmc=NA


cn3<-list()
lfomcmc_set_a<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set_a[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  head(lfomcmc_set_a[[a]])
  nsim<-length(unique( lfomcmc_set_a[[a]]$iteration))
  lfomcmc_set_a[[a]]=lfomcmc_set_a[[a]][c(27,8,18,9,24,15,21,12)] #reorder estimation models

  sc3=apply(lfomcmc_set_a[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set_a[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
       "trendLinearProd1"="dynamic.a", 
        "trendLinearProd1.3"="dynamic.a", 
       "trendLinearProd2"="dynamic.a", 
       "trendLinearProd5"="dynamic.a",
       "trendLinearProd7"="dynamic.a",
       "trendLinearProd10"="dynamic.a",
       "regimeProd1"="regime.a",
       "regimeProd1.3"="regime.a",
       "regimeProd2"="regime.a",      
       "regimeProd5"="regime.a",      
       "regimeProd7"="regime.a",
       "regimeProd10"="regime.a"
      ) 
conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab"     
       ))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM


conf_matrix$OM2<-dplyr::case_match(conf_matrix$OM,
  "trendLinearProd1"~"trend log(a) 1.3 -> 0.03", 
  "trendLinearProd1.3"~"trend log(a) 1.3 -> 0.30",
       "trendLinearProd2"~"trend log(a) 1.3 -> 0.69", 
       "trendLinearProd5"~"trend log(a) 1.3 -> 1.61",
       "trendLinearProd7"~"trend log(a) 1.3 -> 1.95",
       "trendLinearProd10"~"trend log(a) 1.3 -> 2.30",
       "regimeProd1"~"regime log(a) 1.3 -> 0.03", 
       "regimeProd1.3"~"regime log(a) 1.3 -> 0.30",
       "regimeProd2"~"regime log(a) 1.3 -> 0.69",      
       "regimeProd5"~"regime log(a) 1.3 -> 1.61",      
       "regimeProd7"~"regime log(a) 1.3 -> 1.95",
       "regimeProd10"~"regime log(a) 1.3 -> 2.30")

conf_matrix$OM2<-factor(conf_matrix$OM2, levels=c( "trend log(a) 1.3 -> 0.03", 
 "trend log(a) 1.3 -> 0.30",
  "trend log(a) 1.3 -> 0.69", 
  "trend log(a) 1.3 -> 1.61",
  "trend log(a) 1.3 -> 1.95",
  "trend log(a) 1.3 -> 2.30",
  "regime log(a) 1.3 -> 0.03", 
  "regime log(a) 1.3 -> 0.30",
  "regime log(a) 1.3 -> 0.69",      
  "regime log(a) 1.3 -> 1.61",      
   "regime log(a) 1.3 -> 1.95",
   "regime log(a) 1.3 -> 2.30"))




pmclfo_sensa=ggplot(data =  conf_matrix, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle(expression("LFO Sensitivity"~log(alpha)))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo_sensa
ggsave("figures/LFO_MCMC_sensa.png",
 plot=pmclfo_sensa,  width = 9,height = 7)

