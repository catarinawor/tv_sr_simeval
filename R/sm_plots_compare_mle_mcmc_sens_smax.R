#======================================================
#sensitivity smax scenarios case full confusion matrices
#Catarina Wor
#November 2022
#============================================


library("ggpubr")
library(gridExtra)
library(ggplot2)

mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

#========================================================================================================

#========================================================================================================
#sensitivity smax scenario
#read in data

restmb_smax<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")

aic_smax=subset(restmb_smax, parameter=='AIC'&method=='MLE')
bic_smax=subset(restmb_smax, parameter=='BIC'&method=='MLE')

aic_smax$mode[aic_smax$convergence>0]<-Inf
bic_smax$mode[aic_smax$convergence>0]<-Inf

scn<-factor(unique(aic_smax$scenario), levels=c(
  "trendLinearSmax025", 
  "trendLinearSmax050",
  "trendLinearSmax150", 
  "trendLinearSmax200", 
  "trendLinearSmax300",
  "regimeSmax025",      
  "regimeSmax050",      
  "regimeSmax150",      
  "regimeSmax200",      
  "regimeSmax300" 
 ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab")
##Confusion matrices
conf_matrix_smax<-expand.grid(EM=EM,OM=scn)
conf_matrix_smax$w_AIC=NA
conf_matrix_smax$BIC=NA

cn1<-list()
cn2<-list()

#summarize model selection
aic_set=list()
bic_set=list()

o=0
for(a in seq_along(scn)){

   #AIC
  aica<-subset(aic_smax,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_smax,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
    head( bic_set[[a]])

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_smax$w_AIC[myseq]<-cn1[[a]]
  conf_matrix_smax$BIC[myseq]<-cn2[[a]]
  o=max(myseq)

}


conf_matrix_smax$eqem_om <- dplyr::recode(conf_matrix_smax$OM, 
      "trendLinearSmax025"="dynamic.b", 
       "trendLinearSmax050"="dynamic.b",
       "trendLinearSmax150"="dynamic.b", 
       "trendLinearSmax200"="dynamic.b", 
       "trendLinearSmax300"="dynamic.b",
       "regimeSmax025"="regime.b",      
       "regimeSmax050"="regime.b",      
       "regimeSmax150"="regime.b",      
       "regimeSmax200"="regime.b",      
       "regimeSmax300"="regime.b"         
      )   

conf_matrix_smax$eqem_om<-factor(conf_matrix_smax$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab"        
       ))
conf_matrix_smax$diag<-conf_matrix_smax$eqem_om==conf_matrix_smax$EM


conf_matrix_smax$OM2<-dplyr::case_match(conf_matrix_smax$OM,
  "trendLinearSmax025" ~ "trend 25% Smax",
  "trendLinearSmax050" ~ "trend 50% Smax",
  "trendLinearSmax150" ~ "trend 150% Smax",
  "trendLinearSmax200" ~ "trend 200% Smax",
  "trendLinearSmax300" ~ "trend 300% Smax",
  "regimeSmax025" ~ "regime 25% Smax",
  "regimeSmax050" ~ "regime 50% Smax",
  "regimeSmax150" ~ "regime 150% Smax",
  "regimeSmax200" ~ "regime 200% Smax",
  "regimeSmax300"~ "regime 300% Smax")

conf_matrix_a$OM2<-factor(conf_matrix_a$OM2, levels=c( "trend 25% Smax", 
 "trend 50% Smax",
  "trend 150% Smax", 
  "trend 200% Smax",
  "trend 300% Smax",
  "regime 25% Smax",
  "regime 50% Smax", 
  "regime 150% Smax",
  "regime 200% Smax",      
  "regime 300% Smax"))


pa_aic_smax=ggplot(data =  conf_matrix_smax, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle(expression("AIC Sensitivity"~S[max]))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_smax, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_aic_smax
ggsave("figures/AIC_MLE_senssmax.png",
 plot=pa_aic_smax,  width = 9,height = 7)


pa_bic_smax=ggplot(data =  conf_matrix_smax, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle(expression("BIC Sensitivity"~S[max]))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_smax, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_bic_smax
ggsave("figures/BIC_MLE_senssmax.png",
 plot=pa_bic_smax,  width = 9,height = 7)


