#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================


library("ggpubr")
library(gridExtra)
library(ggplot2)
library(cowplot)
library(dplyr)



mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)



#========================================================================================================

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


conf_matrix_a$difflog_a_ch<- dplyr::case_match(conf_matrix_a$OM,
  "trendLinearProd1"~"decrease 1.27", 
  "trendLinearProd1.3"~"decrease 1.00",
  "trendLinearProd2"~"decrease 0.61",
  "trendLinearProd5"~"increase 0.31", 
  "trendLinearProd7"~"increase 0.65", 
   "trendLinearProd10"~"increase 1.0", 
   "regimeProd1"~"decrease 1.27",  
   "regimeProd1.3"~"decrease 1.00",        
   "regimeProd2"~"decrease 0.61",      
   "regimeProd5"~"increase 0.31",     
   "regimeProd7"~"increase 0.65",      
   "regimeProd10"~"increase 1.0")   

conf_matrix_a$direction<- dplyr::case_match(conf_matrix_a$OM,
  "trendLinearProd1"~"decrease", 
  "trendLinearProd1.3"~"decrease",
  "trendLinearProd2"~"decrease",
  "trendLinearProd5"~"increase", 
  "trendLinearProd7"~"increase", 
   "trendLinearProd10"~"increase", 
   "regimeProd1"~"decrease",  
   "regimeProd1.3"~"decrease",        
   "regimeProd2"~"decrease",      
   "regimeProd5"~"increase",     
   "regimeProd7"~"increase",      
   "regimeProd10"~"increase")   



conf_matrix_a$difflog_a<- dplyr::case_match(conf_matrix_a$OM,
  "trendLinearProd1"~1.27, 
  "trendLinearProd1.3"~1.00,
  "trendLinearProd2"~0.61,
  "trendLinearProd5"~0.31, 
  "trendLinearProd7"~0.65, 
   "trendLinearProd10"~1.0, 
   "regimeProd1"~1.27,  
   "regimeProd1.3"~1.00,    
   "regimeProd2"~0.61,      
   "regimeProd5"~0.31,     
   "regimeProd7"~0.65,      
   "regimeProd10"~ 1.0)   


conf_matrix_a$type<- dplyr::case_match(conf_matrix_a$OM,
  "trendLinearProd1"~"decrease trend", 
  "trendLinearProd1.3"~"decrease trend",
  "trendLinearProd2"~"decrease trend",
  "trendLinearProd5"~"increase trend", 
  "trendLinearProd7"~"increase trend", 
   "trendLinearProd10"~"increase trend", 
   "regimeProd1"~"decrease regime", 
   "regimeProd1.3"~"decrease regime",         
   "regimeProd2"~"decrease regime",      
   "regimeProd5"~"increase regime",     
   "regimeProd7"~"increase regime",      
   "regimeProd10"~"increase regime") 



conf_matrix_right_a<-conf_matrix_a[conf_matrix_a$EM=="dynamic.a",]

lineAIC_a<-ggplot(conf_matrix_right_a)+
geom_point(aes(x=difflog_a,y=w_AIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=difflog_a,y=w_AIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with AICc")+
xlab(expression("difference"~"in"~log(alpha)))+
coord_cartesian(ylim=c(0.0,0.9))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineAIC_a



lineBIC_a<-ggplot(conf_matrix_right_a)+
geom_point(aes(x=difflog_a,y=BIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=difflog_a,y=BIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with BIC")+
xlab(expression("difference"~"in"~log(alpha)))+
coord_cartesian(ylim=c(0.0,0.88))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineBIC_a


#other colour scheme
lineAIC_a_br<-ggplot(conf_matrix_right_a)+
geom_point(aes(x=difflog_a,y=w_AIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=difflog_a,y=w_AIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#B2182B","#B2182B","#2166AC",  "#2166AC"))+
scale_linetype_manual(values = c(2,1,2,1))+
ylab("proportion of correct model assignment with AICc")+
xlab(expression("difference"~"in"~log(alpha)))+
coord_cartesian(ylim=c(0.0,0.9))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineAIC_a_br
ggsave("figures/lineAICsensa_br.png",
 plot=lineAIC_a_br,width = 10,height = 6)


lineBIC_a_br<-ggplot(conf_matrix_right_a)+
geom_point(aes(x=difflog_a,y=BIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=difflog_a,y=BIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#B2182B","#B2182B","#2166AC",  "#2166AC"))+
scale_linetype_manual(values = c(2,1,2,1))+
ylab("proportion of correct model assignment with BIC")+
xlab(expression("difference"~"in"~log(alpha)))+
coord_cartesian(ylim=c(0.0,0.88))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineBIC_a_br
ggsave("figures/lineBICsensa_br.png",
 plot=lineBIC_a_br,width = 10,height = 6)



#========================================================================================================
#smax

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
    "trendLinearSmax025"~"trend Smax * 0.25",
       "trendLinearSmax050"~"trend Smax * 0.50",
       "trendLinearSmax150"~"trend Smax * 1.50", 
       "trendLinearSmax200"~"trend Smax * 2.00", 
       "trendLinearSmax300"~"trend Smax * 3.00",
       "regimeSmax025"~"regime Smax * 0.25",      
       "regimeSmax050"~"regime Smax * 0.50",      
       "regimeSmax150"~"regime Smax * 1.50",      
       "regimeSmax200"~"regime Smax * 2.00",      
       "regimeSmax300"~"regime Smax * 3.00" )


conf_matrix_smax$OM2<-factor(conf_matrix_smax$OM2, levels=c("trend Smax * 0.25",
       "trend Smax * 0.50",
       "trend Smax * 1.50", 
       "trend Smax * 2.00", 
       "trend Smax * 3.00",
       "regime Smax * 0.25",      
       "regime Smax * 0.50",     
       "regime Smax * 1.50",      
       "regime Smax * 2.00",      
       "regime Smax * 3.00"   ))



#calculate these 
conf_matrix_smax$diffsmax_ch<- dplyr::case_match(conf_matrix_smax$OM,
  "trendLinearSmax025"~"decrease x0.25",
   "trendLinearSmax050"~"decrease x0.5",
   "trendLinearSmax150"~"increase x1.5",
 "trendLinearSmax200"~"increase x2.0", 
 "trendLinearSmax300"~"increase x3.0", 
 "regimeSmax025"~"decrease x0.25",     
 "regimeSmax050"~"decrease x0.5",      
 "regimeSmax150"~"increase x1.5",      
 "regimeSmax200"~"increase x2.0",      
 "regimeSmax300" ~"increase x3.0")   




conf_matrix_smax$diffsmax<- dplyr::case_match(conf_matrix_smax$OM,
  "trendLinearSmax025"~4,
   "trendLinearSmax050"~2,
   "trendLinearSmax150"~1.5,
 "trendLinearSmax200"~2.0, 
 "trendLinearSmax300"~3.0, 
 "regimeSmax025"~4,     
 "regimeSmax050"~2,      
 "regimeSmax150"~1.5,      
 "regimeSmax200"~2.0,      
 "regimeSmax300" ~3.0)



conf_matrix_smax$type<- dplyr::case_match(conf_matrix_smax$OM,
 "trendLinearSmax025"~"decrease trend",
   "trendLinearSmax050"~"decrease trend",
   "trendLinearSmax150"~"increase trend",
 "trendLinearSmax200"~"increase trend", 
 "trendLinearSmax300"~"increase trend", 
 "regimeSmax025"~"decrease regime",     
 "regimeSmax050"~"decrease regime",      
 "regimeSmax150"~"increase regime",      
 "regimeSmax200"~"increase regime",      
 "regimeSmax300" ~"increase regime")   

 
conf_matrix_right_smax<-conf_matrix_smax[conf_matrix_smax$EM=="dynamic.b",]

lineAICsmax<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=w_AIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=diffsmax,y=w_AIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("     ")+
xlab(expression("multiplication"~"factor"~"for"~S[max]))+
coord_cartesian(ylim=c(0.1,0.88))+
mytheme +
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineAICsmax




lineAICsmax_wlabel<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=w_AIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=diffsmax,y=w_AIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with AICc")+
xlab(expression("multiplication"~"factor"~"for"~S[max]))+
coord_cartesian(ylim=c(0.1,0.88))+
mytheme +
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))
lineAICsmax_wlabel

ggsave("figures/lineAICsenssmax.png",
 plot=lineAICsmax_wlabel,width = 10,height = 6)






lineBICsmax<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=BIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=diffsmax,y=BIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("      ")+
xlab(expression("multiplication"~"factor"~"for"~S[max]))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineBICsmax





lineBICsmax_wlabel<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=BIC,color=type),size=4,show.legend = F)+
geom_line(aes(x=diffsmax,y=BIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with BIC")+
xlab(expression("multiplication"~"factor"~"for"~S[max]))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))
lineBICsmax_wlabel

ggsave("figures/lineBICsenssmax.png",
 plot=lineBICsmax_wlabel,width = 10,height = 6)



#----------------------------------
#need to run plots_cluster_confmat_sensa before this works

linesAIC_alpha_smax <- ggarrange(lineAIC_a, lineAICsmax,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesAIC_alpha_smax

ggsave("figures/lineAICsens_alpha_smax.png",
 plot=linesAIC_alpha_smax, width = 13,height = 6)


linesBIC_alpha_smax <- ggarrange(lineBIC_a, lineBICsmax,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesBIC_alpha_smax

ggsave("figures/lineBICsens_alpha_smax.png",
 plot=linesBIC_alpha_smax, width = 13,height = 6)


#====================================================================
#sigma

res_sigmed<-readRDS(file = "outs/simest/sigmamed_sensitivity/res_sigmed.rds")

res_sigmed<-res_sigmed[res_sigmed$convergence==0,]

res_sigmed$sigscenario<-"medium"
res_sigmed$scenario<-case_match(
res_sigmed$scenario,
 "sigmamed_stationary"~"stationary",
 "sigmamed_decLinearProd"~"decLinearProd",        
 "sigmamed_regimeProd"~ "regimeProd",            
 "sigmamed_sineProd"~ "sineProd",             
 "sigmamed_regimeCap"~ "regimeCap",             
 "sigmamed_decLinearCap"~ "decLinearCap",        
 "sigmamed_regimeProdCap"~ "regimeProdCap",
 "sigmamed_shiftCap"~"shiftCap",
 "sigmamed_decLinearProdshiftCap"~"decLinearProdshiftCap")

aic=subset(res_sigmed,parameter=='AIC'&method=='MLE')
bic=subset(res_sigmed,parameter=='BIC'&method=='MLE')

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
conf_matrix_sigmed<-expand.grid(EM=EM,OM=scn)
conf_matrix_sigmed$w_AIC=NA
conf_matrix_sigmed$BIC=NA

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
  conf_matrix_sigmed$w_AIC[myseq]<-cn1[[a]]
  conf_matrix_sigmed$BIC[myseq]<-cn2[[a]]
  
  o=max(myseq)

}

conf_matrix_sigmed$eqem_om <- dplyr::recode(conf_matrix_sigmed$OM, 
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
conf_matrix_sigmed$eqem_om<-factor( 
  conf_matrix_sigmed$eqem_om,levels=c( "stationary",
 "autocorr", "dynamic.a", "regime.a", "dynamic.b", 
 "regime.b", "dynamic.ab", "regime.ab" ))      

conf_matrix_sigmed$diag<-conf_matrix_sigmed$eqem_om==conf_matrix_sigmed$EM

conf_matrix_sigmed$scentype <-
dplyr::case_match(as.character(conf_matrix_sigmed$OM), 
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

conf_matrix_sigmed$scendesc<-
  dplyr::case_match(as.character(conf_matrix_sigmed$OM), 
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
 
conf_matrix_sigmed$fullname <- apply( 
  conf_matrix_sigmed[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )

conf_matrix_sigmed$fullname<-factor(conf_matrix_sigmed$fullname, 
  levels=c(" stationary",
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


#========================================================================================================
#base case - sigma=0.6
resbase<-readRDS(file = "outs/simest/generic/resbase.rds")

resbase<-resbase[resbase$convergence==0,]

resbase$sigscenario<-"high"

resbase<-resbase[resbase$scenario%in%c("stationary",
                                       "decLinearProd",
                                       "sineProd",
                                       "regimeProd",
                                       "decLinearCap",         
                                       "regimeCap",
                                       "shiftCap",
                                       "decLinearProdshiftCap",
                                       "regimeProdCap"),]

aicbase=subset(resbase,parameter=='AIC'&method=='MLE')
bicbase=subset(resbase,parameter=='BIC'&method=='MLE')

##Confusion matrices
conf_matrix_base<-expand.grid(EM=EM,OM=scn)
conf_matrix_base$w_AIC=NA
conf_matrix_base$BIC=NA


cn1_base<-list()
cn2_base<-list()


aic_set_base<-list()
bic_set_base<-list()
o=0

for(a in seq_along(scn)){

  #AIC
  aica_base<-subset(aicbase,scenario==scn[a])
  aic_set_base[[a]]=tidyr::spread(aica_base[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set_base[[a]]$iteration))
  aic_set_base[[a]]=aic_set_base[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(aic_set[[a]])
  sc1_base=apply(aic_set_base[[a]],1,which.min)
  
  cn1_base[[a]]=summary(factor(sc1_base,levels=seq(1:ncol(aic_set_base[[a]]))))/naicsim

  bica_base=subset(bicbase,scenario==scn[a])
  bic_set_base[[a]]=tidyr::spread(bica_base[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set_base[[a]]$iteration))
  bic_set_base[[a]]=bic_set_base[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models

  sc2_base=apply(bic_set_base[[a]],1,which.min)
 
  cn2_base[[a]]=summary(factor(sc2_base,levels=seq(1:ncol(bic_set_base[[a]]))))/nbicsim

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_base$w_AIC[myseq]<-cn1_base[[a]]
  conf_matrix_base$BIC[myseq]<-cn2_base[[a]]
  o=max(myseq)

}

conf_matrix_base$eqem_om <- dplyr::recode(conf_matrix_base$OM, 
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
conf_matrix_base$eqem_om<-factor( conf_matrix_base$eqem_om,levels=c( "stationary",
 "autocorr", "dynamic.a", "regime.a", "dynamic.b", "regime.b", "dynamic.ab", "regime.ab"
  ))      

conf_matrix_base$diag<-conf_matrix_base$eqem_om==conf_matrix_base$EM

conf_matrix_base$scentype <-
dplyr::case_match(as.character(conf_matrix_base$OM), 
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

                                
                                 
                      
conf_matrix_base$scendesc<-
  dplyr::case_match(as.character(conf_matrix_base$OM), 
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
 
conf_matrix_base$fullname <- apply( 
  conf_matrix_base[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )

conf_matrix_base$fullname<-factor(conf_matrix_base$fullname, 
  levels=c(" stationary",
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



#======================================================================================================================
#aggregate plot

conf_matrix_sigmed$sigma<-0.3
conf_matrix_base$sigma<-0.6

names(conf_matrix_sigmed)
names(conf_matrix_base)

conf_matrix_sigcomp<-rbind(conf_matrix_sigmed,conf_matrix_base)
conf_matrix_sigcomp$OM

conf_matrix_sigcomp<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM%in%c("decLinearProd","regimeProd","decLinearCap","shiftCap"), ]

conf_matrix_sigcomp$OM<-as.character(conf_matrix_sigcomp$OM)
conf_matrix_a$EM<-as.character(conf_matrix_a$EM)


conf_matrix[conf_matrix$OM=="regimeProd",]
conf_matrix_base[conf_matrix_base$OM=="regimeProd",]

conf_matrix_sigcomp_sc1<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="decLinearProd"&conf_matrix_sigcomp$EM=="dynamic.a",]
conf_matrix_sigcomp_sc2<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="regimeProd"&conf_matrix_sigcomp$EM%in%c("dynamic.a"),]
conf_matrix_sigcomp_sc3<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="decLinearCap"&conf_matrix_sigcomp$EM=="dynamic.b",]
conf_matrix_sigcomp_sc4<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="shiftCap"&conf_matrix_sigcomp$EM%in%c("dynamic.b"),]
#conf_matrix_sigcomp_sc5<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="decLinearProdshiftCap"&conf_matrix_sigcomp$EM=="dynamic.ab",]



conf_matrix_sigcomp_right<-rbind(conf_matrix_sigcomp_sc1,
  conf_matrix_sigcomp_sc2,
  conf_matrix_sigcomp_sc3,
  conf_matrix_sigcomp_sc4)



conf_matrix_sigcomp_right$OM2<-case_match(
  conf_matrix_sigcomp_right$OM,
  "decLinearProd" ~"decline log(a)",
  "regimeProd" ~ "shift up log(a)",   
  "decLinearCap" ~"decline Smax",  
  "shiftCap"  ~"shift down Smax")
conf_matrix_sigcomp_right$OM2<-factor(conf_matrix_sigcomp_right$OM2,levels=c("shift up log(a)","decline log(a)",
  "shift down Smax", "decline Smax"))


lineAIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_point(aes(x=sigma,y=w_AIC,color=OM2),size=4,show.legend = F)+
geom_line(aes(x=sigma,y=w_AIC,color=OM2),linewidth=1.2)+
ylab("proportion of correct model assignment with AICc")+
xlab(expression(paste("value of", ~sigma)))+
scale_color_viridis_d(begin=.1, end=.8,option = "A") +
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineAIC_sig

barAIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_bar(stat="identity",aes(x=as.factor(sigma),y=w_AIC,fill=OM2), position=position_dodge())+
ylab("proportion of correct model assignment with AICc")+
xlab(expression(paste("value of", ~sigma)))+
scale_fill_viridis_d(begin=.1, end=.8,option = "B", labels = c(expression(shift~up~log(alpha)) , 
                                                              expression(decline~log(alpha)),
                                                              expression(shift~down~S[max]),
                                                              expression(decline~S[max]))) +
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(1.5, "line"))+
guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
barAIC_sig



lineBIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_point(aes(x=sigma,y=BIC,color=OM2),size=4,show.legend = F)+
geom_line(aes(x=sigma,y=BIC,color=OM2),linewidth=1.2)+
ylab("proportion of correct model assignment with BIC")+
xlab(expression(paste("value of", ~sigma)))+
scale_color_viridis_d(begin=.1, end=.8,option = "B") +
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(1.5, "line"))+
guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineBIC_sig


barBIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_bar(stat="identity",aes(x=as.factor(sigma),y=BIC,fill=OM2), position=position_dodge())+
ylab("proportion of correct model assignment with BIC")+
xlab(expression(paste("value of", ~sigma)))+
scale_fill_viridis_d(begin=.1, end=.9,option = "B", labels = c(expression(shift~up~log(alpha)) , 
                                                              expression(decline~log(alpha)),
                                                              expression(shift~down~S[max]),
                                                              expression(decline~S[max]))) +
mytheme +
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(1.5, "line")) +
guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
barBIC_sig


#----------------------------------
#need to run plots_cluster_confmat_sensa before this works
#require that you run confmat_sens_a first
linesAIC_alpha_smax <- ggarrange(lineAIC_a, lineAICsmax, #lineAIC_sig,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesAIC_alpha_smax




plot_grid(linesAIC_alpha_smax,lineAIC_sig,  rel_widths = c(2, 1) )


plot_grid(linesAIC_alpha_smax,barAIC_sig,  rel_widths = c(2, 1) )


ggsave("figures/lineAICsens_alpha_smax.png",
 plot=linesAIC_alpha_smax, width = 13,height = 6)




linesBIC_alpha_smax <- ggarrange(lineBIC_a, lineBICsmax, 
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesBIC_alpha_smax


linesBIC_alpha_smax_sig<-plot_grid(linesBIC_alpha_smax,lineBIC_sig,  rel_widths = c(2, 1) )
linesBIC_alpha_smax_sig
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineBICsens_alpha_smax_sig.png",
 plot=linesBIC_alpha_smax_sig, width = 18,height = 6)



linesBIC_alpha_smax_sigbar<-plot_grid(linesBIC_alpha_smax,barBIC_sig,  rel_widths = c(2, 1) )
linesBIC_alpha_smax_sigbar
ggsave("figures/lineBICsens_alpha_smax_sigbar.png",
 plot=linesBIC_alpha_smax_sigbar, width = 18,height = 6)


#blue and red 

linesBIC_alpha_smax_br <- ggarrange(lineBIC_a_br, lineBICsmax_br, 
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesBIC_alpha_smax_br
linesBIC_alpha_smax_sigbar_br<-plot_grid(linesBIC_alpha_smax_br,barBIC_sig,  rel_widths = c(2, 1) )
linesBIC_alpha_smax_sigbar_br


ggsave("figures/lineBICsens_alpha_smax_sigbar_br.png",
 plot=linesBIC_alpha_smax_sigbar_br, width = 18,height = 6)






#---------------------------------------------------------------------------------------------
#stan



#===================================================================================
reslfo<-readRDS(file = "outs/simest/sigmamed_sensitivity/resstanloo_sigmed.rds")
head(reslfo)


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
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_med/LFO_MCMC_sigmed.png",
 plot=pmclfo, width = 9,height = 7)


#lfo with 3 last yrs average
reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA


cn4<-list()
lfomcmc_set3<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3[[a]]$iteration))
  lfomcmc_set3[[a]]=lfomcmc_set3[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc4=apply(lfomcmc_set3[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_3yrs.png", plot=pmclfo3)




#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5[[a]]$iteration))
  lfomcmc_set5[[a]]=lfomcmc_set5[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_5yrs.png", plot=pmclfo5)




















#LFO
reslfo<-readRDS(file = "outs/simest/Smax_sensitivity/resstanloo_smax.rds")
head(reslfo)

#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
#    "hmma_last5","hmmb_last5","hmmab_last5"),]

#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last5"="rwa",
#      "rwb_last5"="rwb",
#      "rwab_last5"="rwab",
#    "hmma_last5"="hmma",
#    "hmmb_last5"="hmmb",
#    "hmmab_last5"="hmmab"
#      ) 


#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
#    "hmma_last3","hmmb_last3","hmmab_last3"),]


#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last3"="rwa",
#      "rwb_last3"="rwb",
#      "rwab_last3"="rwab",
#    "hmma_last3"="hmma",
#    "hmmb_last3"="hmmb",
#    "hmmab_last3"="hmmab"
#      ) 


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
conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab" 
         
        
       ))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM







pmclfo_senssmax=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("sens smax LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo_senssmax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/LFO_MCMC_senssmax.png",
 plot=pmclfo_senssmax)


#lfo with 3 last yrs average


reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA


cn4<-list()
lfomcmc_set3_a<-list()


o=0
for(a in seq_along(scn)){
 
  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3_a[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3_a[[a]]$iteration))
  lfomcmc_set3_a[[a]]=lfomcmc_set3_a[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  #head(lfomcmc_set3_a[[a]])

  sc4=apply(lfomcmc_set3_a[[a]],1,which.max)
  cn4[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3_a[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn4[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo3_senssmax=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("sens smax LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3_senssmax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/LFO_MCMC_3yrs_smax.png",
 plot=pmclfo3_senssmax)




#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5_a<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5_a[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5_a[[a]]$iteration))
  lfomcmc_set5_a[[a]]=lfomcmc_set5_a[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5_a[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5_a[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo5_senssmax=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("sens smax LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5_senssmax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/LFO_MCMC_5yrs_smax.png",
 plot=pmclfo5_senssmax)





#========================================================
#old code
#========================================================

simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_smax<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")
#res_smax95<-readRDS(file = "outs/simest/res_gamma_alpha/res_smax_95.rds")


ressmax<-rbind(res_smax)

#res_a <- res_a[res_a$convergence==0,]

aic_smax=subset(ressmax, parameter=='AIC'&method=='MLE')
bic_smax=subset(ressmax, parameter=='BIC'&method=='MLE')
lfo_smax=subset(ressmax, parameter=='LFO'&method=='MLE')


lfo_smax<-lfo_smax[lfo_smax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_smax[is.na(lfo_smax$est),]<--Inf

lfo_smax<-lfo_smax[lfo_smax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_smax[is.na(lfo_smax$est),]<--Inf
aic_smax$est[aic_smax$convergence>0]<-Inf
bic_smax$est[aic_smax$convergence>0]<-Inf

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
    "regimeSmax300"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  dim(aica)
  unique(aica$scenario)


   #AIC
  aica<-subset(aic_smax,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_smax,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_smax,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
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

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p


reslfo<-readRDS(file = "outs/simest/generic/resstanloo.rds")
head(reslfo)

#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
#    "hmma_last5","hmmb_last5","hmmab_last5"),]

#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last5"="rwa",
#      "rwb_last5"="rwb",
#      "rwab_last5"="rwab",
#    "hmma_last5"="hmma",
#    "hmmb_last5"="hmmb",
#    "hmmab_last5"="hmmab"
#      ) 


#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
#    "hmma_last3","hmmb_last3","hmmab_last3"),]


#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last3"="rwa",
#      "rwb_last3"="rwb",
#      "rwab_last3"="rwab",
#    "hmma_last3"="hmma",
#    "hmmb_last3"="hmmb",
#    "hmmab_last3"="hmmab"
#      ) 

scn<-factor(unique(reslfo$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "regimeProd",
  "shiftProd", 
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
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$LFOmcmc=NA


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
  conf_matrix$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
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
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM






pmclfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC.png", plot=pmclfo)


#lfo with 3 last yrs average


reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA


cn4<-list()
lfomcmc_set3<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3[[a]]$iteration))
  lfomcmc_set3[[a]]=lfomcmc_set3[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc4=apply(lfomcmc_set3[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_3yrs.png", plot=pmclfo3)




#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5[[a]]$iteration))
  lfomcmc_set5[[a]]=lfomcmc_set5[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/sensa_LFO_MCMC_5yrs.png", plot=pmclfo5)











#========================================================================================================
#sensitivity smax scenario double alpha
#read in data
simPar <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_smaxda<-readRDS(file = "outs/simest/Smax_sensitivity_doublealpha/res_smaxda.rds")
#res_smaxda56<-readRDS(file = "outs/simest/res_gamma_alpha/res_smaxda_56.rds")


ressmaxda<-rbind(res_smaxda)

#res_a <- res_a[res_a$convergence==0,]

aic_smaxda=subset(ressmaxda, parameter=='AIC'&method=='MLE')
bic_smaxda=subset(ressmaxda, parameter=='BIC'&method=='MLE')
lfo_smaxda=subset(ressmaxda, parameter=='LFO'&method=='MLE')


lfo_smaxda<-lfo_smaxda[lfo_smaxda$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_smaxda[is.na(lfo_smaxda$est),]<--Inf

lfo_smaxda<-lfo_smaxda[lfo_smaxda$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_smaxda[is.na(lfo_smaxda$est),]<--Inf
aic_smaxda$est[aic_smaxda$convergence>0]<-Inf
bic_smaxda$est[aic_smaxda$convergence>0]<-Inf

scn<-factor(unique(aic_smaxda$scenario), levels=c(
    "trendLinearSmax025_da",
    "trendLinearSmax050_da", 
    "trendLinearSmax150_da", 
    "trendLinearSmax200_da", 
    "trendLinearSmax300_da", 
    "regimeSmax025_da",     
    "regimeSmax050_da",      
    "regimeSmax150_da",      
    "regimeSmax200_da",      
    "regimeSmax300_da"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  dim(aica)
  unique(aica$scenario)


   #AIC
  aica<-subset(aic_smaxda,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_smaxda,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_smaxda,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "trendLinearSmax025_da"="dynamic.b",
    "trendLinearSmax050_da"="dynamic.b", 
    "trendLinearSmax150_da"="dynamic.b", 
    "trendLinearSmax200_da"="dynamic.b", 
    "trendLinearSmax300_da"="dynamic.b", 
    "regimeSmax025_da"="regime.b",     
    "regimeSmax050_da"="regime.b",      
    "regimeSmax150_da"="regime.b",      
    "regimeSmax200_da"="regime.b",      
    "regimeSmax300_da"="regime.b"
      )   

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens smax double alpha")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens smax  double alpha")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax  double alpha")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p




#========================================================================================================
#sigma low sensitivity 
#read in data
simPar <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_siglow<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow.rds")

res_siglow36<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow_36.rds")
res_siglow84<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow_84.rds")

ressiglow<-rbind(res_siglow,res_siglow36,res_siglow84)

#res_a <- res_a[res_a$convergence==0,]

aic_siglow=subset(ressiglow, parameter=='AIC'&method=='MLE')
bic_siglow=subset(ressiglow, parameter=='BIC'&method=='MLE')
lfo_siglow=subset(ressiglow, parameter=='LFO'&method=='MLE')


lfo_siglow<-lfo_siglow[lfo_siglow$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_siglow[is.na(lfo_siglow$est),]<--Inf

lfo_siglow<-lfo_siglow[lfo_siglow$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_siglow[is.na(lfo_siglow$est),]<--Inf
aic_siglow$est[aic_siglow$convergence>0]<-Inf
bic_siglow$est[aic_siglow$convergence>0]<-Inf

scn<-factor(unique(aic_siglow$scenario), levels=c(
     "sigmalow_stationary",
     "sigmalow_decLinearProd",        
     "sigmalow_regimeProd",            
     "sigmalow_sineProd",             
     "sigmalow_regimeCap",             
     "sigmalow_decLinearCap",         
     "sigmalow_regimeProdCap",         
     "sigmalow_shiftCap",             
     "sigmalow_decLinearProdshiftCap"
  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()

unique(aica$model)
o=0

aica[7282, ]
aica[8002, ]

for(a in seq_along(scn)){

   #AIC
  aica<-subset(aic_siglow,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_siglow,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_siglow,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "sigmalow_stationary"="stationary",
     "sigmalow_decLinearProd"="dynamic.a",        
     "sigmalow_regimeProd"="regime.a",            
     "sigmalow_sineProd"="dynamic.a",             
     "sigmalow_regimeCap"="regime.b",             
     "sigmalow_decLinearCap"="dynamic.b",         
     "sigmalow_regimeProdCap"="regime.ab",         
     "sigmalow_shiftCap"="regime.b",             
     "sigmalow_decLinearProdshiftCap"="dynamic.ab")
    

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens sigma low")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens sigma low")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax sigma low")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



#========================================================================================================
#sigma med sensitivity 
#read in data
simPar <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_sigmed<-readRDS(file = "outs/simest/sigmamed_sensitivity/res_sigmed.rds")

#res_sigmed73<-readRDS(file = "outs/simest/res_gamma_alpha/res_sigmed_73.rds")


ressigmed<-rbind(res_sigmed)

#res_a <- res_a[res_a$convergence==0,]

aic_sigmed=subset(ressigmed, parameter=='AIC'&method=='MLE')
bic_sigmed=subset(ressigmed, parameter=='BIC'&method=='MLE')
lfo_sigmed=subset(ressigmed, parameter=='LFO'&method=='MLE')


lfo_sigmed<-lfo_sigmed[lfo_sigmed$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_sigmed[is.na(lfo_sigmed$est),]<--Inf

lfo_sigmed<-lfo_sigmed[lfo_sigmed$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_sigmed[is.na(lfo_sigmed$est),]<--Inf
aic_sigmed$est[aic_sigmed$convergence>0]<-Inf
bic_sigmed$est[aic_sigmed$convergence>0]<-Inf

scn<-factor(unique(aic_sigmed$scenario), levels=c(
     "sigmamed_stationary",
     "sigmamed_decLinearProd",        
     "sigmamed_regimeProd",            
     "sigmamed_sineProd",             
     "sigmamed_regimeCap",             
     "sigmamed_decLinearCap",         
     "sigmamed_regimeProdCap",         
     "sigmamed_shiftCap",             
     "sigmamed_decLinearProdshiftCap"
  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

   #AIC
  aica<-subset(aic_sigmed,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_sigmed,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_sigmed,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "sigmamed_stationary"="stationary",
     "sigmamed_decLinearProd"="dynamic.a",        
     "sigmamed_regimeProd"="regime.a",            
     "sigmamed_sineProd"="dynamic.a",             
     "sigmamed_regimeCap"="regime.b",             
     "sigmamed_decLinearCap"="dynamic.b",         
     "sigmamed_regimeProdCap"="regime.ab",         
     "sigmamed_shiftCap"="regime.b",             
     "sigmamed_decLinearProdshiftCap"="dynamic.ab")
    

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens sigma med")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens sigma med")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax sigma med")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p




