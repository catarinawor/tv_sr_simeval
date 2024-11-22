#====================================================
#routines to visualize smulation estimation results
#Catarina Wor
#November 2024
#====================================================

library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)

source("R/utils.R")

mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)


source("R/read_er_data.R")

#========================================================================================================
df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", 
                                      "convergence","conv_warning","pbias","bias"))

df$variable<-factor(df$variable,levels=c("median","mode", "sim"))

df$scenario<-factor(df$scenario,levels=c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "incLinearProd_lowERLowError", 
                                         "highERHighError",               
                                         "highERLowError",               
                                         "lowERHighError",                
                                         "lowERLowError",                
                                          "ShiftERHighError",              
                                          "ShiftERLowError",               
                                          "trendERHighError",              
                                          "trendERLowError" )  ) 

df$ERtrend<-case_match(df$scenario,"decLinearProd_highERLowError"~"highER",
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

df$ERerror<-case_match(df$scenario,"decLinearProd_highERLowError"~"LowError",
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

df$dynamics<-case_match(df$scenario,"decLinearProd_highERLowError"~"decLinear",
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

df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"  ))
 
df$model2<-case_match(df$model,
    "simple"~"stationary",
    "autocorr"~"autocorr", 
    "rwa"~ "rw.a",
    "hmma" ~"hmm.a",
    "rwb" ~"rw.b",
    "hmmb" ~"hmm.b", 
    "rwab"~"rw.ab",
    "hmmab"~ "hmm.ab" )

df$model2<-factor(df$model2,levels=c("stationary",
                                   "autocorr", 
                                   "rw.a",
                                   "hmm.a",
                                   "rw.b",
                                   "hmm.b",
                                   "rw.ab",
                                   "hmm.ab"  ))

summarydf  <- df %>%
   group_by(scenario,parameter,
    method,model,by,variable,dynamics,
    ERerror,ERtrend) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))

summarydf <- reshape2::dcast(data=summarydf,  
    scenario + parameter + method + model + by + variable + dynamics + ERerror + ERtrend  ~prob, 
    value.var= "qs",fun.aggregate=mean)

summarydf$model2<-case_match(summarydf$model,
    "simple"~"stationary",
    "autocorr"~"autocorr", 
    "rwa"~ "rw.a",
    "hmma" ~"hmm.a",
    "rwb" ~"rw.b",
    "hmmb" ~"hmm.b", 
    "rwab"~"rw.ab",
    "hmmab"~ "hmm.ab" )

summarydf$model2<-factor(summarydf$model2, levels=c("stationary","autocorr",
    "rw.a","hmm.a","rw.b","hmm.b","rw.ab", "hmm.ab"))




#=======================================================================================
#stats --- main paper graphs
#meancv


resparam$ERtrend<-case_match(resparam$scenario,"decLinearProd_highERLowError"~"highER",
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



resparam$ERerror<-case_match(resparam$scenario,"decLinearProd_highERLowError"~"LowError",
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


resparam$dynamics<-case_match(resparam$scenario,"decLinearProd_highERLowError"~"decLinear",
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



cvdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by,
                               ERtrend=resparam$ERtrend,
                               ERerror=resparam$ERerror,
                               dynamics=resparam$dynamics), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})


cvdfsmsy<-cvdf[cvdf$parameter=="smsy"&cvdf$method=="MLE"&
              cvdf$scenario%in%c("highERLowError", 
                        "lowERLowError", 
                        "decLinearProd_highERLowError",
                        "incLinearProd_highERLowError", 
                        "decLinearProd_lowERLowError", 
                        "incLinearProd_lowERLowError")&
               cvdf$model%in%c("autocorr", "rwa"),]





meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               method=cvdf$method,
                               model=cvdf$model,
                               ERtrend=cvdf$ERtrend,
                               ERerror=cvdf$ERerror,
                               dynamics=cvdf$dynamics
                              ), function(x){mean(x,na.rm=T)})




meancvdf_redux<-meancvdf[meancvdf$parameter=="smsy"&
meancvdf$scenario%in%c("highERLowError", 
                        "lowERLowError", 
                        "decLinearProd_highERLowError",
                        "incLinearProd_highERLowError", 
                        "decLinearProd_lowERLowError", 
                        "incLinearProd_lowERLowError")&
meancvdf$model%in%c("autocorr", "rwa")&
meancvdf$method=="MLE",]

meancvdf_redux$meancv<-round(meancvdf_redux$x,2)


meancvdf_redux$dynamics2<-case_match(meancvdf_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")

meancvdf_redux$model2<-case_match(meancvdf_redux$model,
    "autocorr"~"autocorr", 
    "rwa"~"rw.a")

meancvdf_redux$model2<-factor(meancvdf_redux$model2, levels=c("autocorr", "rw.a"))

#median
#absolute median percent bias
abspbiasdf<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by,
                               ERtrend=resparam$ERtrend,
                               ERerror=resparam$ERerror,
                               dynamics=resparam$dynamics ), function(x){median(abs(x),na.rm=T)})

meanabspbiasdf<-aggregate(abspbiasdf$x, list(parameter=abspbiasdf$parameter,
                               scenario=abspbiasdf$scenario,
                               method=abspbiasdf$method,
                               model=abspbiasdf$model,
                               ERtrend=abspbiasdf$ERtrend,
                               ERerror=abspbiasdf$ERerror,
                               dynamics=abspbiasdf$dynamics ), function(x){mean(x,na.rm=T)})

meanabspbiasdf$model2<-case_match(meanabspbiasdf$model,
    "simple"~"stationary",
    "autocorr"~"autocorr", 
    "rwa"~ "rw.a",
    "hmma" ~"hmm.a",
    "rwb" ~"rw.b",
    "hmmb" ~"hmm.b", 
    "rwab"~"rw.ab",
    "hmmab"~ "hmm.ab" )


meanabspbiasdf_redux<-meanabspbiasdf[meanabspbiasdf$parameter=="smsy"&
meanabspbiasdf$scenario%in%c("highERLowError", 
                        "lowERLowError", 
                        "decLinearProd_highERLowError",
                        "incLinearProd_highERLowError", 
                        "decLinearProd_lowERLowError", 
                        "incLinearProd_lowERLowError")&
meanabspbiasdf$model%in%c("autocorr", "rwa")&
meanabspbiasdf$method=="MLE",]


meanabspbiasdf_redux$meanabspbias<-round(meanabspbiasdf_redux$x,2)

meanabspbiasdf_redux$model2<-factor(meanabspbiasdf_redux$model2, levels=c("autocorr", "rw.a"))


meanabspbiasdf_redux$dynamics2<-case_match(meanabspbiasdf_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")


#-----------------------------------------------------------------------------------
#summary plots 

summarydf_smsy_sim_redux<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERLowError", 
                                                        "lowERLowError", 
                                                        "decLinearProd_highERLowError",
                                                        "incLinearProd_highERLowError", 
                                                        "decLinearProd_lowERLowError", 
                                                        "incLinearProd_lowERLowError")&
                                summarydf$model%in%c("autocorr", "rwa")&
                                summarydf$method=="MLE",]

summarydf_smsy_sim_redux$ERtrend<-factor(summarydf_smsy_sim_redux$ERtrend,
    levels=c("highER", "lowER"))

summarydf_smsy_sim_redux$dynamics2<-case_match(summarydf_smsy_sim_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")


summarydf_smsy_redux<-summarydf[summarydf$variable=="mode"&
                                summarydf$parameter=="smsy"&
                                summarydf$scenario%in%c("highERLowError", 
                                                        "lowERLowError",
                                                        "decLinearProd_highERLowError",
                                                        "incLinearProd_highERLowError", 
                                                        "decLinearProd_lowERLowError", 
                                                        "incLinearProd_lowERLowError")&
                                    summarydf$model%in%c("autocorr", "rwa")&
                                    summarydf$method=="MLE",]


summarydf_smsy_redux$ERtrend<-factor(summarydf_smsy_redux$ERtrend,
    levels=c("highER",  "lowER"))

summarydf_smsy_redux$dynamics2<-case_match(summarydf_smsy_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")


meancvdf_redux$labels<-paste("CV = ",meancvdf_redux$meancv)




df_smsy_sim<- df[df$parameter%in%c("smsy")&df$variable=="sim",]
df_smsy_est<- df[df$parameter%in%c("smsy")&
                 df$variable=="mode"&
                 df$scenario%in%c("highERLowError", 
                                  "lowERLowError",
                                  "decLinearProd_highERLowError",
                                  "incLinearProd_highERLowError", 
                                  "decLinearProd_lowERLowError", 
                                  "incLinearProd_lowERLowError")&
                 df$model%in%c("autocorr", "rwa")&
                 df$method=="MLE",]


df_smsy_est$dynamics2<-case_match(df_smsy_est$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")


psmsy_highERscn_hair_cv<-ggplot(df_smsy_est) + 
geom_line(aes(x=by-54,y=value/10000, color=model2, group=iteration), alpha=.1)+
guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1)))+
geom_line(data=summarydf_smsy_sim_redux,aes(x=by-54,y=median/10000),color="black", alpha=.8,linewidth=1.2)+
geom_line(data=summarydf_smsy_redux,aes(x=by-54,y= median/10000),color="darkred", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(2,18))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s")))) + 
xlab("year") +
facet_grid(dynamics2+ERtrend ~model, scales="free_y")+
theme(
    legend.text = element_text(size=18),legend.key.width = unit(2, 'cm'))+
geom_text(data = meancvdf_redux, aes(x=-Inf,y=Inf,hjust=0,
                vjust=1.7,label=labels), size=4.5)
psmsy_highERscn_hair_cv

df_smsy_est_redux<- df[df$parameter=="smsy"&
                    df$variable=="mode"&
                    df$scenario%in%c("highERLowError",  
                                     #"ShiftERLowError", 
                                     "lowERLowError",
                                     "decLinearProd_highERLowError",
                                     "incLinearProd_highERLowError",
                                     #"decLinearProd_ShiftERLowError", 
                                     #"incLinearProd_ShiftERLowError",
                                     "decLinearProd_lowERLowError", 
                                     "incLinearProd_lowERLowError")&
                    df$model%in%c("autocorr", "rwa")&
                    df$method=="MLE",]
df_smsy_est_redux$ERtrend<-factor(df_smsy_est_redux$ERtrend,levels=c("highER", "lowER"))


df_smsy_est_redux$dynamics2<-case_match(df_smsy_est_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")


#boxplots
meanabspbiasdf_redux$labels<-paste("|%bias| =",meanabspbiasdf_redux$meanabspbias)


psmsy_erscn_boxplot_abs<-ggplot(df_smsy_est_redux) + 
geom_boxplot(aes(x=model2,y=abs(bias)/10000, fill=model2),outlier.shape = NA, width=0.8)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(dynamics2 ~ERtrend, scales="free_y")+
 ylab(expression(absolute~bias~"in"~ S[MSY]~("10,000s"))) +
 coord_cartesian(ylim = c(0,5.5))+ 
   xlab("estimation model") +
 geom_text(data = meanabspbiasdf_redux, aes(x=model2,y=Inf,hjust=.5,
                vjust=1.7,label=labels), size=4.5)+
mytheme
psmsy_erscn_boxplot_abs

multi.page.abs.smsy_erbxp <- ggarrange(psmsy_highERscn_hair_cv, psmsy_erscn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")



multi.page.abs.smsy_erbxp_title<-annotate_figure(multi.page.abs.smsy_erbxp, 
    top = text_grob(expression("Sensitivity to exploitation history and time-varying"~log(alpha)), 
               face = "bold" , size = 14))

multi.page.abs.smsy_erbxp_title


ggsave("figures/bias_trends_boxplot_smsy_ERscn.png",
    plot=multi.page.abs.smsy_erbxp_title, width = 15,height = 9, bg = "white")
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_boxplot_smsy_ERscn.png",
#    plot=multi.page.abs.smsy_erbxp_title, width = 15,height = 9, bg = "white")





#=======================================================================================================
#summary plots 
#