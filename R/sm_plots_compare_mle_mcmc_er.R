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

summarydf_alpha_sim1<-summarydf[summarydf$parameter=="logalpha"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("decLinearProd_highERLowError",
                                  "decLinearProd_ShiftERLowError", 
                                  "incLinearProd_highERLowError",  
                                  "incLinearProd_ShiftERLowError",
                                  "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                                ]

summarydf_alpha1<-summarydf[summarydf$parameter=="logalpha"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                            ]

scenlab1<-c("highER","ShiftER","lowER", "highER", "ShiftER", "lowER")
names(scenlab1) <- c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "incLinearProd_lowERLowError")


er_alpha1<-ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by-54,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.3,2.4))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab1))
er_alpha1
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha_ERtrendloga.png",
#    plot=er_alpha1, width = 15,height = 8)
ggsave("figures/compareMCMC_MLE_alpha_ERtrendloga.png",
    plot=er_alpha1, width = 15,height = 8)


summarydf_alpha_sim2<-summarydf[summarydf$parameter=="logalpha"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError",
                                          "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                                ]

summarydf_alpha2<-summarydf[summarydf$parameter=="logalpha"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError",
                                          "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                            ]

                                      
scenlab2<-c("highER","lowER", "ShiftER","trendER","highER","lowER", "ShiftER","trendER")
names(scenlab2) <- c( "highERHighError",                              
                      "lowERHighError",                           
                       "ShiftERHighError",                                                                   
                       "trendERHighError",
                       "highERLowError", 
                      "lowERLowError",  
                       "ShiftERLowError",      
                       "trendERLowError")

er_alpha2<-ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by-54,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.3,2.4))+ 
mytheme + 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(ERerror+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab2))
er_alpha2
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha_ERstationary.png",
#    plot=er_alpha2, width = 15,height = 12)
ggsave("figures/compareMCMC_MLE_alpha_ERstationary.png",
    plot=er_alpha2, width = 15,height = 12)

#=======================================================
#smax estimates
summarydf_smax_sim1<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("decLinearProd_highERLowError",
                                  "decLinearProd_ShiftERLowError", 
                                  "incLinearProd_highERLowError",  
                                  "incLinearProd_ShiftERLowError",
                                  "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                                ]

summarydf_smax1<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                            ]

scenlab1<-c("highER","ShiftER","lowER", "highER", "ShiftER", "lowER")
names(scenlab1) <- c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "incLinearProd_lowERLowError")

er_smax1<-ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by-54,y=median/1000,ymin = lower/1000, ymax = upper/1000,  col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by-54,y= median/1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(paste(S[max]~"(1000s)"))) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000)/1000)+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab1))
er_smax1
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax_ERtrendloga.png",
#    plot=er_smax1, width = 15,height = 8)
ggsave("figures/compareMCMC_MLE_smax_ERtrendloga.png",
    plot=er_smax1, width = 15,height = 8)

summarydf_smax_sim2<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError",
                                          "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                                ]

summarydf_smax2<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError",
                                          "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                            ]


er_smax2<-ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by-54,y=median/1000,ymin = lower/1000, ymax = upper/1000, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by-54,y=median/1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(paste(S[max]~"(1000s)"))) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000)/1000)+ 
facet_grid(ERerror+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smax2
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax_ERstationary.png",
#    plot=er_smax2, width = 15,height = 12)
ggsave("figures/compareMCMC_MLE_smax_ERstationary.png",
    plot=er_smax2, width = 15,height = 12)



#=======================================================
#smsy estimates



summarydf_smsy_sim1<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("decLinearProd_highERLowError",
                                  "decLinearProd_ShiftERLowError", 
                                  "incLinearProd_highERLowError",  
                                  "incLinearProd_ShiftERLowError",
                                  "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                                ]

summarydf_smsy1<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                            ]



er_smsy1<-ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by-54, y=median/1000, ymin = lower/1000, 
    ymax = upper/1000, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by-54,y= median/1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(paste(S[MSY]~"(1000s)"))) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000)/1000)+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab1))
er_smsy1
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy_ERtrendloga.png",
#    plot=er_smsy1, width = 15,height = 8)
ggsave("figures/compareMCMC_MLE_smsy_ERtrendloga.png",
    plot=er_smsy1, width = 15,height = 8)


summarydf_smsy_sim2<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError",
                                          "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                                ]

summarydf_smsy2<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError",
                                          "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                            ]

er_smsy2<-ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by-54,y=median/1000,ymin = lower/1000, ymax = upper/1000,
  col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by-54,y= median/1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(paste(S[MSY]~"(1000s)"))) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000)/1000)+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smsy2
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy_ERstationary.png",
#    plot=er_smsy2, width = 15,height = 12)
ggsave("figures/compareMCMC_MLE_smsy_ERstationary.png",
    plot=er_smsy2, width = 15,height = 12)

