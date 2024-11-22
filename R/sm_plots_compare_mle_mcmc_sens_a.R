#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
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

source("R/read_sensa_data.R")
#=================================================================================================================


df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","conv_warning","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))

df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"))


summarydf  <- df %>%
   group_by(scenario,parameter,
    method,model,by,variable, col) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))

summarydf$scenario<-factor(summarydf$scenario,levels=c( "trendLinearProd1" ,
                                                        "trendLinearProd1.3",
                                                        "trendLinearProd2", 
                                                        "trendLinearProd5", 
                                                        "trendLinearProd7",
                                                        "trendLinearProd10" ,    
                                                        "regimeProd1",
                                                        "regimeProd1.3",  
                                                        "regimeProd2" ,     
                                                        "regimeProd5",     
                                                        "regimeProd7",
                                                        "regimeProd10" ))

summarydf$magnitude_ch<-"a 3.7 -> 1.03"
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd1")]<-"a 3.7 -> 1.03"

summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd1.3")]<-"a 3.7 -> 1.35"    
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd1.3")]<-"a 3.7 -> 1.35"    

summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd10")]<-"a 3.7 -> 10"  
summarydf$magnitude_ch[summarydf$scenario%in%c( "regimeProd10")]<-"a 3.7 -> 10"  

summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd2" )]<-"a 3.7 -> 2"   
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd2")]<-"a 3.7 -> 2"   

summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd5")]<-"a 3.7 -> 5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd5")]<-"a 3.7 -> 5" 

summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd7")]<-"a 3.7 -> 7" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd7")]<-"a 3.7 -> 7"               

summarydf$magnitude_ch<-factor(summarydf$magnitude_ch,levels=c("a 3.7 -> 1.03",
       "a 3.7 -> 1.35", 
       "a 3.7 -> 2",
       "a 3.7 -> 5",
       "a 3.7 -> 7", 
       "a 3.7 -> 10"))

summarydf$trendtype<-case_match(summarydf$scenario, "trendLinearProd1" ~"Trend",
                                                        "trendLinearProd1.3"~"Trend",
                                                        "trendLinearProd2"~"Trend", 
                                                        "trendLinearProd5"~"Trend", 
                                                        "trendLinearProd7"~"Trend",
                                                        "trendLinearProd10"~"Trend" ,    
                                                        "regimeProd1"~"Shift",
                                                        "regimeProd1.3"~"Shift",  
                                                        "regimeProd2"~"Shift",     
                                                        "regimeProd5"~"Shift",     
                                                        "regimeProd7"~"Shift",
                                                        "regimeProd10"~"Shift" )



summarydf <- reshape2::dcast(data=summarydf,  scenario + parameter + method + model + by +
 variable  + magnitude_ch + trendtype ~prob, 
    value.var= "qs",fun.aggregate=mean)

summarydf_alpha_sim<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="sim",]

summarydf_alpha<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="mode",]


p_sensa_alpha <- ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by-54,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim,aes(x=by-54,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(trendtype+magnitude_ch~model, scales="free_y")
p_sensa_alpha

ggsave("figures/compareMCMC_MLE_sensa_alpha.png",
    plot=p_sensa_alpha, width = 15,height = 18 )



p_sensa_alpha_zoom <- ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by-54,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim,aes(x=by-54,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.1,2.5))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(trendtype+magnitude_ch~model, scales="free_y")
p_sensa_alpha_zoom
ggsave("figures/compareMCMC_MLE_sensa_alpha_zoom.png",
    plot=p_sensa_alpha_zoom, width = 15,height = 18 )


#smax
summarydf_smax_sim<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="sim",]

summarydf_smax<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="mode",]


p_sensa_smax<-ggplot() + 
geom_pointrange(data=summarydf_smax,aes(x=by-54, y=median/1000 ,ymin =lower /1000,
 ymax = upper/1000, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim,aes(x=by-54,y= median/1000),color="black", 
    alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,370000)/1000)+ 
mytheme+ 
ylab(expression(S[max]~"(1000s)")) +
xlab("year") +
facet_grid(trendtype+magnitude_ch~model, scales="free_y")
p_sensa_smax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smax.png", 
    plot=p_sensa_smax, width = 15,height = 18)



#=======================================================================
#smsy
summarydf_smsy_sim<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="sim",]

summarydf_smsy<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="mode",]


p_sensa_smsy<-ggplot() + 
geom_pointrange(data=summarydf_smsy,aes(x=by-54,y=median/1000 ,ymin =lower/1000 , ymax = upper/1000, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim,aes(x=by-54,y= median/1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+ 
ylab(expression(S[MSY]~"(1000s)")) +
xlab("year") +
facet_grid(trendtype+magnitude_ch~model, scales="free_y")
p_sensa_smsy
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy.png", 
#    plot=p_sensa_smsy, width = 15,height = 18)
ggsave("figures/compareMCMC_MLE_sensa_smsy.png", 
    plot=p_sensa_smsy, width = 15,height = 18)


