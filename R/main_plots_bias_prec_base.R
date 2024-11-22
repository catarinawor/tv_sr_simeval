#============================================
#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================

library(ggplot2)
library(gridExtra)
library(dplyr)
source("R/utils.R")
#source("code/cluster_func_plots.R")

library("ggpubr")

mytheme = list(
    theme_classic(13)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"), legend.title = element_blank(),
            legend.position="bottom", strip.text = element_text( size=13),
            axis.text=element_text(face="bold",size=13),axis.title = element_text(face="bold",size=13),
            legend.text=element_text(size=13),
            plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

#read in data
source("R/read_base_data.R")

#========================================================================================================
#base case
#read in data

df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", 
                                      "convergence","conv_warning","pbias","bias"))


df$scenario<-factor(df$scenario,levels=c("stationary",
                                        "autocorr",
                                        "sigmaShift",
                                        "decLinearProd",
                                        "sineProd", 
                                        "regimeProd", 
                                        "shiftProd",                       
                                        "decLinearCap",
                                        "regimeCap",
                                        "shiftCap",                      
                                        "regimeProdCap",         
                                        "decLinearProdshiftCap"  ))

df$scentype<-dplyr::case_match(df$scenario, 
      "stationary"~"stationary",
      "autocorr"~"autocorrelation",
      "sigmaShift"~"sigma", 
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

df$scentype<-factor(df$scentype, levels=c("stationary","autocorrelation","sigma","log(alpha)",
    "Smax","both" )  ) 

df$scendesc<-dplyr::case_match(df$scenario, 
      "stationary"~".",
      "autocorr"~".",
      "sigmaShift"~"shift up", 
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

df$scendesc<-factor(df$scendesc, levels=c(".","shift up","linear decline",
    "sine trend","shift up + down","shift down + up","shift down",
    "regime both","trend & shift") ) 
     
df$scentrend<-dplyr::case_match(df$scenario, 
      "stationary"~"stationary",
      "autocorr"~"autocorr",
      "sigmaShift"~"sigma shift", 
      "decLinearProd"~"trend",
      "sineProd"~"trend",
      "regimeProd"~"shift",
      "shiftProd"~"shift",
      "decLinearCap"~"trend",
      "regimeCap"~"shift",
      "shiftCap"~"shift", 
      "regimeProdCap"~"shift",
      "decLinearProdshiftCap"~"trend & shift"
      )   


df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"  ))


df$model2<-dplyr::case_match(df$model, 
     "simple"~"stationary",
     "autocorr"~"autocorr", 
     "rwa"~"rw.a",
     "hmma"~"hmm.a",
     "rwb"~"rw.b",
     "hmmb"~"hmm.b",
     "rwab"~"rw.ab",
     "hmmab"~"hmm.ab" 
      )   


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
    method,model,model2,by,variable,scencode,scentype,scentrend,scendesc) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))

summarydf <- reshape2::dcast(data=summarydf,  
    scenario + parameter + method + model +model2 + by + variable + scencode + scentype +scentrend +scendesc~ prob, 
    value.var= "qs",fun.aggregate=mean)

#==================================================================================
#mean CV estimates
cvdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})

meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               method=cvdf$method,
                               model=cvdf$model
                              ), function(x){mean(x,na.rm=T)})

meancvdf$model2<-case_match(meancvdf$model,
    "simple"~"stationary",
    "autocorr"~"autocorr", 
    "rwa"~ "rw.a",
    "hmma" ~"hmm.a",
    "rwb" ~"rw.b",
    "hmmb" ~"hmm.b", 
    "rwab"~"rw.ab",
    "hmmab"~ "hmm.ab" )

meancvdf$scendesc<-dplyr::case_match(meancvdf$scenario, 
      "stationary"~".",
      "autocorr"~".",
      "sigmaShift"~"shift up", 
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

meancvdf$scendesc <- factor(meancvdf$scendesc, levels=c(".","shift up",
    "linear decline",
    "sine trend",
    "shift up + down",
    "shift down + up",
    "shift down",
    "regime both",
    "trend & shift") ) 
     
meancvdf$scentype<-dplyr::case_match(meancvdf$scenario, 
       "stationary"~"stationary",
      "autocorr"~"autocorrelation",
      "sigmaShift"~"sigma", 
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

meancvdf$scentype<-factor(meancvdf$scentype, levels=c("stationary","autocorrelation",
    "sigma","log(alpha)","Smax","both" )  ) 

meancvdf$labels<-round(meancvdf$x,2)

#mean absolute percent bias
abspbiasdf <- aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){median(abs(x),na.rm=T)})

meanabspbiasdf <- aggregate(abspbiasdf$x, list(parameter=abspbiasdf$parameter,
                               scenario=abspbiasdf$scenario,
                               method=abspbiasdf$method,
                               model=abspbiasdf$model
                              ), function(x){mean(x,na.rm=T)})
         
meanabspbiasdf$model2<-case_match(meanabspbiasdf$model,
    "simple"~"stationary",
    "autocorr"~"autocorr", 
    "rwa"~ "rw.a",
    "hmma" ~"hmm.a",
    "rwb" ~"rw.b",
    "hmmb" ~"hmm.b", 
    "rwab"~"rw.ab",
    "hmmab"~ "hmm.ab" )

#===================================================================================
#=========================================================================================================\
# smsy-focusedredux plots
# alpha varies
meancvdf_redux_alpha<-meancvdf[meancvdf$parameter=="smsy"&
meancvdf$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
meancvdf$model%in%c("autocorr", "rwa", "hmma")&
meancvdf$method=="MLE",]


meancvdf_redux_alpha$meancv<-round(meancvdf_redux_alpha$x,2)
meancvdf_redux_alpha$scenario2<-case_match(meancvdf_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

meancvdf_redux_alpha$model2<-factor(meancvdf_redux_alpha$model2, levels=c("autocorr", "rw.a", "hmm.a"))

meancvdf_redux_alpha$scenario2<-factor(meancvdf_redux_alpha$scenario2,levels=c("stationary",
    "linear decline" ,"sine fluctuation","shift decline" ))

#absolute median percent bias
abspbiasdf<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){median(abs(x),na.rm=T)})

meanabspbiasdf<-aggregate(abspbiasdf$x, list(parameter=abspbiasdf$parameter,
                               scenario=abspbiasdf$scenario,
                               method=abspbiasdf$method,
                               model=abspbiasdf$model
                              ), function(x){mean(x,na.rm=T)})

meanabspbiasdf$model2<-case_match(meanabspbiasdf$model,
     "simple"~"stationary",
    "autocorr"~"autocorr", 
     "rwa"~ "rw.a",
     "hmma" ~"hmm.a",
      "rwb" ~"rw.b",
      "hmmb" ~"hmm.b", 
       "rwab"~"rw.ab",
        "hmmab"~ "hmm.ab" )

meanabspbiasdf_redux_alpha<-meanabspbiasdf[meanabspbiasdf$parameter=="smsy"&
meanabspbiasdf$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
meanabspbiasdf$model%in%c("autocorr", "rwa", "hmma")&
meanabspbiasdf$method=="MLE",]

meanabspbiasdf_redux_alpha$meanabspbias<-round(meanabspbiasdf_redux_alpha$x,2)
meanabspbiasdf_redux_alpha$scenario2<-case_match(meanabspbiasdf_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")
meanabspbiasdf_redux_alpha$model2<-factor(meanabspbiasdf_redux_alpha$model2, levels=c("autocorr", "rw.a", "hmm.a"))

meanabspbiasdf_redux_alpha$scenario2<-factor(meanabspbiasdf_redux_alpha$scenario2,levels=c("stationary",
    "linear decline" ,"sine fluctuation","shift decline" ))

#estimates
df_smsy_sim<- df[df$parameter%in%c("smsy")&df$variable=="sim",]
df_smsy_est<- df[df$parameter%in%c("smsy")&df$variable=="mode",]


summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model2=df_smsy_est$model2,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)

summarydf_smsy_sim<-aggregate(df_smsy_sim$value,by=list(scenario=df_smsy_sim$scenario, 
    method=df_smsy_sim$method, 
    model2=df_smsy_sim$model2,
    by=df_smsy_sim$by ),
    function(x) {unique(x)})

summarydf_smsy_sim_redux_alpha<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
summarydf_smsy_sim$model2%in%c("autocorr", "rw.a", "hmm.a")&summarydf_smsy_sim$method=="MLE",]

summarydf_smsy_sim_redux_alpha$scenario2<-case_match(summarydf_smsy_sim_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

summarydf_smsy_sim_redux_alpha$scenario2<-factor(summarydf_smsy_sim_redux_alpha$scenario2, 
    levels=c("stationary",
     "linear decline" ,
     "sine fluctuation",
     "shift decline" ))

summarydf_smsy_redux_alpha<-summarydf_smsy[summarydf_smsy$scenario%in%c("autocorr", "decLinearProd", "shiftProd", "sineProd" )&
  summarydf_smsy$model2%in%c("autocorr", "rw.a", "hmm.a") &summarydf_smsy$method=="MLE",]


summarydf_smsy_redux_alpha$scenario2<-case_match(summarydf_smsy_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

summarydf_smsy_redux_alpha$scenario2<-factor(summarydf_smsy_redux_alpha$scenario2, 
    levels=c("stationary",
     "linear decline",
     "sine fluctuation",
     "shift decline" ))

meancvdf_redux_alpha$labels<-paste("CV = ",meancvdf_redux_alpha$meancv)

#line plot
df_smsy_est_redux_alpha<-df_smsy_est[
           df_smsy_est$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
           df_smsy_est$model2%in%c("autocorr", "rw.a", "hmm.a")&
           df_smsy_est$method=="MLE",]

df_smsy_est_redux_alpha$scenario2<-case_match(df_smsy_est_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

df_smsy_est_redux_alpha$scenario2<-factor(df_smsy_est_redux_alpha$scenario2, 
    levels=c("stationary",
     "linear decline",
     "sine fluctuation",
     "shift decline" ))

df_smsy_sim_redux_alpha<-df_smsy_sim[
           df_smsy_sim$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
           df_smsy_sim$model2%in%c("autocorr", "rw.a", "hmm.a")&
           df_smsy_sim$method=="MLE",]

df_smsy_est_redux_alpha[df_smsy_est_redux_alpha$iteration==1&
           df_smsy_est_redux_alpha$model=="autocorr",]


psmsy_alphascn_hair<-ggplot(df_smsy_est_redux_alpha) + 
geom_line(aes(x=by-54,y=value/10000, color=model2, group=iteration), alpha=.1)+
guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1)))+
geom_line(data=summarydf_smsy_sim_redux_alpha,aes(x=by-54,y= x/10000),color="black", alpha=.8,linewidth=1.2)+
geom_line(data=summarydf_smsy_redux_alpha,aes(x=by-54,y= x.50./10000),color="darkred", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
 facet_grid(scenario2~model2, scales="free_y")+
coord_cartesian(ylim = c(2,19))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s")))) +
xlab("year") +
mytheme+
theme( strip.text.y = element_blank(),strip.text.x = element_blank(),
    legend.text = element_text(size=18),legend.key.width = unit(2, 'cm'))+
geom_text(data = meancvdf_redux_alpha, aes(x=-Inf,y=Inf,hjust=-0.05,
                vjust=1.7,label=labels), size=4.5)
psmsy_alphascn_hair


df_smsy_est_redux_alpha<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]


df_smsy_est_redux_alpha$scenario2<-case_match(df_smsy_est_redux_alpha$scenario,
    "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")


df_smsy_est_redux_alpha$scenario2<-factor(df_smsy_est_redux_alpha$scenario2, levels=c("stationary",
     "linear decline",
     "sine fluctuation",
     "shift decline" ))

meanabspbiasdf_redux_alpha$labels<-paste("|%bias| = ",meanabspbiasdf_redux_alpha$meanabspbias)

psmsy_alphascn_boxplot_abs<-ggplot(df_smsy_est_redux_alpha) + 
geom_boxplot(aes(x=model2,y=abs(bias)/10000, fill=model2), outlier.shape = NA, width=0.8)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(absolute~bias~"in"~ S[MSY]~("10,000s"))) +
  coord_cartesian(ylim = c(0,5))+ 
   xlab("estimation model") +
 geom_text(data = meanabspbiasdf_redux_alpha, aes(x=model2,y=Inf,hjust=0.5,
                vjust=1.7,label=labels), size=4.5)+
mytheme
psmsy_alphascn_boxplot_abs

#hair plot version


multi.page.abs.smsy_alphascn.hairbx <- ggarrange(psmsy_alphascn_hair, psmsy_alphascn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")


multi.page.abs.smsy_alphascn_title.hairbx<-annotate_figure(multi.page.abs.smsy_alphascn.hairbx, 
    top = text_grob(expression("Stationary and time-varying"~log(alpha)), 
               face = "bold" , size = 14))
multi.page.abs.smsy_alphascn_title.hairbx

ggsave("figures/bias_trends_smsy_basealpha_scn_hairbx.png",
    plot=multi.page.abs.smsy_alphascn_title.hairbx,width = 14,height = 8, bg = "white")
#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_smsy_basealpha_scn_hairbx.png",
#    plot=multi.page.abs.smsy_alphascn_title.hairbx,width = 14,height = 8, bg = "white")



#==================================================================
# Bias when beta or both parameters change

meancvdf_redux_smaxboth<-meancvdf[meancvdf$parameter=="smsy"&
meancvdf$scenario%in%c("decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap")&
meancvdf$model%in%c("autocorr","rwb", "rwab", "hmmb", "hmmab")&
meancvdf$method=="MLE",]

meancvdf_redux_smaxboth$meancv<-round(meancvdf_redux_smaxboth$x,2)
meancvdf_redux_smaxboth$scenario2<-case_match(meancvdf_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")
meancvdf_redux_smaxboth$model2<-factor(meancvdf_redux_smaxboth$model2, levels=c("autocorr","rw.b", "rw.ab","hmm.b", "hmm.ab"))

meancvdf_redux_smaxboth$scenario2<-factor(meancvdf_redux_smaxboth$scenario2,levels=c("linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both" ))


#estimates

summarydf_smsy_sim_redux_smaxboth<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("decLinearCap", 
      "shiftCap","regimeProdCap","decLinearProdshiftCap")&
         summarydf_smsy_sim$model2%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab")&
         summarydf_smsy_sim$method=="MLE",]

summarydf_smsy_sim_redux_smaxboth$scenario2<-case_match(summarydf_smsy_sim_redux_smaxboth$scenario,
    "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")


summarydf_smsy_sim_redux_smaxboth$scenario2<-factor(summarydf_smsy_sim_redux_smaxboth$scenario2, levels=c(
    "linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))


summarydf_smsy_redux_smaxboth<-summarydf_smsy[summarydf_smsy$scenario%in%c( "decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap" )&
  summarydf_smsy$model%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab") &summarydf_smsy$method=="MLE",]
#"autocorr"

summarydf_smsy_redux_smaxboth$scenario2<-case_match(summarydf_smsy_redux_smaxboth$scenario,
    "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")

summarydf_smsy_redux_smaxboth$scenario2<-factor(summarydf_smsy_redux_smaxboth$scenario2, levels=c(
    "linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))


meancvdf_redux_smaxboth$labels<-paste("CV = ",meancvdf_redux_smaxboth$meancv)



#hair plot
df_smsy_est_redux_smaxboth<-df_smsy_est[df_smsy_est$scenario%in%c("decLinearCap", 
      "shiftCap","regimeProdCap","decLinearProdshiftCap")&
df_smsy_est$model2%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab")&
df_smsy_est$method=="MLE",]


df_smsy_est_redux_smaxboth$scenario2<-case_match(df_smsy_est_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")


df_smsy_est_redux_smaxboth$scenario2<-factor(df_smsy_est_redux_smaxboth$scenario2, 
    levels=c("linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))


df_smsy_sim_redux_smaxboth<-df_smsy_sim[df_smsy_sim$scenario%in%c("decLinearCap", 
      "shiftCap","regimeProdCap","decLinearProdshiftCap")&
df_smsy_sim$model2%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab")&
df_smsy_sim$method=="MLE",]


psmsy_smaxbothscn_hair<-ggplot(df_smsy_est_redux_smaxboth) + 
geom_line(aes(x=by-54,y=value/10000, color=model2, group=iteration), alpha=.1)+
guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1)))+
geom_line(data=summarydf_smsy_sim_redux_smaxboth,aes(x=by-54,y= x/10000),color="black", alpha=.8,linewidth=1.2)+
geom_line(data=summarydf_smsy_redux_smaxboth,aes(x=by-54,y= x.50./10000),color="darkred", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
 facet_grid(scenario2~model2, scales="free_y")+
coord_cartesian(ylim = c(2,15))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s")))) +
xlab("year") +
theme(strip.text.y = element_blank(),strip.text.x = element_blank(),
    legend.text = element_text(size=18),legend.key.width = unit(2, 'cm'))+
geom_text(data = meancvdf_redux_smaxboth, aes(x=-Inf,y=Inf,hjust=-0.05,
                vjust=1.7,label=labels), size=4.5)
psmsy_smaxbothscn_hair



#calculate pbias

meanabspbiasdf_redux_smaxboth<-meanabspbiasdf[meanabspbiasdf$parameter=="smsy"&
meanabspbiasdf$scenario%in%c("decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap")&
meanabspbiasdf$model%in%c("autocorr","rwb", "rwab","hmmb", "hmmab")&
meanabspbiasdf$method=="MLE",]


meanabspbiasdf_redux_smaxboth$meanabspbias<-round(meanabspbiasdf_redux_smaxboth$x,2)
meanabspbiasdf_redux_smaxboth$scenario2<-case_match(meanabspbiasdf_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")

meanabspbiasdf_redux_smaxboth$model2<-factor(meanabspbiasdf_redux_smaxboth$model2, 
    levels=c("autocorr","rw.b", "hmm.b", "rw.ab", "hmm.ab"))

meanabspbiasdf_redux_smaxboth$scenario2<-factor(meanabspbiasdf_redux_smaxboth$scenario2,
    levels=c("linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))

df_smsy_est_redux_smaxboth<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap")&
df$model%in%c("autocorr","rwb", "rwab","hmmb", "hmmab")&
df$method=="MLE",]
#"autocorr"

df_smsy_est_redux_smaxboth$scenario2<-case_match(df_smsy_est_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")

df_smsy_est_redux_smaxboth$scenario2<-factor(df_smsy_est_redux_smaxboth$scenario2, levels=c(
    "linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))






meanabspbiasdf_redux_smaxboth$labels<-paste("|%bias| =",meanabspbiasdf_redux_smaxboth$meanabspbias)

psmsy_smaxscn_boxplot_abs<-ggplot(df_smsy_est_redux_smaxboth) + 
geom_boxplot(aes(x=model2,y=abs(bias)/10000, fill=model2), outlier.shape = NA, width=0.8)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(absolute~bias~"in"~ S[MSY]~("10,000s"))) +
  xlab("estimation model") +
  coord_cartesian(ylim = c(0,5))+ 
 geom_text(data = meanabspbiasdf_redux_smaxboth, aes(x=model2,y=Inf,hjust=.6,
                vjust=1.7,label=labels), size=4.5)+
mytheme
psmsy_smaxscn_boxplot_abs

#hair plot version
multi.page.abs.smsy_smaxbothscn.hairbx <- ggarrange(psmsy_smaxbothscn_hair, psmsy_smaxscn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")

multi.page.abs.smsy_smaxbothscn_title.hairbx<-annotate_figure(multi.page.abs.smsy_smaxbothscn.hairbx, 
    top = text_grob(expression("Time-varying"~S[max]~"or both parameters"), 
               face = "bold" , size = 14))
multi.page.abs.smsy_smaxbothscn_title.hairbx

#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_boxplot_smsy_smaxbothscn_scn_hairbx.png",
#    plot=multi.page.abs.smsy_smaxbothscn_title.hairbx,width = 22,height = 8, bg = "white")
ggsave("figures/bias_trends_boxplot_smsy_smaxbothscn_scn_hairbx.png",
    plot=multi.page.abs.smsy_smaxbothscn_title.hairbx,width = 22,height = 8, bg = "white")


