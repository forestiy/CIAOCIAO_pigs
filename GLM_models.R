rm(list=ls()) 

#Loading data
load("FDF.Rdata");load("AB_Classes_DF_long.Rdata")


#loading packages
library(haven)
library(tibble)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(SciViews)
library(data.table)
library(car)
library(ggpubr)
library(Hmisc)
#library(broom)
library(broom.mixed)
library(lme4)
library(lmerTest)
library(forcats)
library(epiDisplay)

#loading few functions
lnn<-function(x){ln(x+((x^2) + 1)^(1/2))}
un_ln <- function (x) {(exp(2*(x))-1)/(2*exp((x)))}
`%notin%` <- Negate(`%in%`)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#____________________START________LMER_BW_for_SOWS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ALLRESULTS<-list()

AIC_lmer<-NULL
AIC_val<-NULL
mod_test<-NULL
model_bw_Weaners_full<-NULL
model_bw_Weaners_R1<-NULL
f_lmer<-NULL
AIC_val_intr<-NULL
md_intr<-NULL
AIC_intr<-NULL
#Sows lmer
#######

FDFs<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),]

DS<-colnames(FDF[which(grepl("__Bin",colnames(FDF),fixed=T) & grepl("Ze|Zu",colnames(FDF),fixed=F))])
DS<-DS[which(DS %notin% c("Zeugen_NotmentioningABclass__Bin",                        
                          "Zeugen_Not_mentioningagegroup__Bin",   
                          "Zeugen__DIGESTIETRACTUS____Bin",
                          "Zeugen__HUIDenWONDINFECTIES____Bin",
                          "Zeugen__RESPIRATIETRACTUS____Bin" ,
                          "ZeugenGeneral_Indications__Bin" ,
                          "Zuigendebiggen__RESPIRATIETRACTUS____Bin",
                          "Zuigendebiggen__SPECIFIEKESYSTEMISCHEAANDOENINGEN____Bin"))]

f_lmer<-as.formula(paste0(paste0("DDDA_sows_piglets_tr_"," ~ "), paste0(DS,collapse=" * "),'+',
                          "  Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#



model_bw_Sows_full<-lmer(f_lmer,
                         REML   = FALSE, 
                         control= lmerControl(optimizer="bobyqa"),
                         data   = FDFs)

DSintr<-rownames(coef(summary(model_bw_Sows_full)))[-c(1,6,10:16)]#,10
DSintr<-DSintr[c(4,1:3,5:length(DSintr))]

mod_test_lmer<-NULL
AIC_val_lmer<-NULL
for (u in 1:length(DSintr)){
  f_lmer<- as.formula(paste0(paste0("DDDA_sows_piglets_tr_"," ~ "), paste0( DSintr[-u],collapse=" + "),
                             "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
  
  model_bw_Sows_R1<-lmer(f_lmer,
                         REML   = FALSE, 
                         control= lmerControl(optimizer="bobyqa"),
                         data   = FDFs)
  
  
  mod_test_lmer<- anova(model_bw_Sows_R1,model_bw_Sows_full)
  
  
  
  AIC_val_lmer<-append(AIC_val_lmer,AIC(model_bw_Sows_R1))
  
  
  
  
}
md_lmer<- which(AIC_val_lmer==min(AIC_val_lmer, na.rm = T))
AIC_lmer<- (AIC_val_lmer[md_lmer] - AIC(model_bw_Sows_full))


while(AIC_lmer<0 ){#& abs(AIC_lmer)>=1.9
  
  ##exclude interaction when main effect is off
  if((grepl(":",DSintr[md_lmer], fixed=T)==F)){
    md_lmer<-c(md_lmer,which(grepl(":",DSintr, fixed=T) & grepl(DSintr[md_lmer],DSintr, fixed=T)))
  }else{md_lmer<-md_lmer}
  
  print(paste("OUT_lmer->",DSintr[md_lmer]))
  
  
  model_bw_Sows_full<-update(model_bw_Sows_full,formula=as.formula(paste0(paste0("DDDA_sows_piglets_tr_"," ~ "), paste0( DSintr[-md_lmer],collapse=" + "),
                                                                          "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)")))#
  
  
  
  DSintr<-DSintr[-md_lmer]
  AIC_val<-c()
  
  for (u in 1:length(DSintr)){
    f_lmer<- as.formula(paste0(paste0("DDDA_sows_piglets_tr_"," ~ "), paste0( DSintr[-u],collapse=" + "),
                               "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
    model_bw_Sows_R1<-update(model_bw_Sows_R1,formula=f_lmer)
    
    mod_test<- anova(model_bw_Sows_R1,model_bw_Sows_full)
    
    AIC_val<-append(AIC_val,AIC(model_bw_Sows_R1))
    
  }
  md_lmer<- which(AIC_val==min(AIC_val, na.rm = T))
  AIC_lmer<- (AIC_val[md_lmer] - AIC(model_bw_Sows_full))
  
  #print(paste("OUT_lmer2->",DSintr[md_lmer]))
  
}
######
model_bw_Sows_full<-update(model_bw_Sows_full, REML   = T)
summary(model_bw_Sows_full)
#library(MuMIn)
#AICc(model_bw_Sows_full)

kl<-data.frame(coef(summary(model_bw_Sows_full)))
kl$names<-rownames(kl)

kl<-cbind(kl,confint.merMod(model_bw_Sows_full, method="profile")[-c(1:3),])
rownames(kl)<-NULL

kl<-kl[,c(6,1,7:8,2:5)]
print(kl)
kl$Estim_Scale<-un_ln(kl$Estimate)
kl$Estim_lo_Scale<-un_ln(kl$`2.5 %`)
kl$Estim_hi_Scale<-un_ln(kl$`97.5 %`)
ALLRESULTS[[length(ALLRESULTS)+1]]<-kl

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#_________________________END___LMER_BW_for_SOWS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________START______LMER_BW_for_WEANERS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#Bw selection for diseases and AMU in weaners
#Weaners lmer
#############
md_lmer<-NULL
AIC_lmer<-NULL
AIC_val<-NULL
mod_test<-NULL
model_bw_Weaners_full<-NULL
model_bw_Weaners_R1<-NULL
f_lmer<-NULL

FDFw3<-FDF[-which(FDF$Bed_score=="-None-"),]

DS<-colnames(FDF[which(grepl("__Bin",colnames(FDF),fixed=T) & grepl("Ges|Wea",colnames(FDF),fixed=F))])

DS<-DS[which(DS %notin% c("Gespeendebiggen__HUIDenWONDINFECTIES____Bin",
                          "Gespeendebiggen__SPECIFIEKESYSTEMISCHEAANDOENINGEN____Bin",
                          "GespeendebiggenGeneral_Indications__Bin" ,
                          "Weaners_Not_mentioningagegroup__Bin",
                          "Weaners_NotmentioningABclass__Bin" ))]


f_lmer<-as.formula(paste0(paste0("DDDA_weaners_tr_"," ~ "), paste0(DS,collapse=" * "),'*',
                          "Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))



model_bw_Weaners_full<-lmer(f_lmer,
                            REML   = FALSE, 
                            control= lmerControl(optimizer="bobyqa"),
                            data   = FDFw3)

DSintr<-rownames(coef(summary(model_bw_Weaners_full)))[-c(1,6,10:28)]

mod_test_lmer<-NULL
AIC_val_lmer<-NULL
for (u in 1:length(DSintr)){
  f_lmer<- as.formula(paste0(paste0("DDDA_weaners_tr_"," ~ "), paste0( DSintr[-u],collapse=" + "),
                             " + Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
  
  model_bw_Weaners_R1<-lmer(f_lmer,
                            REML   = FALSE, 
                            control= lmerControl(optimizer="bobyqa"),
                            data   = FDFw3)
  
  
  mod_test_lmer<- anova(model_bw_Weaners_R1,model_bw_Weaners_full)
  
  
  
  AIC_val_lmer<-append(AIC_val_lmer,mod_test_lmer$AIC[1])
  
  
  
  
}
md_lmer<- which(AIC_val_lmer==min(AIC_val_lmer, na.rm = T))
AIC_lmer<- (AIC_val_lmer[md_lmer] - AIC(model_bw_Weaners_full))


while(AIC_lmer<0 ){
  
  ##exclude interaction when main effect is off
  if((grepl(":",DSintr[md_lmer], fixed=T)==F)){
    md_lmer<-c(md_lmer,which(grepl(":",DSintr, fixed=T) & grepl(DSintr[md_lmer],DSintr, fixed=T)))
  }else{md_lmer<-md_lmer}
  
  print(paste("OUT_lmer->",DSintr[md_lmer]))
  
  
  model_bw_Weaners_full<-update(model_bw_Weaners_full,formula=as.formula(paste0(paste0("DDDA_weaners_tr_"," ~ "), paste0( DSintr[-md_lmer],collapse=" + "),
                                                                                "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)")))#
  
  
  
  DSintr<-DSintr[-md_lmer]
  AIC_val<-c()
  
  for (u in 1:length(DSintr)){
    f_lmer<- as.formula(paste0(paste0("DDDA_weaners_tr_"," ~ "), paste0( DSintr[-u],collapse=" + "),
                               " + Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
    model_bw_Weaners_R1<-update(model_bw_Weaners_R1,formula=f_lmer)
    
    mod_test<- anova(model_bw_Weaners_R1,model_bw_Weaners_full)
    
    AIC_val<-append(AIC_val,mod_test$AIC[1])
    
  }
  md_lmer<- which(AIC_val==min(AIC_val, na.rm = T))
  AIC_lmer<- (AIC_val[md_lmer] - AIC(model_bw_Weaners_full))
  
  #print(paste("OUT_lmer->",DSintr[md_lmer]))
  
}
model_bw_Weaners_full<-update(model_bw_Weaners_full, REML   = T)

summary(model_bw_Weaners_full)



kl<-data.frame(coef(summary(model_bw_Weaners_full)))
kl$names<-rownames(kl)

kl<-cbind(kl,confint.merMod(model_bw_Weaners_full, method="profile")[-c(1:3),])
rownames(kl)<-NULL

kl<-kl[,c(6,1,7:8,2:5)]
kl$Estim_Scale<-un_ln(kl$Estimate)
kl$Estim_lo_Scale<-un_ln(kl$`2.5 %`)
kl$Estim_hi_Scale<-un_ln(kl$`97.5 %`)
# 
# kl$mean_resp<-ifelse(kl$names=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin:Gespeendebiggen__RESPIRATIETRACTUS____Bin",(sum(kl$Estimate[c(3,4,7)])),
#                      (kl$Estimate))
# 
# 
# kl$low_resp<-ifelse(kl$names=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin:Gespeendebiggen__RESPIRATIETRACTUS____Bin",(sum(kl$`2.5 %`[c(3,4,7)])),
#                     (kl$`2.5 %`))
# 
# kl$high_resp<-ifelse(kl$names=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin:Gespeendebiggen__RESPIRATIETRACTUS____Bin",(sum(kl$`97.5 %`[c(3,4,7)])),
#                      (kl$`97.5 %`))

print(kl)


##########
ALLRESULTS[[length(ALLRESULTS)+1]]<-kl
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________END______LMER_BW_for_WEANERS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________START______BIN_BW_for_SOWS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#####Create set with all ABs and etiology
############
gg_AB<-AB_Classes_DF_long[-which(AB_Classes_DF_long$Etiology_All_sowsUP=="Zero use"),]

#testing how much is unaccounted by class very low 0.222 DDDA
test_AB<-data.frame(AB_Classes_DF_long %>%
                      dplyr::group_by(FarmID_1_) %>%
                      dplyr::summarise(DDDA_sows_class=sum(DDDA_sows_class), H9_ii_Whole_DDDA_sows_499_=min(H9_ii_Whole_DDDA_sows_499_)))

test_AB$diff<-test_AB$H9_ii_Whole_DDDA_sows_499_-test_AB$DDDA_sows_class

#sum(test_AB$diff)
test_AB$diff[which(test_AB$diff<0)]<-0
test_AB$diff[which(test_AB$diff<10^-8)]<-0
test_AB$Etiology<-"Unaccounted by AB class"


D<-data.frame(gg_AB %>% 
                dplyr::group_by(Antibiotic_Classes, Etiology_All_sowsUP) %>%
                dplyr::summarise(count = n() ,DDDA_sows=sum(DDDA_sows_class)))
#D
D$Etiology_All_sowsUP<-ifelse(D$Etiology_All_sowsUP=="Zeugen__LOCOMOTIE/ZENUWSTELSEL__||Zeugen__LOCOMOTIE/ZENUWSTELSEL__","Zeugen__LOCOMOTIE/ZENUWSTELSEL__",
                              ifelse(D$Etiology_All_sowsUP=="Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__||Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__",
                                     "Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__",
                                     ifelse(D$Etiology_All_sowsUP=="Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__||Zeugen__LOCOMOTIE/ZENUWSTELSEL__",
                                            "Zeugen__LOCOMOTIE/ZENUWSTELSEL__||Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__",
                                            ifelse(D$Etiology_All_sowsUP=="Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__||Zuigende biggen__DIGESTIETRACTUS__",
                                                   "Zuigende biggen__DIGESTIETRACTUS__||Zuigende biggen__LOCOMOTIE/ZENUWSTELSEL__",
                                                   D$Etiology_All_sowsUP))))

re<-D%>%
  dplyr::group_by(Antibiotic_Classes) %>%
  dplyr::summarise(count = n() ,DDDA_S=sum(DDDA_sows))
#re
cats<-rev(re$Antibiotic_Classes[order(re$DDDA_S)])

D$Antibiotic_Classes<-factor(D$Antibiotic_Classes, levels=cats, ordered=TRUE)

#View(df2_s)
df2_s <-  D %>% 
  # group by Factor_1
  dplyr::group_by(Antibiotic_Classes) %>% 
  # within each group, sort rows according to value of Indep_var
  dplyr::arrange(desc(DDDA_sows)) %>% 
  # create new column with rank (or row number) within each group
  dplyr::mutate(order_in_group=row_number()) %>% 
  # transform rank values into a factor
  dplyr::mutate(order_in_group=as.factor(order_in_group)) %>%
  # remove grouping
  dplyr::ungroup()

#########

library(gtools)
library(glmmTMB)

#####Binomial

AMU_class_s<-colnames(FDF[which(grepl("AMU_Ind_sows",colnames(FDF),fixed=T))])
#i<-"PenicillinsAMU_Ind_sows"
#FDFs
DS<-NULL
#i<-"OtherAMU_Ind_sows"
ABclasses_res_sows<-NULL
for(i in AMU_class_s){
  
  
  DSa<-pull(df2_s[which(df2_s[,1]==substring(i,1,nchar(i)-12)),2])
  DSa<-gsub(" ","",gsub("/","",gsub("-","", unique(unlist(strsplit(DSa,split="||", fixed=T))))))
  DSa<- unique(ifelse(DSa=="NotmentioningABclass","Zeugen_Indtr",ifelse(DSa=="Not_mentioningagegroup","Zeugen_Indtr",DSa)))
  
  DS<-colnames(FDF[which(grepl("__Bin",colnames(FDF),fixed=T) & grepl("Ze|Zu",colnames(FDF),fixed=F))])
  
  DS<-DS[which(DS %notin% c("Zeugen_NotmentioningABclass__Bin",                        
                            "Zeugen_Not_mentioningagegroup__Bin",   
                            "Zeugen__DIGESTIETRACTUS____Bin",
                            "Zeugen__HUIDenWONDINFECTIES____Bin",
                            "Zeugen__RESPIRATIETRACTUS____Bin" ,
                            "ZeugenGeneral_Indications__Bin" ,
                            "Zuigendebiggen__RESPIRATIETRACTUS____Bin",
                            "Zuigendebiggen__SPECIFIEKESYSTEMISCHEAANDOENINGEN____Bin") & 
                 substring(DS,1,nchar(DS)-5) %in% DSa)]
  print(paste("-------------------",i,"-----------------------"))
  
  #FDF$Gamma_outcome<-FDF[,i]
  FDFs$Binomial_outcome<-ifelse(FDFs[,i]>0,1,0)
  
  if(length(table(FDFs$Binomial_outcome))>1 & min(table(FDFs$Binomial_outcome))>15 ){
    
    
    f_g_in<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0(DS,collapse=" + "),
                               "+ Animals_Total_tr  + (1|Vet_practise_name_3_/Vet_name_2_)"))#
    
    
    #First full model to see that it is running =============================================================================================================
    
    ZIG_model_full<-glmer(f_g_in,
                          family=binomial(link="logit"),
                          control= glmerControl(optimizer="bobyqa"),
                          data   = FDFs)
    
    
    if(is.na(AIC(ZIG_model_full))==F ){
      
      if(length(DS[-which(DS=="Zeugen_Indtr__Bin")])>1){
        jjj<-expand.grid(DS[-which(DS=="Zeugen_Indtr__Bin")],DS[-which(DS=="Zeugen_Indtr__Bin")])
        jjj<-jjj[-which(jjj$Var1==jjj$Var2),]
        
        jjj<-data.frame(t(apply(jjj, 1, 
                                FUN=function(x) sort(x, decreasing=TRUE))))
        
        jjj<-jjj[-which(duplicated(paste0(jjj[,1],jjj[,2]))),]
        for(m in 1:nrow(jjj)){
          jjj$inter[m]<-paste0(c(jjj[m,1],jjj[m,2]), collapse=":")
        }
        DSintrZIG<-c(DS,jjj$inter)
      }else{
        DSintrZIG<-DS
      }
      
      
      
      
      #DSintrZIG<-DSintrZIG[-c(18:29)]
      f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0(DSintrZIG,collapse=" + "),
                              "+ Animals_Total_tr+ (1|Vet_practise_name_3_/Vet_name_2_)"))# 
      
      
      
      #First full model to start selection =============================================================================================================
      
      ZIG_model_full<-glmer(f_g,
                            family=binomial(link="logit"),
                            control= glmerControl(optimizer="bobyqa"),
                            data   = FDFs)
      
      
      
      AIC_val_g<-c()
      
      AIC(ZIG_model_full)
      
      while(is.na(AIC(ZIG_model_full))){
        print("---------------INSIDE of 1stwhile----------")
        
        
        j_id<-c()
        for (j in rev(which(grepl(":",DSintrZIG,fixed=T)))){#29:6
          
          j_id<-append(j_id,j)
          print(j_id)
          f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0(DSintrZIG[-j_id],collapse=" + "),
                                  "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
          
          
          
          #Second model for selection in case interactions are NA  =============================================================================================================
          
          ZIG_model_full<-glmer(f_g,
                                family=binomial(link="logit"),
                                control= glmerControl(optimizer="bobyqa"),
                                data   = FDFs)
          
          
          
          print(AIC(ZIG_model_full))
          if(is.na(AIC(ZIG_model_full))==F){
            print("Founded before :) limit of interactions")
            break
          }
        }
        
        if(j==min(which(grepl(":",DSintrZIG,fixed=T)))){
          break
        }
      }
      print("---------------OUT of 1stwhile----------")
      #print(AIC(ZIG_model_full))
      
      DSintrZIG <- rownames(coef(summary(ZIG_model_full)))[-c(1,which(rownames(coef(summary(ZIG_model_full)))=="Animals_Total_tr"))]
      
      print(DSintrZIG)
      
      
      for (u in 1:length(DSintrZIG)){
        f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                                "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
        
        
        ZIG_model_R1<-glmer(f_g,
                            family=binomial(link="logit"),
                            control= glmerControl(optimizer="bobyqa"),
                            data   = FDFs)
        
        mod_test_g<- anova(ZIG_model_R1,ZIG_model_full)
        
        
        
        AIC_val_g<-append(AIC_val_g,mod_test_g$AIC[1])
        
        
        
        
      }
      md_g<- which(AIC_val_g==min(AIC_val_g, na.rm = T))
      AICmd_g<- (AIC_val_g[md_g] - AIC(ZIG_model_full))
      #print(paste("OUT_zig->",DSintrZIG[md_g]))
      
      
      while(AICmd_g<=0 | is.na(AICmd_g)){ #| (AICmd_g>0 & AICmd_g<1.95)
        print("---------------INSIDE of 2nd_while----------")
        
        ##exclude interaction when main effect is off
        if((grepl(":",DSintrZIG[md_g], fixed=T)==F)){
          md_g<-c(md_g,which(grepl(":",DSintrZIG, fixed=T) & grepl(DSintrZIG[md_g],DSintrZIG, fixed=T)))
        }else{md_g<-md_g}
        
        print(paste("OUT_bin->",DSintrZIG[md_g]))
        
        ZIG_model_full<-update(ZIG_model_full,formula=as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0( DSintrZIG[-md_g],collapse=" + "),
                                                                        "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)")))#
        
        
        
        DSintrZIG<-DSintrZIG[-md_g]
        AIC_val<-c()
        
        for (u in 1:length(DSintrZIG)){
          f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                                  "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
          ZIG_model_R1<-update(ZIG_model_R1,formula=f_g)
          
          mod_test<- anova(ZIG_model_R1,ZIG_model_full)
          
          AIC_val<-append(AIC_val,mod_test$AIC[1])
          
        }
        
        md_g<- which(AIC_val==min(AIC_val, na.rm = T))
        AICmd_g<- (AIC_val[md_g] - AIC(ZIG_model_full))
        print(paste("Remove because of lower AIC->",AICmd_g,names(md_g)))
        
        if (any(coef(summary(ZIG_model_full))[-c(1,nrow(coef(summary(ZIG_model_full)))),2]>150) & AICmd_g>0){
          print(AICmd_g)
          AICmd_g<- -1
          print(AICmd_g)
          md_g<- which.max(coef(summary(ZIG_model_full))[-c(1,nrow(coef(summary(ZIG_model_full)))),2])
          print(paste("Remove because of infalted SE->",names(md_g)))
        }
        
      }
      
      #print(ZIG_model_full$call$formula)
      print(summary(ZIG_model_full))
      
      print("End of selection")
      
      FDFs$Binomial_outcome<-NULL
      inf_<-NULL
      tryCatch({
        print(paste(i,"CI_before"))
        tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
        tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="profile",exponentiate=F,effects="fixed"))
      }, error=function(e){cat("*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&ERROR with profile, calculating boot:",conditionMessage(e), "\n")
        tidy_obj<-NULL
        tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
        
      })
      print(tidy_obj)
      inf_<-cbind(AB=i,Preds=rownames(coef(summary(ZIG_model_full)))[-1],
                  Coef=coef(summary(ZIG_model_full))[-1,1],
                  OR=tidy_obj[-1,"estimate"],
                  lower=tidy_obj[-1,"conf.low"],
                  higher=tidy_obj[-1,"conf.high"],
                  Pval=tidy_obj[-1,"p.value"])
      rownames(inf_)<-NULL
      inf_<-data.frame(inf_)
      
      
      if(is.null(inf_$Pval)==F){
        inf_$Bonf<-ifelse(as.numeric(inf_$Pval)<(0.05/ifelse(nrow(inf_)<=1,1,(nrow(inf_)-1))),"Sign","Not")
      }
      inf_$AIC<-AIC(ZIG_model_full)
      if(is.null(inf_$Pval)){
        inf_$Preds<-"None"
        inf_$Coef<-"None"
        inf_$OR<-"None"
        inf_$lower<-"None"
        inf_$higher<-"None"
        inf_$Pval<-"None"
        inf_$Bonf<-"None"
        inf_$AIC<-"None"
        #inf_<-inf_[c(1,3:7,2)]
      }
      
      ABclasses_res_sows<-rbind(ABclasses_res_sows,inf_)
      
      print(inf_)
    }else{print("----Undefined 1st full model----")
      FDFs$Binomial_outcome<-NULL
      
    }
  }else{print("----Not contrast in AMU ind AB----")
    FDFs$Binomial_outcome<-NULL}
}

######
ALLRESULTS[[length(ALLRESULTS)+1]]<-ABclasses_res_sows


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________END______BIN_BW_for_SOWS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________START______ZIG_BW_for_SOWS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#####Sows_ZIGAMMA

AMU_class_s<-colnames(FDF[which(grepl("AMU_Ind_sows",colnames(FDF),fixed=T))])

DS<-NULL

ABclasses_res_sows<-NULL
for(i in AMU_class_s){
  
  
  DSa<-pull(df2_s[which(df2_s[,1]==substring(i,1,nchar(i)-12)),2])
  DSa<-gsub(" ","",gsub("/","",gsub("-","", unique(unlist(strsplit(DSa,split="||", fixed=T))))))
  DSa<- unique(ifelse(DSa=="NotmentioningABclass","Zeugen_Indtr",ifelse(DSa=="Not_mentioningagegroup","Zeugen_Indtr",DSa)))
  
  DS<-colnames(FDF[which(grepl("__Bin",colnames(FDF),fixed=T) & grepl("Ze|Zu",colnames(FDF),fixed=F))])
  
  DS<-DS[which(DS %notin% c("Zeugen_NotmentioningABclass__Bin",                        
                            "Zeugen_Not_mentioningagegroup__Bin",   
                            "Zeugen__DIGESTIETRACTUS____Bin",
                            "Zeugen__HUIDenWONDINFECTIES____Bin",
                            "Zeugen__RESPIRATIETRACTUS____Bin" ,
                            "ZeugenGeneral_Indications__Bin" ,
                            "Zuigendebiggen__RESPIRATIETRACTUS____Bin",
                            "Zuigendebiggen__SPECIFIEKESYSTEMISCHEAANDOENINGEN____Bin") & 
                 substring(DS,1,nchar(DS)-5) %in% DSa)]
  print(paste("-------------------",i,"-----------------------"))
  
  FDFs$Gamma_outcome<-FDFs[,i]
  f_g_in<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0(DS,collapse=" + "),
                             "+ Animals_Total_tr  + (1|Vet_practise_name_3_/Vet_name_2_)"))#
  
  
  
  
  
  ZIG_model_full<-glmmTMB(f_g_in,
                          family=ziGamma(link="log"),
                          ziformula=~1,
                          data= FDFs)
  if(is.na(AIC(ZIG_model_full))==F ){
    
    if(length(DS[-which(DS=="Zeugen_Indtr__Bin")])>1){
      jjj<-expand.grid(DS[-which(DS=="Zeugen_Indtr__Bin")],DS[-which(DS=="Zeugen_Indtr__Bin")])
      jjj<-jjj[-which(jjj$Var1==jjj$Var2),]
      
      jjj<-data.frame(t(apply(jjj, 1, 
                              FUN=function(x) sort(x, decreasing=TRUE))))
      
      jjj<-jjj[-which(duplicated(paste0(jjj[,1],jjj[,2]))),]
      for(m in 1:nrow(jjj)){
        jjj$inter[m]<-paste0(c(jjj[m,1],jjj[m,2]), collapse=":")
      }
      DSintrZIG<-c(DS,jjj$inter)
    }else{
      DSintrZIG<-DS
    }
    
    
    
    
    #DSintrZIG<-DSintrZIG[-c(18:29)]
    f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0(DSintrZIG,collapse=" + "),
                            "+ Animals_Total_tr+ (1|Vet_practise_name_3_/Vet_name_2_)"))# 
    
    
    
    
    
    ZIG_model_full<-glmmTMB(f_g,
                            family=ziGamma(link="log"),
                            ziformula=~1,
                            data= FDFs)
    
    AIC_val_g<-c()
    
    AIC(ZIG_model_full)
    
    while(is.na(AIC(ZIG_model_full))){
      
      print("---------------Inside 1st of while----------")
      
      j_id<-c()
      for (j in rev(which(grepl(":",DSintrZIG,fixed=T)))){#29:6
        
        j_id<-append(j_id,j)
        print(j_id)
        f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0(DSintrZIG[-j_id],collapse=" + "),
                                "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
        
        
        
        
        
        ZIG_model_full<-glmmTMB(f_g,
                                family=ziGamma(link="log"),
                                ziformula=~1,
                                data= FDFs)
        
        print(AIC(ZIG_model_full))
        if(is.na(AIC(ZIG_model_full))==F){
          print("Founded before :) limit of interactions")
          break
        }
      }
      
      if(j==min(which(grepl(":",DSintrZIG,fixed=T)))){
        break
      }
    }
    print("---------------OUT 1st of while----------")
    #print(AIC(ZIG_model_full))
    
    DSintrZIG <- rownames(coef(summary(ZIG_model_full))$cond)[-c(1,which(rownames(coef(summary(ZIG_model_full))$cond)=="Animals_Total_tr"))]
    
    print(DSintrZIG)
    
    
    for (u in 1:length(DSintrZIG)){
      f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                              "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
      
      ZIG_model_R1<-glmmTMB(f_g,
                            family=ziGamma(link="log"),
                            ziformula=~1,
                            data= FDFs)
      
      
      mod_test_g<- anova(ZIG_model_R1,ZIG_model_full)
      
      
      
      AIC_val_g<-append(AIC_val_g,mod_test_g$AIC[1])
      
      
      
      
    }
    md_g<- which(AIC_val_g==min(AIC_val_g, na.rm = T))
    AICmd_g<- (AIC_val_g[md_g] - AIC(ZIG_model_full))
    #print(paste("OUT_zig->",DSintrZIG[md_g]))
    
    
    while(AICmd_g<=0 | is.na(AICmd_g)){ #| (AICmd_g>0 & AICmd_g<1.95)
      
      ##exclude interaction when main effect is off
      if((grepl(":",DSintrZIG[md_g], fixed=T)==F)){
        md_g<-c(md_g,which(grepl(":",DSintrZIG, fixed=T) & grepl(DSintrZIG[md_g],DSintrZIG, fixed=T)))
      }else{md_g<-md_g}
      
      print(paste("OUT_zig->",DSintrZIG[md_g]))
      
      ZIG_model_full<-update(ZIG_model_full,formula=as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0( DSintrZIG[-md_g],collapse=" + "),
                                                                      "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)")))#
      
      
      
      DSintrZIG<-DSintrZIG[-md_g]
      AIC_val<-c()
      if(length(DSintrZIG)!=0){
        for (u in 1:length(DSintrZIG)){
          f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                                  "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
          ZIG_model_R1<-update(ZIG_model_R1,formula=f_g)
          
          mod_test<- anova(ZIG_model_R1,ZIG_model_full)
          
          AIC_val<-append(AIC_val,mod_test$AIC[1])
          
        }
        md_g<- which(AIC_val==min(AIC_val, na.rm = T))
        AICmd_g<- (AIC_val[md_g] - AIC(ZIG_model_full))
      }else{
        AICmd_g<-1
      }
      
    }
    
    #print(ZIG_model_full$call$formula)
    print(summary(ZIG_model_full))
    
    print("End of 1st selection")
    
    FDFs$Gamma_outcome<-NULL
    inf_<-NULL
    
    tryCatch({
      print(paste(AMU_class_s[i],"CI_before"))
      tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
      tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="profile",exponentiate=F,effects="fixed"))
    }, error=function(e){cat("*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&ERROR with profile, calculating boot:",conditionMessage(e), "\n")
      tidy_obj<-NULL
      tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
    })
    
    inf_<-cbind(AB=i,Preds=rownames(coef(summary(ZIG_model_full))$cond)[-1],
                Coef=coef(summary(ZIG_model_full))$cond[-1,1],
                OR=tidy_obj[-c(1,nrow(tidy_obj)),"estimate"],
                lower=tidy_obj[-c(1,nrow(tidy_obj)),"conf.low"],
                higher=tidy_obj[-c(1,nrow(tidy_obj)),"conf.high"],
                Pval=tidy_obj[-c(1,nrow(tidy_obj)),"p.value"])
    
    # inf_<-cbind(AB=i,Preds=rownames(coef(summary(ZIG_model_full))$cond)[-1],coef=coef(summary(ZIG_model_full))$cond[-1,1],
    #             Ser=coef(summary(ZIG_model_full))$cond[-1,2],
    #             Pval=coef(summary(ZIG_model_full))$cond[-1,4])
    rownames(inf_)<-NULL
    inf_<-data.frame(inf_)
    
    if(is.null(inf_$Pval)==F){
      inf_$Bonf<-ifelse(as.numeric(inf_$Pval)<(0.05/ifelse(nrow(inf_)<=1,1,(nrow(inf_)-1))),"Sign","Not")
    }
    inf_$AIC<-AIC(ZIG_model_full)
    if(is.null(inf_$Pval)){
      inf_$Preds<-"None"
      inf_$Coef<-"None"
      inf_$OR<-"None"
      inf_$lower<-"None"
      inf_$higher<-"None"
      inf_$Pval<-"None"
      inf_$Bonf<-"None"
      inf_$AIC<-"None"
      #inf_<-inf_[c(1,3:7,2)]
    }
    
    ABclasses_res_sows<-rbind(ABclasses_res_sows,inf_)
    
    
  }else{print("----Undefined 1st full model----")
    FDFs$Gamma_outcome<-NULL
    
  }
}

ALLRESULTS[[length(ALLRESULTS)+1]]<-ABclasses_res_sows

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________END______ZIG_BW_for_SOWS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________START______BIN_BW_for_WEANERS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

###Weaners_________________________________________________________________________________________________________________________________________________________________
#####Create set with all ABs and etiology

#####
gg_AB<-AB_Classes_DF_long[-which(AB_Classes_DF_long$Etiology_All_weanersUP =="Zero use"),]

#testing how much is unaccounted by class low 3.26 DDDA
test_AB<-data.frame(AB_Classes_DF_long %>%
                      dplyr::group_by(FarmID_1_) %>%
                      dplyr::summarise(DDDA_weaner_class=sum(DDDA_weaner_class), H9_ii_Whole_DDDA_wp_500_=min(H9_ii_Whole_DDDA_wp_500_)))

test_AB$diff<-test_AB$H9_ii_Whole_DDDA_wp_500_-test_AB$DDDA_weaner_class

#sum(test_AB$diff)
test_AB$diff[which(test_AB$diff<0)]<-0
test_AB$diff[which(test_AB$diff<10^-8)]<-0
test_AB$Etiology<-"Unaccounted by AB class"



D<-data.frame(gg_AB %>% 
                dplyr::group_by(Antibiotic_Classes, Etiology_All_weanersUP) %>%#Prevalence_weaners_class
                dplyr::summarise(count = n() ,DDDA_weaners=sum(DDDA_weaner_class)))

D$Etiology_All_weanersUP<-ifelse(D$Etiology_All_weanersUP=="Gespeende biggen__HUID- en WONDINFECTIES__||Gespeende biggen__LOCOMOTIE/ZENUWSTELSEL__","Gespeende biggen__LOCOMOTIE/ZENUWSTELSEL__||Gespeende biggen__HUID- en WONDINFECTIES__",
                                 ifelse(D$Etiology_All_weanersUP=="Gespeende biggen__LOCOMOTIE/ZENUWSTELSEL__||Gespeende biggen__RESPIRATIETRACTUS__","Gespeende biggen__RESPIRATIETRACTUS__||Gespeende biggen__LOCOMOTIE/ZENUWSTELSEL__",
                                        D$Etiology_All_weanersUP))

unique(D$Etiology_All_weanersUP)


re<-D%>%
  dplyr::group_by(Antibiotic_Classes) %>%
  dplyr::summarise(count = n() ,DDDA_W=sum(DDDA_weaners))
#re
cats<-rev(re$Antibiotic_Classes[order(re$DDDA_W)])

#sum(df2$DDDA_weaners)
#sum(FDF$H9_ii_Whole_DDDA_wp_500_, na.rm=TRUE)
D$Antibiotic_Classes<-factor(D$Antibiotic_Classes, levels=cats, ordered=TRUE)

df2_w <-  D %>% 
  # group by Factor_1
  dplyr::group_by(Antibiotic_Classes) %>% 
  # within each group, sort rows according to value of Indep_var
  dplyr::arrange(desc(DDDA_weaners)) %>% 
  # create new column with rank (or row number) within each group
  dplyr::mutate(order_in_group=row_number()) %>% 
  # transform rank values into a factor
  dplyr::mutate(order_in_group=as.factor(order_in_group)) %>%
  # remove grouping
  dplyr::ungroup()


########



#####Binomial
AMU_class_w<-colnames(FDF[which(grepl("AMU_Ind_weaners",colnames(FDF),fixed=T))])

#FDFw3
DS<-NULL
ABclasses_res_weaners<-NULL
for(i in AMU_class_w){
  
  
  
  DSa<-pull(df2_w[which(df2_w[,1]==substring(i,1,nchar(i)-15)),2])
  DSa<-gsub(" ","",gsub("/","",gsub("-","", unique(unlist(strsplit(DSa,split="||", fixed=T))))))
  DSa<- unique(ifelse(DSa=="NotmentioningABclass","Weaners_Indtr",ifelse(DSa=="Not_mentioningagegroup","Weaners_Indtr",DSa)))
  
  DS<-colnames(FDF[which(grepl("__Bin",colnames(FDF),fixed=T) & grepl("Ges|Wea",colnames(FDF),fixed=F))])
  
  DS<-DS[which(DS %notin% c("Gespeendebiggen__HUIDenWONDINFECTIES____Bin",
                            "Gespeendebiggen__SPECIFIEKESYSTEMISCHEAANDOENINGEN____Bin",
                            "GespeendebiggenGeneral_Indications__Bin" ,
                            "Weaners_Not_mentioningagegroup__Bin",
                            "Weaners_NotmentioningABclass__Bin" ) & 
                 substring(DS,1,nchar(DS)-5) %in% DSa)]
  
  print(paste("-------------------",i,"-----------------------"))
  
  #FDF$Gamma_outcome<-FDF[,i]
  FDFw3$Binomial_outcome<-ifelse(FDFw3[,i]>0,1,0)
  
  if(length(table(FDFw3$Binomial_outcome))>1 & min(table(FDFw3$Binomial_outcome))>15 ){
    
    f_g_in<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0(DS,collapse=" + "),
                               " + Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#+ Animals_Total_tr
    
    
    
    
    #First full model to see that it is running =============================================================================================================
    
    ZIG_model_full<-glmer(f_g_in,
                          family=binomial(link="logit"),
                          control= glmerControl(optimizer="bobyqa"),
                          data   = FDFw3)
    
    
    if(is.na(AIC(ZIG_model_full))==F ){
      
      if(length(DS[-which(DS=="Weaners_Indtr__Bin")])>1){
        jjj<-expand.grid(DS[-which(DS=="Weaners_Indtr__Bin")],DS[-which(DS=="Weaners_Indtr__Bin")])
        jjj<-jjj[-which(jjj$Var1==jjj$Var2),]
        
        jjj<-data.frame(t(apply(jjj, 1, 
                                FUN=function(x) sort(x, decreasing=TRUE))))
        
        jjj<-jjj[-which(duplicated(paste0(jjj[,1],jjj[,2]))),]
        for(m in 1:nrow(jjj)){
          jjj$inter[m]<-paste0(c(jjj[m,1],jjj[m,2]), collapse=":")
        }
        DSintrZIG<-c(DS,jjj$inter)
      }else{
        DSintrZIG<-DS
      }
      
      
      
      #DSintrZIG<-DSintrZIG[-c(18:29)]
      f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0(DSintrZIG,collapse=" + "),
                              "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))# 
      
      
      
      
      #First full model to start selection =============================================================================================================
      
      ZIG_model_full<-glmer(f_g,
                            family=binomial(link="logit"),
                            control= glmerControl(optimizer="bobyqa"),
                            data   = FDFw3)
      
      AIC_val_g<-c()
      
      AIC(ZIG_model_full)
      
      while(is.na(AIC(ZIG_model_full))){
        
        
        j_id<-c()
        for (j in rev(which(grepl(":",DSintrZIG,fixed=T)))){#
          
          j_id<-append(j_id,j)
          print(j_id)
          f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0(DSintrZIG[-j_id],collapse=" + "),
                                  "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))# 
          
          
          
          
          #Second model for selection in case interactions are NA  =============================================================================================================
          
          ZIG_model_full<-glmer(f_g,
                                family=binomial(link="logit"),
                                control= glmerControl(optimizer="bobyqa"),
                                data   = FDFw3)
          
          print(AIC(ZIG_model_full))
          if(is.na(AIC(ZIG_model_full))==F){
            break
          }
        }
        
        if(j==min(which(grepl(":",DSintrZIG,fixed=T)))){
          break
        }
      }
      print("---------------OUT of while----------")
      #print(AIC(ZIG_model_full))
      #*
      #DSintrZIG <- rownames(coef(summary(ZIG_model_full))$cond)[-1]
      DSintrZIG <- rownames(coef(summary(ZIG_model_full)))[-c(1,which(rownames(coef(summary(ZIG_model_full)))=="Animals_Total_tr"))]
      #*
      
      print(DSintrZIG)
      
      #THE SELECTION time   =============================================================================================================
      for (u in 1:length(DSintrZIG)){
        f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                                "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))# 
        
        
        ZIG_model_R1<-glmer(f_g,
                            family=binomial(link="logit"),
                            control= glmerControl(optimizer="bobyqa"),
                            data   = FDFw3)
        
        
        mod_test_g<- anova(ZIG_model_R1,ZIG_model_full)
        
        
        
        AIC_val_g<-append(AIC_val_g,mod_test_g$AIC[1])
        
        
        
        
      }
      md_g<- which(AIC_val_g==min(AIC_val_g, na.rm = T))
      AICmd_g<- (AIC_val_g[md_g] - AIC(ZIG_model_full))
      #print(paste("OUT_zig->",DSintrZIG[md_g]))
      
      
      while(AICmd_g<=0 | is.na(AICmd_g)){ #| (AICmd_g>0 & AICmd_g<1.95)
        
        ##exclude interaction when main effect is off
        if((grepl(":",DSintrZIG[md_g], fixed=T)==F)){
          md_g<-c(md_g,which(grepl(":",DSintrZIG, fixed=T) & grepl(DSintrZIG[md_g],DSintrZIG, fixed=T)))
        }else{md_g<-md_g}
        
        print(paste("OUT_bin->",DSintrZIG[md_g]))
        
        ZIG_model_full<-update(ZIG_model_full,formula=as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0( DSintrZIG[-md_g],collapse=" + "),
                                                                        "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)")))#
        
        
        
        DSintrZIG<-DSintrZIG[-md_g]
        AIC_val<-c()
        
        for (u in 1:length(DSintrZIG)){
          f_g<- as.formula(paste0(paste0("Binomial_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                                  "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
          ZIG_model_R1<-update(ZIG_model_R1,formula=f_g)
          
          mod_test<- anova(ZIG_model_R1,ZIG_model_full)
          
          AIC_val<-append(AIC_val,mod_test$AIC[1])
          
        }
        md_g<- which(AIC_val==min(AIC_val, na.rm = T))
        AICmd_g<- (AIC_val[md_g] - AIC(ZIG_model_full))
        print(paste("Remove because of lower AIC->",AICmd_g,names(md_g)))
        
        if (any(coef(summary(ZIG_model_full))[-c(1,nrow(coef(summary(ZIG_model_full)))),2]>150) & AICmd_g>0){
          print(AICmd_g)
          AICmd_g<- -1
          print(AICmd_g)
          md_g<- which.max(coef(summary(ZIG_model_full))[-c(1,nrow(coef(summary(ZIG_model_full)))),2])
          print(paste("Remove because of infalted SE->",names(md_g)))
        }
        
      }
      
      #print(ZIG_model_full$call$formula)
      print(summary(ZIG_model_full))
      
      print("End of selection")
      
      FDFw3$Binomial_outcome<-NULL
      inf_<-NULL
      
      tryCatch({
        print(paste(i,"CI_before"))
        tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
        tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="profile",exponentiate=F,effects="fixed"))
      }, error=function(e){cat("*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&ERROR with profile, calculating boot:",conditionMessage(e), "\n")
        tidy_obj<-NULL
        tidy_obj<-data.frame(broom.mixed::tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
      })
      
      inf_<-cbind(AB=i,Preds=rownames(coef(summary(ZIG_model_full)))[-1],
                  Coef=coef(summary(ZIG_model_full))[-1,1],
                  OR=tidy_obj[-1,"estimate"],
                  lower=tidy_obj[-1,"conf.low"],
                  higher=tidy_obj[-1,"conf.high"],
                  Pval=tidy_obj[-1,"p.value"])
      
      
      rownames(inf_)<-NULL
      inf_<-data.frame(inf_)
      
      if(is.null(inf_$Pval)==F){
        inf_$Bonf<-ifelse(as.numeric(inf_$Pval)<(0.05/ifelse(nrow(inf_)<=1,1,(nrow(inf_)-1))),"Sign","Not")
      }
      inf_$AIC<-AIC(ZIG_model_full)
      if(is.null(inf_$Pval)){
        inf_$Preds<-"None"
        inf_$Coef<-"None"
        inf_$OR<-"None"
        inf_$lower<-"None"
        inf_$higher<-"None"
        inf_$Pval<-"None"
        inf_$Bonf<-"None"
        inf_$AIC<-"None"
        #inf_<-inf_[c(1,3:7,2)]
      }
      
      ABclasses_res_weaners<-rbind(ABclasses_res_weaners,inf_)
      
      
    }else{print("----Undefined 1st full model----")
      FDFw3$Binomial_outcome<-NULL
      
    }
  }else{print("----Not contrast in AMU ind AB----")
    FDFw3$Binomial_outcome<-NULL}
}

ALLRESULTS[[length(ALLRESULTS)+1]]<-ABclasses_res_weaners

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________END______BIN_BW_for_WEANERS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________START______ZIG_BW_for_WEANERS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#####ZIGAMMA
AMU_class_w<-colnames(FDF[which(grepl("AMU_Ind_weaners",colnames(FDF),fixed=T))])


DS<-NULL

ABclasses_res_weaners<-NULL
for(i in AMU_class_w){
  
  
  
  DSa<-pull(df2_w[which(df2_w[,1]==substring(i,1,nchar(i)-15)),2])
  DSa<-gsub(" ","",gsub("/","",gsub("-","", unique(unlist(strsplit(DSa,split="||", fixed=T))))))
  DSa<- unique(ifelse(DSa=="NotmentioningABclass","Weaners_Indtr",ifelse(DSa=="Not_mentioningagegroup","Weaners_Indtr",DSa)))
  
  DS<-colnames(FDF[which(grepl("__Bin",colnames(FDF),fixed=T) & grepl("Ges|Wea",colnames(FDF),fixed=F))])
  
  DS<-DS[which(DS %notin% c("Gespeendebiggen__HUIDenWONDINFECTIES____Bin",
                            "Gespeendebiggen__SPECIFIEKESYSTEMISCHEAANDOENINGEN____Bin",
                            "GespeendebiggenGeneral_Indications__Bin" ,
                            "Weaners_Not_mentioningagegroup__Bin",
                            "Weaners_NotmentioningABclass__Bin" ) & 
                 substring(DS,1,nchar(DS)-5) %in% DSa)]
  
  print(paste("-------------------",i,"-----------------------"))
  
  FDFw3$Gamma_outcome<-FDFw3[,i]
  
  
  f_g_in<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0(DS,collapse=" + "),
                             " + Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))
  
  
  
  
  #First full model to see that it is running =============================================================================================================
  ZIG_model_full<-glmmTMB(f_g_in,
                          family=ziGamma(link="log"),
                          ziformula=~1,
                          data= FDFw3)
  
  
  
  if(is.na(AIC(ZIG_model_full))==F ){
    
    if(length(DS[-which(DS=="Weaners_Indtr__Bin")])>1){
      jjj<-expand.grid(DS[-which(DS=="Weaners_Indtr__Bin")],DS[-which(DS=="Weaners_Indtr__Bin")])
      jjj<-jjj[-which(jjj$Var1==jjj$Var2),]
      
      jjj<-data.frame(t(apply(jjj, 1, 
                              FUN=function(x) sort(x, decreasing=TRUE))))
      
      jjj<-jjj[-which(duplicated(paste0(jjj[,1],jjj[,2]))),]
      for(m in 1:nrow(jjj)){
        jjj$inter[m]<-paste0(c(jjj[m,1],jjj[m,2]), collapse=":")
      }
      DSintrZIG<-c(DS,jjj$inter)
    }else{
      DSintrZIG<-DS
    }
    
    
    
    #DSintrZIG<-DSintrZIG[-c(18:29)]
    f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0(DSintrZIG,collapse=" + "),
                            "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))# 
    
    
    
    
    #First full model to start selection =============================================================================================================
    ZIG_model_full<-glmmTMB(f_g,
                            family=ziGamma(link="log"),
                            ziformula=~1,
                            data= FDFw3)
    
    
    AIC_val_g<-c()
    
    AIC(ZIG_model_full)
    
    while(is.na(AIC(ZIG_model_full))){
      print("---------------Insude 1st of while----------")
      
      
      j_id<-c()
      for (j in rev(which(grepl(":",DSintrZIG,fixed=T)))){#
        
        j_id<-append(j_id,j)
        print(j_id)
        f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0(DSintrZIG[-j_id],collapse=" + "),
                                "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))# 
        
        
        
        
        ZIG_model_full<-glmmTMB(f_g,
                                family=ziGamma(link="log"),
                                ziformula=~1,
                                data= FDFw3)
        
        
        print(AIC(ZIG_model_full))
        if(is.na(AIC(ZIG_model_full))==F){
          break
        }
      }
      
      if(j==min(which(grepl(":",DSintrZIG,fixed=T)))){
        break
      }
    }
    print("---------------OUT of 1st while----------")
    #print(AIC(ZIG_model_full))
    
    
    DSintrZIG <- rownames(coef(summary(ZIG_model_full))$cond)[-c(1,which(rownames(coef(summary(ZIG_model_full))$cond)=="Animals_Total_tr"))]
    
    
    
    print(DSintrZIG)
    
    #THE SELECTION time   =============================================================================================================
    for (u in 1:length(DSintrZIG)){
      f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                              "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))# 
      
      
      ZIG_model_R1<-glmmTMB(f_g,
                            family=ziGamma(link="log"),
                            ziformula=~1,
                            data= FDFw3)
      
      
      
      mod_test_g<- anova(ZIG_model_R1,ZIG_model_full)
      
      
      
      AIC_val_g<-append(AIC_val_g,mod_test_g$AIC[1])
      
      
      
      
    }
    md_g<- which(AIC_val_g==min(AIC_val_g, na.rm = T))
    AICmd_g<- (AIC_val_g[md_g] - AIC(ZIG_model_full))
    #print(paste("OUT_zig->",DSintrZIG[md_g]))
    
    
    while(AICmd_g<0 | is.na(AICmd_g)){ #| (AICmd_g>0 & AICmd_g<1.95)
      
      ##exclude interaction when main effect is off
      if((grepl(":",DSintrZIG[md_g], fixed=T)==F)){
        md_g<-c(md_g,which(grepl(":",DSintrZIG, fixed=T) & grepl(DSintrZIG[md_g],DSintrZIG, fixed=T)))
      }else{md_g<-md_g}
      
      print(paste("OUT_zig->",DSintrZIG[md_g]))
      
      ZIG_model_full<-update(ZIG_model_full,formula=as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0( DSintrZIG[-md_g],collapse=" + "),
                                                                      "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)")))#
      
      
      
      DSintrZIG<-DSintrZIG[-md_g]
      AIC_val<-c()
      
      for (u in 1:length(DSintrZIG)){
        f_g<- as.formula(paste0(paste0("Gamma_outcome"," ~ "), paste0( DSintrZIG[-u],collapse=" + "),
                                "+ Animals_Total_tr + (1|Vet_practise_name_3_/Vet_name_2_)"))#
        ZIG_model_R1<-update(ZIG_model_R1,formula=f_g)
        
        mod_test<- anova(ZIG_model_R1,ZIG_model_full)
        
        AIC_val<-append(AIC_val,mod_test$AIC[1])
        
      }
      md_g<- which(AIC_val==min(AIC_val, na.rm = T))
      AICmd_g<- (AIC_val[md_g] - AIC(ZIG_model_full))
      
      
    }
    
    #print(ZIG_model_full$call$formula)
    print(summary(ZIG_model_full))
    
    print("End of selection")
    
    FDFw3$Gamma_outcome<-NULL
    inf_<-NULL
    
    tryCatch({
      print(paste(AMU_class_s[i],"CI_before"))
      tidy_obj<-data.frame(tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
      tidy_obj<-data.frame(tidy(ZIG_model_full,conf.int=TRUE,conf.method="profile",exponentiate=F,effects="fixed"))
    }, error=function(e){cat("*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&ERROR with profile, calculating boot:",conditionMessage(e), "\n")
      tidy_obj<-NULL
      tidy_obj<-data.frame(tidy(ZIG_model_full,conf.int=TRUE,conf.method="Wald",exponentiate=F,effects="fixed"))
    })
    
    inf_<-cbind(AB=i,Preds=rownames(coef(summary(ZIG_model_full))$cond)[-1],
                Coef=coef(summary(ZIG_model_full))$cond[-1,1],
                OR=tidy_obj[-c(1,nrow(tidy_obj)),"estimate"],
                lower=tidy_obj[-c(1,nrow(tidy_obj)),"conf.low"],
                higher=tidy_obj[-c(1,nrow(tidy_obj)),"conf.high"],
                Pval=tidy_obj[-c(1,nrow(tidy_obj)),"p.value"])
    
    
    rownames(inf_)<-NULL
    inf_<-data.frame(inf_)
    
    if(is.null(inf_$Pval)==F){
      inf_$Bonf<-ifelse(as.numeric(inf_$Pval)<(0.05/ifelse(nrow(inf_)<=1,1,(nrow(inf_)-1))),"Sign","Not")
    }
    inf_$AIC<-AIC(ZIG_model_full)
    if(is.null(inf_$Pval)){
      inf_$Preds<-"None"
      inf_$Coef<-"None"
      inf_$OR<-"None"
      inf_$lower<-"None"
      inf_$higher<-"None"
      inf_$Pval<-"None"
      inf_$Bonf<-"None"
      inf_$AIC<-"None"
    }
    
    ABclasses_res_weaners<-rbind(ABclasses_res_weaners,inf_)
    
    
  }else{print("----Undefined 1st full model----")
    FDFw3$Gamma_outcome<-NULL
    
  }
  
}

ALLRESULTS[[length(ALLRESULTS)+1]]<-ABclasses_res_weaners


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#______________________END______ZIG_BW_for_WEANERS________________________________
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#Ploting the each result

#%%%%%%%%%%___Plot_models___%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ALLRESULTS[[1]]$Bonf<-ifelse(ALLRESULTS[[1]]$Pr...t..<(0.05/(nrow(ALLRESULTS[[1]])-1)),"Yes","No")

ALLRESULTS[[1]]$names<-c("(Intercept)"                                 
                         ,"Individual Treatments in sows/sucklings"                           
                         ,"Musculoskeletal/Neurological in sows (Group Treatments)"        
                         ,"Digestive in sucklings (Group Treatments)"      
                         ,"Musculoskeletal/Neurological in sucklings (Group Treatments)"
                         ,"Number of animals (ln transformed)"  )
S_lmer_plot<-ggplot(ALLRESULTS[[1]][-1,], aes(reorder(names, Estim_Scale), Estim_Scale)) + 
  geom_hline(yintercept=0, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estim_lo_Scale, ymax=Estim_hi_Scale,lty=factor(ALLRESULTS[[1]][-1,"Bonf"])), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[1]][-1,"Bonf"])),size=3) +
  ylim(-0.1, 2.5)+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance",# shape= "Relevance",
       size="Gini importance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))




ALLRESULTS[[2]]$Bonf<-ifelse(ALLRESULTS[[2]]$Pr...t..<(0.05/(nrow(ALLRESULTS[[2]])-1)),"Yes","No")


ALLRESULTS[[2]]$names<-c("(Intercept)"                                 
                         ,"Digestive in weaners (Group Treatments)"                           
                         ,"  Musculoskeletal/Neurological in weaners (Group Treatments)"        
                         ,"Respiratory in weaners (Group Treatments)"      
                         ,"Individual Treatments in weaners"
                         ,"Number of animals (ln transformed)"  )

W_lmer_plot<-ggplot(ALLRESULTS[[2]][-1,], aes(reorder(names, Estim_Scale), Estim_Scale)) + 
  geom_hline(yintercept=0, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=Estim_lo_Scale, ymax=Estim_hi_Scale,lty=factor(ALLRESULTS[[2]][-1,"Bonf"])), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[2]][-1,"Bonf"])),size=3) +
  ylim(-0.1, 2.5)+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance",# shape= "Relevance",
       size="Gini importance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))

#Two lmer models for Sows and weaners
library(cowplot)
plot_grid(S_lmer_plot, W_lmer_plot, labels = c("(A)","(B)"), ncol=1)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
glimpse(ALLRESULTS[[3]])
ALLRESULTS[[3]][,c(3:7)]<-sapply(ALLRESULTS[[3]][,c(3:7)],as.numeric)

ALLRESULTS[[3]]$Preds<-c("Musculoskeletal/Neurological in sucklings (Group Treatments)",
                         "Individual Treatments in sows/sucklings",
                         "Number of animals (ln transformed)",
                         "Individual Treatments in sows/sucklings",
                         "Number of animals (ln transformed)",
                         "        Musculoskeletal/Neurological in sows (Group Treatments)",
                         "Individual Treatments in sows/sucklings",
                         "Number of animals (ln transformed)",
                         "        Musculoskeletal/Neurological in sows (Group Treatments)",
                         "Individual Treatments in sows/sucklings",
                         "Number of animals (ln transformed)",
                         "                                   Digestive in sucklings (Group Treatments)",
                         "Individual Treatments in sows/sucklings",
                         "Number of animals (ln transformed)")

inter<-0
S_Mac_bin_plot<-ggplot(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="Macrolides_lincosamidesAMU_Ind_sows"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="Macrolides_lincosamidesAMU_Ind_sows"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  #ggtitle("(Macrolides/Lincosamides)")+
  geom_point(aes(lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="Macrolides_lincosamidesAMU_Ind_sows"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-1,5.5),breaks  =  seq(-1, 5.5,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance",# shape= "Relevance",
       size="Gini importance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        #plot.title = element_text(face="bold", size=15,hjust=-1.28,vjust=2),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))

S_Pen_bin_plot<-ggplot(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="PenicillinsAMU_Ind_sows"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="PenicillinsAMU_Ind_sows"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  #ggtitle("(Penicillins)")+
  geom_point(aes(lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="PenicillinsAMU_Ind_sows"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-1,5.5),breaks  =  seq(-1, 5.5,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance",# shape= "Relevance",
       size="Gini importance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        #plot.title = element_text(face="bold", size=15,hjust=-0.89,vjust=2),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))

S_Tetra_bin_plot<-ggplot(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="TetracyclinesAMU_Ind_sows"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="TetracyclinesAMU_Ind_sows"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  #ggtitle("(Tetracyclines)")+
  geom_point(aes(lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="TetracyclinesAMU_Ind_sows"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-1,5.5),breaks  =  seq(-1, 5.5,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance",# shape= "Relevance",
       size="Gini importance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        #plot.title = element_text(face="bold", size=15,hjust=-0.95,vjust=2),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))


S_Trim_bin_plot<-ggplot(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="Trimethoprim_SulfonamidesAMU_Ind_sows"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="Trimethoprim_SulfonamidesAMU_Ind_sows"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  #ggtitle("(Trimethoprim/Sulfonamides)")+
  geom_point(aes(lty=factor(ALLRESULTS[[3]][which(ALLRESULTS[[3]]$AB=="Trimethoprim_SulfonamidesAMU_Ind_sows"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-1,5.5),breaks  =  seq(-1, 5.5,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance",# shape= "Relevance",
       size="Gini importance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        #plot.title = element_text(face="bold", size=15,hjust=-1.42,vjust=2),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))



#All BIN  models for sows

plot_grid(S_Mac_bin_plot,S_Pen_bin_plot,S_Tetra_bin_plot, S_Trim_bin_plot,
          hjust = 0, vjust = 1,
          labels = c("(BIN: Macrolides-Lincosamides)","(BIN: Penicillins)","(BIN: Tetracyclines)","(BIN: Trimethoprim-Sulphanamides)"),
          ncol=1)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
glimpse(ALLRESULTS[[5]])
ALLRESULTS[[5]][,c(3:7)]<-sapply(ALLRESULTS[[5]][,c(3:7)],as.numeric)


ALLRESULTS[[5]]$Preds<-c("Digestive in weaners (Group Treatments)",
                         "Respiratory in weaners (Group Treatments)",
                         "Individual Treatments in weaners",
                         "Number of animals (ln transformed)",
                         "           Digestive * Respiratory in weaners (Group Treatments)",
                         "Individual Treatments in weaners",
                         "Number of animals (ln transformed)",
                         "Musculoskeletal/Neurological in weaners (Group Treatments)",
                         "Individual Treatments in weaners",
                         "Number of animals (ln transformed)",
                         "                               Respiratory in weaners (Group Treatments)",
                         "Individual Treatments in weaners",
                         "Number of animals (ln transformed)",
                         "Digestive in weaners (Group Treatments)",
                         "                               Respiratory in weaners (Group Treatments)",
                         "Individual Treatments in weaners",
                         "Number of animals (ln transformed)")

W_Mac_bin_plot<-ggplot(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-4,9),breaks  =  seq(-4, 9,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))

W_Pen_bin_plot<-ggplot(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="PenicillinsAMU_Ind_weaners"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="PenicillinsAMU_Ind_weaners"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="PenicillinsAMU_Ind_weaners"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-4,9),breaks  =  seq(-4, 9,by=1))+
  scale_linetype_manual(values=c(1,2))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))

W_Tetra_bin_plot<-ggplot(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="TetracyclinesAMU_Ind_weaners"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="TetracyclinesAMU_Ind_weaners"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="TetracyclinesAMU_Ind_weaners"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-4,9),breaks  =  seq(-4, 9,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))


W_Trim_bin_plot<-ggplot(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="Trimethoprim_SulfonamidesAMU_Ind_weaners"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="Trimethoprim_SulfonamidesAMU_Ind_weaners"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[5]][which(ALLRESULTS[[5]]$AB=="Trimethoprim_SulfonamidesAMU_Ind_weaners"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-4,9),breaks  =  seq(-4, 9,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))



ALLRESULTS[[6]][,c(3:7)]<-sapply(ALLRESULTS[[6]][,c(3:7)],as.numeric)


ALLRESULTS[[6]][which(ALLRESULTS[[6]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),]$Preds<-c("Respiratory in weaners",
                                                                                               "                                            Number of animals (ln transformed)")


W_Macro_ZIG_plot<-ggplot(ALLRESULTS[[6]][which(ALLRESULTS[[6]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),], aes(reorder(Preds, OR),OR)) + 
  geom_hline(yintercept=inter, lty=3, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower, ymax=higher,lty=factor(ALLRESULTS[[6]][which(ALLRESULTS[[6]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),]$Bonf)), 
                lwd=1.5, width=0) +coord_flip()+
  geom_point(aes(lty=factor(ALLRESULTS[[6]][which(ALLRESULTS[[6]]$AB=="Macrolides_lincosamidesAMU_Ind_weaners"),]$Bonf)),size=3) +
  scale_y_continuous(limits=c(-4,9),breaks  =  seq(-4, 9,by=1))+
  scale_linetype_manual(values=c(2,1))+
  labs(x="Aetiology",y="b-coefficient",lty="Bonferonni adj. significance")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.key.width = unit(2, 'cm'),
        legend.text = element_text(face="bold", size=15),
        legend.title = element_text(face="bold", size=15),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=15),
        axis.text.y=element_text(face="bold", size=15),
        axis.title.y = element_text(vjust=2,face="bold", size=15),
        axis.title.x = element_text(vjust=0,face="bold", size=15),
        panel.grid = element_line(colour = "#e3e3e3"))




#All BIN and ZIG models for weaners
plot_grid(W_Mac_bin_plot,W_Macro_ZIG_plot,W_Pen_bin_plot,W_Tetra_bin_plot, W_Trim_bin_plot,
          hjust = 0, vjust = 1,
          labels = c("(BIN: Macrolides-Lincosamides)","(ZIG: Macrolides-Lincosamides)","(BIN: Penicillins)","(BIN: Tetracyclines)","(BIN: Trimethoprim-Sulphanamides)"), ncol=1)
