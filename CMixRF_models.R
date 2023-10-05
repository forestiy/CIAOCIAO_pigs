rm(list=ls()) 

library(party)
library(partykit)
library(tibble)
library(lme4)
library(lmerTest)
library(stats)
library(SciViews)
library(Rmisc)
library(dplyr)

`%notin%` <- Negate(`%in%`)

load("finalDF.Rdata"); 

FDF<-finalDF




secondary_questions_in<-c("E3ii_Easily_recognised_170_"
                          #,"PCV2_Zuig_Rep_Vac_Bin"#"PCV2_ZuigVac_Bin"
                          #,"Myco_Zuig_Rep_Vac_Bin"#"Myco_ZuigVac_Bin"
                          #,"PRRS_Zuig_Rep_Vac_Bin"#"PRRS_ZuigVac_Bin"
                          ,"Ecoli_Sow_Rep_Vac_Bin"#"Ecoli_SowVac_Bin"
                          #,"Myco_Sow_Rep_Vac_Bin"#"Myco_SowVac_Bin"
                          #,"Ery_Sow_Rep_Vac_Bin" #"Ery_SowVac_Bin"
                          ,"Ssuis_Sow_Rep_Vac_Bin"#"Ssuis_SowVac_Bin"
                          #,"Rota_Sow_Rep_Vac_Bin"#"Rota_SowVac_Bin"
                          ,"PPV_Sow_Rep_Vac_Bin")# "PPV_SowVac_Bin"

diseases<-c(#"Zeugen_Indtr__Bin",
  "Zeugen__LOCOMOTIEZENUWSTELSEL____Bin",
  "Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",
  "Gespeendebiggen__DIGESTIETRACTUS____Bin",
  "Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",
  "Gespeendebiggen__RESPIRATIETRACTUS____Bin"#,
  # "Weaners_Indtr__Bin",
  # "H9_ii_Whole_DDDA_sows_499_",
  # "H9_ii_Whole_DDDA_wp_500_",
  # "H9_ii_Whole_DDDA_fp_501_",
  # "Macrolides_lincosamidesAMU_Ind_sows_Bin",
  # "Macrolides_lincosamidesAMU_Ind_weaners_Bin",
  # "Macrolides_lincosamidesAMU_Ind_weaners",
  # "PenicillinsAMU_Ind_sows_Bin",
  # "PenicillinsAMU_Ind_weaners_Bin",
  # "TetracyclinesAMU_Ind_sows_Bin",
  # "TetracyclinesAMU_Ind_weaners_Bin",
  # "Trimethoprim_SulfonamidesAMU_Ind_sows_Bin",
  # "Trimethoprim_SulfonamidesAMU_Ind_weaners_Bin"
)



Vacsdf<-data.frame(VACC=c(rep("Ecoli_SowVac_Bin",2),rep("Clostr_SowVac_Bin",1),rep("Ssuis_SowVac_Bin",5),
                          #rep("Ssuis_SowVacREP_BIN2",5),
                          rep("Rota_SowVac_Bin",3),rep("PPV_SowVac_Bin",1),#rep("Ery_SowVac_Bin",3),
                          rep("Influ_SowVac_Bin",4),rep("Myco_SowVac_Bin",4),rep("PRRS_SowVac_Bin",4),
                          rep("PCV2_ZuigVac_Bin",3),rep("Myco_ZuigVac_Bin",3),rep("PRRS_ZuigVac_Bin",4),
                          rep("PCV2_WeanVac_Bin",3)),
                   DIS=c("H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_",
                         "H9_ii_Whole_DDDA_wp_500_",
                         "Zeugen__LOCOMOTIEZENUWSTELSEL____Bin","Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin","Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_",
                         #"Zeugen__LOCOMOTIEZENUWSTELSEL____Bin","Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin","Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_",
                         "Gespeendebiggen__DIGESTIETRACTUS____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_",
                         "H9_ii_Whole_DDDA_sows_499_",
                         #"H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__RESPIRATIETRACTUS____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__RESPIRATIETRACTUS____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__RESPIRATIETRACTUS____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__DIGESTIETRACTUS____Bin","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__RESPIRATIETRACTUS____Bin","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__RESPIRATIETRACTUS____Bin","H9_ii_Whole_DDDA_sows_499_","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_",
                         "Gespeendebiggen__DIGESTIETRACTUS____Bin","H9_ii_Whole_DDDA_wp_500_","H9_ii_Whole_DDDA_fp_501_"))

treez<-500
rfMaxit_B<- 1 #200
rfMaxit<- 1 #100
rfMaxit0<-1 #3
B_reps = 0
activate_bootstrap=F
#glimpse(FDF)


#Functions needed
MixRF_mod<-    function (Y, X, random, data, initialRandomEffects = 0, ErrorTolerance = 0.001, 
                         MaxIterations = 100) 
{
  Target = Y
  ContinueCondition = TRUE
  iterations <- 0
  AdjustedTarg <<- Target - initialRandomEffects
  oldLogLik <- -Inf
  while (ContinueCondition) {
    iterations <- iterations + 1
    #rf = randomForest(X, AdjustedTarget, keep.inbag=TRUE, proximity=TRUE)
    #resi = Target - rf$predicted
    ncol = length(strsplit(X,split="[+]")[[1]])
    mtry = max(floor(ncol/3), 1)
    f1 = as.formula(paste0('AdjustedTarg ~ ', X))
    rf =party::cforest(f1, data=data, control = cforest_unbiased(mtry = mtry, ntree = treez))#, weights = weights,
    resi = Target - predict(rf)
    
    f0 = as.formula(paste0("resi ~ -1 + ", random))
    lmefit <- lmer(f0, data = data)
    newLogLik <- as.numeric(logLik(lmefit))
    ContinueCondition <- (abs(newLogLik - oldLogLik) > ErrorTolerance & 
                            iterations < MaxIterations)
    oldLogLik <- newLogLik
    AllEffects <- predict(lmefit)
    AdjustedTarg <<- Target - AllEffects
  }
  result <- list(forest = rf, MixedModel = lmefit, RandomEffects = ranef(lmefit), 
                 IterationsUsed = iterations)
  return(result)
}





MixRFb <- function(Y, x, random, data, initialRandomEffects=0,
                   ErrorTolerance=0.001, MaxIterations=200, NTREES=treez,
                   ErrorTolerance0=0.001, MaxIterations0=15, verbose=FALSE) {
  
  # Condition that indicates the loop has not converged or run out of iterations
  ContinueCondition0 <- TRUE
  iterations0 = 0
  
  # Get initial values
  
  mu = rep(mean(as.numeric(as.character((Y)))),length(Y))
  eta = log(mu/(1-mu))
  y = eta + (as.numeric(as.character((Y)))-mu)/(mu*(1-mu))
  weights = mu*(1-mu)
  
  AdjustedTarget <- y - initialRandomEffects
  
  f1 = as.formula(paste0('AdjustedTarget ~ ', x))
  f0 = as.formula(paste0('resi ~ -1 + ', random))
  
  # mimic randomForest's mtry
  ncol = length(strsplit(x,split="[+]")[[1]])
  mtry = if (!is.null(y) && !is.factor(Y))
    max(floor(ncol/3), 1) else floor(sqrt(ncol))
  
  oldLogLik = oldEta = -Inf
  
  # PQL
  while(ContinueCondition0) {
    
    iterations0 <- iterations0 + 1
    
    iterations = 0
    ContinueCondition <- TRUE
    
    # random forest + lmer
    while(ContinueCondition) {
      
      iterations <- iterations + 1
      
      # random Forest
      data$AdjustedTarget = AdjustedTarget
      rf = party::cforest(f1, data=data, weights = weights, control = cforest_unbiased(mtry = mtry, ntree = treez))
      
      # y - X*beta (out-of-bag prediction)
      pred = predict(rf, OOB = TRUE)
      resi = y - pred
      
      ## Estimate New Random Effects and Errors using lmer
      lmefit <- lmer(f0, data=data, weights=weights)
      
      # check convergence
      LogLik <- as.numeric(logLik(lmefit))
      
      ContinueCondition <- (abs(LogLik-oldLogLik)>ErrorTolerance & iterations < MaxIterations)
      oldLogLik <- LogLik
      
      # Extract (the only) random effects part (Zb) to make the new adjusted target
      AllEffects <- predict(lmefit)
      
      #  y-Zb
      AdjustedTarget <- y - AllEffects
      
      # monitor the change the of logLikelihood
      if(verbose) print(c(iterations0,iterations,LogLik))
    }
    
    eta = pred + AllEffects
    mu = 1/(1+exp(-eta))
    y = eta + (as.numeric(as.character((Y)))-mu)/(mu*(1-mu))
    AdjustedTarget <- y - AllEffects
    weights = as.vector(mu*(1-mu))
    
    print(c(iter = iterations0, maxEtaChange=max(abs(eta-oldEta))))
    
    ContinueCondition0 <- (max(abs(eta-oldEta))>ErrorTolerance0 & iterations0 < MaxIterations0)
    oldEta <- eta
  }
  
  result <- list(forest=rf, MixedModel=lmefit, RandomEffects=ranef(lmefit),
                 IterationsUsed=iterations0)
  
  return(result)
}

########################
#as.vector(which(sapply(FDF, class)=="numeric"|sapply(FDF, class)=="factor"))
# predict the link transformed response (eta)
############
predict.MixRF <- function(object, newdata, EstimateRE=TRUE){
  
  forestPrediction <- predict(object$forest,newdata=newdata,OOB=T)
  
  # If not estimate random effects, just use the forest for prediction.
  if(!EstimateRE){
    return(forestPrediction)
  }
  
  RandomEffects <- predict(object$MixedModel, newdata=newdata, allow.new.levels=T)
  
  completePrediction = forestPrediction + RandomEffects
  
  return(completePrediction)
}
#################

#taildocking is none!
#Find the correct variables for each disease
Sep_d<-data.frame(Vars=colnames(FDF),
                  Cat=c("Gen",rep(NA,length(FDF)-1)))
#####
Sep_d$Cat[which(colnames(FDF)=="FarmID_1_"):which(colnames(FDF)=="Vet_practise_name_3_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="A8i_Farm_vocation_20_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="H9_ii_Whole_DDDA_sows_499_")]<-"Sow"
Sep_d$Cat[which(colnames(FDF)=="H9_ii_Whole_DDDA_wp_500_")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="H9_ii_Whole_DDDA_fp_501_")]<-"Fat"
Sep_d$Cat[which(colnames(FDF)=="DDDA_sows_piglets_tr_"):which(colnames(FDF)=="calc_fat")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Animals_Total")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="Animals_Total_tr")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Zeugen__LOCOMOTIEZENUWSTELSEL____Bin"):which(colnames(FDF)=="Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="AmphenicolsAMU_Ind_sows"):which(colnames(FDF)=="Trimethoprim_SulfonamidesAMU_Ind_sows")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Gespeendebiggen__DIGESTIETRACTUS____Bin"):which(colnames(FDF)=="Zeugen_Indtr__Bin")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="Macrolides_lincosamidesAMU_Ind_weaners"):which(colnames(FDF)=="Trimethoprim_SulfonamidesAMU_Ind_weaners_Bin")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="A3_Employees_status_6_"):which(colnames(FDF)=="A5i_Personnel_num_8_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="A5ii_FTE_num_9_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="A7iib_Cattle_14_")]<-"Sow|Wean"
Sep_d$Cat[which(colnames(FDF)=="A7iic_SheepGoats_16_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="A7i_Other_animals_19_")]<-"Gen"


Sep_d$Cat[which(colnames(FDF)=="A21_Lactating_piglets_losses_24_"):which(colnames(FDF)=="A24_Sows_losses_30_")]<-"NONE"
Sep_d$Cat[c(which(colnames(FDF)=="A25_Live_births_per_litter_32_"),
            which(colnames(FDF)=="A27_Weaned_per_sow_per_year_36_"),
            which(colnames(FDF)=="A28_Weaned_per_sow_per_litter_38_"))]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="A14_Delivery_freq_33_"):which(colnames(FDF)=="A15_UBNs_num_35_")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="A27_Weaned_per_sow_per_year_36_"):which(colnames(FDF)=="A19_Sows_removal_freq_43_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="A31_Open_sow_period_44_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="A20_Worp_Index_45_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="A34_Ventilation_check_48_"):which(colnames(FDF)=="A36_Resp_Areas2_112_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="B1_Loading_bay_Hygiene_Score_113_")]<-"NONE"


Sep_d$Cat[which(colnames(FDF)=="C1_Water_source_127_"):which(colnames(FDF)=="C2_Water_acidification_128_")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="C8iii_Organic_acids_146_"):which(colnames(FDF)=="C8vi_Myc_binder_149_")]<-"Wean"

Sep_d$Cat[which(colnames(FDF)=="D4ii_Design_Score1_156_"):which(colnames(FDF)=="E2_Cleaned_equip_168_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="E3i_Stored_individually_169_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="E3ii_Easily_recognised_170_")]<-"Gen"#Leave one for the code to ork in the if
Sep_d$Cat[which(colnames(FDF)=="E4_Clean_equip_before_arrival_171_"):which(colnames(FDF)=="E5_Disinf_equip_before_arrival_172_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="E6_Foot_baths_score_173_")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="F2_Fence_score_217_")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="G2_Free_sow_220_"):which(colnames(FDF)=="G3_AllinAllout_farrow_221_")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="G4_Anorexia_length_222_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="G5_Slow_growers_mix_223_"):which(colnames(FDF)=="G9_Piglet_nest_score_245_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="G11_Aggression_score_247_")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="G14_Inspect_sick_250_")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="H5_Health_status_origin_check_257_")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="H6iii_DewormingS_260_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="H6iv_ScabiesS_261_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="H6v_Drying_powderS_262_"):which(colnames(FDF)=="H6vii_SoapS_264_")]<-"Sow"

Sep_d$Cat[which(colnames(FDF)=="H6iii_DewormingW_267_")]<-"Wean"

Sep_d$Cat[which(colnames(FDF)=="H6iv_ScabiesW_268_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="H6v_Drying_powderW_269_")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="H6vi_Disinfect_powderW_270_")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="H6vii_SoapW_271_")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="I4ii_Drying_days_506_")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="I5_Cleaning_score_507_"):which(colnames(FDF)=="A1_Province_4_PercDen")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="Educational_level_A6c")]<-"None"
Sep_d$Cat[which(colnames(FDF)=="num_diff_species_")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="Weaners_density_c_A32")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Weaners_density_Overbin_A32")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Weaners_density_Underbin_A32")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Weaners_density_Equalbin_A32")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="CGH_with_sows_A33"):which(colnames(FDF)=="CGH_with_sowsngilts_A33")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="CGH_with_sowsORgilts_A33")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="CGH_with_sowsORgiltsnoFL_A33")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="CGH_with_vleesvarkens_A33")]<-"Wean|Fat"
Sep_d$Cat[which(colnames(FDF)=="CGH_general_A33")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="Loading_bay_bin_B1")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="Incom_Anim_Bin"):which(colnames(FDF)=="Incom_Gilts")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Water_check_c")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Water_check_More_C3")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="Feed_unlimited_C4i")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Diseased_visits_end_D1")]<-"Sow|Wean"

Sep_d$Cat[which(colnames(FDF)=="Diseased_visits_start_D1"):which(colnames(FDF)=="Age_visits_JO")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="Swine_free_contact_24h_D3")]<-"Sow|Fat"
Sep_d$Cat[which(colnames(FDF)=="Foot_baths_bin_E6")]<-"Gen"


Sep_d$Cat[which(colnames(FDF)=="Foot_baths_c_E6")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Fence_score_c_F2")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Companion_animals_cont_with_dier_F3"):which(colnames(FDF)=="Companion_animals_only_out_F3")]<-"Gen"


Sep_d$Cat[which(colnames(FDF)=="Weekly_system_1_G1")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Weekly_system_12_G1"):which(colnames(FDF)=="Weekly_system_123_G1")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="Weekly_system_345_G1"):which(colnames(FDF)=="Weekly_system_45_G1")]<-"NONE"


Sep_d$Cat[which(colnames(FDF)=="CF_plan_bin_G7i")]<-"Sow|Wean"
Sep_d$Cat[which(colnames(FDF)=="CF_plan_c_G7i")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Sick_bay_absence_G12"):which(colnames(FDF)=="Sick_bay_other_comp_or_building_G12")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="H6ii_Alt_treat_reasonS_259_estrus_bin")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Cleaning_frq_DCc_I4i")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Air_Score"):which(colnames(FDF)=="Env_measures_T_bin")]<-"Wean"
#Sep_d$Cat[which(colnames(FDF)=="Env_measures_TRV_bin")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Env_Score"):which(colnames(FDF)=="Water_equip_Score")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Bed_presence")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Bed_score")]<-"None"
Sep_d$Cat[which(colnames(FDF)=="Dir_cont")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Floor_heat")]<-"Wean"

Sep_d$Cat[which(colnames(FDF)=="Floor_Grid_bin")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Floor_Grid")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Beton_PGrid")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Plastic_PGrid")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Beton_partialgrid_floor"):which(colnames(FDF)=="Plastic_grid_floor")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Plastic_partialgrid_floor")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Auto_vent")]<-"Wean"

Sep_d$Cat[which(colnames(FDF)=="Arrival_animals_Bin")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Arrival_animals_Gelten")]<-"Gen"

Sep_d$Cat[which(colnames(FDF)=="Isolat_room_Gilts")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Quarant_room_Gilts")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Pen_room_Gilts")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Castration")]<-"Sow"
Sep_d$Cat[which(colnames(FDF)=="Tail_docking")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Teeth_cut")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Aprt_dis_Freq")]<-"Sow"
Sep_d$Cat[which(colnames(FDF)=="Aprt_dis_AB")]<-"NONE"


Sep_d$Cat[which(colnames(FDF)=="Aprt_dis_None")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Aprt_dis_Else")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Ecoli_SowVac_Bin"):which(colnames(FDF)=="PRRS_ZuigVac_Bin")]<-"Vacs"

Sep_d$Cat[which(colnames(FDF)=="PCV2_WeanVac_Bin")]<-"Vacs"

Sep_d$Cat[c(which(colnames(FDF)=="Kapot_Needle_Zeu_Bin"),which(colnames(FDF)=="hokken_Needle_Zeu_Bin"),which(colnames(FDF)=="nest_Needle_Zeu_Bin"))]<-"Sow"
Sep_d$Cat[c(which(colnames(FDF)=="Kapot_Needle_Zuig_Bin"),which(colnames(FDF)=="hokken_Needle_Zuig_Bin"),
            which(colnames(FDF)=="nest_Needle_Zuig_Bin"), which(colnames(FDF)=="Naaldvrij_Needle_Zuig_Bin"))]<-"Suck"
Sep_d$Cat[c(which(colnames(FDF)=="hokken_Needle_Gesp_Bin"))]<-"Wean"

Sep_d$Cat[which(colnames(FDF)=="Kapot_Needle_ZuigZeu_bin"):which(colnames(FDF)=="Naaldvrij_Needle_ZuigZeu_bin")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="C5_Number_of_feeds_135_"):which(colnames(FDF)=="Fosf")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="E3i_ii_comb_CONT")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="E3i_ii_comb_BIN")]<-"NONE"
Sep_d$Cat[c(which(colnames(FDF)=="Ssuis_SowVacREP_BIN2"))]<-"Vacs"
Sep_d$Cat[which(colnames(FDF)=="BiologischBeterleven3"):which(colnames(FDF)=="Conventioneel")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="ConventioneelBeterleven1")]<-"NONE"



Sep_d$Cat[which(colnames(FDF)=="havovwombo24")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="hbowo")]<-"Gen"
Sep_d$Cat[which(colnames(FDF)=="vbomavovmbombo1")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="Geslotenbedrijfzeugenbiggenenvleesvarkens"):which(colnames(FDF)=="Vermeerderingzeugenenbiggentotca23kg")]<-"Sow|Wean"

Sep_d$Cat[which(colnames(FDF)=="Aanvoervananderbedrijf")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Eigenopfok")]<-"NONE"

Sep_d$Cat[which(colnames(FDF)=="KIalleenzoekbeeraangevoerdvananderbedrijf"):which(colnames(FDF)=="KIalleenzoekbeeruiteigenopfok")]<-"Sow|Wean"

Sep_d$Cat[which(colnames(FDF)=="Brijvoer"):which(colnames(FDF)=="Kruimel")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="twee_weekssysteem"):which(colnames(FDF)=="Wekelijkssysteem")]<-"NONE"
Sep_d$Cat[which(colnames(FDF)=="Pergeslacht"):which(colnames(FDF)=="Pertoom")]<-"Wean"
Sep_d$Cat[which(colnames(FDF)=="Pr_Fat")]<-"Sow|Wean"

########
Sep_d
any(is.na(Sep_d))



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!_____________________Here starts the selection algorithm _____________________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!______________________________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#run the code 20 times so we will select the ones that appears 20 times out of 20
TEN_TIMES_RES<-list()
TEN_TIMES_res_adj<-list()

REPET<-2
for (BS in 1:REPET){
  print(paste("BS->",BS))
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  #First step perform BinomRF
  
  #initially, need to select correct observations and candidate variables for each outcome
  
  Cand_Vars_listSQ<-list()
  Vars_Sl_list_All<-list()
  for (t in 1:length(diseases)){
    
    if (grepl("Gespeendebiggen__DIGESTIETRACTUS____Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__DIGESTIETRACTUS____Bin")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_digestive_if")
      
    }else if  (grepl("Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_locomotive_if")
      
    }else if  (grepl("Gespeendebiggen__RESPIRATIETRACTUS____Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__RESPIRATIETRACTUS____Bin")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_resp_if")
      
    }else if  (grepl("H9_ii_Whole_DDDA_fp_501_",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Fat",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_fp_501_")])
      FDF_rep<-FDF[which(grepl("leesvark",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("fat_DDDA_if")
      
    }else if(grepl("H9_ii_Whole_DDDA_sows_499_",diseases[t]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_sows_499_")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("Sow_DDDA_if")
      
    }else if(grepl("Macrolides_lincosamidesAMU_Ind_sows_Bin|PenicillinsAMU_Ind_sows_Bin|TetracyclinesAMU_Ind_sows_Bin|Trimethoprim_SulfonamidesAMU_Ind_sows_Bin",diseases[t]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      #vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_sows_499_")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("Sow_AB_Classes_if")
      
    }
    else if  (grepl("H9_ii_Whole_DDDA_wp_500_",diseases[t]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_wp_500_")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_DDDA_if")
      
      
    }else if  (grepl("Macrolides_lincosamidesAMU_Ind_weaners_Bin|Macrolides_lincosamidesAMU_Ind_weaners|PenicillinsAMU_Ind_weaners_Bin|TetracyclinesAMU_Ind_weaners_Bin|Trimethoprim_SulfonamidesAMU_Ind_weaners_Bin",diseases[t]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      #vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_wp_500_")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_AB_Classes_if")
      
      
    }
    else if  (grepl("Zeugen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Zeugen__LOCOMOTIEZENUWSTELSEL____Bin")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("sows_locomotion_if")
      
    }else if  (grepl("Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(Sep_d$Cat=="Gen" | Sep_d$Cat=="Sow"|Sep_d$Cat=="Suck"),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("suckling_locomotion_if")
      
    }else if  (grepl("Zeugen_Indtr__Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(Sep_d$Cat=="Gen" | Sep_d$Cat=="Sow" |Sep_d$Cat=="Suck"),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("Sows_ind_if")
      
    }else if  (grepl("Weaners_Indtr__Bin",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_ind_if")
    }else if  (grepl("Dis_wean",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("Dis_wean_if")
    }
    else {print("Errror!!!!")
      
    }
    
    
    
    
    Cand_Vars<-c(colnames(FDF_rep[which(colnames(FDF_rep)=="A3_Employees_status_6_"):length(FDF_rep)]))
    
    Cand_Vars_list<-list()
    
    # apply the averge BinomialRF method by splitting the variables per 10% increase
    for ( w in seq(1/10,1/1,0.1)){
      Cand_Vars_list[[length(Cand_Vars_list)+1]]<-sample(Cand_Vars,round(w*length(Cand_Vars)), replace = F)
    }
    Cand_Vars_listSQ[[length(Cand_Vars_listSQ)+1]]<-Cand_Vars_list
    Vars_Sl<-list()
    
    for( o in 1:length(Cand_Vars_list)){
      #fitting the model depending if it is a binary or continuous outcome
      if(length(table(FDF_rep[,diseases[t]]))==2){
        tr<-MixRFb(Y= FDF_rep[,diseases[t]],
                   x= paste(c(Cand_Vars_list[[o]],"Animals_Total"),collapse="+"),
                   random= "(1 |Vet_practise_name_3_/Vet_name_2_)", 
                   data=FDF_rep, 
                   initialRandomEffects = 0, 
                   ErrorTolerance = 0.001,
                   NTREES=treez,
                   MaxIterations = rfMaxit_B, 
                   ErrorTolerance0 = 0.001, 
                   MaxIterations0 = rfMaxit0,
                   verbose = F)
      }else if(length(table(FDF_rep[,diseases[t]]))>2){
        
        if(any(is.na(FDF_rep[,diseases[t]]))){
          print("missings in response")
          FDF_rep<-FDF_rep[-which(is.na(FDF_rep[,diseases[t]])),]
        }
        
        tr <- MixRF_mod(Y = FDF_rep[,diseases[t]],
                        X = paste(c(Cand_Vars_list[[o]],"Animals_Total"),collapse="+"),
                        random = "(1 |Vet_practise_name_3_/Vet_name_2_)",
                        data = FDF_rep,
                        initialRandomEffects = 0,
                        ErrorTolerance = 0.01,
                        MaxIterations = rfMaxit)
      }
      
      #find the root node variables
      root_vars<-c()
      internal_1vars<-c()
      internal_2vars<-c()
      for ( i in 1:treez){
        out <- party:::prettytree(tr$forest@ensemble[[i]], names(tr$forest@data@get("input")))
        
        root_vars<-append(root_vars,ifelse(is.null( out$psplit$variableName),NA, out$psplit$variableName))
        internal_1vars<-append(internal_1vars,ifelse(is.null(out$left$psplit$variableName),NA,out$left$psplit$variableName))
        internal_2vars<-append(internal_2vars,ifelse(is.null(out$right$psplit$variableName),NA,out$right$psplit$variableName))
        
      }
      
      
      
      
      P<-length(FDF_rep[which(colnames(FDF_rep)=="A3_Employees_status_6_"):length(FDF_rep)])+1 #why + 1 here its because of animals total add
      #P<-length(VARR)+1
      
      m=max(floor(P/3),1)
      #m=floor(sqrt(P))
      
      store<-c()
      #store_int<-c()
      #store_Kway<-c()
      #Proot_Kway<-c()
      for (k in 1:1){#if more than main effects k>1 then it is 1:k e.g. 1:3 three way interaction
        for (i in 1:m){
          
          if(k==1){
            store<-append(store,(P-i)/(P-(i-1)))
            
            #store_int<-append(store_int,((P-1)-i)/((P-1)-(i-1)))#this is for two way interaction will see how to...
            
          }
          #store_Kway<-append(store_Kway,((P-k)-i)/((P-k)-(i-1)))
        }
        #Proot_Kway<-append(Proot_Kway,(1-prod(store_Kway))*(1/m))
        #store_Kway<-NULL
      }
      Proot<-(1-prod(store))*(1/m)
      
      
      
      gs<-data.frame(sort(table(root_vars),decreasing=T))
      
      colnames(gs)[1]<-"Vars"
      
      
      #find significant variables
      
      gs$significance <- sapply(gs$Freq, function(zzz) stats::binom.test(zzz,n = treez, p = Proot, alternative = "greater")$p.value)
      gs$adjSignificance <- stats::p.adjust(gs$significance, method = "bonferroni")
      gs$Sign<-ifelse(gs$adjSignificance<0.05,T,F)
      gs$Effect<-"Main"
      
      
      
      aqw<-gs[which(gs$Sign==T),]
      
      aqw$DS<-diseases[t]
      Vars_Sl[[length(Vars_Sl)+1]]<-aqw
      
      print(paste(o,"from",length(Cand_Vars_list),"Disease->",diseases[t],"BS->",BS))
    }
    
    
    Vars_Sl_list_All[[length(Vars_Sl_list_All)+1]]<-Vars_Sl
    
    print(paste("Disease->",t,diseases[t],"BS->",BS))
  }
  
  Vars_Sl_list_All
  
  
  
  
  #Step2 for completing the binomialRF modeling average method 
  All_res<-list()
  All_res_adj<-list()
  for( f in 1:length(diseases)){#:2 or number of diseases
    
    results<-do.call(rbind,Vars_Sl_list_All[[f]])
    results$Vars<-as.character(results$Vars)
    
    
    totalcountsALL<-NULL
    
    for (i in 1:length(Cand_Vars_listSQ[[f]])){
      totalcounts<-c()
      
      nnams<-names(sort(table(results$Vars),decreasing=T))
      
      for( r in 1:length(nnams)){
        
        if(grepl(" : ",nnams[r],fixed=T)){
          totalcounts<-append(totalcounts,  ifelse(sum(unlist(strsplit(nnams[r]," : ",fixed = T))%in% c(Cand_Vars_listSQ[[f]][[i]],"Animals_Total"))==2,1,0))
          
        }else{
          totalcounts<-append(totalcounts,nnams[r] %in% Cand_Vars_listSQ[[f]][[i]])
        }
        
      }
      totalcountsALL<-cbind(totalcountsALL,as.numeric(totalcounts))
    }
    
    
    totalcountsALL<-data.frame(totalcountsALL)
    totalcountsALL$sum<-rowSums(data.frame(totalcountsALL))
    resultss<-sort(table(results$Vars),decreasing=T)/totalcountsALL$sum
    
    if(length(which(names(resultss)=="Animals_Total"))==1){
      resultss<-resultss[-which(names(resultss)=="Animals_Total")]
    }
    else {
      resultss<-resultss
    }
    
    All_res[[length(All_res)+1]]<-resultss
    
    All_res_adj[[length(All_res_adj)+1]]<-cbind(data.frame(sort(resultss[which(resultss>0.5)],decreasing=T)),diseases[f])
    
  }
  
  
  #Step3 calcualte effect sizes for selected variables in each outcome
  
  #Necessary function which calculate the effect size based on partial dependence method from the rUtilities package
  Eff_Size_Bin_Manual<- function(RFMOD,YVAR,XVAR,DATA){
    n <- nrow(DATA)
    x.pt <- sort(unique(DATA[, XVAR]))
    y.pt <- c()
    w <- rep(1, n)
    for (i in seq(along = x.pt)) {
      
      x.data <- DATA
      x.data[, XVAR] <- rep(x.pt[i], n)
      
      
      if (length(table(DATA[, YVAR]))==2){
        y.pt[i] <- stats::weighted.mean( 1/(1+exp(-predict.MixRF(RFMOD, x.data, EstimateRE=TRUE))),# type="prob")[,2],
                                         w, na.rm = TRUE)
      }else{
        y.pt[i] <- stats::weighted.mean( predict.MixRF(RFMOD, x.data, EstimateRE=TRUE),
                                         w, na.rm = TRUE)
      }
    }
    y.pt<<-y.pt
    x.pt<<-x.pt
    # if(length(table(DATA[, XVAR]))==2){
    #   effect.size <- y.pt[2]-y.pt[1]
    #   a <- effect.size
    #   names(a) <- paste("Effect size for",XVAR, sep = " ")
    #   return(a)
    # } else if (length(table(DATA[, XVAR]))>=3){
    #w <- table(y.pt)
    #ww<-rep(1/w,w)
    #we<-rep(1,length(x.pt))#here on text says explanatory variable in code they are placing y?
    we<-as.numeric(table(x.pt))
    effect.size <- stats::lm(y.pt ~ x.pt , data = DATA, weights = we)
    a <- stats::coefficients(effect.size)[2]
    names(a) <- paste("Effect size for",XVAR, sep = " ")
    return(a)
    #}
  }
  #Eff_Size_Bin_Manual(RFMOD=RF_model,XVAR=as.character(All_res_adj[[1]]$Var1[10]),DATA=FDF)
  
  #need to select correct observations and candidate variables for each outcome
  ES_BOOT_list<-list()
  EF_list<-list()
  for(k in 1:length(diseases)){
    
    #####
    if (grepl("Gespeendebiggen__DIGESTIETRACTUS____Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__DIGESTIETRACTUS____Bin")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_digestive_if")
      
    }else if  (grepl("Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_locomotive_if")
      
    }else if  (grepl("Gespeendebiggen__RESPIRATIETRACTUS____Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__RESPIRATIETRACTUS____Bin")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_resp_if")
      
    }else if  (grepl("H9_ii_Whole_DDDA_fp_501_",diseases[k],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Fat",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_fp_501_")])
      FDF_rep<-FDF[which(grepl("leesvark",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("fat_DDDA_if")
      
    }else if(grepl("H9_ii_Whole_DDDA_sows_499_",diseases[k]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_sows_499_")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("Sow_DDDA_if")
      
    }else if(grepl("Macrolides_lincosamidesAMU_Ind_sows_Bin|PenicillinsAMU_Ind_sows_Bin|TetracyclinesAMU_Ind_sows_Bin|Trimethoprim_SulfonamidesAMU_Ind_sows_Bin",diseases[k]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      #vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_sows_499_")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("Sow_AB_Classes_if")
      
    }
    else if  (grepl("H9_ii_Whole_DDDA_wp_500_",diseases[k]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_wp_500_")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_DDDA_if")
      
      
    }else if  (grepl("Macrolides_lincosamidesAMU_Ind_weaners_Bin|Macrolides_lincosamidesAMU_Ind_weaners|PenicillinsAMU_Ind_weaners_Bin|TetracyclinesAMU_Ind_weaners_Bin|Trimethoprim_SulfonamidesAMU_Ind_weaners_Bin",diseases[k]) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      #vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_wp_500_")])
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_AB_Classes_if")
      
      
    }else if  (grepl("Zeugen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Zeugen__LOCOMOTIEZENUWSTELSEL____Bin")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("sows_locomotion_if")
      
    }else if  (grepl("Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(Sep_d$Cat=="Gen" | Sep_d$Cat=="Sow"|Sep_d$Cat=="Suck"),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")])
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("suckling_locomotion_if")
      
    }else if  (grepl("Zeugen_Indtr__Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(Sep_d$Cat=="Gen" | Sep_d$Cat=="Sow"|Sep_d$Cat=="Suck"),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
      
      print("Sows_ind_if")
      
    }else if  (grepl("Weaners_Indtr__Bin",diseases[k],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("weaner_ind_if")
    }else if  (grepl("Dis_wean",diseases[t],fixed=T) ){
      vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
      
      vars<-vars[-which(vars %in% secondary_questions_in)]
      FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
      
      print("Dis_wean_if")
    }else {print("Errror!!!!")
      
    }
    
    #####
    
    
    
    if(length(table(FDF_rep[,diseases[k]]))==2){
      FDF_qw<-FDF_rep
      RF_model<-MixRFb(Y= FDF_qw[,diseases[k]],
                       x= paste(c(as.character(All_res_adj[[k]]$Var1),"Animals_Total"),collapse="+"),
                       random= "(1 |Vet_practise_name_3_/Vet_name_2_)", 
                       data=FDF_qw, 
                       initialRandomEffects = 0, 
                       ErrorTolerance = 0.001,
                       NTREES=treez,
                       MaxIterations = rfMaxit_B, 
                       ErrorTolerance0 = 0.001, 
                       MaxIterations0 = rfMaxit0,
                       verbose = F)
      
    }else if(length(table(FDF_rep[,diseases[k]]))>2){
      
      if(any(is.na(FDF_rep[,diseases[k]]))){
        print("missings in response")
        FDF_qw<-FDF_rep[-which(is.na(FDF_rep[,diseases[k]])),]
      }else{FDF_qw<-FDF_rep}
      
      RF_model <- MixRF_mod(Y = FDF_qw[,diseases[k]],
                            X = paste(c(as.character(All_res_adj[[k]]$Var1),"Animals_Total"),collapse="+"),
                            random = "(1 |Vet_practise_name_3_/Vet_name_2_)",
                            data = FDF_qw,
                            initialRandomEffects = 0,
                            ErrorTolerance = 0.01,
                            MaxIterations = rfMaxit)
    }
    
    
    effect_sizes<-NULL
    ES_BOOT_DF<-NULL
    for(i in 1:length(All_res_adj[[k]]$Var1)){
      EFs<-cbind(as.character(All_res_adj[[k]]$Var1[i]),
                 NA,
                 as.numeric(Eff_Size_Bin_Manual(RFMOD=RF_model,YVAR=as.character(All_res_adj[[k]]$`diseases[f]`[i]),XVAR=as.character(All_res_adj[[k]]$Var1[i]),DATA=FDF_qw)),
                 NA,
                 diseases[k])
      

      effect_sizes<-rbind(effect_sizes,EFs)#
      print(c(i,diseases[k]))#
    }
    
    
    
    colnames(effect_sizes)<-c("PREDS","upper_ES","mean_ES","lower_ES","RESP")
    EF_list[[length(EF_list)+1]]<-effect_sizes
    
    
  }
 
  TEN_TIMES_RES[[length(TEN_TIMES_RES)+1]]<-EF_list
  TEN_TIMES_res_adj[[length(TEN_TIMES_res_adj)+1]]<-All_res_adj
  
  
}


#&&&&&&&&&&&&&&&&&&&This is the full file for Diseases_EF_pigs analysis &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#Finding the variables that appeared across 20 iterations
library(dplyr)
library(stringr)

for(i in 1:length(TEN_TIMES_RES)){
  
  TEN_TIMES_RES[[i]]<-data.frame(do.call("rbind", TEN_TIMES_RES[[i]]))
}
EF_listDF1 <- data.frame(do.call("rbind", TEN_TIMES_RES))

EF_listDF1$Comb<-paste0(EF_listDF1$PREDS,"IOI",EF_listDF1$RESP)

factors<-names(which(table(EF_listDF1$Comb)==REPET)) #this is to be 20 as the REPET defined in the beginning of the script


ved<-str_split(factors,"IOI",simplify = F)

ved_df<-data.frame()
for (i in 1:length(ved)){
  hjh<- matrix(ved[[i]],ncol=2,byrow=F)
  
  ved_df<-rbind(ved_df,hjh)
}

#EF_listDF1

ffp<-ved_df[order(ved_df$V2),]
ffp$comb<-paste0(ved_df$V1,"IOI",ved_df$V2)


EF_listDF1<-EF_listDF1[which(EF_listDF1$Comb %in% ffp$comb),]

EF_listDF1<-aggregate(as.numeric(EF_listDF1$mean_ES), by=list(PREDS=EF_listDF1$PREDS,RESP=EF_listDF1$RESP,Category=EF_listDF1$Comb), FUN=mean)

names(EF_listDF1)[4]<-"mean_ES"
EF_listDF1$Category<-NULL
#EF_listDF1[order(EF_listDF1$RESP),]


EF_listDF1[,c(3)]<-sapply(EF_listDF1[,c(3)],as.numeric)
EF_listDF1<-data.frame(EF_listDF1%>%group_by(PREDS)%>%dplyr::mutate(count=n()))

library(hutils)

for (i in 1:nrow(EF_listDF1)){
  
  EF_listDF1$Sign[i]<-ifelse(all_same_sign(EF_listDF1$mean_ES[which(EF_listDF1$PREDS==EF_listDF1$PREDS[i])]),"Yes","No")
}

EF_listDF1[order(EF_listDF1$RESP,(EF_listDF1$mean_ES),decreasing=T),]


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#&&&&&&&&&&&&&&&&&&&&&&&&&&&_____________Here is the 95%CI_bootstrap_______&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



EF_listDF2 <- EF_listDF1

colnames(EF_listDF1)[1]<-"Var1"
colnames(EF_listDF1)[2]<-'diseases[f]'
All_res_adj<-split(EF_listDF1,f=EF_listDF1$`diseases[f]`)

diseases_names<-names(All_res_adj)[which(names(All_res_adj) %in%diseases)]

names(All_res_adj)<-NULL

diseases<-diseases_names


B_reps = 2
activate_bootstrap=T



#Step4 calculate the bootstrap, as before initially we need to select the correct obs and variables for each outcome
ES_BOOT_list<-list()
EF_list_St4<-list()
for(k in 1:length(diseases)){
  
  #####
  if (grepl("Gespeendebiggen__DIGESTIETRACTUS____Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__DIGESTIETRACTUS____Bin")])
    FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
    
    print("weaner_digestive_if")
    
  }else if  (grepl("Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")])
    FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
    
    print("weaner_locomotive_if")
    
  }else if  (grepl("Gespeendebiggen__RESPIRATIETRACTUS____Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Gespeendebiggen__RESPIRATIETRACTUS____Bin")])
    FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
    
    print("weaner_resp_if")
    
  }else if  (grepl("H9_ii_Whole_DDDA_fp_501_",diseases[k],fixed=T) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Fat",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_fp_501_")])
    FDF_rep<-FDF[which(grepl("leesvark",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
    
    print("fat_DDDA_if")
    
  }else if(grepl("H9_ii_Whole_DDDA_sows_499_",diseases[k]) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_sows_499_")])
    FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
    
    print("Sow_DDDA_if")
    
  }else if(grepl("Macrolides_lincosamidesAMU_Ind_sows_Bin|PenicillinsAMU_Ind_sows_Bin|TetracyclinesAMU_Ind_sows_Bin|Trimethoprim_SulfonamidesAMU_Ind_sows_Bin",diseases[k]) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    #vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_sows_499_")])
    FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
    
    print("Sow_AB_Classes_if")
    
  }
  else if  (grepl("H9_ii_Whole_DDDA_wp_500_",diseases[k]) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_wp_500_")])
    FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
    
    print("weaner_DDDA_if")
    
    
  }else if  (grepl("Macrolides_lincosamidesAMU_Ind_weaners_Bin|Macrolides_lincosamidesAMU_Ind_weaners|PenicillinsAMU_Ind_weaners_Bin|TetracyclinesAMU_Ind_weaners_Bin|Trimethoprim_SulfonamidesAMU_Ind_weaners_Bin",diseases[k]) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    #vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="H9_ii_Whole_DDDA_wp_500_")])
    FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
    
    print("weaner_AB_Classes_if")
    
    
  }
  else if  (grepl("Zeugen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Sow",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Zeugen__LOCOMOTIEZENUWSTELSEL____Bin")])
    FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
    
    print("sows_locomotion_if")
    
  }else if  (grepl("Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(Sep_d$Cat=="Gen" | Sep_d$Cat=="Sow"|Sep_d$Cat=="Suck"),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    vars<-c(vars,Vacsdf$VACC[which(Vacsdf$DIS=="Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin")])
    FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
    
    print("suckling_locomotion_if")
    
  }else if  (grepl("Zeugen_Indtr__Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(Sep_d$Cat=="Gen" | Sep_d$Cat=="Sow" |Sep_d$Cat=="Suck"),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    FDF_rep<-FDF[which(grepl("eugen",FDF$A8i_Farm_vocation_20_,fixed=T)),vars]
    
    print("Sows_ind_if")
    
  }else if  (grepl("Weaners_Indtr__Bin",diseases[k],fixed=T) ){
    vars<-Sep_d[which(grepl("Gen",Sep_d$Cat,fixed=T) | grepl("Wean",Sep_d$Cat,fixed=T)),1]
    
    vars<-vars[-which(vars %in% secondary_questions_in)]
    FDF_rep<-FDF[-which(is.na(FDF$Bed_score)),vars]
    
    print("weaner_ind_if")
  }else {print("Errror!!!!")
    
  }
  
  #####
  
  Cand_Vars<-EF_listDF2$PREDS[which(EF_listDF2$RESP==diseases[k])]
  
  if(length(table(FDF_rep[,diseases[k]]))==2){
    FDF_qw<-FDF_rep
    RF_model<-MixRFb(Y= FDF_qw[,diseases[k]],
                     x= paste(c(Cand_Vars,"Animals_Total"),collapse="+"),
                     random= "(1 |Vet_practise_name_3_/Vet_name_2_)", 
                     data=FDF_qw, 
                     initialRandomEffects = 0, 
                     ErrorTolerance = 0.001,
                     NTREES=treez,
                     MaxIterations = rfMaxit_B, 
                     ErrorTolerance0 = 0.001, 
                     MaxIterations0 = rfMaxit0,
                     verbose = F)
    
  }else if(length(table(FDF_rep[,diseases[k]]))>2){
    
    if(any(is.na(FDF_rep[,diseases[k]]))){
      print("missings in response")
      FDF_qw<-FDF_rep[-which(is.na(FDF_rep[,diseases[k]])),]
    }else{FDF_qw<-FDF_rep}
    
    RF_model <- MixRF_mod(Y = FDF_qw[,diseases[k]],
                          X = paste(c(Cand_Vars,"Animals_Total"),collapse="+"),
                          random = "(1 |Vet_practise_name_3_/Vet_name_2_)",
                          data = FDF_qw,
                          initialRandomEffects = 0,
                          ErrorTolerance = 0.01,
                          MaxIterations = rfMaxit)
  }
  
  
  effect_sizes<-NULL
  ES_BOOT_DF<-NULL
  for(i in 1:length(All_res_adj[[k]]$Var1)){
    EFs<-cbind(as.character(All_res_adj[[k]]$Var1[i]),
               NA,
               as.numeric(Eff_Size_Bin_Manual(RFMOD=RF_model,YVAR=as.character(All_res_adj[[k]]$`diseases[f]`[i]),XVAR=as.character(All_res_adj[[k]]$Var1[i]),DATA=FDF_qw)),
               NA,
               diseases[k])
    
    if(activate_bootstrap==T){
      #999 #Number of bootstrap iterations
      n = nrow(FDF_qw)
      ES_BOOT <- vector()
      print("Inside_Boot")
      for(db in 1:B_reps) {
        boot.samples <- FDF_qw[sample(1:nrow(FDF_qw), n, replace = TRUE),]
        
        if(length(table(FDF_rep[,diseases[k]]))==2){
          RF_model<-MixRFb(Y= boot.samples[,diseases[k]],
                           x= paste(c(Cand_Vars,"Animals_Total"),collapse="+"),
                           random= "(1 |Vet_practise_name_3_/Vet_name_2_)", 
                           data=boot.samples, 
                           initialRandomEffects = 0, 
                           ErrorTolerance = 0.001,
                           NTREES=treez,
                           MaxIterations = rfMaxit_B, 
                           ErrorTolerance0 = 0.001, 
                           MaxIterations0 = rfMaxit0,
                           verbose = F)
          
        }else if(length(table(FDF_rep[,diseases[k]]))>2){
          RF_model <- MixRF_mod(Y = boot.samples[,diseases[k]],
                                X = paste(c(Cand_Vars,"Animals_Total"),collapse="+"),
                                random = "(1 |Vet_practise_name_3_/Vet_name_2_)",
                                data = boot.samples,
                                initialRandomEffects = 0,
                                ErrorTolerance = 0.01,
                                MaxIterations = rfMaxit)
        }
        
        ES_BOOT <- append(ES_BOOT,  as.numeric(Eff_Size_Bin_Manual(RFMOD=RF_model,YVAR=as.character(All_res_adj[[k]]$`diseases[f]`[i]),XVAR=as.character(All_res_adj[[k]]$Var1[i]),DATA=boot.samples)))
        
        print(c(i,db,diseases[k]))
      }
      ES_BOOT_DF<-cbind(ES_BOOT_DF,as.numeric(c(EFs[3],ES_BOOT)))
      ES_TOT<-t(data.frame(CI(as.numeric(c(EFs[3],ES_BOOT)),ci=(1-0.05)^(1/ncol(ES_BOOT_DF)))))
      rownames(ES_TOT)<-NULL
      EFs<-cbind(EFs[1],ES_TOT,EFs[5])
      
    }
    
    effect_sizes<-rbind(effect_sizes,EFs)#
    print(c(i,diseases[k]))#
  }
  
  if(activate_bootstrap==T){
    ES_BOOT_list[[length(ES_BOOT_list)+1]]<-ES_BOOT_DF
  }
  
  colnames(effect_sizes)<-c("PREDS","upper_ES","mean_ES","lower_ES","RESP")
  EF_list_St4[[length(EF_list_St4)+1]]<-effect_sizes
  
  
}
EF_list_St4
ES_BOOT_list

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#Plotign the results
EF_list_St4<-data.frame(do.call("rbind", EF_list_St4))
EF_listDF1<-EF_list_St4

EF_listDF1[,c(2:4)]<-sapply(EF_listDF1[,c(2:4)],as.numeric)
EF_listDF1<-EF_listDF1[-grep("Indtr__Bin",EF_listDF1$RESP,fixed=T),]
EF_listDF1<-data.frame(EF_listDF1%>%group_by(PREDS)%>%dplyr::mutate(count=n()))

library(hutils)

for (i in 1:nrow(EF_listDF1)){
  
  EF_listDF1$Sign[i]<-ifelse(all_same_sign(EF_listDF1$mean_ES[which(EF_listDF1$PREDS==EF_listDF1$PREDS[i])]),"Yes","No")
}
EF_listDF1$CI_cross<-ifelse(EF_listDF1$upper_ES * EF_listDF1$lower_ES >0,"ok","crosses")

#()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()






EF_listDF1$RESP<-ifelse(EF_listDF1$RESP=="Gespeendebiggen__RESPIRATIETRACTUS____Bin","Respiratory in weaners",
                 ifelse(EF_listDF1$RESP=="Zeugen__LOCOMOTIEZENUWSTELSEL____Bin","Musculoskeletal/Neurological in sows",
                 ifelse(EF_listDF1$RESP=="Zuigendebiggen__LOCOMOTIEZENUWSTELSEL____Bin","Musculoskeletal/Neurological in sucklings",
                 ifelse(EF_listDF1$RESP=="Gespeendebiggen__DIGESTIETRACTUS____Bin","Digestive in weaners",
                 ifelse(EF_listDF1$RESP=="Gespeendebiggen__LOCOMOTIEZENUWSTELSEL____Bin","Musculoskeletal/Neurological in weaners",
                 ifelse(EF_listDF1$RESP=="Zeugen_Indtr__Bin","Individual Treatments in sows",
                 ifelse(EF_listDF1$RESP=="Weaners_Indtr__Bin","Individual Treatments in weaners",
                           "WRONG")))))))

# #############_______Names

# # ################
EF_listDF1$PREDS<-as.character(EF_listDF1$PREDS)
EF_listDF1$PREDS<-ifelse(grepl("Naaldvrij_Needle_Zuig_Bin",EF_listDF1$PREDS,fixed=T), "Use of needle-free vaccines in sucklings (Ref:No)",
                  ifelse(grepl("G2_Free_sow_220_",EF_listDF1$PREDS,fixed=T), "Free sow during lactation (Ref:No)",
                  ifelse(grepl("Isolat_room",EF_listDF1$PREDS,fixed=T), "Incoming animals are placed in an isolated room within the main production building (Ref:No)",
                  ifelse(grepl("Overijsel",EF_listDF1$PREDS,fixed=T), "Low pig density area (Ref: High)",
                   ifelse(grepl("H6vi_Disinfect_powderS_263_",EF_listDF1$PREDS,fixed=T), "Disinfecting powder is used in sows (Ref:No)",
                   ifelse(grepl("Overijssel",EF_listDF1$PREDS,fixed=T), "!!Overijssel!!",
                   ifelse(grepl("A36_Sep_Areas2_111_",EF_listDF1$PREDS,fixed=T), "There is clear separation between clean and dirty areas in the indoors of the farm (Ref:No)",
                   ifelse(grepl("Age_visits_JO",EF_listDF1$PREDS,fixed=T), "Daily workflow is from young to older pigs (Ref:No)",
                   ifelse(grepl("H6vii_SoapS_264_",EF_listDF1$PREDS,fixed=T), "Sows are washed with water and soap (Ref:No)",
                   ifelse(grepl("E2_Cleaned_equip_168_",EF_listDF1$PREDS,fixed=T), "Equipment from other farms is always cleaned and disinfected first (Ref:No)",
                   ifelse(grepl("A15_UBNs_num_35_",EF_listDF1$PREDS,fixed=T), "Number of farms supplying pigs to the farm within the year of 2019",
                   ifelse(grepl("CGH_with_sows_A33",EF_listDF1$PREDS,fixed=T), "At least one building contains weaners with sows (Ref:No)",
                   ifelse(grepl("Arrival_animals_Gelten",EF_listDF1$PREDS,fixed=T), "Gilts are introduced externally (Ref:No)",
                   ifelse(grepl("Tail_docking",EF_listDF1$PREDS,fixed=T), "Tail docking of sucklings (Ref:No)",
                   ifelse(grepl("Conventioneel",EF_listDF1$PREDS,fixed=T), "The farm is conventional (Ref:No)",
                   ifelse(grepl("Companion_animals_only_out_F3",EF_listDF1$PREDS,fixed=T), "Companion animals have access only in the outdoor area of the farm (Ref:No)",
                   ifelse(grepl("H5_Health_status_origin_check_257_",EF_listDF1$PREDS,fixed=T), "Origin farm(s) of incoming animals are of same of higher health status (Ref:No)",
                          ifelse(grepl("Korrel",EF_listDF1$PREDS,fixed=T), "Feed is in form of pellet (Ref:No)",
                          ifelse(grepl("E3i_Stored_individually_169_",EF_listDF1$PREDS,fixed=T), "The equipment is stored and used specifically for each barn separately (Ref:No)",
                          ifelse(grepl("G11_Aggression_score_247_",EF_listDF1$PREDS,fixed=T), "Aggressive behavior score between weaned piglets",
                          ifelse(grepl("E4_Clean_equip_before_arrival_171_",EF_listDF1$PREDS,fixed=T), "The equipment is thoroughly cleaned before each new arrival of incoming pigs (Ref:No)",
                          ifelse(grepl("C8iv_OrgnAc_Storage_147_",EF_listDF1$PREDS,fixed=T), "Organic acids are used in feed storage facilities (Ref:No)",
                          ifelse(grepl("E1_Use_on_other_farms_167_",EF_listDF1$PREDS,fixed=T), "Equipment of the farm is also used by external farms (Ref:No)",
                          ifelse(grepl("Kruimel",EF_listDF1$PREDS,fixed=T), "Feed is in form of crumble (Ref:No)",
                          ifelse(grepl("BiologischBeterleven3",EF_listDF1$PREDS,fixed=T), "The farm is biological/Beter_Leven_3 (Ref:No)",
                           ifelse(grepl("Companion_animals_cont_with_dier_F3",EF_listDF1$PREDS,fixed=T), "Companion animals come in contact with the production animals (Ref:No)",
                           ifelse(grepl("Ssuis_SowVac_Bin",EF_listDF1$PREDS,fixed=T), "Streptococcus Suis vaccination in sows (Ref:No)",
                           ifelse(grepl("Floor_Grid",EF_listDF1$PREDS,fixed=T), "Type of grid on weaner’s floor is fully slatted (Ref:Fully solid)",
                           ifelse(grepl("Quarant_room_Gilts",EF_listDF1$PREDS,fixed=T), "Incoming gilts are quarantined (Ref:No)",
                           ifelse(grepl("Aprt_dis_per_cycle",EF_listDF1$PREDS,fixed=T), "Wound inflicting apparatus for sucklings is dinsfected per cycle (Ref:No)",
                           ifelse(grepl("A36_Resp_Areas1_110_",EF_listDF1$PREDS,fixed=T), "Score of complying with measures of clean and dirty outdoor areas of farm",
                           ifelse(grepl("Cleaning_frq_DCc_I4i",EF_listDF1$PREDS,fixed=T), "Stables of weaners and sucklings are cleaned and disinfected per cycle (Ref:No)",
                          ifelse(grepl("CGH_with_gilts_A33",EF_listDF1$PREDS,fixed=T), "At least one builidng contains weaners with gilts  (Ref:No)",
                           ifelse(grepl("D4ii_Design_Score1_156_",EF_listDF1$PREDS,fixed=T), "Score for the infrastructure of the farm's hygiene lock",
                           ifelse(grepl("Bed_score",EF_listDF1$PREDS,fixed=T), "Score of bedding conditions",
                         ifelse(grepl("PRRS_ZuigVac_Bin",EF_listDF1$PREDS,fixed=T), "PRRS vaccination in sucklings (Ref:No)",
                         ifelse(grepl("Myco_ZuigVac_Bin",EF_listDF1$PREDS,fixed=T), "Mycoplasma hyopneumoniae vaccination in sucklings (Ref:No)",
                         ifelse(grepl("Build_Ins_years",EF_listDF1$PREDS,fixed=T), "Age of the interior of buildings with weaners",
                         ifelse(grepl("KIalleenzoekbeeraangevoerdvananderbedrijf",EF_listDF1$PREDS,fixed=T), "Only search boars present and supplied from external farm  (Ref:No)",
                         ifelse(grepl("A4_Buildings_num_7_",EF_listDF1$PREDS,fixed=T), "Number of total buildings on the farm",
                         ifelse(grepl("A14_Delivery_freq_33_",EF_listDF1$PREDS,fixed=T), "Delivery frequency of pigs from farm",
                         ifelse(grepl("CFib",EF_listDF1$PREDS,fixed=T), "Crude fiber levels in weaners feed",
                         ifelse(grepl("Perkraamkamer",EF_listDF1$PREDS,fixed=T), "Weaners are organized per nursery unit (Ref:No)",
                         ifelse(grepl("Sick_bay_other_comp_or_building_G12",EF_listDF1$PREDS,fixed=T), "Sick bay is present in separate compartment within the weaners' stable or a different building  (Ref:No)",
                         ifelse(grepl("A36_Sep_Areas1_109_",EF_listDF1$PREDS,fixed=T), "There is clear separation between clean and dirty areas in the outdoors of the farm (Ref:No)",
                         ifelse(grepl("A26_Lactation_length_34_",EF_listDF1$PREDS,fixed=T), "Lactation period length (d)",
                         ifelse(grepl("I4ii_Drying_days_506_",EF_listDF1$PREDS,fixed=T), "Drying period length after cleaning and before new animals arrive (d)",
                         ifelse(grepl("CGH_with_sowsngilts_A33",EF_listDF1$PREDS,fixed=T), "At least one building contains weaners with sows/gilts (Ref:No)",
                         ifelse(grepl("Loading_bay_bin_B1",EF_listDF1$PREDS,fixed=T), "Presence of loading bay (Ref:No)",
                                EF_listDF1$PREDS)))))))))))))))))))))))))))))))))))))))))))))))))
EF_listDF1$PREDS<-ifelse(grepl("KIalleenzoekbeeruiteigenopfok",EF_listDF1$PREDS,fixed=T), "Only search boars present and supplied from own production (Ref:No)",
                  ifelse(grepl("J1_Consistency_score_509_",EF_listDF1$PREDS,fixed=T), "Score for quality of workflow",
                  ifelse(grepl("Eigenopfok",EF_listDF1$PREDS,fixed=T), "Gilts originate from own production",
                  ifelse(grepl("Pen_room_Gilts",EF_listDF1$PREDS,fixed=T), "Incoming animals are placed in an isolated pens within the main production building (Ref:No)",
                  ifelse(grepl("Plastic_grid_floor",EF_listDF1$PREDS,fixed=T), "Plastic slatted floor is present (Ref:No)",
                  ifelse(grepl("Sick_bay_Apart_hok_or_none_G12",EF_listDF1$PREDS,fixed=T), "Sick bay is either absent or just one neighbouring pen (Ref:No)",
                  ifelse(grepl("hbowo",EF_listDF1$PREDS,fixed=T), "Tertiary education by farmer (Ref:No)",
                  ifelse(grepl("CGH_general_A33",EF_listDF1$PREDS,fixed=T), "At least one builidng contains weaners with other age groups (Ref:No)",
                  ifelse(grepl("E6_Foot_baths_score_173_",EF_listDF1$PREDS,fixed=T), "Score for use of foot baths",
                  ifelse(grepl("Brijvoer",EF_listDF1$PREDS,fixed=T), "Feed is in liquid form (Ref:No)",
                  ifelse(grepl("J2_Invest_score_510_",EF_listDF1$PREDS,fixed=T), "Farmer's investment motivation score",
                  ifelse(grepl("G14_Inspect_sick_250_",EF_listDF1$PREDS,fixed=T), "Inspection of sick animals is taking place at least twice a day with immediate isolation (Ref:No)",
                  ifelse(grepl("Fosf",EF_listDF1$PREDS,fixed=T), "Phosphorus levels in feed",
                  ifelse(grepl("I6_Disinfecting_score_508_",EF_listDF1$PREDS,fixed=T), "Score of disinfection procedure",
                  ifelse(grepl("C5_Number_of_feeds_135_",EF_listDF1$PREDS,fixed=T), "Number of different feeds used throughout the weaning period",
                  ifelse(grepl("Diseased_visits_start_D1",EF_listDF1$PREDS,fixed=T), "Diseased animals are visited in the beginning of the round (Ref:No)",
                  ifelse(grepl("CF_plan_bin_G7i",EF_listDF1$PREDS,fixed=T), "Cross fostering is taking place with a specific protocol (Ref:No)",
                  ifelse(grepl("Educational_level_A6c",EF_listDF1$PREDS,fixed=T), "Having higher educational level (continuous)",
                  ifelse(grepl("A5i_Personnel_num_8_",EF_listDF1$PREDS,fixed=T), "Number of workers on the farm",
                  ifelse(grepl("Build_Out_years",EF_listDF1$PREDS,fixed=T), "Age of the exterior of buildings with weaners",
                  ifelse(grepl("G5_Slow_growers_mix_223_",EF_listDF1$PREDS,fixed=T), "There is mixing of slow growers",
                  ifelse(grepl("A7iib_Cattle_14_",EF_listDF1$PREDS,fixed=T), "Cattle are present on the farm (Ref:No)",
                  ifelse(grepl("Weekly_system_123_G1",EF_listDF1$PREDS,fixed=T), "The farrowing rhythm is every 1, 2 or 3 weeks (Ref:4 or 5)",
                  ifelse(grepl("H6v_Drying_powderS_262_",EF_listDF1$PREDS,fixed=T), "Drying powders are used in sows (Ref:No)",
                  ifelse(grepl("C1_Water_source_127_",EF_listDF1$PREDS,fixed=T), "Water is supplied from public network (Ref:Private source)",
                  ifelse(grepl("hokken_Needle_Zeu_Bin",EF_listDF1$PREDS,fixed=T), "Needles are changed per pen in sows (Ref:No)",
                  ifelse(grepl("B1_Loading_bay_Hygiene_Score_113_",EF_listDF1$PREDS,fixed=T), "Hygiene score of loading bay",
                  ifelse(grepl("G3_AllinAllout_farrow_221_",EF_listDF1$PREDS,fixed=T), "AIAO system is used during farrowing (Ref:No)",
                  ifelse(grepl("Water_check_More_C3",EF_listDF1$PREDS,fixed=T), "Water is checked more than once per year (Ref:No)",
                  ifelse(grepl("Weekly_system_12_G1",EF_listDF1$PREDS,fixed=T), "The farrowing rhythm is every 1 or 2 weeks (Ref:3, 4 or 5)",
                  ifelse(grepl("Aprt_dis_Else",EF_listDF1$PREDS,fixed=T), "Non-antibiotic disinfectants for wound inflicting apparatus are used",
                  ifelse(grepl("Diseased_visits_noroder_D1",EF_listDF1$PREDS,fixed=T), "No order is followed when visiting diseased animals",
                  ifelse(grepl("E5_Disinf_equip_before_arrival_172_",EF_listDF1$PREDS,fixed=T), "Equipment is disinfected when arrives on farm (Ref:No)",
                  ifelse(grepl("D4i_Hyg_lock2_158_",EF_listDF1$PREDS,fixed=T), "A hygiene lock is present for each building (Ref:No)",
                  ifelse(grepl("Ecoli_SowVac_Bin",EF_listDF1$PREDS,fixed=T), "E. coli vaccination in sows (Ref:No)",
                  ifelse(grepl("Aprt_dis_None",EF_listDF1$PREDS,fixed=T), "No disinfection for wound inflicting apparatus is applied",
                  ifelse(grepl("Ery_SowVac_Bin",EF_listDF1$PREDS,fixed=T), "Erysipelas vaccination in sows (Ref:No)",
                  ifelse(grepl("A1_Province_4_PercDen",EF_listDF1$PREDS,fixed=T), "Pig density of province where the farm is located (heads/Ut.Agr.Ar. in ha)",
                  ifelse(grepl("Myco_SowVac_Bin",EF_listDF1$PREDS,fixed=T), "Mycoplasma vaccination in sows (Ref:No)",
                  ifelse(grepl("A31_Open_sow_period_44_",EF_listDF1$PREDS,fixed=T), "Open sow period length (d)",
                  ifelse(grepl("Weaners_density_Overbin_A32",EF_listDF1$PREDS,fixed=T), "Density in weaners is above requirements respective to the type of production (Ref:No)",
                  ifelse(grepl("Aprt_dis_Freq",EF_listDF1$PREDS,fixed=T), "Frequency of cleaning wound inflicting apparatuses for sucklings",
                  ifelse(grepl("Rota_SowVac_Bin",EF_listDF1$PREDS,fixed=T), "Rotavirus vaccination in sows (Ref:No)",
                  ifelse(grepl("CGH_with_sowsORgilts_A33",EF_listDF1$PREDS,fixed=T), "At least one builidng contains weaners with sows or gilts (Ref:No)",
                  ifelse(grepl("Diseased_visits_end_D1",EF_listDF1$PREDS,fixed=T), "Diseased animals are visited in the end of the round (Ref:No)",
                                EF_listDF1$PREDS)))))))))))))))))))))))))))))))))))))))))))))


EF_listDF1$PREDS<-ifelse(grepl("Auto_vent",EF_listDF1$PREDS,fixed=T), "In weaners buildings there is automatic ventilation (Ref:No)",
                  ifelse(grepl("C8vi_Myc_binder_149_",EF_listDF1$PREDS,fixed=T), "Frequency of using mycotoxins binders in feed for weaners",
                  ifelse(grepl("Swine_free_contact_24h_D3",EF_listDF1$PREDS,fixed=T), "Visitors of the farm have no contact with other pigs for at least 24h (Ref:No)",
                  ifelse(grepl("D4iii_Use_Score1_157_",EF_listDF1$PREDS,fixed=T), "Score for properly using the hygiene lock of the whole farm",
                  ifelse(grepl("Bed_presence",EF_listDF1$PREDS,fixed=T), "Using of bedding material in weaners buildings",
                  ifelse(grepl("Foot_baths_bin_E6",EF_listDF1$PREDS,fixed=T), "Presence of disinfection foot baths on the farm",
                         EF_listDF1$PREDS))))))
################
EF_listDF1$PREDS<-make.unique(EF_listDF1$PREDS, sep = " |")
EF_listDF1$PREDS<-as.character(EF_listDF1$PREDS)

EF_listDF1$PREDS<-ifelse(grepl(" |1",EF_listDF1$PREDS,fixed=T),paste0(" ",substr(EF_listDF1$PREDS,0,(nchar(EF_listDF1$PREDS)-3))),
                         ifelse(grepl(" |2",EF_listDF1$PREDS,fixed=T),paste0("  ",substr(EF_listDF1$PREDS,0,(nchar(EF_listDF1$PREDS)-3))),
                                EF_listDF1$PREDS ))

EF_listDF1<-EF_listDF1[order(EF_listDF1$RESP,EF_listDF1$mean_ES,decreasing=T),]
EF_listDF1$PREDS<-factor(EF_listDF1$PREDS,levels=c(EF_listDF1$PREDS[order(EF_listDF1$RESP,EF_listDF1$mean_ES,decreasing=T)]))
EF_listDF1$RESP<-factor(EF_listDF1$RESP,levels=unique(EF_listDF1$RESP[order(EF_listDF1$RESP,EF_listDF1$mean_ES,decreasing=T)]))

glimpse(EF_listDF1)


palette <- c(  "#BEE942" ,
               "#E195A3", 
               "#C43EE9", 
               "#74E9B1",
               "#7CA7DB", 
               "#E0CA55", 
               "#888D92") 

sizz<-12



library(ggplot2)
library(cowplot)


#Diseases Weaners plot

Wean_dis_plot<-ggplot(EF_listDF1[which(grepl("weaners",EF_listDF1$RESP,fixed=T)),], aes(PREDS, mean_ES)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower_ES, ymax=upper_ES,col=RESP), 
                lwd=0.75, width=0) +
  coord_flip()+
 geom_point(aes(col=RESP,size=factor(count,levels=c("1","2","3","4")), shape=factor(Sign,levels=c("No","Yes")))) +
  scale_colour_manual(values=c( "dark green",
                                "dark green",#Locomotion weaners
                                "dark blue",#Resp_weaners
                                "dark red"),breaks=c("Digestive in weaners", "", "Musculoskeletal/Neurological in weaners",
                                                     "Respiratory in weaners"))+
  scale_shape_manual(values=c(17,19))+
  scale_size_manual(values=c(2,4,6,8))+
  labs(y="Effect Size",x="On-farm practices",colour="Disease aetiology for group treatments",size="Number of diseases it was selected for",shape="Uniform direction")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.text = element_text(face="bold", size=sizz),
        legend.title = element_text(face="bold", size=sizz),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=sizz),
        axis.text.y=element_text(face="bold", size=sizz),
        axis.title.y = element_text(vjust=2,face="bold", size=sizz),
        axis.title.x = element_text(vjust=0,face="bold", size=sizz),
        panel.grid = element_line(colour = "#e3e3e3"))+
  guides(colour = guide_legend(override.aes = list(size=5),order = 1),
         size=guide_legend(order = 2),
         shape = guide_legend(override.aes = list(size=3),order = 3))
Wean_dis_plot




Sows_dis_plot<-ggplot(EF_listDF1[which(grepl("sows|suck",EF_listDF1$RESP)),], aes(PREDS, mean_ES)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lower_ES, ymax=upper_ES,col=RESP), 
                lwd=0.75, width=0) +
  coord_flip()+
  geom_point(aes(col=RESP,size=factor(count,levels=c("1","2")), shape=factor(Sign,levels=c("No","Yes")))) +
  scale_colour_manual(values=c("#E195A3", #nd weaners
                               "#C43EE9"),breaks=c( "Musculoskeletal/Neurological in sows", "Musculoskeletal/Neurological in sucklings"),
                      labels=c( "Musculoskeletal/Neurological in sows", "Musculoskeletal/Neurological in sucklings"))+
  scale_shape_manual(values=c(17,19))+
  scale_size_manual(values=c(2,4))+
  labs(x="On-farm practises",y="Effect Size",colour="Disease aetiology for group treatments",size="Number of diseases it was selected for",shape="Uniform direction")+
  theme(legend.position="right",panel.background = element_rect(fill = 'white'),
        legend.text = element_text(face="bold", size=sizz),
        legend.title = element_text(face="bold", size=sizz),
        axis.text.x=element_text(angle=0,hjust=0.3,vjust=1,face="bold", size=sizz),
        axis.text.y=element_text(face="bold", size=sizz),
        axis.title.y = element_text(vjust=2,face="bold", size=sizz),
        axis.title.x = element_text(vjust=0,face="bold", size=sizz),
        panel.grid = element_line(colour = "#e3e3e3"))+
  guides(colour = guide_legend(override.aes = list(size=5),order = 1),
         size=guide_legend(order = 2),
         shape = guide_legend(override.aes = list(size=3),order = 3))
# guides(colour = "none",
#        size="none",
#        shape = "none")

Sows_dis_plot
