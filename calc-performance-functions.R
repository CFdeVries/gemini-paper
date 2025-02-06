calc_performance = function(recall, cancer, R2_read = 0, R3_read = 0, XR_read = 0){
  
  # this function calculates the performance of a given workflow in terms of its
  ##  recall rate (RR), cancer detection rate (CDR), sensitivity (sens),
  ##  specificity (spec), positive predictive value (PPV), and workload
  ##  it then prints each metric using the perf_print() function, combines them
  ##  and returns these metrics in a dataframe (output)
  
  # expected input: for the workflow, whether a case is recalled (recall), has a
  ##  breast cancer detected (cancer), whether Reader 2 and Reader 3 read the
  ##  case (R2_read/R3_read; optional), and whether whether the case would
  ##  have been additionally read (XR_read; optional)
  
  # TP - true positives; TN - true negatives
  # FP - false positives; FN - false negatives
  TP <- recall == 1 & cancer == 1
  TN <- recall == 0 & cancer == 0
  FP <- recall == 1 & cancer == 0
  FN <- recall == 0 & cancer == 1
  
  sum_TP <- sum(TP)
  sum_TN <- sum(TN)
  sum_FP <- sum(FP)
  sum_FN <- sum(FN)
  
  # Double Reading (DR) workload (with arbitration)
  DR_workload <- 22342
  
  # RR <- sum(recall == 1)/sum(recall == 1 | recall == 0) 
  RR <- perf_print(binom.confint(sum(recall == 1), sum(recall == 1 | recall == 0), conf.level = 0.95, methods = "wilson"))
  
  # CDR <- sum(cancer[recall == 1])/sum(cancer == 1 | cancer == 0)
  CDR <- perf_print(binom.confint(sum(cancer[recall == 1]), sum(cancer == 1 | cancer == 0), conf.level = 0.95, methods = "wilson"), scale = "permil")
  
  # sens <- sum(TP)/(sum(TP) + sum(FN))
  sens <- perf_print(binom.confint(sum_TP, sum_TP + sum_FN, conf.level = 0.95, methods = "wilson"))
  
  # spec <- sum(TN)/(sum(TN) + sum(FP))
  spec <- perf_print(binom.confint(sum_TN, sum_TN + sum_FP, conf.level = 0.95, methods = "wilson"))
  
  # PPV <- sum(TP)/(sum(TP) + sum(FP))
  PPV <- perf_print(binom.confint(sum_TP, sum_TP + sum_FP, conf.level = 0.95, methods = "wilson"))
  
  workload <- length(recall) + sum(R2_read) + sum(R3_read) + sum(XR_read)
  
  workload <- (1-(workload)/DR_workload)*100
  
  workload <- sprintf("%.1f%%", round(workload, 1))
  
  output <- data.frame(CDR, RR, sens, spec, PPV, workload)
  
  colnames(output)[6] <- "workload reduction"
  
  return(output)
}

calc_rel_change = function(recall, recall_DR, cancer){
  
  # this function calculates the performance of a given workflow in terms of its
  ##  relative change compared to routine double reading (DR) for the following
  ##  metrics: recall rate (RR), cancer detection rate (CDR), sensitivity (sens),
  ##  specificity (spec), positive predictive value (PPV), and workload (WL)
  ##  
  ##  it returns these metrics in a dataframe (table)
  
  # expected input: for the workflow, whether a case is recalled (recall), has a
  ##  breast cancer detected (cancer), whether Reader 2 and Reader 3 read the
  ##  case (R2_read/R3_read; optional), and whether whether the case would
  ##  have been additionally read (XR_read; optional)
  
  # other expected inputs: where a case is recalled by routine double reading
  ##  (DR) and whether the third reader read it as part of routine double
  ##  (R3_read_DR)
  
  # TP - true positives; TN - true negatives
  # FP - false positives; FN - false negatives
  TP <- recall == 1 & cancer == 1
  TN <- recall == 0 & cancer == 0
  FP <- recall == 1 & cancer == 0
  FN <- recall == 0 & cancer == 1
  
  TP_DR <- recall_DR == 1 & cancer == 1
  TN_DR <- recall_DR == 0 & cancer == 0
  FP_DR <- recall_DR == 1 & cancer == 0
  FN_DR <- recall_DR == 0 & cancer == 1
  
  CDR <- sum(cancer[recall == 1])/sum(cancer == 1 | cancer == 0)
  CDR_DR <- sum(cancer[recall_DR == 1])/sum(cancer == 1 | cancer == 0)
  
  RR <- sum(recall == 1)/sum(recall == 1 | recall == 0) 
  RR_DR <- sum(recall_DR == 1)/sum(recall_DR == 1 | recall_DR == 0) 
  
  sens <- sum(TP)/(sum(TP) + sum(FN))
  sens_DR <- sum(TP_DR)/(sum(TP_DR) + sum(FN_DR))
  
  spec <- sum(TN)/(sum(TN) + sum(FP))
  spec_DR <- sum(TN_DR)/(sum(TN_DR) + sum(FP_DR))
  
  PPV <- sum(TP)/(sum(TP) + sum(FP))
  PPV_DR <- sum(TP_DR)/(sum(TP_DR) + sum(FP_DR))
  
  table <- data.frame((CDR-CDR_DR)/CDR_DR*100,
                      (RR-RR_DR)/RR_DR*100,
                      (sens-sens_DR)/sens_DR*100,
                      (spec-spec_DR)/spec_DR*100,
                      (PPV-PPV_DR)/PPV_DR*100)
  
  colnames(table) <- c("CDR", "RR", "sens", "spec", "PPV")
  
  return(table)
}

calc_SA_WF = function(df, OPs, scenario){
  
  # this function returns the performance for the standalone (SA) workflows,
  ##  i.e. triage (triage or triage negatives) or AI-Additional Read
  
  # it calls the calc_performance() function, which calculates the recall rate,
  ##  cancer detection rate, sensitivity, specificity, positive predictive value,
  ##  and workload for a given workflow, and returns a datafrane (table) with
  ##  these metrics
  
  # expected input: a dataframe (df) that contains, for the workflow (WF),
  ##  whether a case is recalled (df$"name of the workflow"), has a breast
  ##  cancer detected (df$cancer), whether Reader 2 and Reader 3 read the case
  ##  (R2/R3_read_"name of triage workflow"), and whether whether the case would
  ##  have been additionally read (only for AI-Additional Read workflow)
  ##
  ##  naming convention for workflow name -> e.g. Tneg_Mia_OP2, number is
  ##  OP for the workflow
  ##  Tneg -> triage negatives; alternatives are Triage and
  ##  XR -> AI-Additional Read
  
  # other expected inputs: the operating points for the workflow (OPs) and the
  ##  scenario ("Tneg", "Triage" or "XR")
  
  table <- data.frame()
  
  for(i in OPs){
    if(scenario == "XR"){
      XR_read_col <- paste0("XR_read_Mia_OP", i)
      R2_reads_col = "R2_read"
      R3_reads_col = "R3_read"
      workflow <- paste0("DR_XR_A", i)
    } else {
      workflow <- paste0(scenario, "_Mia_OP", i)
      XR_read_col <- NA
      R2_reads_col = paste0("R2_read_", workflow)
      R3_reads_col = paste0("R3_read_", workflow)
      }

    performance <- calc_performance(df[[workflow]], df$cancer, df[[R2_reads_col]], df[[R3_reads_col]], df[[XR_read_col]])
    
    table <- rbind(table, performance)
  }
  
  table <- cbind(paste0("OP", as.character(OPs)), table)
  colnames(table)[1] <- ""
  
  return(table)
}

calc_combo_WF = function(df, OPs_triage, OPs_AR, triage_scenario){
  
  # this function returns the performance for the combination workflows,
  ##  combining triage with AI-Additional Read
  
  # it calls the calc_performance() function, which calculates the recall rate,
  ##  cancer detection rate, sensitivity, specificity, positive predictive value,
  ##  and workload for a given workflow, and returns a table with these metrics
  
  # expected input: a dataframe (df) that contains, for the workflow (WF),
  ##  whether a case is recalled (df$"name of the workflow"), has a breast
  ##  cancer detected (df$cancer), whether Reader 2 and Reader 3 read the case
  ##  (R2/R3_read_"name of triage workflow"), whether the AI would have
  ##  recalled the case if used standalone (SA) at its AI-Additional Read
  ##  operating point (OP), and whether the case would have been additionally
  ##  read (XR_read)
  ##
  ##  naming convention for workflow name -> e.g. Tneg_XR_A1_A2, first number is
  ##  OP for the triage workflow, second number is OP for AI-Additional Read
  ##  Tneg -> triage negatives; alternative is Triage. XR refers to
  ##  AI-Additional Read
  
  # other expected inputs: the operating points for triage (OPs_triage) and
  ##  AI-Additional Read (OPs_AR) and the triage scenario (triage_scenario,
  ##  either "Tneg" or "Triage")

  table <- data.frame()
    
  WF_OPs = c()
  for(i in OPs_triage){
    for(j in OPs_AR){
      triage = paste0(triage_scenario, "_Mia_OP", i)
      
      workflow = paste0(triage_scenario, "_XR_A", i, "_A", j)
      
      R2_reads_col = paste0("R2_read_", triage)
      R3_reads_col = paste0("R3_read_", triage)
      
      Mia_SA = paste0("Mia_OP", j)
      
      XR_read = if_else(df[[triage]] == 0 & df[[Mia_SA]] == 1, 1, 0)
      
      performance <- calc_performance(df[[workflow]], df$cancer, df[[R2_reads_col]], df[[R3_reads_col]], XR_read)
      
      table <- rbind(table, performance)
      
      WF_OPs <- append(WF_OPs, paste0("OP", i, "+OP",j))
    }
  }
  table <- cbind(WF_OPs, table)
  colnames(table)[1] <- ""
  return(table)
}


perf_print = function(value, scale){
  
  # this function is used by the calc_performance() function to format and print
  # this formatting includes rounding, converting to percentage or per thousand,
  ##  including numerator & denominator, and confidence interval (CI)
  
  # expected input: the output from binom.confint() and scale (percent or per
  ##  thousand)
  
  if(missing(scale)){
    scale = "percent"
  }
  
  if(scale == "percent"){
    sprintf("%.1f%% [%i / %i], CI: %.1f-%.1f",
            round(value$mean*100, 1),
            value$x,
            value$n,
            value$lower*100,
            value$upper*100
    )
  } else if(scale == "permil"){
    sprintf("%.1f per thousand [%i / %i], CI: %.1f-%.1f",
            round(value$mean*1000, 1),
            value$x,
            value$n,
            value$lower*1000,
            value$upper*1000
    )
  }
}