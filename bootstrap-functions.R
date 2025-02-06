calc_sens_boot = function(data, indices){
  
  # this function produces a statistic to be bootstrapped
  # it returns the following statistic: the relative difference between the
  ##  sensitivity (sens) of the AI workflow (WF) and standard double
  ##  reading (DR)

  # expected input: a dataframe (data) with columns WF recall, DR recall, and
  ##  cancer diagnoses. 1 indicates Yes, 0 indicates No
  
  # the function is used as an input for boot() from the 'boot' package
  
  dt <- data[indices,]
  recall_WF <- dt[,1]
  recall_DR <- dt[,2]
  cancer <- dt[,3]
  
  # TP - true positives; TN - true negatives
  # FP - false positives; FN - false negatives
  TP_WF <- recall_WF == 1 & cancer == 1
  FN_WF <- recall_WF == 0 & cancer == 1
  
  TP_DR <- recall_DR == 1 & cancer == 1
  FN_DR <- recall_DR == 0 & cancer == 1
  
  sens_WF <- sum(TP_WF)/(sum(TP_WF) + sum(FN_WF))
  sens_DR <- sum(TP_DR)/(sum(TP_DR) + sum(FN_DR))
  
  ratio <- sens_WF/sens_DR
  
  return(ratio)
}

calc_spec_boot = function(data, indices){
  
  # this function produces a statistic to be bootstrapped
  # it returns the following statistic: the relative difference between the
  ##  specificity (spec) of the AI workflow (WF) and standard double
  ##  reading (DR)
  
  # expected input: a dataframe (data) with columns WF recall, DR recall, and
  ##  cancer diagnoses. 1 indicates Yes, 0 indicates No
  
  # the function is used as an input for boot() from the 'boot' package
  
  dt <- data[indices,]
  recall_WF <- dt[,1]
  recall_DR <- dt[,2]
  cancer <- dt[,3]
  
  # TP - true positives; TN - true negatives
  # FP - false positives; FN - false negatives
  TN_WF <- recall_WF == 0 & cancer == 0
  FP_WF <- recall_WF == 1 & cancer == 0
  
  TN_DR <- recall_DR == 0 & cancer == 0
  FP_DR <- recall_DR == 1 & cancer == 0
  
  spec_WF <- sum(TN_WF)/(sum(TN_WF) + sum(FP_WF))
  spec_DR <- sum(TN_DR)/(sum(TN_DR) + sum(FP_DR))
  
  ratio <- spec_WF/spec_DR
  
  return(ratio)
}


calc_RR_boot = function(data, indices){
  
  # this function produces a statistic to be bootstrapped
  # it returns the following statistic: the relative difference between the
  ##  recall rate (RR) of the AI workflow (WF) and standard double reading (DR)
  
  # expected input: a dataframe (data) with columns WF recall, DR recall.
  # 1 indicates Yes, 0 indicates No
  
  # the function is used as an input for boot() from the 'boot' package
  
  dt <- data[indices,]
  recall_WF <- dt[,1]
  recall_DR <- dt[,2]
  
  # N: number of cases in dataset
  N <- nrow(data)
  
  RR_WF <- sum(recall_WF == 1)/N
  RR_DR <- sum(recall_DR == 1)/N
  
  ratio <- RR_WF/RR_DR
  
  return(ratio)
}


calc_CDR_boot = function(data, indices){
  
  # this function produces a statistic to be bootstrapped
  # it returns the following statistic: the relative difference between the
  ##  cancer detection rate (CDR) of the AI workflow (WF) and standard double
  ##  reading (DR)
  
  # expected input: a dataframe (data) with columns WF recall, DR recall, and
  ##  cancer diagnoses. 1 indicates Yes, 0 indicates No
  
  # the function is used as an input for boot() from the 'boot' package
  
  dt <- data[indices,]
  recall_WF <- dt[,1]
  recall_DR <- dt[,2]
  cancer <- dt[,3]
  
  # N: number of cases in dataset
  N <- nrow(data)
  
  CDR_WF <- sum(cancer[recall_WF == 1])/N
  CDR_DR <- sum(cancer[recall_DR == 1])/N
  
  ratio <- CDR_WF/CDR_DR
  
  return(ratio)
}


calc_PPV_boot = function(data, indices){
  
  # this function produces a statistic to be bootstrapped
  # it returns the following statistic: the relative difference between the
  ##  postie predictive value (PPV) of the AI workflow (WF) and standard double
  ##  reading (DR)
  
  # expected input: a dataframe (data) with columns WF recall, DR recall, and
  ##  cancer diagnoses. 1 indicates Yes, 0 indicates No
  
  # the function is used as an input for boot() from the 'boot' package
  
  dt <- data[indices,]
  recall_WF <- dt[,1]
  recall_DR <- dt[,2]
  cancer <- dt[,3]
  
  # TP - true positives; TN - true negatives
  # FP - false positives; FN - false negatives
  TP_WF <- recall_WF == 1 & cancer == 1
  FP_WF <- recall_WF == 1 & cancer == 0
  
  TP_DR <- recall_DR == 1 & cancer == 1
  FP_DR <- recall_DR == 1 & cancer == 0
  
  PPV_WF <- sum(TP_WF)/(sum(TP_WF) + sum(FP_WF))
  PPV_DR <- sum(TP_DR)/(sum(TP_DR) + sum(FP_DR))
  
  ratio <- PPV_WF/PPV_DR
  
  return(ratio)
}

boot_CI = function(reps, workflow, metric, CI_range){
  
  # this function generates and returns two-sided confidence intervals (CI)
  # based on the repetitions produced by boot()
  
  # expected input: repetitions produced by boot(),
  ##  name of the workflow (character string),
  ##  name of the metric compared using bootstrapping (character string),
  ##  and range (e.g. 90, or a vector with multiple ranges)
  
  rows <- data.frame()
  
  for(conf in CI_range){
    CI <- boot.ci(reps, conf = conf, type = c("perc"))
    
    if(!is_empty(CI)){
      row <- data.frame(workflow, metric, conf, CI[[4]][4], CI[[4]][5])
    } else {
      row <- data.frame(workflow, metric, conf, "NA", "NA")
    }
    rows <- rbind(rows, row)
  }
  return(rows)
}