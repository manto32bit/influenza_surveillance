###Created for influenza surveillance project###
### By Channy###


#determine the influenza outbreak
#output a vector of indicators for influenza outbreak (0 or 1)
cal_epidemic = function(df,response,quantile_cut){
  d_var <-rbind(df)[response]
  epidemic.pre = epidemic = rep(0,nrow(d_var))
  for (i in 1:nrow(d_var)){
    if (i <= 53){d_var_yr <- d_var[1:53,1]
    }else{
      d_var_yr <- d_var[(i - 52):i,1]
    }
    c <- quantile(d_var_yr[d_var_yr!=0], quantile_cut, na.rm = T)
    if (!is.na(d_var[i,1])){
      if (d_var[i,1] > c) epidemic.pre[i] = 1
    }
  }
  for (j in 3:nrow(d_var)){
    if ((epidemic.pre[j] == 1) & (epidemic.pre[j-1] == 1) & (epidemic.pre[j-2] == 1)) {
      epidemic[j - 1] = epidemic[j - 2] = epidemic[j]  = 1
    }
  }
  df$epidemic = epidemic
  return(epidemic)
}

####helper functions####

# calculate the mean of cases in previous m weeks
# m in c(2:12) will be tested
moving_mean = function(metric,m){
  test = NULL
  for (i in 1:nrow(df)){
    if (i<=m){
      test[[i]] = NA
    }else{
      test[[i]] = mean(metric[(i-m):i])
    }
  }
  return(unlist(test))
}
# calculate the sd of cases in previous m weeks
# m in c(2:12) will be tested
moving_sd = function(metric,m){
  test = NULL
  for (i in 1:nrow(df)){
    if (i<=m){
      test[[i]] = NA
    }else{
      test[[i]] = sd(metric[(i-m):i])
    }
  }
  return(unlist(test))
}
# calcluate the zm based on mean and sd in previous m weeks
cal_zm = function(metric,m){
  d1 = moving_mean(metric,m=2)
  d2 = moving_sd(metric,m=2)
  return((metric-d1)/d2)
}


#the vector of predicted influenza outbreak will be output

####Non-risk adjusted model####
# CUSUM
# k in [0.5,1,1.5,2] will be tested
cal_unadjusted_cusum = function(metric,m,k){
  zm = cal_zm(metric,m=m)
  upper_cumulative_sum = NULL
  non_na_zm = zm[-seq(1:m)]
  for (i in 1:length(non_na_zm)){
    if (i==1){
      single = max(0,(non_na_zm[i]-k))
    }else{
      single = max(0,(non_na_zm[i]-k+upper_cumulative_sum[i-1]))
    }
    upper_cumulative_sum = append(upper_cumulative_sum,single)
  }
  return(upper_cumulative_sum)
}
# EWMA
# k in [0.5,1,1.5,2] will be tested
# lambda in [0.6,0.7,0.8,0.9] will be tested
cal_unadjusted_ewma = function(metric,m,k,lambda){
  zm = cal_zm(metric,m=m)
  upper_cumulative_sum = NULL
  non_na_zm = zm[-seq(1:m)]
  for (i in 1:length(non_na_zm)){
    if (i==1){
      single = max(0,(non_na_zm[i]-k))
    }else{
      new_zm = non_na_zm[i]*lambda + (1-lambda)*non_na_zm[i-1]
      single = max(0,(new_zm-k+upper_cumulative_sum[i-1]))
    }
    upper_cumulative_sum = append(upper_cumulative_sum,single)
  }
  upper_cumulative_sum = append(rep(NA,m),upper_cumulative_sum)
  return(upper_cumulative_sum)
}
#serfling regresion model
cal_serfling_regressiion = function(df,metric_name,upper_threshold = 0.95){
  new_df = df %>% mutate(t = row_number(),
                         theta = 2*t/52,
                         sin_2 = sinpi(theta),
                         cos_2 = cospi(theta))
  equation = paste(metric_name,"~t + sin_2 + cos_2",sep = "")
  serfling_reg = lm(equation, data = new_df, na.action = na.exclude)
  pred_y = predict(serfling_reg,data=new_df)
  output = pred_y>quantile(pred_y,upper_threshold=0.95)
  return(as.numeric(output))
}


####RA-models #####
# these models will adjust the meteorological factors#

#### calculate the kt for risk-adjusted models ####
# outbreak: the indicator of influenza outbreak (1 or 0)
# OR0: odd ratio under null hypothesis
# ORA: odd ratio under alternative hypothesis
# pred_t: predictions generated from the logit regression with adjustment of weather covariates
cal_kt = function(outbreak,OR0,ORA,pred_t){
  kt = NULL
  for (i in 1:length(outbreak)){
    if(i==0){
      kt_single = log((1-pred_t[i]+OR0*pred_t[i])/(1-pred_t[i]+ORA*pred_t[i])*ORA/OR0)
    }else{
      kt_single = log((1-pred_t[i]+OR0*pred_t[i])/(1-pred_t[i]+ORA*pred_t[i]))
    }
  kt = append(kt,kt_single)
  }
  return(kt)
}

#### RA-CUSUM ####
cal_adjusted_cusum = function(metric,m,outbreak,OR0,ORA,pred_t){
  kt = cal_kt(outbreak,OR0,ORA,pred_t)
  zm = cal_zm(metric,m=m)
  upper_cumulative_sum = NULL
  non_na_zm = zm[-seq(1:m)]
  for (i in 1:length(non_na_zm)){
    if (i==1){
      single = max(0,(non_na_zm[i]-kt[i]))
    }else{
      single = max(0,(non_na_zm[i]-kt[i]+upper_cumulative_sum[i-1]))
    }
    upper_cumulative_sum = append(upper_cumulative_sum,single)
  }
  return(upper_cumulative_sum)
}

#### RA-EWMA ####
cal_adjusted_ewma = function(metric,m,lambda,outbreak,OR0,ORA,pred_t){
  kt = cal_kt(outbreak,OR0,ORA,pred_t)
  zm = cal_zm(metric,m=m)
  upper_cumulative_sum = NULL
  non_na_zm = zm[-seq(1:m)]
  for (i in 1:length(non_na_zm)){
    if (i==1){
      single = max(0,(non_na_zm[i]-kt[i]))
    }else{
      new_zm = non_na_zm[i]*lambda + (1-lambda)*non_na_zm[i-1]
      single = max(0,(new_zm-kt[i]+upper_cumulative_sum[i-1]))
    }
    upper_cumulative_sum = append(upper_cumulative_sum,single)
  }
  upper_cumulative_sum = append(rep(NA,m),upper_cumulative_sum)
  return(upper_cumulative_sum)
}

#### RA-serfling ####
cal_adjusted_serfling_regressiion = function(df,metric_name,covarite_names,upper_threshold = 0.95){
  new_df = df %>% mutate(t = row_number(),
                         theta = 2*t/52,
                         sin_2 = sinpi(theta),
                         cos_2 = cospi(theta))
  equation = paste(metric_name,"~t + sin_2 + cos_2 + ",
                   paste(covarite_names,collapse = " + "),sep = "")
  serfling_reg = lm(equation, data = new_df, na.action = na.exclude)
  pred_y = predict(serfling_reg,data=new_df)
  output = pred_y>quantile(pred_y,upper_threshold=0.95)
  return(as.numeric(output))
}


