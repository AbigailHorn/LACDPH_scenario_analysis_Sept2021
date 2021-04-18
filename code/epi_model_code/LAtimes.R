## data extracted from LA Times count-level cases, death, hospital, and in-icu data 
## https://github.com/datadesk/california-coronavirus-data


latimes_readin <- function(){
  
  # load COVID county-level hospital data
  hospital_readin <- as.data.frame(data.table::fread("https://raw.githubusercontent.com/datadesk/california-coronavirus-data/master/cdph-hospital-patient-county-totals.csv",
                                                     stringsAsFactors = TRUE) )
  hospital_readin <- hospital_readin %>% filter(county=="Los Angeles")
  hospital_readin$Htot <- hospital_readin$positive_patients + hospital_readin$suspected_patients
  hospital_readin$Q <- hospital_readin$icu_positive_patients + hospital_readin$icu_suspected_patients
  # hospital_readin$Htot <- hospital_readin$positive_patients 
  # hospital_readin$Q <- hospital_readin$icu_positive_patients
  hospital <- select(hospital_readin, c(date, Htot, Q))
  hospital$date <- as.Date(hospital$date)
  
  # latest_data_dates <- latest_data
  # latest_data_dates$date <- as.Date("2020-03-01") + 0:(nrow(latest_data)-1)
  # 
  # data_joined <- left_join(latest_data_dates, hospital, by="date", keep=FALSE)
  
  
  # load COVID county-level case and death data
  case_readin <- as.data.frame(data.table::fread("https://raw.githubusercontent.com/datadesk/california-coronavirus-data/master/latimes-county-totals.csv",
                                                 stringsAsFactors = TRUE) )
  case_readin <- case_readin %>% filter(county=="Los Angeles")
  case_readin$Idetectcum <- case_readin$confirmed_cases
  case_readin$I_detect_new <- case_readin$new_confirmed_cases
  case_readin$D <- case_readin$deaths
  case_readin$D_new <- case_readin$new_deaths
  cases <- select(case_readin, c(date,Idetectcum,I_detect_new,D,D_new))
  cases$date <- as.Date(cases$date)
  cases <- cases %>% filter(date > "2020-02-29") %>% arrange(date)
  
  la_data <- left_join(cases, hospital, by="date", keep=FALSE)
  #la_data[is.na(la_data)] <- 0
  
  la_data$I_detect_new = zoo::rollmean(la_data$I_detect_new, k = 7, fill = NA, align = 'right') %>% round(digits=0)
  la_data$D_new = zoo::rollmean(la_data$D_new, k = 7, fill = NA, align = 'right') %>% round(digits=0)
  
  return(la_data)
}

#la_data <- latimes_readin()

latimes_age_readin <- function(){
  ages_readin <- as.data.frame(data.table::fread("https://raw.githubusercontent.com/datadesk/california-coronavirus-data/master/cdph-age.csv", stringsAsFactors = TRUE) )
  
  ages_readin <- unique(ages_readin)
  
  ages_daily = ages_readin %>%
    group_by(age) %>%
    arrange(date) %>%  # first sort by day
    mutate(Diff_day = date - lag(date),  # Difference in time (just in case there are gaps)
           new_cases = confirmed_cases_total - lag(confirmed_cases_total),
           new_deaths = deaths_total - lag(deaths_total)) %>% # Difference in case between days 
    arrange(age)
  
  age_0.18 = ages_readin %>% filter(age %in% c("0-4","5-17")) %>% group_by(date) %>% mutate(case_pct = sum(confirmed_cases_percent), death_pct = sum(deaths_percent)) %>% filter(age=="0-4")
  age_0.18$age.strata = "age_0.18"
  age_0.18$age = NULL
  age_0.18$confirmed_cases_percent = NULL
  age_0.18$deaths_percent = NULL
  
  age_19.49 = ages_readin %>% filter(age %in% c("18-34","35-49")) %>% group_by(date) %>% mutate(case_pct = sum(confirmed_cases_percent), death_pct = sum(deaths_percent)) %>% filter(age=="18-34")
  age_19.49$age.strata = "age_19.49"
  age_19.49$age = NULL
  age_19.49$confirmed_cases_percent = NULL
  age_19.49$deaths_percent = NULL
  
  age_50.64 = ages_readin %>% filter(age %in% c("50-59","60-64")) %>% group_by(date) %>% mutate(case_pct = sum(confirmed_cases_percent), death_pct = sum(deaths_percent)) %>% filter(age=="50-59")
  age_50.64$age.strata = "age_50.64"
  age_50.64$age = NULL
  age_50.64$confirmed_cases_percent = NULL
  age_50.64$deaths_percent = NULL
  
  age_65.79 = ages_readin %>% filter(age %in% c("65-69","70-74","75-79")) %>% group_by(date) %>% mutate(case_pct = sum(confirmed_cases_percent), death_pct = sum(deaths_percent)) %>% filter(age=="65-69")
  age_65.79$age.strata = "age_65.79"
  age_65.79$age = NULL
  age_65.79$confirmed_cases_percent = NULL
  age_65.79$deaths_percent = NULL
  
  age_80. = ages_readin %>% filter(age %in% "80+") %>% group_by(date) %>% mutate(case_pct = sum(confirmed_cases_percent), death_pct = sum(deaths_percent)) %>% filter(age=="80+")
  age_80.$age.strata = "age_80."
  age_80.$age = NULL
  age_80.$confirmed_cases_percent = NULL
  age_80.$deaths_percent = NULL
  
  ages_pct <- rbind(age_0.18, age_19.49, age_50.64, age_65.79, age_80.) %>% arrange(date) %>% as.data.frame()
  ages_pct$age.strata <- as.factor(ages_pct$age.strata)
  ages_pct$date <- as.Date(ages_pct$date)
  return(ages_pct)
}

#ages_pct <- latimes_age_readin()

#ages_pct_test <- ages_pct %>% group_by(date) %>% mutate(test.tot = sum(case_pct), test.tot.death = sum(death_pct))




