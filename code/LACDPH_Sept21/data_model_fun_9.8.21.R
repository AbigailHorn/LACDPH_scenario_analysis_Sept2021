

###################################################################################################
## COVID INPUT DATA
# latest_data: cumulative counts for "Htotcum","D","Vcum","Idetectcum","H_new","D_new"
# no_obs: number of observation days
## Read and process the data

latest_covid_data <- function(truncate=0){
  cum_file <- sort(dir_ls(data.dir, regexp = "cum_counts21_"), decreasing = TRUE)[1]
  latest_data = t(read.csv(cum_file, sep=",",stringsAsFactors = FALSE))
  colnames<-c("Htotcum","D","Vcum","Idetectcum","H_new","D_new","I_detect_new","V_new")
  nvars <- ncol(latest_data)
  colnames(latest_data) <- colnames
  latest_data <- as.data.frame(latest_data)
  latest_data <- latest_data[-1,]
  latest_data[1:nvars] <- sapply(latest_data[1:nvars],as.character)
  latest_data[1:nvars] <- sapply(latest_data[1:nvars],as.numeric)
  latest_data <- latest_data %>% dplyr::select(-V_new)
  no_obs <- nrow(latest_data)
  
  latest_data <- latest_data[c(1:(no_obs-truncate)),]
  
  ## Change date to number
  # step <- c(1:no_obs)
  # data <- cbind(step,latest_data)
  return(latest_data)
}


########################################################################################
## PLOTTING INPUT
########################################################################################

### IF PLOTTING ALL VARIABLES
all.variables <- c("S",
                   "I_detect_new",
                   "I",
                   "Idetectcum",
                   "Itot",
                   "Itotcum",
                   "H_new",
                   "Htot",
                   "Htotcum",
                   "Q",
                   "Qcum",
                   "V",
                   "Vcum",
                   "D_new",
                   "D",
                   "R"
)

vars.plus.R <- all.variables

### IF ONLY PLOTTING VARIABLES AGAINST DATA
only.vars.with.data <- c(
  "I_detect_new",
  "Idetectcum",
  "H_new",
  "Htotcum",
  "Vcum",
  "D",
  "D_new"
)

########################################################################################
## HELPER FUNCTIONS

## CONFIDENCE INTERVALS
posterior.CI <- function(posterior.var, round.by=3){
  median = quantile(posterior.var, c(.5), na.rm=TRUE)
  low_95 = quantile(posterior.var, c(.025), na.rm=TRUE)
  low_50 = quantile(posterior.var, c(.25), na.rm=TRUE)
  mean = mean(posterior.var)
  up_50 = quantile(posterior.var, c(.75), na.rm=TRUE)
  up_95 = quantile(posterior.var, c(.975), na.rm=TRUE)
  posterior.CI <- as.data.frame(cbind(low_95,low_50,median,mean,up_50,up_95))
  posterior.CI <- round(posterior.CI, digits=round.by)
  return(posterior.CI)
}

get.CFR.IFR.by.date <- function(traj.CI, CFR.or.IFR, date.in, round.by=4){
  if (CFR.or.IFR=="CFR"){
    state.name.in="CFRobs"
  }
  if (CFR.or.IFR=="IFR"){
    state.name.in="CFRactual"
  }
  posterior.CI.out <- traj.CI %>% filter(date %in% as.Date(date.in))
  posterior.CI.out <- posterior.CI.out %>% filter(state.name==state.name.in) %>% select(-c(state.name,N)) %>% mutate_if(is.numeric, round, digits=round.by)
  return(posterior.CI.out)
}

# FORMAT AS "mean(low_95,up_95)"
posterior.CI.FORMAT <- function(var.CI,use.mean=1){
  var.CI <- as.data.frame(var.CI)
  low_95 <- var.CI$low_95
  mean <- var.CI$mean
  up_95 <- var.CI$up_95
  median <- var.CI$median
  if (use.mean==1){
    var.95.CI <- paste0(mean, " (", low_95, ",", up_95,")")
  }
  if (use.mean==0){
    var.95.CI <- paste0(median, " (", low_95, ",", up_95,")")
  }
  return(var.95.CI)
}

var.format <- function(var.CI,use.mean.select){
  n.idx <- nrow(var.CI)
  var.format = as.data.frame(matrix(nrow=n.idx, ncol=1))
  for (i in 1:n.idx){
    var.format[i,] <- posterior.CI.FORMAT(var.CI[i,],use.mean=use.mean.select)
  }
  return(var.format)
}


########################################################################################
## SPECIFYING EPIDEMIC MODEL TO BE SIMULATED AND SCENARIOS
########################################################################################

replicates = 20
iter=10
time.steps=500
scenario.IDs = c("Scenario_EveryoneVac_Efficacy.85", "Scenario_NooneVac_Efficacy.85")

#traj.raw <- correlated.param.SIM(scenario.IDs=scenario.IDs,replicates=replicates,iter=iter,time.steps=time.steps)


correlated.param.SIM <- function(scenario.IDs,replicates,iter,time.steps) {
  
  all.scenarios = vector("list", length(scenario.IDs))
  
  for (this.scenario in 1:length(scenario.IDs)){
    
    print(paste0("Scenario: ", scenario.IDs[this.scenario]))
    
    ##########################################################################################
    ## Read in variables from csv
    p_tab <- read.csv(path(data.dir, paste0(scenario.IDs[this.scenario],".csv")), sep=",", header=TRUE)
    start_date <- p_tab$Value[which(p_tab$Var=="start_date")] %>% as.Date()
    #p_tab$Value <-gsub(",","",p_tab$Value)
    p_tab$Value <- as.numeric(as.character(p_tab$Value))
    p_tab$Lower <- as.numeric(as.character(p_tab$Lower))
    p_tab$Upper <- as.numeric(as.character(p_tab$Upper))
    
    ########## Fixed parameters
    everyone_vac <- p_tab$Value[which(p_tab$Var=="everyone_vac")]
    N_LA <- p_tab$Value[which(p_tab$Var=="N_LA")]
    N_eligible <- p_tab$Value[which(p_tab$Var=="N_eligible")]
    
    ########## Variables to be sampled over replicates
    ## Define functions to generate replicates
    cut01 <- function(x) pmax(pmin(x, 1), 0)
    getRandomNorm <- function(rep, value) {
      replicatedVar <- 
        rnorm(rep,
              mean=p_tab$Value[which(p_tab$Var==value)],
              sd=calcStDevFromLowerUpper(mean=p_tab$Value[which(p_tab$Var==value)],
                                         lower=p_tab$Lower[which(p_tab$Var==value)],
                                         upper=p_tab$Upper[which(p_tab$Var==value)]))
      return(replicatedVar)
    }
    calcStDevFromLowerUpper <- function(mean=mean, lower=lower, upper=upper, CI=.95) {
      z.stat <- qnorm((CI + (1-CI)/2))
      se.u <- (upper-mean)/z.stat
      se.l <- (mean-lower)/z.stat
      mean(c(se.u, se.l))
    }
    ## Get variables and replicates
    vac_efficacy <- cut01(getRandomNorm(replicates, "vac_efficacy"))
    N_vac <- getRandomNorm(replicates, "N_vac")
    isolate_prop <- cut01(getRandomNorm(replicates, "isolate_prop"))
    #isolate_prop[is.na(isolate_prop)] <- 0
    vac_prop <- N_vac/N_LA                                        # Proportion of LAC population vaccinated
    
    I_ini <- round(getRandomNorm(replicates, "I_ini"))
    E_ini <- round(getRandomNorm(replicates, "E_ini"))
    ## Get S_ini based on scenario 
    if(everyone_vac == 0){
      S_ini = N_LA*(1-isolate_prop)             # Susceptible are all except those that self-isolate.  
    }  else if(everyone_vac == .5){
      S_ini = N_LA*(1-isolate_prop)             # Everyone vaccinated is still susceptible, but at reduced transmission rate
    }  else if(everyone_vac == 1){
      S_ini = N_LA - N_eligible*vac_efficacy    # Everyone vaccinated is protected. Susceptible are all <12 yrs and those for whom vaccine was not effective.  
    }
    
    Alpha <- cut01(getRandomNorm(replicates, "Alpha"))
    Kappa <- cut01(getRandomNorm(replicates, "Kappa"))
    Delta <- cut01(getRandomNorm(replicates, "Delta"))
    r <- cut01(getRandomNorm(replicates, "r"))
    
    ##########################################################################################
    ########## Get Rt over time
    
    ## Evaluate Rt from csv file
    fn_t_readin_path <- path(data.dir, "fn_t_readin.csv")
    fn_t_readin = as.data.frame(read.csv(fn_t_readin_path, sep=",",stringsAsFactors = FALSE))
    Rt_dates <- as.Date(fn_t_readin$t)
    Rt_t <- round(as.numeric(Rt_dates - start_date) + start_time)
    Rt_y <- fn_t_readin$Rt
    
    ## Get Rt in vaxed and unvaxed populations using weighting equations
    Rt_uv_v <- function(Rt_eff, vac_efficacy, vac_prop, everyone_vac){
      w_v = vac_efficacy*vac_prop
      w_uv = 1-w_v
      Rt = Rt_eff * 1/w_uv  # For LA overall, adjusted to get Rt based on behavior alone, not including size of susceptible population
      Rt_v = round(Rt/(5-4*w_v),2)   # For vaccinated
      Rt_uv = round(5*Rt_v,2)        # For unvaccinated
      if (everyone_vac==1) return(Rt_v) else return(Rt_uv)
      if(everyone_vac == 0){
        return(Rt_uv)            # No one is vaccinated
      }  else if(everyone_vac == .5){
        return(Rt_v)             # Everyone vaxed but vaccinated can get infected
      }  else if(everyone_vac == 1){
        return(Rt_uv)            # Everyone vaxed and protected, except for those that were not included 
      }
    }
    
    ## Apply function to get Rt(t,y) over replicates of vac_efficacy and vac_prop
    Rt_vac=sapply(Rt_y,Rt_uv_v,vac_efficacy=vac_efficacy,vac_prop=vac_prop,everyone_vac=everyone_vac)               # Get Rt(t,y) over iterations - use this version if want to simulate over multiple Rt values
    #Rt_vac = mapply(Rt_uv_v,Rt_eff=Rt_y,vac_efficacy=vac_efficacy, vac_prop=vac_prop)   # Get Rt(t,y) - use this version if want one Rt_y for each Rt_t
    
    ########################################################################
    ########## Iterate over replicates and stochastic model realizations
    
    TEST.out <- vector("list", replicates)
    for (idx in 1:replicates) {
      Rt_sim = Rt_vac[idx,]
      r_sim = r[idx]
      Alpha_sim = Alpha[idx]
      Kappa_sim = Kappa[idx]
      Delta_sim = Delta[idx]
      S_ini_sim = S_ini[idx]
      E_ini_sim = E_ini[idx]
      I_ini_sim = I_ini[idx]
      
      ## COMPILE
      x <- seihqdr_generator(Rt_t=Rt_t, Rt_y=Rt_sim, r=r_sim, Alpha=Alpha_sim, Kappa=Kappa_sim, Delta=Delta_sim, S_ini=S_ini_sim, E_ini=E_ini_sim, I_ini=I_ini_sim, A_ini=3*I_ini_sim, R_ini=N_LA-S_ini_sim)
      
      ## SIMULATE
      TEST<-as.data.frame(plyr::rdply(iter, x$run(0:time.steps),.id="iter"))
      #TEST <- cbind(data.frame(date = TEST$step), TEST)
      #TEST$step <- NULL
      
      ## BIND INCLUDING OFFSETING OBSERVED DATA BY START DATE
      TEST.out[[idx]] <- cbind(data.frame(date=TEST$step, replicate = idx), TEST)
    }
    
    ## ADD TO DATAFRAME OVER ALL PARAMETER VALUES
    TEST.out <- do.call(rbind, TEST.out)
    TEST.out$step <- NULL
    all.scenarios[[this.scenario]] <- cbind(data.frame(scenario=scenario.IDs[this.scenario]), TEST.out)
  }
  
  scenarios.out = do.call(rbind, all.scenarios)
  return(scenarios.out)
}


########################################################################################
## GETTING MODEL OUTPUT + SUMMARY STATISTICS FUNCTION
########################################################################################
# num.to.sample <- 20
# ABC.out.mat <- ABC_out$param[1:num.to.sample,]
# par.vec.length <- num.to.sample
# iter <- 10
# time.steps <- 300
# vars.to.plot <- vars.plus.R

#traj.0 = model.output.to.plot.SIM(scenario.IDs, replicates=replicates,iter=iter,time.steps=time.steps)

model.output.to.plot.SIM <- function(scenario.IDs, replicates, iter, time.steps, start_date="2021-06-01") {
  
  ## Raw model output over each replicate and iteration
  traj.raw <- correlated.param.SIM(scenario.IDs, replicates=replicates,iter=iter,time.steps=time.steps)
  
  ## Add CFR and IFR
  traj <- dplyr::mutate(traj.raw, Itot=I+A, CFRobs=(D/Idetectcum), CFRactual=(D/(Itotcum)) )
  #traj <-  dplyr::select(traj,c(1:4,CFRobs,CFRactual,vars.to.plot))#,I_tot_new))
  rm(traj.raw)
  
  print("Starting CI calc")
  
  ## Get summary statistics
  df.traj <- reshape2::melt(traj, measure.vars = c(5:ncol(traj)), variable.name = "state.name")
  df.traj_dt <- as.data.table(df.traj)
  
  traj.CI <- df.traj_dt[, list(
    N=.N,
    mean = mean(value),
    median = quantile(value, c(.5),na.rm=TRUE),
    low_95 = quantile(value, c(.025),na.rm=TRUE),
    up_95 = quantile(value, c(.975),na.rm=TRUE),
    up_50 = quantile(value,.75,na.rm=TRUE),
    low_50 = quantile(value,.25,na.rm=TRUE)),
    by = c("date", "scenario", "state.name")]
  traj.CI <- as.data.frame(traj.CI)
  rm(traj)
  
  ## Align dates
  init.date = as.Date(start_date)
  traj.CI[["date"]] <- traj.CI[["date"]] + init.date
  
  return(traj.CI)
  
}


########################################################################################
########################################################################################
## FUNCTIONS FOR ABC
########################################################################################
########################################################################################

###################################################################################################
## "SUMMARY STATISTICS":
## The cumulative number of cases at all (trusted) time points
###################################################################################################

sum.stats.SIMTEST <- function(data){
  
  #print(no_obs)
  #no_obs <- nrow(data)
  #last_date <- as.numeric(as.Date(last_date_data) - as.Date("2020-03-01"))
  
  # Which values of variables to consider
  I.trust.n <- c(10:no_obs)  # The first 9 days of illness cases are unreliable/unavailable
  #H.trust.n <- c(17:last_date)  # The first 16 days of hospitalizations are unreliable/unavailable
  #V.trust.n <- c(19:last_date)  # The first 18 days of ventilation are unreliable/unavailable
  D.trust.n <- c(18:no_obs)  # The first 17 days of mortality are unreliable/unavailable
  #Hnew.trust.n <- c(19:last_date) # The first 18 days of new hospitalizations are unreliable/unavailable
  Dnew.trust.n <- c(28:no_obs) # The first 28 days of new deaths are unreliable/unavailable
  HQ.trust.n <- c(29:no_obs)
  
  ss.Icum <- data$Idetectcum[I.trust.n]
  ss.I <- data$I_detect_new[I.trust.n]
  #ss.H <- data$Htotcum[H.trust.n]
  #ss.V <- data$Vcum[V.trust.n]
  ss.D <- data$D[D.trust.n]
  #ss.Hnew <- data$H_new[Hnew.trust.n]
  ss.Dnew <- data$D_new[Dnew.trust.n]
  ss.Htot <- data$Htot[HQ.trust.n]
  ss.Q <- data$Q[HQ.trust.n]
  
  # print(length(ss.Icum))
  # print(length(ss.I))
  # print(length(ss.D))
  # print(length(ss.Dnew))
  # print(length(ss.Htot))
  # print(length(ss.Q))
  
  # Which variables to consider
  
  summarystats = c(ss.I, ss.Icum, ss.Htot, ss.Q, ss.Dnew, ss.D)
  
  return(summarystats)
}


###################################################################################################
## SIMULATION MODEL FUNCTION TO COMPUTE FOR ABC ALGORITHM
## A function implementing the model to be simulated
## It must take as arguments a vector of model parameter values par
## and it must return a vector of summary statistics (compartmental model variables from simulation)
###################################################################################################

model.1sim.stats.no.R <- function(par){
  
  R0 <- par[1]
  r1 <- par[2]
  start_time <- par[3]
  R0_redux1 <- par[4]
  Delta1 <- par[5]
  Alpha1 <- par[6]
  Kappa1 <- par[7]
  p_V <- par[8]
  R0_redux2 <- par[9]
  Delta2 <- par[10]
  Alpha2 <- par[11]
  Kappa2 <- par[12]
  r2 <- par[13]
  R0_redux3 <- par[14]
  
  ### BRING IN BETA_T ALPHA_T KAPPA_T DELTA_T FUNCTIONS
  fn_t_readin_code <- path(code.paper.dir, "fn_t_readin_code_FULL.R")
  source(fn_t_readin_code, local=TRUE)
  
  # length.B <- length(Beta_y)
  # Beta_y[length.B] <- Beta_y[length.B]*0.9
  # Beta_y[length.B] <- Beta_y[length.B-1]*0.9
  
  # print(paste0("Dimensions of Beta_y",(Beta_y)))
  # print(paste0("Dimensions of Beta_t",(Beta_t)))
  # print(paste0("Dimensions of r_y",(r_y)))
  # print(paste0("no_obs",(no_obs)))
  
  ### GENERATE SIMULATION
  x <- seihqdr_generator(Alpha_t=Alpha_t, Alpha_y=Alpha_y, Kappa_t=Kappa_t, Kappa_y=Kappa_y, Delta_t=Delta_t, Delta_y=Delta_y, Beta_t=Beta_t, Beta_y=Beta_y, r_t=r_t, r_y=r_y, S_ini=1e7, E_ini=10, p_QV=p_V)
  st <- start_time
  one_sim <- as.data.frame(x$run(0:(st+no_obs))[(st+1):(st+no_obs),])
  
  # print(head(one_sim))
  # print(tail(one_sim))
  #print(no_obs)
  
  ### SUMMARY STATISTICS COMPUTED ON MODEL OUTPUT:
  summarymodel <- sum.stats.SIMTEST(one_sim)
  
  #print(length(summarymodel))
  
  return(summarymodel)
}


