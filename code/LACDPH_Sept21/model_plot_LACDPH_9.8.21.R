

########################################################################################
## SPECIFYING EPIDEMIC MODEL TO BE SIMULATED AND SCENARIOS
########################################################################################

# replicates = 20
# iter=10
# time.steps=500
# scenario.IDs = c("Scenario_EveryoneVac_Efficacy.85", "Scenario_NooneVac_Efficacy.85")

#traj.raw <- correlated.param.SIM(scenario.IDs=scenario.IDs,replicates=replicates,iter=iter,time.steps=time.steps)

correlated.param.SIM <- function(scenario.IDs,replicates,iter,time.steps) {
  
  all.scenarios = vector("list", length(scenario.IDs))
  
  for (this.scenario in 1:length(scenario.IDs)){
    
    print(paste0("Scenario: ", scenario.IDs[this.scenario]))
    
    ##########################################################################################
    ## Read in variables from csv
    p_tab <- read.csv(path(data.dir.update, paste0(scenario.IDs[this.scenario],".csv")), sep=",", header=TRUE)
    start_date <- p_tab$Value[which(p_tab$Var=="start_date")] %>% as.Date()
    #p_tab$Value <-gsub(",","",p_tab$Value)
    p_tab$Value <- as.numeric(as.character(p_tab$Value))
    p_tab$Lower <- as.numeric(as.character(p_tab$Lower))
    p_tab$Upper <- as.numeric(as.character(p_tab$Upper))
    
    ########## Fixed parameters
    everyone_vac <- p_tab$Value[which(p_tab$Var=="everyone_vac")]
    N_LA <- p_tab$Value[which(p_tab$Var=="N_LA")]
    N_eligible <- p_tab$Value[which(p_tab$Var=="N_eligible")]
    prop_eligible <- N_eligible/N_LA
    
    Rt_factor <- p_tab$Value[which(p_tab$Var=="Rt_factor")]
    
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
    Rt_t <- round(as.numeric(Rt_dates - start_date))
    Rt_y <- fn_t_readin$Rt * Rt_factor  # Modified by Rt_factor increase or decrease
    
    #print(Rt_y/(1-vac_prop))
    
    ## Get Rt in vaxed and unvaxed populations using weighting equations
    Rt_uv_v <- function(Rt_eff, vac_efficacy, vac_prop, everyone_vac, prop_eligible){
      w_v = rep(.525,100)  #vac_prop #*vac_efficacy #rep(.55,100) 
      w_uv = 1-w_v
      Rt = Rt_eff * 1/w_uv  # For LA overall, adjusted to get Rt based on behavior alone, not including size of susceptible population
      Rt_v = round(Rt/(5-4*w_v),2)  # round(Rt/(5-4*vac_prop),2)  # For vaccinated
      Rt_uv = round(5*Rt_v,2)        # For unvaccinated
      if(everyone_vac == 0){
        return(Rt_uv)            # No one is vaccinated
      }  else if(everyone_vac == .5){
        return(prop_eligible*Rt_v + (1-prop_eligible)*Rt_uv)   # Weighted average of Rt based on prop_eligible+vaccinated vs. prop_ineligible  
      }  else if(everyone_vac == 1){
        return(Rt_uv)            # Everyone vaxed and "successful" is fully protected, all "unsuccessful" + ineligible are fully unprotected 
      }
    }
    
    ## Apply function to get Rt(t,y) over replicates of vac_efficacy and vac_prop
    Rt_vac=sapply(Rt_y,Rt_uv_v,vac_efficacy=vac_efficacy,vac_prop=vac_prop,everyone_vac=everyone_vac,prop_eligible=prop_eligible)               # Get Rt(t,y) over iterations - use this version if want to simulate over multiple Rt values
    #Rt_vac = mapply(Rt_uv_v,Rt_eff=Rt_y,vac_efficacy=vac_efficacy, vac_prop=vac_prop)   # Get Rt(t,y) - use this version if want one Rt_y for each Rt_t
    
    print(Rt_vac[1,])

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
  traj.CI[["date"]] <- traj.CI[["date"]] + as.Date(start_date)
  
  return(traj.CI)
  
}

########################################################################################
## PLOTTING ALL FACETED FUNCTION
########################################################################################

# data.in = la_data
# time.steps.plot = as.numeric(as.Date("2022-03-31") - as.Date(start_date))
# date.offset.4plot = 15
# 
# plot.model.data.all(traj.CI=traj.plot, data.in=data.in, time.steps.plot=time.steps.plot, vars.to.plot=vars.to.plot)

plot.model.data.all <- function(traj.CI, data.in, time.steps.plot, vars.to.plot) {
  
  ## Select only more recent dates
  startDatePlot <- start_date
  endDatePlot <- startDatePlot + time.steps.plot # - 40
  traj.CI <- traj.CI %>% dplyr::filter(date < endDatePlot) #dplyr::filter(date > startDatePlot-1) %>% dplyr::filter(date < endDatePlot)
  
  ## Filter to variables selected to plot
  traj.CI <- traj.CI %>% dplyr::filter(state.name==c(vars.to.plot))
  
  if(!is.null(data.in)){
    data <- data.in %>% dplyr::filter(date > startDatePlot-1)
    data <- reshape2::melt(data, measure.vars = c(2:ncol(data)), variable.name = "state.name")
    dataCum = data %>% filter(state.name %in% c("D", "Idetectcum")) %>% group_by(state.name) %>% mutate(value = value - dplyr::first(value)) %>% as.data.frame()
    dataNew = data %>% filter(state.name %in% c("D", "Idetectcum") == FALSE) %>% as.data.frame()
    data=as.data.frame(rbind(dataCum,dataNew))
    
    # rebase0 <- function(data,thisVar){
    #   data_thisVar = data %>% dplyr::filter(state.name %in% thisVar) %>% mutate(value = value-min(value))
    #   return(rbind(data_thisVar,data))
    # }
    # data = rebase0(data, "D")
    # data = rebase0(data,"Idetectcum")
  }
  
  ## PLOTTING
  #traj.CI.line <- reshape2::melt(traj.CI[c("date", "state.name", "mean", "median")], id.vars = c("date", "state.name"))
  traj.CI.line <- reshape2::melt(traj.CI[c("date", "state.name", "median")], id.vars = c("date", "state.name"))
  traj.CI.area <- reshape2::melt(traj.CI[c("date", "state.name", "low_95", "low_50", "up_50", "up_95")], id.vars = c("date", "state.name"))
  traj.CI.area$type <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj.CI.area$CI <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj.CI.area$variable <- NULL
  traj.CI.area <- reshape2::dcast(traj.CI.area, "date+state.name+CI~type")
  
  p <- ggplot(transform(traj.CI.area, state.name = factor(state.name, levels=vars.to.plot)))
  
  longnames <- c("Susceptible",
                 "New Obs. Infected",
                 "Current Obs. Infected",
                 "Cum. Obs. Infected",
                 "Current Tot. Infected",
                 "Cum. Tot. Infected",
                 "New in Hospital",
                 "Current in Hospital",
                 "Cum. in Hospital",
                 "Current in ICU",
                 "Cum. in ICU",
                 "Current Ventilation",
                 "Cum. Ventilation",
                 "New Deaths",
                 "Cum. Deaths",
                 "Recovered")
  
  names(longnames) <-  c("S",
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
  
  p <- p + facet_wrap(~state.name, labeller = labeller(state.name = longnames), scales = "free_y")
  p <- p + geom_ribbon(data = traj.CI.area, aes_string(x = "date", ymin = "low", ymax = "up", alpha = "CI", fill = "state.name"),show.legend = c(fill=FALSE))
  p <- p + geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable", colour = "state.name"), size = 1, show.legend = c(colour=FALSE))
  
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.20, "50" = 0.50), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  
  
  ## ADD DATA
  if(!is.null(data.in)){
    p <- p + geom_point(data = data, aes_string(x = "date", y = "value"), size = .5, colour = "black")
  }
  
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot-31,endDatePlot)), date_breaks = "1 month" , date_labels = "%d-%b-%y")
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 12, face = "bold"))
  p <- p + ylab("Numbers in Compartments") + xlab(NULL)
  p <- p + scale_y_continuous(labels = scales::comma)
  p
}


