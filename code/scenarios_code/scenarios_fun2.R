

# checking.output <- correlated.param.SCENARIOS.out[[1]]
# 
# all.scenarios <- do.call(rbind,correlated.param.SCENARIOS.out)
# 
# geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable", colour = "state.name"))
# 
# bp <- ggplot(all.scenarios %>% filter(state.name %in% c("D","CFRobs","H_new")), aes(x=date, y=median)) + 
#   geom_line(scales='free')+
#   facet_grid(state.name ~ scenario.id, scales='free')
# 
# view(checking.output)
# 
# weighted.avg.scenarios.overall <- 
#   weighted.avg.protect.SCENARIOS(X.mat=X.mat, freq.PREV.q = freq.PREV.q, 
#                                  freq.LAC.obs.age=freq.LAC.obs.age, logit.SEIR.est=logit.SEIR.est,
#                                  psi.mat=psi.mat,percent.to.remove=percent.to.remove,
#                                  factor.to.remove=factor.to.remove)
# 
# protect.scenarios.mat <- matrix(NA, nrow=3, ncol=2)
# colnames(protect.scenarios.mat) <- c("index","protect.name")
# protect.scenarios.mat[,"index"] <- c(1:3)
# protect.scenarios.mat[,"protect.name"] <- rownames(weighted.avg.scenarios.overall)
# protect.scenarios.mat <- as.data.frame(protect.scenarios.mat)
# 
# NPI.scenarios.mat <- matrix(NA, nrow=3, ncol=3)
# colnames(NPI.scenarios.mat) <- c("index","NPI.name","Beta_y_vals")
# NPI.scenarios.mat[,"index"] <- c(1:3)
# NPI.scenarios.mat[,"NPI.name"] <- c("obs","moderate","nothing")
# NPI.scenarios.mat[,"Beta_y_vals"] <- c("Beta_y","Beta_y_moderate","Beta_y_nothing")
# NPI.scenarios.mat <- as.data.frame(NPI.scenarios.mat)
# 
# correlated.param.SCENARIOS.out <- correlated.param.SCENARIOS(ABC.out.mat=ABC.out.mat[c(1:5),],iter=5,time.steps=200,weighted.avg.scenarios=weighted.avg.scenarios.overall,protect.scenarios.mat=protect.scenarios.mat, NPI.scenarios.mat=NPI.scenarios.mat)
# 

########################################################################################
## Get Alpha, Kappa, Delta with populations protected
########################################################################################

weighted.avg.protect.SCENARIOS <- function(X.mat, freq.PREV.q, freq.LAC.obs.age, logit.SEIR.est, psi.mat, percent.to.remove, factor.to.remove){
  
  n.dates <- ncol(freq.LAC.obs.age)
  
  #percent.to.remove <- c(1,0.5) # Remove 100%
  #factor.to.remove <- 4  # Age.65.
  weighted.avg <- vector("list", n.dates)
  weighted.avg.scenarios <- vector("list",length(percent.to.remove))
  weighted.avg.scenarios.overall <- NULL
  
  for (remove in 1:length(percent.to.remove)){
    for (t in 1:n.dates){
      percent.remove <- percent.to.remove[remove]
      time = t
      name.date <- colnames(freq.LAC.obs.age)[t]
      
      freq.I <- freq.ILL(X.mat, freq.PREV.q, freq.LAC.obs.age, time = time)
      
      Pr.H.I.q <- get.Pr.q(model=1, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.I[[2]],  psi.mat)
      freq.H <- get.population.freq(prob.q.vec.IN = Pr.H.I.q, freq.q.IN = freq.I[[1]], X.mat)
      
      Pr.Q.H.q <- get.Pr.q(model=2, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.H[[2]],  psi.mat)
      freq.Q <- get.population.freq(prob.q.vec.IN = Pr.Q.H.q, freq.q.IN = freq.H[[1]], X.mat)
      
      Pr.D.Q.q <- get.Pr.q(model=3, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.Q[[2]],  psi.mat)
      freq.D <- get.population.freq(prob.q.vec.IN = Pr.D.Q.q, freq.q.IN = freq.Q[[1]], X.mat)
      
      ##
      freq.I.filter <- freq.I[[1]]*(1-percent.remove*X.mat[,factor.to.remove])
      weighted.avg.H.filter <- t(Pr.H.I.q) %*% freq.I.filter
      
      freq.H.filter <- (get.population.freq(Pr.H.I.q, freq.I.filter, X.mat))[[1]]
      weighted.avg.Q.filter <- t(Pr.Q.H.q) %*% freq.H.filter
      
      freq.Q.filter <- (get.population.freq(Pr.Q.H.q, freq.H.filter, X.mat))[[1]]
      weighted.avg.D.filter <- t(Pr.D.Q.q) %*% freq.Q.filter
      
      weighted.avg.t <- cbind( weighted.avg.H.filter, weighted.avg.Q.filter, weighted.avg.D.filter )
      colnames(weighted.avg.t) <- c("Alpha.","Kappa.","Delta.")
      colnames(weighted.avg.t) <- paste0( colnames(weighted.avg.t), name.date)
      weighted.avg[[t]] <- weighted.avg.t
      ##
    }
    weighted.avg.scenarios[[remove]] <- do.call(cbind, weighted.avg)
  }
  weighted.avg.scenarios.overall <- do.call(rbind,weighted.avg.scenarios)
  #weighted.avg.scenarios.overall <- as.data.frame(weighted.avg.scenarios)
  rownames(weighted.avg.scenarios.overall) <- paste0("Protect.",percent.to.remove*100)
  
  return(weighted.avg.scenarios.overall)
  
}


########################################################################################
## SPECIFYING EPIDEMIC MODEL TO BE SIMULATED AND SCENARIOS
########################################################################################

correlated.param.SCENARIOS <- function(ABC.out.mat,iter,time.steps,weighted.avg.scenarios,protect.scenarios.mat, NPI.scenarios.mat) {
  
  ##########################################################################################
  ## Read in csv files with Beta_t, Alpha_t, Kappa_t, Delta_t
  
  fn_t_readin_path <- path(data.dir, "fn_t_readin.csv")
  fn_t_readin = as.data.frame(read.csv(fn_t_readin_path, sep=",",stringsAsFactors = FALSE))
  
  ## Get r_t
  r_t_dates <- as.Date(fn_t_readin$r_t)
  r_t_dates <- na.omit(r_t_dates)
  r_t_dates[1] <- r_t_dates[1]-start_time
  r_t <- round(as.numeric(r_t_dates - as.Date("2020-03-01")) + start_time)
  
  ## Get Beta_t
  start_time <- 45
  Beta_t_dates <- as.Date(fn_t_readin$Beta_t)
  Beta_t_dates[1] <- Beta_t_dates[1]-start_time
  Beta_t <- round(as.numeric(Beta_t_dates - as.Date("2020-03-01")) + start_time)
  
  ## Get Alpha_t
  alpha_t_readin_path <- path(data.dir, "alpha_t_readin.csv")
  alpha_t_readin = as.data.frame(read.csv(alpha_t_readin_path, sep=",",stringsAsFactors = FALSE))
  
  Alpha_t_dates <- as.Date(alpha_t_readin$Alpha_t)
  Alpha_t_dates[1] <- Alpha_t_dates[1]-start_time
  Alpha_t <- round(as.numeric(Alpha_t_dates - as.Date("2020-03-01")) + start_time)
  Kappa_t <- Alpha_t
  Delta_t <- Alpha_t

  ##########################################################################################
  ## INITIALIZE SCENARIOS.out
  SCENARIOS.out <- vector("list", 9)
  scenario.idx <- 1
  
  for (protect.idx in 1:3){
    
    protect.scenario = protect.idx
    protect.name <- protect.scenarios.mat[protect.idx, "protect.name"]
    
    for (NPI.idx in 1:3){
      
      print(paste0("Starting scenario ", scenario.idx))
      
      NPI.scenario = NPI.idx
      NPI.name <- NPI.scenarios.mat[NPI.idx, "NPI.name"]
      Beta_y_vals <- NPI.scenarios.mat[NPI.idx, "Beta_y_vals"]
      
      TEST.out <- vector("list", nrow(ABC.out.mat))
      
      for (idx in 1:nrow(ABC.out.mat)){
        
        R0 <- ABC.out.mat[idx,1]
        r1 <- ABC.out.mat[idx,2]
        start_time <- round(ABC.out.mat[idx,3])
        R0_redux1 <- ABC.out.mat[idx,4]
        Delta1 <- ABC.out.mat[idx,5]
        Alpha1 <- ABC.out.mat[idx,6]
        Kappa1 <- ABC.out.mat[idx,7]
        p_V <- ABC.out.mat[idx,8]
        R0_redux2 <- ABC.out.mat[idx,9]
        Delta2 <- ABC.out.mat[idx,10]
        Alpha2 <- ABC.out.mat[idx,11]
        Kappa2 <- ABC.out.mat[idx,12]
        r2 <- ABC.out.mat[idx,13]
        
        ##########################################################################################
        ## Alpha Kappa Delta t1 and t2 if not using the "Protect.0" scenario
        
        if (protect.name!="Protect.0"){
         
          Alpha1 <- weighted.avg.scenarios[protect.scenario,"Alpha.t1"]
          Alpha2 <- weighted.avg.scenarios[protect.scenario,"Alpha.t2"]
          
          Kappa1 <- weighted.avg.scenarios[protect.scenario,"Kappa.t1"]
          Kappa2 <- weighted.avg.scenarios[protect.scenario,"Kappa.t2"]
          
          Delta1 <- weighted.avg.scenarios[protect.scenario,"Delta.t1"]
          Delta2 <- weighted.avg.scenarios[protect.scenario,"Delta.t2"]
          
        }
        
        ################################
        ## r_y
        r_y_chr <- fn_t_readin$r_y
        assign("r1",r1)
        assign("r2", r1)
        
        r_y <- as.vector(length(r_t))
        for (z in 1:length(r_t)){
          r_y[z] = get(r_y_chr[z])
        }
        
        
        ## Alpha_t
        
        Alpha_y_chr <- alpha_t_readin$Alpha_y
        assign("Alpha1",Alpha1)
        assign("Alpha2", Alpha2)
        
        Alpha_y <- as.vector(length(Alpha_t))
        for (z in 1:length(Alpha_t)){
          Alpha_y[z] = get(Alpha_y_chr[z])
        }
        
        ## Kappa_t
        
        Kappa_y_chr <- alpha_t_readin$Kappa_y
        assign("Kappa1",Kappa1)
        assign("Kappa2", Kappa2)
        
        Kappa_y <- as.vector(length(Alpha_t))
        for (z in 1:length(Alpha_t)){
          Kappa_y[z] = get(Kappa_y_chr[z])
        }
        
        ## Delta_t
        
        Delta_y_chr <- alpha_t_readin$Delta_y
        assign("Delta1",Delta1)
        assign("Delta2", Delta2)
        
        Delta_y <- as.vector(length(Alpha_t))
        for (z in 1:length(Alpha_t)){
          Delta_y[z] = get(Delta_y_chr[z])
        }
        
        ##########################################################################################
        ## Beta_y
        
        # if (NPI.scenario == 1){ Beta_y_vals <- "Beta_y" }
        # if (NPI.scenario == 2) { Beta_y_vals <- "Beta_y_moderate" }
        # if (NPI.scenario == 3) {Beta_y_vals <- "Beta_y_nothing" }
        
        Beta_y_vals <- as.character(NPI.scenarios.mat[NPI.scenario,"Beta_y_vals"])
        
        mu_y_chr <- fn_t_readin[,Beta_y_vals]
        assign("mu.0",1)
        assign("mu.1", R0_redux1)
        assign("mu.2", R0_redux2)
        
        mu_y <- as.vector(length(Beta_t))
        for (z in 1:length(Beta_t)){
          mu_y[z] = get(mu_y_chr[z])
        }
        R0_y <- R0*mu_y
        
        ## Get Beta_y as a function of R0, R0_redux, r, and Alpha
        
        Br.function <- function(R0.in, r.in, Alpha.in){
          d_IH <- 10   #days between illness onset and hospitalization
          d_IR <- 7    #days between illness onset and recovery (hospitalization not required)
          Br <- R0.in * ( 1 / ( (r.in/ ((Alpha.in/d_IH) + ((1-Alpha.in)/d_IR)))  + (1-r.in)*d_IR ))
          return(Br)
        }
        
        Beta_y<- c(
          Br.function(R0.in<-R0_y[1], r.in<-r1, Alpha.in<-Alpha1) ,
          Br.function(R0.in<-R0_y[2], r.in<-r1, Alpha.in<-Alpha1),
          Br.function(R0.in<-R0_y[3], r.in<-r1, Alpha.in<-Alpha1),
          Br.function(R0.in<-R0_y[4], r.in<-r1, Alpha.in<-Alpha1),
          Br.function(R0.in<-R0_y[5], r.in<-r1, Alpha.in<-Alpha1),
          Br.function(R0.in<-R0_y[6], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[7], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[8], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[9], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[10], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[11], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[12], r.in<-r2, Alpha.in<-Alpha2),
          Br.function(R0.in<-R0_y[13], r.in<-r2, Alpha.in<-Alpha2)
        )
        
        ##########################################################################################
        ## RUN THE MODEL
        
        ## COMPILE 
        x <- seihqdr_generator(Alpha_t=Alpha_t, Alpha_y=Alpha_y, Kappa_t=Kappa_t, Kappa_y=Kappa_y, Delta_t=Delta_t, Delta_y=Delta_y, Beta_t=Beta_t, Beta_y=Beta_y, r_t=r_t, r_y=r_y, S_ini=1e7, E_ini=10, p_QV=p_V)
        
        ## SIMULATE
        TEST<-as.data.frame(plyr::rdply(iter, x$run(0:time.steps),.id="iter"))
        
        ## BIND INCLUDING OFFSETING OBSERVED DATA BY START DATE
        TEST.out[[idx]] <- cbind(data.frame(par.id = idx, date = -start_time+TEST$step), TEST)
        
        ## BIND INCLUDING OFFSETING OBSERVED DATA BY START DATE
        TEST.out[[idx]] <- cbind(data.frame(scenario.id = paste0(scenario.idx, "_", NPI.name, "_", protect.name ) , protect.id = protect.name, NPI.id = NPI.name, par.id = idx, date = -start_time+TEST$step), TEST)
        #TEST.out[[idx]] <- cbind(data.frame(scenario.id = scenario.idx, protect.id = protect.name, NPI.id = NPI.name, par.id = idx, date = -start_time+TEST$step), TEST)
        
      }  # end over idx
      
      ## Add to a dataframe over all idx
      TEST.out <- do.call(rbind, TEST.out)
      
      ## Get CI for scenario
      SCENARIO.CI <- scenarios.get.CI(TEST.out = TEST.out)
        
      ## Put TEST.out dataframe into SCENARIOS.out
      # SCENARIOS.out[[scenario.idx]] <- TEST.out#SCENARIO.CI
      SCENARIOS.out[[scenario.idx]] <- SCENARIO.CI
      
      
      ##
      scenario.idx <- scenario.idx + 1 
      
      rm(TEST.out)
      
    } # end over y NPI.scenario
  } # end over x protect.scenario protect.idx
  
  # Save memory

  # Output
  return(SCENARIOS.out)
}


########################################################################################
## FUNCTION TO GET CI FOR EACH SCENARIO (embedded in correlated.param.SCENARIOS)
########################################################################################

scenarios.get.CI <- function(TEST.out) {
  
  library(data.table)
  init.date.data="2020-03-01"
  
  ### Add CFR and IFR to the list (EXTRA STEP NOW THAT THIS IS BEING USED ALSO FOR summary_table)
  traj <- dplyr::mutate(TEST.out, Itot=I+A, CFRobs=(D/Idetectcum), CFRactual=(D/(Itotcum)) )
  traj <-  dplyr::select(traj,c(1:7,CFRobs,CFRactual,vars.plus.R))
  ###
  
  ## TO SAVE MEMORY
  rm(TEST.out)
  
  print("Starting CI calc")
  
  ### MELTING AND APPLYING SUMMARY STAT FUNCTIONS
  df.traj <- reshape2::melt(traj, measure.vars = c(8:ncol(traj)), variable.name = "state.name")
  df.traj_dt <- as.data.table(df.traj)
  
  traj.CI <- df.traj_dt[, list(
    N=.N,
    mean = mean(value),
    median = quantile(value, c(.5),na.rm=TRUE),
    low_95 = quantile(value, c(.025),na.rm=TRUE),
    up_95 = quantile(value, c(.975),na.rm=TRUE),
    up_50 = quantile(value,.75,na.rm=TRUE),
    low_50 = quantile(value,.25,na.rm=TRUE)),
    by = c("date", "state.name","scenario.id","protect.id","NPI.id")]
  traj.CI <- as.data.frame(traj.CI)
  
  ## TO ALIGN DATES: MODEL
  init.date = init.date.data #"2020-01-03"
  init.date <- as.Date(init.date)
  traj.CI[["date"]] <- traj.CI[["date"]] + init.date
  
  return(traj.CI)
  
}

