
get.PAR.tables <- function(ABC_out, traj.CI, date.in){
  
  #############################################
  #############################################
  #############################################
  ## PARAMETER INPUTS
  
  ABC.par <- ABC_out$param
  
  ## Confidence intervals
  posterior.CI <- function(posterior.var, round.by=4){
    median = quantile(posterior.var, c(.5), na.rm=TRUE)
    low_95 = quantile(posterior.var, c(.025), na.rm=TRUE)
    low_50 = quantile(posterior.var, c(.25), na.rm=TRUE)
    mean = mean(posterior.var)
    up_50 = quantile(posterior.var, c(.75), na.rm=TRUE)
    up_95 = quantile(posterior.var, c(.975), na.rm=TRUE)
    posterior.CI <- as.data.frame(cbind(low_95,low_50,mean,median,up_50,up_95))
    posterior.CI <- round(posterior.CI, digits=round.by)
    return(posterior.CI)
  }
  
  ## Put fn(t) in format for plots
  format.4.plot <- function(fn_t, fn_y_chr, fn.posterior.CI, fn.name){
    fn_y <- as.data.frame(matrix(nrow=length(fn_t), ncol=ncol(fn.posterior.CI) ))
    for (i in 1:length(fn_t)){
      fn_y[i,] = get(fn_y_chr[i])
    }
    colnames(fn_y) <- colnames(fn.posterior.CI)
    
    fn_plot <- as.vector(fn_y)
    rownames(fn_plot) <- 1:length(fn_t)
    fn_plot$date <- fn_t
    fn_plot$state.name <- rep(fn.name, length(fn_t))
    
    return(fn_plot)
  }
  
  #################################
  # GET VARIABLES AND APPLY CI
  ABC.par.CI <- apply(ABC_out$param, MARGIN=2, FUN=posterior.CI)
  
  #############################################
  ## Alpha(t) Kappa(t) Delta(t)
  
  Alpha1.CI <- ABC.par.CI[[6]]
  Alpha2.CI <- ABC.par.CI[[11]]
  Alpha3.CI <- Alpha2.CI*1.2
  Alpha4.CI <- Alpha2.CI*0.7
  Kappa1.CI <- ABC.par.CI[[7]]
  Kappa2.CI <- ABC.par.CI[[12]]
  Kappa3.CI <- Kappa2.CI*1.25
  Kappa4.CI <- Kappa2.CI*1.1
  Delta1.CI <- ABC.par.CI[[5]]
  Delta2.CI <- ABC.par.CI[[10]]
  Delta3.CI <- Delta2.CI*1.1
  Delta4.CI <- Delta2.CI*1.2
  
  # GET ORDER OF VALUES
  start_time = round(mean(ABC.par[,3]))
  alpha_t_readin_path <- path(data.dir, "alpha_t_readin.csv")
  alpha_t_readin = as.data.frame(read.csv(alpha_t_readin_path, sep=",",stringsAsFactors = FALSE))
  Alpha_t_dates <- as.Date(alpha_t_readin$Alpha_t)
  Alpha_t_dates[1] <- Alpha_t_dates[1]-start_time
  Alpha.t <- Alpha_t_dates
  Alpha.t[length(Alpha.t)] <- endDatePlot-1
  
  # ALPHA
  Alpha.chr <- alpha_t_readin$Alpha_y
  assign("Alpha1",Alpha1.CI)
  assign("Alpha2", Alpha2.CI)
  assign("Alpha3", Alpha3.CI)
  assign("Alpha4", Alpha4.CI)
  Alpha_plot <- format.4.plot(fn_t = Alpha.t, fn_y_chr = Alpha.chr, fn.posterior.CI=Alpha1.CI, fn.name="Alpha_t" )
  
  # KAPPA
  Kappa.chr <- alpha_t_readin$Kappa_y
  assign("Kappa1",Kappa1.CI)
  assign("Kappa2", Kappa2.CI)
  assign("Kappa3", Kappa3.CI)
  assign("Kappa4", Kappa4.CI)
  Kappa_plot <- format.4.plot(fn_t = Alpha.t, fn_y_chr = Kappa.chr, fn.posterior.CI=Kappa1.CI, fn.name="Kappa_t" )
  
  # DELTA
  Delta.chr <- alpha_t_readin$Delta_y
  assign("Delta1",Delta1.CI)
  assign("Delta2", Delta2.CI)
  assign("Delta3", Delta3.CI)
  assign("Delta4", Delta4.CI)
  Delta_plot <- format.4.plot(fn_t = Alpha.t, fn_y_chr = Delta.chr, fn.posterior.CI=Delta1.CI, fn.name="Delta_t" )
  
  #############################################
  ## R(t)
  
  # GET ORDER OF VALUES
  out_R0 <- ABC.par[,1]
  out_R0redux1<- ABC.par[,4]
  out_R0redux2<- ABC.par[,9]
  out_R0redux3<- ABC.par[,14]
  R0_x_redux1 <- out_R0*out_R0redux1
  R0_x_redux2 <- out_R0*out_R0redux2
  R0_x_redux3 <- out_R0*out_R0redux3
  
  # GET QUANTILES FOR VARIABLE
  R0.CI <- posterior.CI(out_R0,4)
  R0.redux1.CI <- posterior.CI(R0_x_redux1,4)
  R0.redux2.CI <- posterior.CI(R0_x_redux2,4)
  R0.redux3.CI <- posterior.CI(R0_x_redux3,4)
  
  # ASSIGN ADDITIONAL R0 VALUES
  R0.redux4.CI <- R0.redux3.CI*1.2
  R0.redux5.CI <- R0.redux3.CI*0.52
  
  # GET ORDER OF VALUES
  fn_t_readin_path <- path(data.dir, "fn_t_readin.csv")
  fn_t_readin = as.data.frame(read.csv(fn_t_readin_path, sep=",",stringsAsFactors = FALSE))
  Beta_t_dates <- as.Date(fn_t_readin$Beta_t)
  Beta_t_dates[1] <- Beta_t_dates[1]-start_time
  Rt.t <- Beta_t_dates
  Rt.t[length(Rt.t)] <- as.Date("2022-01-01")
  
  Rt.chr <- fn_t_readin$Beta_y
  assign("mu.0",R0.CI)
  assign("mu.1", R0.redux1.CI)
  assign("mu.2", R0.redux2.CI)
  assign("mu.3",R0.redux3.CI)
  assign("mu.4",R0.redux4.CI)
  assign("mu.5",R0.redux5.CI)
  
  # PUT IN FORMAT FOR PLOTTING
  Rt_plot <- format.4.plot(fn_t = Rt.t, fn_y_chr = Rt.chr, fn.posterior.CI=R0.CI, fn.name="Rt" )
  Rt_plot$step = as.numeric(Rt_plot$date - as.Date("2020-03-01"))
  
  ## Interpolate to get values for each day
  interpolate.fn <- function(plot.in, col.idx){
    x.vals <- plot.in[,"step"]
    y.vals <- plot.in[,col.idx]
    interpolate.out <- approx(x=x.vals, y=y.vals, method="linear", n=max(x.vals)-min(x.vals)+1)
    return(as.data.frame(interpolate.out))
  }
  cols.to.get <- c(1:6)
  n.cols <- length(cols.to.get)
  Rt_add <- vector("list", length=n.cols)
  for (col.idx in cols.to.get){
    interpolate.out <- interpolate.fn(Rt_plot, col.idx)
    Rt_add[[col.idx]] <- interpolate.out$y
  }
  Rt_add <- do.call(cbind, Rt_add)
  x.vals <- Rt_plot$step
  Rt_add <- cbind(min(x.vals):max(x.vals),Rt_add) %>% as.data.frame()
  colnames(Rt_add) <- c("step",colnames(Rt_plot[1:6]))
  Rt_add$date <- Rt_add$step + as.Date("2020-03-01")
  
  #################################################
  ## Multiply R(t) by fraction of susceptibles to get R(t)_effective
  traj.S <- filter(traj.0, state.name=="S")
  
  ###########################################################################
  ### MULTIPLYING EACH INTERVAL
  
  TEST <- Rt_add
  TEST = TEST %>% mutate(low_95 = low_95* (traj.S$low_95 / 1e7)  )
  TEST = TEST %>% mutate(low_50 = low_50* (traj.S$low_50 / 1e7)  )
  TEST = TEST %>% mutate(median = median* (traj.S$median / 1e7)  )
  TEST = TEST %>% mutate(mean = mean* (traj.S$mean / 1e7) )
  TEST = TEST %>% mutate(up_50 = up_50* (traj.S$up_50 / 1e7)  )
  TEST = TEST %>% mutate(up_95 = up_95* (traj.S$up_95 / 1e7)  )
  TEST = TEST[c(1:min(nrow(traj.S),nrow(TEST))),]
  TEST$state.name = rep("Rt_eff", times = nrow(TEST))
  
  #############################################
  ## r(t)
  
  # GET ORDER OF VALUES
  r1.CI <- ABC.par.CI[[2]]
  r2.CI <- ABC.par.CI[[13]]
  
  # GET ORDER OF VALUES
  fn_t_readin_path <- path(data.dir, "fn_t_readin.csv")
  fn_t_readin = as.data.frame(read.csv(fn_t_readin_path, sep=",",stringsAsFactors = FALSE))
  r_t_dates <- as.Date(fn_t_readin$r_t)
  r_t_dates[1] <- r_t_dates[1]
  r_t_dates <- na.omit(r_t_dates)
  r.t <- r_t_dates
  r.t[length(r.t)] <- endDatePlot-1
  
  r.chr <- fn_t_readin$r_y
  assign("r1",r1.CI)
  assign("r2", r2.CI)
  
  # PUT IN FORMAT FOR PLOTTING
  r_plot <- format.4.plot(fn_t = r.t, fn_y_chr = r.chr, fn.posterior.CI=r1.CI, fn.name="r_t" )
  
  #############################################
  #############################################
  #############################################
  ## FUNCTIONS
  
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
  
  ## Interpolate to get values for each day
  interpolate.fn <- function(plot.in, col.idx){
    x.vals <- plot.in[,"step"]
    y.vals <- plot.in[,col.idx]
    interpolate.out <- approx(x=x.vals, y=y.vals, method="linear", n=max(x.vals)-min(x.vals)+1)
    return(as.data.frame(interpolate.out))
  }
  
  #################################################
  # R(t) TABLE
  Rt.get.dates <- Rt_add %>% filter(date %in% as.Date(date.in)) %>% mutate_if(is.numeric, round, digits=2) 
  Rt.eff.get.dates <- TEST %>% filter(date %in% as.Date(date.in)) %>% mutate_if(is.numeric, round, digits=2)
  
  Rt.format <- var.format(Rt.get.dates, 1)
  rownames(Rt.format) <- date.in
  colnames(Rt.format) <- "R(t)"
  
  Rt.eff.format <- var.format(Rt.eff.get.dates, 1)
  rownames(Rt.eff.format) <- date.in
  colnames(Rt.eff.format) <- "R(t) Effective"
  Rt.table <- cbind(Rt.format,Rt.eff.format)
  colnames(Rt.table) <- c("R(t) mean (95% CI)", "R(t) Effective mean (95% CI)")
  
  #################################################
  # Alpha Kappa Delta TABLE
  
  parameter_expand <- function(PAR_plot){
    cols.to.get <- c(1:6)
    n.cols <- length(cols.to.get)
    PAR_plot$step = as.numeric(PAR_plot$date - as.Date("2020-03-01"))
    PAR_add <- vector("list",length=n.cols)
    for (col.idx in cols.to.get){
      interpolate.out <- interpolate.fn(PAR_plot, col.idx)
      PAR_add[[col.idx]] <- interpolate.out$y
    }
    PAR_add <- do.call(cbind, PAR_add)
    x.vals <- PAR_plot$step
    PAR_add <- cbind(min(x.vals):max(x.vals),PAR_add) %>% as.data.frame()
    colnames(PAR_add) <- c("step",colnames(PAR_plot[1:6]))
    PAR_add$date <- PAR_add$step + as.Date("2020-03-01")
    return(PAR_add)
  }
  
  ## Expand to get value for each date
  Alpha_expand <- parameter_expand(Alpha_plot)
  Kappa_expand <- parameter_expand(Kappa_plot)
  Delta_expand <- parameter_expand(Delta_plot)
  
  ## Select dates
  Alpha.get.dates <- Alpha_expand %>% filter(date %in% as.Date(date.in)) %>% mutate_if(is.numeric, round, digits=3) 
  Kappa.get.dates <- Kappa_expand %>% filter(date %in% as.Date(date.in)) %>% mutate_if(is.numeric, round, digits=3) 
  Delta.get.dates <- Delta_expand %>% filter(date %in% as.Date(date.in)) %>% mutate_if(is.numeric, round, digits=3) 
  
  ## Format
  Alpha.format <- var.format(Alpha.get.dates,1)
  Kappa.format <- var.format(Kappa.get.dates,1)
  Delta.format <- var.format(Delta.get.dates,1)
  AKD.table <- cbind(Alpha.format,Kappa.format,Delta.format)
  rownames(AKD.table) <- date.in
  colnames(AKD.table) <- c("Alpha(t) mean (95% CI)", "Kappa(t) mean (95% CI)", "Delta(t) mean (95% CI)")
  
  
  #################################################
  # r(t) TABLE
  r_expand <- parameter_expand(r_plot)
  r.get.dates <- r_expand %>% filter(date %in% as.Date(date.in)) %>% mutate_if(is.numeric, round, digits=2) 
  r.table <- var.format(r.get.dates,1)
  rownames(r.table) <- date.in
  colnames(r.table) <- "r(t) mean (95% CI)"
  
  
  #################################################
  # CFR IFR TABLE
  
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
  
  CFR.posterior.CI <- get.CFR.IFR.by.date(traj.CI=traj.CI, CFR.or.IFR = "CFR", date.in=date.in, round.by=4)
  IFR.posterior.CI <- get.CFR.IFR.by.date(traj.CI=traj.CI, CFR.or.IFR = "IFR", date.in=date.in, round.by=4)
  
  CFR.format <- var.format(CFR.posterior.CI, use.mean=1)
  rownames(CFR.format) <- date.in
  colnames(CFR.format) <- "CFR"
  IFR.format <- var.format(IFR.posterior.CI, use.mean=1)
  rownames(IFR.format) <- date.in
  colnames(IFR.format) <- "IFR"
  CFR.IFR.table <- cbind(CFR.format,IFR.format)
  colnames(CFR.IFR.table) <- c("CFR mean (95% CI)", "IFR mean (95% CI)")
  
  
  #################################################
  # BRING TOGETHER TO CREATE TABLE
  
  PAR.table <- cbind(Rt.table, r.table, AKD.table)#, CFR.IFR.table)
  
  table.out <- vector(mode="list", length=4)
  table.out[[1]] <- PAR.table
  table.out[[2]] <- r.table
  table.out[[3]] <- AKD.table
  table.out[[4]] <- CFR.IFR.table
  
  return(table.out)
}

