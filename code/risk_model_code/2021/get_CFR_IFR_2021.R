
##############################################################################################################
##############################################################################################################
## GET CFR AND IFR BY PROFILE
##############################################################################################################
##############################################################################################################

##############################################################################################################
## Create dataframe summarizing risk profile info, prevalence of risk profile in general population,
##    and prevalence of risk profile in deceased population
##############################################################################################################

#data.pop.prev <- data.FULL
data.pop.prev <- dataI
colnames(data.pop.prev) <- c("Profile", paste0( "I.prev.", times.dates))
data.D.prev <- dataD
colnames(data.D.prev) <- c("Profile", paste0( "D.prev.", times.dates))
data.prev <- merge(data.pop.prev, data.D.prev, by = "Profile")
data.prev <- cbind(data.prev, freq.PREV.q)
data.prev <- arrange(data.prev, data.prev$Profile)

##############################################################################################################
## Function for risk table with CFR and IFR
##############################################################################################################

get.CFR.IFR.profiles <- function(ABC.out.mat=ABC.out.mat, time.steps=time.steps, iter=iter, data.prev=data.prev, times.dates=times.dates, round.by=round.by) {
  
  ## MODEL OUTPUT TO PLOT
  
  # ABC.out.mat.test = ABC.out.mat[1:100,]
  # model.out <- correlated.param.SIM(ABC.out.mat=ABC.out.mat.test,iter=iter,time.steps=time.steps)
  
  model.out <- correlated.param.SIM(ABC.out.mat=ABC.out.mat,iter=iter,time.steps=time.steps)
  
  # Align dates
  init.date = as.Date("2020-03-01")
  model.out[["date"]] <- model.out[["date"]] + init.date
  model.out.filter <- model.out %>% filter(date %in% times.dates) %>% select(par.id, date, iter, Idetectcum, Itotcum, D)
  
  ########################################
  ## Get CFR and IFR for each iteration
  ########################################
  
  ## Function to multiply each value of I / D by the prevalence of each profile
  prev.multiply <- function(model.out, model.var, data.prev, t) {
    prev.vec <- data.prev[,t+1]
    var.matrix <- matrix(model.out[,model.var], nrow=nrow(model.out), ncol=length(prev.vec), byrow = FALSE)  ## Replicate the variable (I or D) n.profiles times
    out.prev <- sweep(var.matrix, MARGIN=2, prev.vec, `*`)
    return(out.prev)
  }
  
  ## Get CFR and IFR for each date
  n.profiles <- nrow(data.prev)
  CFR.OUT <- vector("list", n.times)
  IFR.OUT <- vector("list", n.times)
  
  for (t in 1:n.times){
    
    date.t <- as.Date(times.dates[t])
    model.out.t <- model.out.filter %>% filter(date==date.t)
    I.prev.t <- prev.multiply(model.out.t, "Idetectcum",data.pop.prev,t=t)
    Itot.prev.t <- prev.multiply(model.out.t, "Itotcum",data.pop.prev,t=t)
    D.prev.t <- prev.multiply(model.out.t, "D",data.D.prev,t=t)
    CFR.t <- D.prev.t / I.prev.t
    IFR.t <- D.prev.t / Itot.prev.t
    
    # Get median, up_95, low_95
    get.CI.FR <- function(FR.t, is.CFR, round.by){
      colnames(FR.t) <- Profile
      FR.t[!is.finite(FR.t)] <- 0
      FR.t <- reshape2::melt(as.data.frame(FR.t), measure.vars = c(1:ncol(FR.t)), variable.name = "Profile")
      FR.dt <- as.data.table(FR.t)
      FR.CI <- FR.dt[, list(
        median = quantile(value, c(.5), na.rm=TRUE),
        low_95 = quantile(value, c(.025), na.rm=TRUE),
        up_95 = quantile(value, c(.975), na.rm=TRUE)),
        by = "Profile"]
      FR.CI <- as.data.frame(FR.CI) %>% mutate_if(is.numeric, round, digits=round.by)
      if (is.CFR==TRUE) FR <- "CFR_" else FR = "IFR_"
      colnames(FR.CI) <- c("Profile", paste0(FR, paste0(c("","low.95_","up.95_"), times.dates[t] )))
      rm(FR.dt)
      rm(FR.t)
      FR.CI$Profile <- NULL
      return(FR.CI)
    }
    CFR.CI.t <- get.CI.FR(CFR.t, is.CFR=TRUE, round.by=round.by)
    IFR.CI.t <- get.CI.FR(IFR.t, is.CFR=FALSE, round.by=round.by)
    
    # Save
    CFR.OUT[[t]] <- CFR.CI.t
    IFR.OUT[[t]] <- IFR.CI.t
  } # end over times i
  
  limit.1 <- function(DF){
    for (i in 1:nrow(DF)){
      for (j in 1:ncol(DF)){
        if (DF[i,j]>1) DF[i,j]=1
      }
    }
    return(DF)
  }
  
  CFR.OUT <- do.call(cbind, CFR.OUT)
  CFR.OUT <- limit.1(CFR.OUT)
  IFR.OUT <- do.call(cbind, IFR.OUT)
  IFR.OUT <- limit.1(IFR.OUT)
  
  # ### SANITY CHECKS: take weighted average of frequency of each profile and CFR, should equal population-average CFR
  # sum(CFR.OUT[,"CFR.2020-05-15"]*freq.OUT[,"freq.I.2020-05-15"])
  # sum(CFR.OUT[,"CFR.t2"]*freq.OUT[,"freq.I.t2"])
  # sum(CFR.OUT[,"CFR.t3"]*freq.OUT[,"freq.I.t3"])
  
  # ########################################
  # ## Create table output
  # ########################################
  # 
  # freq.OUT.I <- freq.OUT[,c(4*(seq.times-1)+1)]
  # table.CFR.IFR <- cbind(CFR.OUT %>% select(c(paste0("CFR.",times.dates))),IFR.OUT %>% select(c(paste0("IFR.",times.dates))))
  # table.CFR.IFR$Profile <- Profile
  # data.freq.CFR.IFR <- data.FULL %>% select("Profile", "Age","BMI","Smoker","Comorbidity","freq.PREV.q") %>% cbind(freq.OUT.I)
  # data.freq.CFR.IFR <- merge(data.freq.CFR.IFR, table.CFR.IFR, by = "Profile")
  # 
  # # FILTER PROFILES TO SHOW IN THE TABLE TO ONLY THOSE > 0 PREVALENCE
  # #data.Pr.CFR.IFR <- filter(data.Pr.CFR.IFR, paste0("freq.I.t1") > 0.000000004)
  # 
  # ########################################
  # ## Add risk group
  # ## Note: Groupings done by CFR on t1
  # ########################################
  # # 
  # # riskprofile <- "Blank"
  # # riskprofile[data.Pr.CFR.IFR$"CFR.t1"<.01] <- "Risk 5"
  # # riskprofile[(data.Pr.CFR.IFR$"CFR.t1"<.04) & (data.Pr.CFR.IFR$"CFR.t1">.01)] <- "Risk 4"
  # # riskprofile[(data.Pr.CFR.IFR$"CFR.t1"<.08) & (data.Pr.CFR.IFR$"CFR.t1">.04)] <- "Risk 3"
  # # riskprofile[(data.Pr.CFR.IFR$"CFR.t1"<.16) & (data.Pr.CFR.IFR$"CFR.t1">.08)] <- "Risk 2"
  # # riskprofile[(data.Pr.CFR.IFR$"CFR.t1">.16)] <- "Risk 1"
  # # riskprofile <- factor(riskprofile, levels = c("Risk 1", "Risk 2", "Risk 3", "Risk 4", "Risk 5" ))
  # # 
  # # data.Pr.CFR.IFR <- tibble::add_column(data.Pr.CFR.IFR, riskprofile, .after = "Profile")
  # 
  # ##### Arrange
  # #data.Pr.CFR.IFR <- arrange(data.Pr.CFR.IFR, desc(data.Pr.CFR.IFR$"CFR.t1"))
  
  data.CFR.IFR <- vector("list", 2)
  data.CFR.IFR[[1]]<-CFR.OUT
  data.CFR.IFR[[2]]<-IFR.OUT
  return(data.CFR.IFR)
  
}

###############################################################################################################################
## CFR and IFR table to print

get.CFR.CI.table <- function(n.times, CFR.OUT, data.FULL){
  
  ## Function to put variables in median (low_95, up_95) format
  var.format <- function(var.CI){
    posterior.CI.FORMAT <- function(var.CI){
      var.CI <- as.data.frame(var.CI)
      low_95 <- var.CI$low_95
      up_95 <- var.CI$up_95
      median <- var.CI$median
      if (median==0) var.95.CI = "NA"
      else var.95.CI <- paste0(median, " (", low_95, ",", up_95,")")
      return(var.95.CI)
    }
    n.idx <- nrow(var.CI)
    var.format = as.data.frame(matrix(nrow=n.idx, ncol=1))
    for (i in 1:n.idx){
      var.format[i,] <- posterior.CI.FORMAT(var.CI[i,])
    }
    return(var.format)
  }
  
  ## Get the table
  CFR.all <- vector("list",n.times)
  for (t in 1:n.times){
    seq = c((3*t-2):(3*t))
    CFR.date <- CFR.OUT[,seq]
    colnames(CFR.date) <- c("median", "low_95","up_95")
    CFR.date.CI <- var.format(CFR.date)
    colnames(CFR.date.CI) <- sapply(colnames(CFR.OUT)[t*3-2], function(x) {str_split(x, "_")[[1]][2]})
    CFR.all[[t]] <- CFR.date.CI
  }
  CFR.all <- do.call(cbind, CFR.all)
  
  
  ## Add descriptive variables
  CFR.all.OUT <- cbind(data.FULL[,1],data.FULL[,c("Age","BMI","Smoker","Comorbidity")],data.FULL[,2],CFR.all)
  colnames(CFR.all.OUT)[1]="Profile"
  colnames(CFR.all.OUT)[6]="Prevalence"
  CFR.all.OUT[CFR.all.OUT[,]==0]="NA"

  return(CFR.all.OUT)
}



