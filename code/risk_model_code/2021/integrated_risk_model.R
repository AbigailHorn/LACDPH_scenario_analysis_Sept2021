


########################################################################################################################################
########################################################################################################################################

##################################################
##### INPUTS: RISK MODEL

### Read in JAM produced conditional effect estimates
# psi.mat = read.csv(path(data.dir, "psi.conditional.effect.estimates.csv"), sep=",", header=TRUE,row.names = 1)  #https://docs.google.com/spreadsheets/d/17QkcmFUCnxxaJPIfMkqWfIBhU2MaxaHGA03EtUBUM7o/edit#gid=0
# rownames(psi.mat) <- colnames(X.mat)

### Profile matrix
#X.mat <- as.matrix(dplyr::select(Pr.H.filter, -1))

### Estimated frequency of each risk profile
#freq.PREV.q <- as.vector(profile.cnt.SPAs[,9])



##################################################
##### INPUTS: OBSERVED LAC DATA ON ILLNESSES
#freq.LAC.obs.age

##################################################
##### INPUTS: EPIDEMIC MODEL

### Read in SEIR estimated Alpha, Kappa, Delta
#risk.probs.POSTERIORS



########################################################################################################################################
########################################################################################################################################
## Function to get estimated frequency of risk profiles in ILLNESS population re-scaled using observed LAC data

freq.ILL <- function(X.mat, freq.PREV.q, freq.LAC.obs.age, time){

  # Number of risk factors
  n.p <- ncol(X.mat)

  X.mat.AGE <- X.mat[,c(1:5)] %>% as.matrix()

  # MARGINALIZE: Get marginal frequency of each risk factor (BY AGE GROUPS ONLY)
  freq.PREV.p <- t(freq.PREV.q) %*% X.mat.AGE

  # Matrix multiply X.mat with marginal vector of risk factor prevalence frequencies FOR AGE ONLY to get the marginal AGE GROUP frequency in PREVALENCE corresponding to each profile
  marginal.PREV.freq.for.q <- X.mat.AGE %*% t(freq.PREV.p)

  # Get LAC OBSERVED marginal frequency of each AGE GROUP in ILLNESSES (adjusting so it is in 4 age groups vs. the 3 observed)
  t <- time
  g.age.p <- freq.LAC.obs.age[t,c(2:6)]

  # Matrix multiply X.mat with LAC OBSERVED marginal frequency of each AGE GROUP in ILLNESSES to get the marginal AGE GROUP frequency in ILLNESSES corresponding to each profile
  marginal.ILL.freq.for.q <- X.mat.AGE %*% t(g.age.p)

  # Dot product individual vector components, all q x 1, to get frequency in ILLNESSES over each risk profile
  freq.ILL.q <- freq.PREV.q * marginal.ILL.freq.for.q / marginal.PREV.freq.for.q

  # MARGINALIZE to get frequency in ILLNESSES over each risk factor
  freq.ILL.p <- t(freq.ILL.q) %*% as.matrix(X.mat)

  freq.ILL <- vector(mode="list", length=2)
  freq.ILL[[1]]<-freq.ILL.q
  freq.ILL[[2]]<-freq.ILL.p

  return(freq.ILL)

}

####################################################################################
## FUNCTION TO GET POPULATION FREQUENCY FUNCTION FOR H, Q, AND D POPULATIONS
## The general idea is to get:
## freq.H.q <- P.H.I.q * freq.ILL.q / (P.H.I.q %*% t(freq.ILL.q))

get.population.freq <- function(prob.q.vec.IN, freq.q.IN, X.mat) {
  X.mat <- as.matrix(X.mat)
  weighted.avg <- as.numeric( t(prob.q.vec.IN) %*% freq.q.IN )
  freq.q.OUT <- prob.q.vec.IN * freq.q.IN / weighted.avg  # FREQUENCY IN POPULATION
  freq.p.OUT <- t(freq.q.OUT) %*% X.mat
  freq.OUT <- vector(mode="list", length=2) #as.list[1,2]
  freq.OUT[[1]] <- freq.q.OUT
  freq.OUT[[2]] <- freq.p.OUT
  freq.OUT[[3]] <- weighted.avg

  return(freq.OUT)
}

####################################################################################
## FUNCTION TO ESTIMATE RISK PROBABILITIES WITH BASELINE FROM EPIDEMIC MODEL

get.Pr.q <- function(model, time, logit.SEIR.est, X.mat, freq.q.IN,  psi.mat){

  X.mat <- as.matrix(X.mat)
  
  # INITIALIZE / EXTRACT
  time <- time
  model <- model
  logit.SEIR.est.m.t <- as.numeric(logit.SEIR.est[time, model])
  logit.SEIR.est.m.t.VEC <- rep( logit.SEIR.est.m.t,  times=nrow(X.mat))
  psi.m.vec <- psi.mat[,model]
  #freq.IN.q <- freq.IN[[2]]
  n.profiles <- nrow(X.mat)
  freq.q.IN.mat <- matrix(rep(freq.q.IN,each=n.profiles),nrow=n.profiles)

  # GET PROBABILITY VECTOR
  Pr.q <- expit( logit.SEIR.est.m.t.VEC + (X.mat - freq.q.IN.mat) %*% psi.m.vec )

  return(Pr.q)

}

####################################################################################
## FUNCTION TO LOOP OVER ALL TIME VALUES

Pr.freq.FUN <- function(X.mat, freq.PREV.q, freq.LAC.obs.age, risk.probs.POSTERIORS, psi.mat){
  
  ## Transform epidemic model probability means for AKD to logit
  logit.SEIR.est <- risk.probs.POSTERIORS
  logit.SEIR.est[] <- lapply(risk.probs.POSTERIORS, function(x) round( log(x/(1-x)), 4 ) )
  
  ## Iterate over dates
  n.dates <- nrow(freq.LAC.obs.age)
  
  Pr.OUT <- vector("list", n.dates)
  freq.OUT <- vector("list", n.dates)
  #Pr.OUT <- matrix(data=NA, nrow=39, ncol=4)
  
  for (t in 1:n.dates){
    time = t
    name.date <- freq.LAC.obs.age$date[t]
    
    freq.I <- freq.ILL(X.mat, freq.PREV.q, freq.LAC.obs.age, time = time)
    
    Pr.H.I.q <- get.Pr.q(model=1, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.I[[2]],  psi.mat)
    freq.H <- get.population.freq(prob.q.vec.IN = Pr.H.I.q, freq.q.IN = freq.I[[1]], X.mat)
    
    Pr.Q.H.q <- get.Pr.q(model=2, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.H[[2]],  psi.mat)
    freq.Q <- get.population.freq(prob.q.vec.IN = Pr.Q.H.q, freq.q.IN = freq.H[[1]], X.mat)
    
    Pr.D.Q.q <- get.Pr.q(model=3, time=time, logit.SEIR.est, X.mat, freq.q.IN = freq.Q[[2]],  psi.mat)
    freq.D <- get.population.freq(prob.q.vec.IN = Pr.D.Q.q, freq.q.IN = freq.Q[[1]], X.mat)
    
    Pr.D.I.q <- Pr.H.I.q*Pr.Q.H.q*Pr.D.Q.q
    
    Pr.OUT.t <- round( cbind(Pr.H.I.q, Pr.Q.H.q, Pr.D.Q.q, Pr.D.I.q), 7)
    colnames(Pr.OUT.t) <- c("P(H|I)_", "P(Q|H)_", "P(D|Q)_", "P(D|I)_")
    colnames(Pr.OUT.t) <- paste0( colnames(Pr.OUT.t), name.date)
    Pr.OUT[[t]] <- Pr.OUT.t
    
    freq.OUT.t <- round ( cbind( freq.I[[1]], freq.H[[1]], freq.Q[[1]], freq.D[[1]] ), 7)
    colnames(freq.OUT.t) <- c("freq.I.", "freq.H.", "freq.Q.", "freq.D.")
    colnames(freq.OUT.t) <- paste0( colnames(freq.OUT.t), name.date)
    freq.OUT[[t]] <- freq.OUT.t
    
  }
  Pr.OUT <- do.call(cbind, Pr.OUT)
  freq.OUT <- do.call(cbind, freq.OUT)
  
  FUN.out <- vector("list", 2)
  FUN.out[[1]]<-Pr.OUT
  FUN.out[[2]]<-freq.OUT
  return(FUN.out)
  
}


##################################################################################################################
##################################################################################################################








