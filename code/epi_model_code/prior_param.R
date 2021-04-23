## Prior parameter distributions

prior.st <- c("unif",44,46)

R0 <- 3.65
prior.R0 <- c("normal",R0,.03)

R0_redux1 <- .21
R0_redux2 <- .33
R0_redux3 <- .55 #.575#.6 #65

prior.R0.redux1 <- c("unif", R0_redux1 - 0.05, R0_redux1 + 0.05)
prior.R0.redux2 <- c("unif", R0_redux2 - 0.05, R0_redux2 + 0.05)
prior.R0.redux3 <- c("unif", R0_redux3 - 0.05, R0_redux3 + 0.05)


# prior.r1 <- c("unif",0.1, 0.75) #c("unif", 0.13, 0.16)
# prior.r2 <- c("unif",0.1, 0.75) #c("unif",0.12, 0.14) #0.30

prior.r1 <- c("unif",0.03, 0.35)
prior.r2 <- c("unif", 0.3, 0.65) #c("unif",0.1, 0.85)


Alpha1 <- 0.375#.155
Kappa1 <- 0.35 #.65
Delta1 <- .6#.63

Alpha2 <- .15 #.09 
Kappa2 <- .2 #.185
Delta2 <- .55 #.575 #.8 

stdev <- 0.003 #0.003
prior.Delta1 <- c("normal",Delta1, stdev)#.001) #0.01)
prior.Alpha1 <- c("normal",Alpha1, stdev)#.002) #0.003)
prior.Kappa1 <- c("normal",Kappa1, stdev)#.002) # 0.03)

# prior.Delta2 <- c("unif",Alpha2-stdev,Alpha2+stdev) #c("normal",Delta2, stdev)#.001)
# prior.Alpha2 <- c("unif",Kappa2-stdev,Kappa2+stdev) # c("normal",Alpha2, stdev)#.001)
# prior.Kappa2 <- c("unif",Delta2-stdev,Delta2+stdev) #c("normal",Kappa2, stdev)#.002)
prior.Delta2 <- c("normal",Delta2, stdev)#.001)
prior.Alpha2 <- c("normal",Alpha2, stdev)#.001)
prior.Kappa2 <- c("normal",Kappa2, stdev)#.002)

p_V <- 0.27 #.28
prior.p_V <- c("normal",p_V, stdev) # c("unif", p_V-0.1, p_V+0.5 ) #c("normal",p_V, 0.08)

prior.par <- list(
  prior.R0,
  prior.r1,
  prior.st,
  prior.R0.redux1,
  prior.Delta1,
  prior.Alpha1,
  prior.Kappa1,
  prior.p_V,
  prior.R0.redux2,
  prior.Delta2, #10
  prior.Alpha2, #11
  prior.Kappa2, #12
  prior.r2,
  prior.R0.redux3)

## TIME-VARYING FUNCTIONS ARE READ IN VIA:
## fn_t_readin_path <- path(data.dir, "fn_t_readin.csv")
