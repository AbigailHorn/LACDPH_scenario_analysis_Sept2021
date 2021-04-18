# COVID19_JAM function

Covid_JAM = function(marginE, p, CorrelationMatrix, N, ridgeTerm = FALSE){
  
  ## Identify the number of categories
  z = as.vector(rep(NA, length(p)))
  
  n0 = N*p
  n1 = N*(1-p)
    
  y0 = -(n1*marginE)/(n0+n1)
  y1 = y0+marginE
  z = n1*y1
  
  ## Compute Covariance matrix based on the correlation matrix
  Dm = 2*p*(1-p)*N
  D_sqrt = diag(sqrt(Dm))
  xTx.scaled = D_sqrt %*% CorrelationMatrix %*% D_sqrt
  
  ## Add a ridge term in case X'X.scaled is singular
  ridgeValue = ifelse(ridgeTerm, min(1,min(diag(xTx.scaled)*.001)), 0)
  xTx.ridge = xTx.scaled + ridgeValue*diag(length(marginE))
  
  ## Get L matrix and zL vector
  L = chol(xTx.ridge)
  zL = solve(t(L))%*%z
  
  condE = as.vector(solve(t(L)%*%L)%*%t(L)%*%zL)
  
  # Get the conditional alpha vectors
  L_sub = L[, c(1,2)]
  tmp.tscore = summary(stats::lm(zL ~ 0 + L_sub))$coef[, 3]
  for(i_L in 3:ncol(L)){
    L_sub = L[, c(1, i_L)]
    tmp.tscore = rbind(tmp.tscore, summary(stats::lm(zL ~ 0 + L_sub))$coef[, 3])
  }
  tscore = rep(NA, ncol(L))
  tscore[1] = mean(tmp.tscore[,1])
  tscore[2:ncol(L)] = tmp.tscore[, 2]
  condSE = abs(condE/tscore)

  return(list(Conditional.betas = condE, Conditional.se = condSE))
}