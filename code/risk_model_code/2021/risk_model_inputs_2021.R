
###############################################################################################################
## Read in calculated conditional log relative risk (psi) for each risk factor (H|I), (V|H), (D|V)
## Data coming from: Guan et al. 2020, Petrilli et al. 2020
psi.mat <- read.csv(path(data.dir, "risk_model_2021/psi.mat.csv"), sep=",", header=TRUE,row.names = 1)

## Get profile matrix X.mat
# X.mat <- read.csv(path(data.dir, "risk_model_2021/X.mat.csv"), sep=",", header=TRUE,row.names = 1)
# X.mat$obese_3BMI.35 <- NULL
# colnames(X.mat)[6:8] <- rownames(psi.mat)[6:8]
# X.mat <- unique(X.mat)
# rownames(X.mat) <- c(1:nrow(X.mat))
# X.mat.clean <- X.mat
# idx=1
# for (i in 1:nrow(X.mat)){
#   if (sum(X.mat[i,c(6:8)])!=0){
#     X.mat.clean[idx,] <- X.mat.clean[i,]
#     idx=idx+1
#   }
# }
# X.mat <- X.mat.clean[1:(idx-1),]
#write.csv(X.mat, file = path(data.dir, "risk_model_2021/X.mat.csv"))

# ## Define other necessary variables
# n_factors <- ncol(X.mat)
# n_profiles <- nrow(X.mat)

## Readin correlation matrix
SIGMA <- read.csv(path(data.dir, "risk_model_2021/SIGMA.csv"), sep=",", header=TRUE,row.names = 1)

## calculating any.comorbidity for LA County
# This requires first calculating the profile.cnt.freq from the JAM.2020 inputs
# profile.cnt.freq.anyComb <- profile.cnt.freq %>% filter(ComorbidityYes==1)
# sum(profile.cnt.freq.anyComb$`LA County`)
# [1] 0.471981

## Read in prevalence of each risk factor p in the population for LA County
## Data coming from:
# [Los Angeles County Health Survey](http://www.publichealth.lacounty.gov/ha/hasurveyintro.htm)
# [California Health Information Survey](http://ask.chis.ucla.edu)
# American Community Survey (via "tidycensus" R package)
pop.prev.2021 <- read.csv(path(data.dir, "risk_model_2021/Pop.prevalence.2021.csv"), sep=",", header=TRUE,row.names = 1) %>% as.data.frame()



