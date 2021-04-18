
###############################################################################################################
## Read in calculated conditional log relative risk (psi) for each risk factor (H|I), (V|H), (D|V)
# Also includes marginal log relative risks
# Data coming from: Ioannou et al. 2020
effects_H = read.table(path(data.update.dir, "Ordinal.Obesity.Effects.Hospitalization.txt"), sep="\t", header=TRUE )
effects_V = read.table(path(data.update.dir, "Ordinal.Obesity.Effects.Ventilation.txt"), sep="\t", header=TRUE )
effects_D = read.table(path(data.update.dir, "Ordinal.Obesity.Effects.Death.txt"), sep="\t", header=TRUE )
#inputs_JAM = read.table(path(data.update.dir, "Ordinal.Comb.marginE.input.txt"), sep=" ", header=TRUE )

## Create psi.mat from effects tables
psi.mat.2 = as.data.frame(effects_H[,1])
psi.mat.2$H.I = effects_H$Conditional.betas
psi.mat.2$V.H = effects_V$Conditional.betas
psi.mat.2$D.V = effects_D$Conditional.betas
rownames(psi.mat.2) <- psi.mat.2[,1]
psi.mat.2$`effects_H[, 1]` <- NULL
age_0.18 <- rep(0, times=3)
age_19.50 <- rep(0, times=3)
obese_0BMI.24 <- rep(0, times=3)
psi.mat.2 <- rbind(age_0.18, age_19.50, psi.mat.2[c(1:3),], obese_0BMI.24, psi.mat.2[c(4:6,11:12),])
rownames(psi.mat.2)[1:2] <- c("age_0.18","age_19.49")
rownames(psi.mat.2)[6:8] <- c("obese_0BMI.24","obese_1BMI25.29","obese_2BMI30.34")
psi.mat = psi.mat.2
#write.csv(psi.mat, file = path(data.dir, "risk_model_2021/psi.mat.csv"))
#view(psi.mat.2)

## Get profile matrix X.mat
X.profiles = read.table(path(data.update.dir, "Ordinal.Obesity.DummyRiskProfile.Hospitalization.txt"), sep=" ", header=TRUE )
X.mat = X.profiles %>% select(-c(Pr, RaceBlack, RaceAsian, RaceHispanic, SexMale))
X.mat = unique(X.mat)
colnames(X.mat) = rownames(psi.mat.2)

## Add profiles for age 0-18
first.row <- X.mat[1,]
n <- 7
replicated <- do.call("rbind", replicate(n, first.row, simplify=FALSE))
X.mat <- rbind(replicated, X.mat)
obese.0 <- c(1, 1, rep(0,6))
obese.1 <- c(0, 0, 1, 1, rep(0,4))
obese.2 <- c(rep(0,4), 1, 1, 0, 0)
obese.3 <- c(rep(0,6), 1, 1)
smoking.vec <- c(1,0,1,0,1,0,1,0)
X.mat$obese_0BMI.24[1:8] <- obese.0
X.mat$obese_1BMI25.29[1:8] <- obese.1
X.mat$obese_2BMI30.34[1:8] <- obese.2
X.mat$obese_3BMI.35[1:8] <- obese.3
X.mat$smoker[1:8] <- smoking.vec
#write.csv(X.mat, file = path(data.dir, "risk_model_2021/X.mat.csv"))

## Define other necessary variables
n_factors <- ncol(X.mat)
n_profiles <- nrow(X.mat)

## Readin correlation matrix
SIGMA_update = read.table(path(data.update.dir, "Dummy.Age.BMI.Corr.Race.txt"), sep=" ", header=TRUE )
SIGMA_update = SIGMA_update[-c(10:14),-c(10:14)]
rownames(SIGMA_update) = rownames(psi.mat.2)
colnames(SIGMA_update) = rownames(psi.mat.2)
SIGMA = SIGMA_update
#write.csv(SIGMA, file = path(data.dir, "risk_model_2021/SIGMA.csv"))

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
#pop.prev.2021 <- read.csv(path(data.update.dir, "Pop.prevalence.2021.csv"), sep=",", header=TRUE,row.names = 1) %>% as.data.frame()



########################################################################################################
## Create RR_table (for paper; not used in analysis)
# Note: This is missing marginal RR CI
RR_table = as.data.frame(effects_H$Risk.Factors)
RR_table$V1 = effects_H$Conditional.CI
RR_table$V2 = effects_V$Conditional.CI
RR_table$V3 = effects_D$Conditional.CI
RR_table$`effects_H$Risk.Factors` = NULL
rownames(RR_table) = effects_H$Risk.Factors
#colnames(RR_table) = c("Conditional RR (H|I)", "Conditional RR (V|H)", "Conditional RR (D|V)")
age_0.18 <- as.data.frame(t(rep(1, times=3)))
age_19.50 <- as.data.frame(t(rep(1, times=3)))
RR_table <- rbind(age_0.18, age_19.50, RR_table)
rownames(RR_table)[1:2] = c("age_0.18", "age_19.50")
colnames(RR_table) = c("Conditional RR (H|I)", "Conditional RR (V|H)", "Conditional RR (D|V)")
#write.csv(RR_table, file = path(data.dir, "risk_model_2021/RR_table.csv"))
#view(RR_table)



