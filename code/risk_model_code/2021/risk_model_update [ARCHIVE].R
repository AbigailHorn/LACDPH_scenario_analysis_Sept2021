
###############################################################################################################
## Read in marginal risk probability P(H|I), P(Q|H), P(D|Q) for each risk profile q, calculated by the JAM model
## Data coming from: Guan et al. 2020, Petrilli et al. 2020
Pr.H <- read.csv(path(data.dir, "Pr.H.csv"), sep=",", header=TRUE,row.names = 1)
Pr.Q <- read.csv(path(data.dir, "Pr.Q.csv"), sep=",", header=TRUE,row.names = 1)
Pr.D <- read.csv(path(data.dir, "Pr.D.csv"), sep=",", header=TRUE,row.names = 1)

## Get profile matrix X.mat
X.mat <- as.matrix(dplyr::select(Pr.H, -1))
n_factors <- ncol(X.mat)
n_profiles <- nrow(X.mat)

#########################################################################################################
## Read in prevalence of each risk factor p in the population for LA County
## Data coming from:
#### [Los Angeles County Health Survey](http://www.publichealth.lacounty.gov/ha/hasurveyintro.htm)
#### [California Health Information Survey](http://ask.chis.ucla.edu)
Pop.prevalence <- read.csv(path(data.dir, "Pop.prevalence.LA.csv"), sep=",", header=TRUE,row.names = 1 )
Pop.prev.matrix<-as.matrix(Pop.prevalence)

#########################################################################################################
## Read in correlation structure, $\Sigma$, between the risk factors p
## Correlation structure is weighted according to the prevalence of each race/ethnicity in LA County
## Calculated using data coming from: The National Health and Nutrition Examination Survey (NHANES)
SIGMA <- read.csv(path(data.dir, "SIGMA.csv"), sep=",", header=TRUE,row.names = 1)


#########################################################################################################
## Read in risk model updates from Lai

## Risk table 2020 (to compare)
effects_H_2020 = read.table("/Users/abigailhorn/Dropbox/0 2019/COVID-19/2020/LACounty_modeling/data/corrmats/Ordinal.Obesity_2020_04_20/Ordinal.Obesity.Effects.Hospitalization.txt", sep=" ",header=TRUE)


###############################################################################################################
## Read in calculated conditional log relative risk (psi) for each risk factor (H|I), (V|H), (D|V)
# Also includes marginal log relative risks
# Data coming from: Ioannou et al. 2020
effects_H = read.table(path(data.update.dir, "Ordinal.Obesity.Effects.Hospitalization.txt"), sep="\t", header=TRUE )
effects_V = read.table(path(data.update.dir, "Ordinal.Obesity.Effects.Ventilation.txt"), sep="\t", header=TRUE )
effects_D = read.table(path(data.update.dir, "Ordinal.Obesity.Effects.Death.txt"), sep="\t", header=TRUE )
#inputs_JAM = read.table(path(data.update.dir, "Ordinal.Comb.marginE.input.txt"), sep=" ", header=TRUE )

## Get profile matrix X.mat
X.profiles = read.table(path(data.update.dir, "Ordinal.Obesity.DummyRiskProfile.Hospitalization.txt"), sep=" ", header=TRUE )
X.mat = X.profiles %>% select(-c(Pr, RaceBlack, RaceAsian, RaceHispanic, SexMale))
X.mat = unique(X.mat)

## Define other necessary variables
n_factors <- ncol(X.mat)
n_profiles <- nrow(X.mat)

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
view(psi.mat.2)

## correlation matrix
SIGMA_update = read.table(path(data.update.dir, "Dummy.Age.BMI.Corr.Race.txt"), sep=" ", header=TRUE )
SIGMA_update = SIGMA_update[-c(10:14),-c(10:14)]
rownames(SIGMA_update) = rownames(psi.mat.2)
colnames(SIGMA_update) = rownames(psi.mat.2)
#write.csv(SIGMA_update, file = path(data.update.dir, "Dummy.Age.BMI.Corr.Race.csv"))

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
pop.prev.2021 <- read.csv(path(data.update.dir, "Pop.prevalence.2021.csv"), sep=",", header=TRUE,row.names = 1) %>% as.data.frame()





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
view(RR_table)



