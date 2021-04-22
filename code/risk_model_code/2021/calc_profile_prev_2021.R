
###############################################################################################################
## CALCULATE FREQUENCY OF EACH RISK PROFILE IN LAC POPULATION
###############################################################################################################

#########################################################################################################
### INPUTS REQUIRED
n_profiles = n_profiles
X.mat = X.mat
Pop.prev.matrix = as.matrix(pop.prev.2021)
SIGMA = SIGMA

#########################################################################################################
### Calculate frequency of each profile in LAC population

### Create sample population using marginal prevalence of each risk factor and weighted corelation matrix between risk profiles
# n_samples <- 100000
# profile.cnt.SPAs <- matrix(0, nrow = n_profiles, ncol = 1)
# sample.pop.df <- as.data.frame(rmvbin(n_samples, margprob=Pop.prev.matrix, sigma=SIGMA))
# colnames(sample.pop.df)<-colnames(Pop.prev.matrix)

## Create popualtion prevalence vector that inflates population prevalence of all other age groups besides 19-49,
## which will result in appropriate marginal distributions in the sampled population
Pop.prev.matrix.inflate = Pop.prev.matrix
Pop.prev.matrix.inflate[,"age_80."] = .1
Pop.prev.matrix.inflate[,"age_65.79"] = .12
Pop.prev.matrix.inflate[,"age_50.64"] = .18
Pop.prev.matrix.inflate[,"age_19.49"] = Pop.prev.matrix.inflate[,"age_19.49"] -.0675
Pop.prev.matrix.inflate[,"obese_0BMI.30"] = .61
Pop.prev.matrix.inflate[,"obese_1BMI30.40"] = .31
Pop.prev.matrix.inflate[,"obese_2BMI40."] = .08

n_samples = 100000
sample.pop.df <- matrix(0, nrow=n_samples, ncol=ncol(X.mat))

## Select only samples with 1 age and 1 obesity category
idx=1
while (idx < n_samples+1){
  current_row = (rmvbin(1, margprob=Pop.prev.matrix.inflate, sigma=SIGMA))
  if (sum(current_row[1:5])==1  & sum(current_row[6:8])==1) {
    sample.pop.df[idx, ] = current_row
    idx = idx+1
  }
}
sample.pop.df <- as.data.frame(sample.pop.df)
colnames(sample.pop.df)<-colnames(Pop.prev.matrix)

## Ensure comorbidity=0 if age_0.18==1
for (i in 1:n_samples){
  if (sample.pop.df[i,"age_0.18"]==1) { sample.pop.df[i,"comorbidity.yes"]=0 }
}

## Checking marginal distributions of sample population = actual population
sum(sample.pop.df[,c(1:5)])
sum(sample.pop.df[,c(6:8)])
sample.pop.df %>% filter(age_0.18==1) %>% nrow()/1e5
sample.pop.df %>% filter(age_19.49==1) %>% nrow()/1e5
sample.pop.df %>% filter(age_50.64==1) %>% nrow()/1e5
sample.pop.df %>% filter(age_65.79==1) %>% nrow()/1e5
sample.pop.df %>% filter(age_80.==1) %>% nrow()/1e5
sample.pop.df %>% filter(obese_0BMI.30==1) %>% nrow()/1e5
sample.pop.df %>% filter(obese_1BMI30.40==1) %>% nrow()/1e5
sample.pop.df %>% filter(obese_2BMI40.==1) %>% nrow()/1e5

### Count population prevalence of profiles (much faster with as.matrix())
profiles <- X.mat
sample.pop.table <- as.matrix(sample.pop.df)

profile.cnt <- c(rep(0,dim(profiles)[1]))
for (i in 1:dim(profiles)[1]) {
  print(i)
  for (j in 1:dim(sample.pop.table)[1]) {
    if(all(profiles[i,]==sample.pop.table[j,]))
    {profile.cnt[i]<-(profile.cnt[i]+1)}
  }
}

### Get frequency distribution over profiles
profile.cnt.freq <- round(profile.cnt/sum(profile.cnt), digits=6)
profile.cnt.freq<-as.data.frame(cbind(profile.cnt.freq,X.mat))
colnames(profile.cnt.freq) <- c("LA County",colnames(X.mat))
write.csv(profile.cnt.freq, file = path(data.dir, "risk_model_2021/profile.cnt.freq.csv"))
write.csv(profile.cnt.freq, file = path(data.dir, "risk_model_2021/profile.cnt.freq.SAVE.csv"))

TEST <- profile.cnt.freq %>% filter(age_80.==1)
sum(TEST$`LA County`)

# formattable(profile.cnt.freq)

#########################################################################################################
### OUTPUTS
# profile.cnt.freq



