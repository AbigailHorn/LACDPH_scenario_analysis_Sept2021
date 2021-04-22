


############################################################################################################
############################################################################################################
## CLEAN AND PROCESS DATA FOR RISK ESTIMATES
############################################################################################################
############################################################################################################

Pr.OUT.print <- apply(as.data.frame(Pr.OUT), MARGIN=2, FUN=round, 4)
freq.OUT.print <- apply(as.data.frame(freq.OUT), MARGIN=2, FUN=round, 4)
freq.PREV.q.print <- round(freq.PREV.q, 6)
i=1:(n.times-1)

# profile.idx <- as.data.frame(c(1:n_profiles))
# colnames(profile.idx) <- "Profile ID"

n.profiles <- nrow(X.mat)
n.times <- length(times.dates) 
seq.times <- 1:n.times

## 1) RENAME VARIABLES #############################################################################

# Change the profile numbers to factors and name the variable
Profile <- factor(seq(1,n.profiles))

## dataI ##
seqI <- 4*(seq.times-1)+1
dataI <- as.data.frame(freq.OUT[,seqI])
colnames(dataI) <- times.dates
dataI <- cbind(Profile, dataI)

dataI.melted <- melt(dataI, id = 'Profile')
dataI.melted$DATE <- dataI.melted$variable
dataI.melted$variable <- NULL
#dataI.melted$DATE <- factor(dataI.melted$DATE, levels =  times.dates)
dataI.melted$prevalence.I <- dataI.melted$value
dataI.melted$value <- NULL

## dataH ##
seqH <- 4*(seq.times-1)+2
dataH <- as.data.frame(freq.OUT[,seqH])
colnames(dataH) <- times.dates
dataH <- cbind(Profile, dataH)

dataH.melted <- melt(dataH, id = 'Profile')
dataH.melted$DATE <- dataH.melted$variable
dataH.melted$variable <- NULL
#dataH.melted$DATE <- factor(dataH.melted$DATE, levels =  times.dates)
dataH.melted$prevalence.H <- dataH.melted$value
dataH.melted$value <- NULL

## dataQ ##
seqQ <- 4*(seq.times-1)+3
dataQ <- as.data.frame(freq.OUT[,seqQ])
colnames(dataQ) <- times.dates
dataQ <- cbind(Profile, dataQ)

dataQ.melted <- melt(dataQ, id = 'Profile')
dataQ.melted$DATE <- dataQ.melted$variable
dataQ.melted$variable <- NULL
#dataQ.melted$DATE <- factor(dataQ.melted$DATE, levels =  times.dates)
dataQ.melted$prevalence.Q <- dataQ.melted$value
dataQ.melted$value <- NULL

## dataD ##
seqD <- 4*(seq.times)
dataD <- as.data.frame(freq.OUT[,seqD])
colnames(dataD) <- times.dates
dataD <- cbind(Profile, dataD)

dataD.melted <- melt(dataD, id = 'Profile')
dataD.melted$DATE <- dataD.melted$variable
dataD.melted$variable <- NULL
#dataD.melted$DATE <- factor(dataD.melted$DATE, levels =  times.dates)
dataD.melted$prevalence.D <- dataD.melted$value
dataD.melted$value <- NULL

##################################################################################################
## 2) DUMMY VARIABLES <- LEVELS ##################################################################

#data.FULL <- as.data.frame(cbind(Profile, X.mat, freq.PREV.q, Pr.OUT, freq.OUT))
data.FULL <- as.data.frame(cbind(Profile, X.mat, freq.PREV.q, Pr.OUT.print))

data.FULL$Age <- "Blank"
data.FULL$Age[data.FULL$age_0.18==1] <- "0-18"
data.FULL$Age[data.FULL$age_19.49==1] <- "19-49"
data.FULL$Age[data.FULL$age_50.64==1] <- "50-64"
data.FULL$Age[data.FULL$age_65.79==1] <- "65-79"
data.FULL$Age[data.FULL$age_80.==1] <- "80+"
data.FULL$age_0.18 <- NULL
data.FULL$age_19.49 <- NULL
data.FULL$age_50.64 <- NULL
data.FULL$age_65.79 <- NULL
data.FULL$age_80. <- NULL
data.FULL$Age <- factor(data.FULL$Age, levels =  c("0-18", "19-49", "50-64", "65-79", "80+"))

data.FULL$BMI <- "Blank"
data.FULL$BMI[data.FULL$obese_0BMI.30==1] <- "BMI<30"
data.FULL$BMI[data.FULL$obese_1BMI30.40==1] <- "30<BMI<40"
data.FULL$BMI[data.FULL$obese_2BMI40.==1] <- "BMI>40"
data.FULL$obese_0BMI.30 <- NULL
data.FULL$obese_1BMI30.40 <- NULL
data.FULL$obese_2BMI40. <- NULL
data.FULL$BMI <- factor(data.FULL$BMI, levels =  c("BMI<30", "30<BMI<40", "BMI>40"))

data.FULL$smoke2 <- "Blank"
data.FULL$smoke2[data.FULL$smoker==1] <- "Smoker"
data.FULL$smoke2[data.FULL$smoker==0] <- "Non Smoker"
data.FULL$smoke2 <- factor(data.FULL$smoke2, levels = c("Smoker", "Non Smoker"))
data.FULL$smoker = NULL
data.FULL$Smoker = data.FULL$smoke2
data.FULL$smoke2 = NULL

data.FULL$Comorbidity <- "Blank"
data.FULL$Comorbidity[data.FULL$comorbidity.yes==1] <- "Comorbidity"
data.FULL$Comorbidity[data.FULL$comorbidity.yes==0] <- "No Comorbidity"
data.FULL$Comorbidity <- factor(data.FULL$Comorbidity, levels = c("Comorbidity", "No Comorbidity"))
data.FULL$comorbidity.yes <- NULL


###############################################################################################################################
## Probabilities table to print

extract.Pr.table <- function(n.times,order.pr){
  i=1:(n.times-1)
  seq.this <- c(1,1+4*i) + 1 + order.pr
  Pr.this <- data.FULL[,seq.this] %>% as.data.frame()
  colnames(Pr.this) <- sapply(colnames(Pr.this)[], function(x) {str_split(x, "_")[[1]][2]})
  Pr.this <- cbind(data.FULL[,1],data.FULL[,c("Age","BMI","Smoker","Comorbidity")],data.FULL[,2],Pr.this)
  colnames(Pr.this)[1]="Profile"
  colnames(Pr.this)[6]="Prevalence"
  return(Pr.this)
}

# Pr.HI <- extract.Pr.table(order.pr = 1)
# Pr.QH <- extract.Pr.table(order.pr = 2)
# Pr.DQ <- extract.Pr.table(order.pr = 3)
# Pr.DI <- extract.Pr.table(order.pr = 4)
# 
# kable(Pr.HI) %>%
#   kable_styling(bootstrap_options = "striped", full_width = F)
# 
# if (print.output==TRUE) {
#   write.csv(Pr.HI, file = path(output.dir/Pr.tables, "Pr.HI.csv"))
#   write.csv(Pr.QH, file = path(output.dir/Pr.tables, "Pr.QH.csv"))
#   write.csv(Pr.DQ, file = path(output.dir/Pr.tables, "Pr.DQ.csv"))
#   write.csv(Pr.DI, file = path(output.dir/Pr.tables, "Pr.DI.csv"))
# }
