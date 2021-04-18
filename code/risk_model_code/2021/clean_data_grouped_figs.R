
################################################################
## Get data in form for plotting grouped figures
################################################################

# profile.prev.data <- risk_table_CFR_FULL %>% select(c(Profile,riskprofile,age, BMI, smoking, comorbidity, freq.PREV.q))
# if (n.times > 1) {
#   for (idx in 2:n.times) {    
#     name <- times[idx]
#     profile.prev.data[,name] <- profile.prev.data$t1 
#   }}

profile.prev.data <- Pr.HI %>% select(c(Profile,Age,BMI,Smoker,Comorbidity,Prevalence))
#colnames(profile.prev.data)[colnames(profile.prev.data)=="Prevalence"] = times.dates[1]
# if (n.times > 1) {
#   for (idx in 2:n.times) {    
#     name <- times.dates[idx]
#     profile.prev.data[,name] <- profile.prev.data$t1 
#   }}

data.melted <- melt(profile.prev.data, id = c("Profile","Age","BMI","Smoker","Comorbidity"))
data.melted$DATE <- data.melted$variable
data.melted$variable <- NULL
data.melted$DATE <- factor(data.melted$DATE, levels = times.dates)
data.melted$prevalence.gen.pop <- data.melted$value
data.melted$value <- NULL

data.melted.all <- vector("list", length=n.times)
for (i in 1:n.times){
  data.melted.current <- data.melted
  data.melted.current$DATE <- times.dates[i]
  data.melted.all[[i]] <- data.melted.current
}
data.melted.all <- do.call(rbind, data.melted.all)
data.melted <- data.melted.all

## Merge with everything else
full.data <- merge(data.melted, dataI.melted, by = c('Profile', 'DATE' ) )
full.data <- merge(full.data, dataH.melted, by = c('Profile', 'DATE' ) )
full.data <- merge(full.data, dataQ.melted, by = c('Profile', 'DATE' ) )
full.data <- merge(full.data, dataD.melted, by = c('Profile', 'DATE' ) )

# ## Save
# full.data.orig <- full.data

## Melt prevalence variable
full.data <- full.data %>% 
  gather(keys, values, prevalence.gen.pop:prevalence.D )
full.data$stage <- full.data$keys
full.data$keys <- NULL 
full.data$stage <- factor(full.data$stage, levels =  c("prevalence.gen.pop", "prevalence.I","prevalence.H", "prevalence.Q", "prevalence.D"))

## Name the time periods depicted
times.names <- times.dates
names(times.names) <- times.dates

################################################################
## Plot function for grouped figures
################################################################



grouped.var.figs <- function(full.data, var, var.name, times.dates){
  
  ## Name the time periods depicted
  times.names <- times.dates
  names(times.names) <- times.dates
  
  ## Enable using aes_string()
  stage.str <- "stage"
  values.str <- "values"
  var.str <- as.character(var)
  
  ## Fig
  grouped.fig <- ggplot(full.data, aes_string(x = stage.str, y = values.str, fill = var)) +
    geom_bar(position="stack", stat = 'identity') + 
    facet_wrap(DATE ~ ., labeller=labeller(DATE = times.names)) +
    labs(title = paste0(var.name, " by stage of disease"), x = NULL, y = "Frequency") +
    scale_x_discrete(labels = c("Prevalence", "Infected", "Hospitalized", "ICU", "Dead")) +
    theme(axis.text.x = element_text(angle = 90))
  grouped.fig
}


