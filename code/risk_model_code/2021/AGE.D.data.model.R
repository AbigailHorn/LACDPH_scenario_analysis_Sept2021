
X.mat.AGE <- X.mat[,c(1:5)] %>% as.matrix()

################################################################################

## Frequency of age in deaths: data
freq.obs.age.D <- ages_pct %>% filter(date %in% closest.dates) %>% select(date,death_pct,age.strata) #%>% unique()
AGE = recode_factor(freq.obs.age.D$age.strata, age_0.18 = "0-18", age_19.49 = "19-49", age_50.64="50-64", age_65.79="65-79", age_80.="80+")
freq.obs.age.D$Age = AGE
freq.obs.age.D$age.strata = NULL
freq.obs.age.D$type <- "data"
freq.obs.age.D$date <- as.character(freq.obs.age.D$date)
freq.obs.age.D$date <- recode(freq.obs.age.D$date, "2020-05-15"="2020-05-15", "2020-08-01"="2020-08-01", "2020-10-14"="2020-10-15", "2020-11-16"="2020-11-15", "2020-12-15"="2020-12-15", "2021-01-13"="2021-01-15", "2021-02-17"="2021-02-15", "2021-03-03"= "2021-03-01")
freq.obs.age.D$date <- as.Date(freq.obs.age.D$date)

## Frequency of age in deaths: model
freq.model.age.D <- vector("list",n.times)
for (i in 2:ncol(dataD)){
  this.date <- as.Date(colnames(dataD)[i])
  freq.D.q <- dataD[,i]
  # MARGINALIZE: Get marginal frequency of each risk factor (BY AGE GROUPS ONLY)
  freq.D.p <- t(freq.D.q) %*% X.mat.AGE
  freq.D.p <- t(freq.D.p) %>% as.data.frame()
  freq.D.p$Age <- c("0-18","19-49", "50-64", "65-79", "80+")
  freq.D.p$date <- this.date
  freq.model.age.D[[i]]<-freq.D.p
}
freq.model.age.D <- do.call(rbind, freq.model.age.D)
freq.model.age.D$death_pct <- freq.model.age.D$V1
freq.model.age.D$V1 <- NULL
freq.model.age.D$type <- "model"
rownames(freq.model.age.D) <- NULL

## Combine
freq.model.age.D <- freq.model.age.D[,c("date","Age","death_pct","type")]
freq.obs.age.D <- freq.obs.age.D[,c("date","Age","death_pct","type")]
freq.age.D <- rbind(freq.model.age.D,freq.obs.age.D)



################################################################
## Plot function for age in deaths
################################################################

plot.freq.age.D <- function(freq.age.D){
  
  times.names <- times.dates
  names(times.names) <- times.dates
  
  p <- ggplot(freq.age.D, aes(x=type, y=death_pct, fill=Age)) + 
    geom_bar(position="stack", stat = 'identity') + 
    facet_wrap(date ~ ., labeller=labeller(date = times.names)) +
    labs(title = paste0("Deaths by age group"), x = NULL, y = "Frequency") +
    scale_x_discrete(labels = c("data","model")) +
    theme(axis.text.x = element_text(angle = 90))
  return(p)
}




