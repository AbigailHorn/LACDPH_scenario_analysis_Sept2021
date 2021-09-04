
#########################################################################################################
#########################################################################################################
## Plot mean (min,max) of each probability over profiles within each age group

## Get AGE.idx as a column

AGE.idx = vector(length=n_profiles)
for(i in 1:n_profiles){
  if (X.mat$age_0.18[i]==1){AGE.idx[i]=1}
  else if (X.mat$age_19.49[i]==1){AGE.idx[i]=2}
  else if (X.mat$age_50.64[i]==1){AGE.idx[i]=3}
  else if (X.mat$age_65.79[i]==1){AGE.idx[i]=4}
  else if (X.mat$age_80.[i]==1){AGE.idx[i]=5}
}
AGE.idx <- as.data.frame(AGE.idx)

## Get min and max

Pr.AGE = cbind(Pr.OUT,AGE.idx)

min.Pr.AGE <- matrix(0,nrow=5,ncol=ncol(Pr.AGE))
max.Pr.AGE <- matrix(0,nrow=5,ncol=ncol(Pr.AGE))
for (i in 1:5){
  min.this.AGE <- Pr.AGE %>% filter(AGE.idx==i) %>% apply(2,min)
  min.Pr.AGE[i,] = min.this.AGE
  max.this.AGE <- Pr.AGE %>% filter(AGE.idx==i) %>% apply(2,max)
  max.Pr.AGE[i,] = max.this.AGE
  i=i+1
}
min.Pr.AGE <- as.data.frame(min.Pr.AGE)
colnames(min.Pr.AGE) <- colnames(Pr.AGE)
min.Pr.AGE$stat <- "low_95"
max.Pr.AGE <- as.data.frame(max.Pr.AGE)
colnames(max.Pr.AGE) <- colnames(Pr.AGE)
max.Pr.AGE$stat <- "up_95"

## Get mean
freq.X <- cbind(freq.OUT, AGE.idx)
#freq.AGE <- matrix(0,nrow=5,ncol=ncol(Pr.AGE))

weighted.Pr <- Pr.OUT*freq.OUT
weighted.Pr <- cbind(weighted.Pr, AGE.idx)
mean.Pr.AGE <- matrix(0,nrow=5,ncol=ncol(Pr.AGE))
for (i in 1:5){
  freq.this.AGE <- freq.X %>% filter(AGE.idx==i) %>% apply(2,sum)
 # freq.AGE[i,] = freq.this.AGE
  mean.this.AGE <- weighted.Pr %>% filter(AGE.idx==i) %>% apply(2,sum)
  mean.Pr.AGE[i,] = mean.this.AGE/freq.this.AGE
  i=i+1
}
mean.Pr.AGE <- as.data.frame(mean.Pr.AGE)
colnames(mean.Pr.AGE) <- colnames(Pr.AGE)
mean.Pr.AGE$stat <- "mean"
mean.Pr.AGE$AGE.idx <- NULL
mean.Pr.AGE$AGE.idx <- 1:5

## Put together and format for plotting
Pr.AGES.stats <- rbind(min.Pr.AGE,mean.Pr.AGE,max.Pr.AGE)
Pr.AGES.stats$AGE.idx <- as.factor(Pr.AGES.stats$AGE.idx)
Pr.AGES <- reshape2::melt(Pr.AGES.stats, id.vars= c("stat", "AGE.idx" ))
Pr.AGES$type <- sapply(Pr.AGES$variable, function(x) {str_split(x, "_")[[1]][1]})
Pr.AGES$date <- sapply(Pr.AGES$variable, function(x) {str_split(x, "_")[[1]][2]})
Pr.AGES$variable = NULL
Pr.AGES$date <- as.Date(Pr.AGES$date)
Pr.AGES$stat <- as.factor(Pr.AGES$stat)
Pr.AGES$type <- as.factor(Pr.AGES$type)
Pr.AGES <- spread(Pr.AGES, stat, value)

## Printing mean, up_95, low_95 for each age group over time
# Ages.names <- colnames(X.mat)[1:5] %>% as.data.frame()
# colnames(Ages.names) <- "age_group"
# Ages.names$AGE.idx <- as.factor(c(1:5))
# illness.probs.ages <- left_join(Ages.names, Pr.AGES, by="AGE.idx")
# illness.probs.ages$AGE.idx <- NULL
# write.csv(illness.probs.ages, file = path("/Users/abigailhorn/Google Drive/UCDavis/illness_progression_probs/ages_illness_probs.csv"))

#########################################################################################################
## PLOT

# Use position_dodge to move overlapped errorbars horizontally
plot.Pr.AGES <- function(Pr.AGES, type.in){
  
  traj.CI <- Pr.AGES %>% filter(type==type.in)
  
  ## Include in plot
  longnames <- c("Ages 0-18", "Ages 19-49", "Ages 50-65", "Ages 65-79", "Ages 80+")
  names(longnames) <- unique(traj.CI$AGE.idx)
  startDatePlot = as.Date("2020-03-01")
  endDatePlot = as.Date("2021-03-15")
  
  p <- ggplot(traj.CI, aes(x=date, y=mean, group=AGE.idx, color=AGE.idx)) + 
    geom_errorbar(aes(ymin=low_95, ymax=up_95), width=.2, 
                  position=position_dodge(30)) +
    geom_line(position=position_dodge(30)) + geom_point(position=position_dodge(30))+
    scale_color_brewer(palette="Paired", labels = longnames)+theme_minimal()
  
  p <- p + theme_bw()
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%b-%y")
 # p <- p + scale_y_continuous(limits = c(0,round((max(traj.CI$up_95)+.1),1)), breaks = seq(from = 0, to = round((max(traj.CI$up_95)+.1),1), by = round((max(traj.CI$up_95)+.1),1)/5))
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 12, face = "bold"))
  p <- p + ylab(type.in) + xlab(NULL) 
  p <- p + theme(legend.title=element_blank())
  p
  
}

plot.HI.AGES <- plot.Pr.AGES(Pr.AGES, "P(H|I)")
plot.QH.AGES <- plot.Pr.AGES(Pr.AGES, "P(Q|H)")
plot.DQ.AGES <- plot.Pr.AGES(Pr.AGES, "P(D|Q)")
plot.DI.AGES <- plot.Pr.AGES(Pr.AGES, "P(D|I)")

#plot.HI.AGES + plot.QH.AGES + plot.DQ.AGES + plot.DI.AGES



