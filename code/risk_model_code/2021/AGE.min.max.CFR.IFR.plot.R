

#########################################################################################################
#########################################################################################################
## Plot mean (min,max) of each probability over profiles within each age group

## Get AGE.idx as a column

plot.CFR.IFR.ages <- function(is.CFR, CFR.OUT, freq.OUT, X.mat){
  
  n_profiles <- nrow(X.mat)
  
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
  #CFR.AGE <- CFR.OUT %>% select(c(paste0("CFR_",times.dates))) %>% cbind(AGE.idx)
  if (is.CFR==1) CFR.AGE <- CFR.OUT %>% select(c(paste0("CFR_",times.dates))) %>% cbind(AGE.idx)
  if (is.CFR==0) CFR.AGE <- CFR.OUT %>% select(c(paste0("IFR_",times.dates))) %>% cbind(AGE.idx)
  CFR.AGE[CFR.AGE[,]==0]=NA
  min.CFR.AGE <- matrix(0,nrow=5,ncol=ncol(CFR.AGE))
  max.CFR.AGE <- matrix(0,nrow=5,ncol=ncol(CFR.AGE))
  for (i in 1:5){
    min.this.AGE <- CFR.AGE %>% filter(AGE.idx==i) %>% apply(2,min,na.rm=TRUE)
    min.CFR.AGE[i,] = min.this.AGE
    max.this.AGE <- CFR.AGE %>% filter(AGE.idx==i) %>% apply(2,max,na.rm=TRUE)
    max.CFR.AGE[i,] = max.this.AGE
    i=i+1
  }
  min.CFR.AGE <- as.data.frame(min.CFR.AGE)
  colnames(min.CFR.AGE) <- colnames(CFR.AGE)
  min.CFR.AGE$stat <- "low_95"
  max.CFR.AGE <- as.data.frame(max.CFR.AGE)
  colnames(max.CFR.AGE) <- colnames(CFR.AGE)
  max.CFR.AGE$stat <- "up_95"
  
  ## Get mean
  freq.OUT.I <- freq.OUT[,c(4*(seq.times-1)+1)]
  freq.X <- cbind(freq.OUT.I, AGE.idx)
  weighted.Pr <- (CFR.AGE %>% select(-"AGE.idx"))*freq.OUT.I
  weighted.Pr <- cbind(weighted.Pr, AGE.idx)
  mean.CFR.AGE <- matrix(0,nrow=5,ncol=ncol(CFR.AGE))
  for (i in 1:5){
    freq.this.AGE <- freq.X %>% filter(AGE.idx==i) %>% apply(2,sum,na.rm=TRUE)
    mean.this.AGE <- weighted.Pr %>% filter(AGE.idx==i) %>% apply(2,sum,na.rm=TRUE)
    mean.CFR.AGE[i,] = mean.this.AGE/freq.this.AGE
    i=i+1
  }
  mean.CFR.AGE <- as.data.frame(mean.CFR.AGE)
  colnames(mean.CFR.AGE) <- colnames(CFR.AGE)
  mean.CFR.AGE$stat <- "mean"
  mean.CFR.AGE$AGE.idx <- NULL
  mean.CFR.AGE$AGE.idx <- 1:5
  
  ## Put together and format for plotting
  CFR.AGES.stats <- rbind(min.CFR.AGE,mean.CFR.AGE,max.CFR.AGE)
  CFR.AGES.stats$AGE.idx <- as.factor(CFR.AGES.stats$AGE.idx)
  CFR.AGES <- reshape2::melt(CFR.AGES.stats, id.vars= c("stat", "AGE.idx" ))
  CFR.AGES$type <- sapply(CFR.AGES$variable, function(x) {str_split(x, "_")[[1]][1]})
  CFR.AGES$date <- sapply(CFR.AGES$variable, function(x) {str_split(x, "_")[[1]][2]})
  CFR.AGES$variable = NULL
  CFR.AGES$date <- as.Date(CFR.AGES$date)
  CFR.AGES$stat <- as.factor(CFR.AGES$stat)
  CFR.AGES$type <- as.factor(CFR.AGES$type)
  CFR.AGES <- spread(CFR.AGES, stat, value)
  
  #########################################################################################################
  ## PLOT
  
  # Use position_dodge to move overlapped errorbars horizontally
  plot.CFR.AGES <- function(CFR.AGES, type.in){
    
    traj.CI <- CFR.AGES %>% filter(type==type.in)
    
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
  
  ## Get plot
  if (is.CFR==1) plot.CFR.AGES <- plot.CFR.AGES(CFR.AGES, "CFR")
  if (is.CFR==0) plot.CFR.AGES <- plot.CFR.AGES(CFR.AGES, "IFR")
  return(plot.CFR.AGES)
}

plot.CFR.NEW = plot.CFR.IFR.ages(is.CFR=1, CFR.OUT=CFR.OUT, freq.OUT=freq.OUT, X.mat=X.mat)
plot.IFR.NEW = plot.CFR.IFR.ages(is.CFR=0, CFR.OUT=IFR.OUT, freq.OUT=freq.OUT, X.mat=X.mat)
plot.IFR.NEW


