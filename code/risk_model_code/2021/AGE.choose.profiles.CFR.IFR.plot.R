
CFR.choose.profile.plot <- function(is.CFR, Profile, freq.PREV.q, CFR.OUT){
  
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
  
  if (is.CFR==1) CFR.profiles <- cbind(Profile,AGE.idx, freq.PREV.q, CFR.OUT)
  if (is.CFR==0) CFR.profiles <- cbind(Profile,AGE.idx, freq.PREV.q, IFR.OUT)
  
  max.prev.CFR.AGE <- matrix(0,nrow=5,ncol=ncol(CFR.profiles)) %>% as.data.frame
  for (i in 1:5){
    CFR.this.AGE <- CFR.profiles %>% filter(AGE.idx==i) %>% as.data.frame() #%>% apply(2, max, na.rm=TRUE)
    max.profile <- CFR.this.AGE[which.max(CFR.this.AGE$freq.PREV.q),] %>% as.data.frame()
    max.prev.CFR.AGE[i,] = max.profile
  }
  colnames(max.prev.CFR.AGE)<-colnames(CFR.profiles)
  
  max.prev.CFR.AGE$freq.PREV.q=NULL
  TEST.melt <- reshape2::melt(max.prev.CFR.AGE, id.vars= c("AGE.idx","Profile"))
  
  for (i in 1:nrow(TEST.melt)){
    if (is.na( sapply(TEST.melt$variable[i], function(x) {str_split(x, "_")[[1]][3]}))) { 
      TEST.melt$stat[i] <- "mean"
      TEST.melt$date[i] <- sapply(TEST.melt$variable[i], function(x) {str_split(x, "_")[[1]][2]})
    }
    if (!is.na( sapply(TEST.melt$variable[i], function(x) {str_split(x, "_")[[1]][3]}))) {
      TEST.melt$stat[i] <- sapply(TEST.melt$variable[i], function(x) {str_split(x, "_")[[1]][2]})
      TEST.melt$date[i] <- sapply(TEST.melt$variable[i], function(x) {str_split(x, "_")[[1]][3]})
    }
  }
  TEST.melt$variable = NULL
  TEST.melt$date <- as.Date(TEST.melt$date)
  TEST.melt$stat <- as.factor(TEST.melt$stat)
  TEST.melt$Profile <- as.factor(TEST.melt$Profile)
  TEST.melt <- spread(TEST.melt, stat, value)
  
  #####################
  ### PLOT
  #####################
  
  longnames <- c("Profile 2: 0-18|BMI<25|~Comorb|~Smoker",
                 "Profile 9: 19-49|BMI<25|~Comorb|~Non Smoker",
                 "Profile 14: 50-64|25<BMI<30|~Comorb|~Smoker",
                 "Profile 43: 65-79|BMI<25|Comorb|~Smoker",
                 "Profile 44: 80+|BMI<25|Comorb|~Smoker"
  )
  
  #longnames <- as.character(unique(TEST.melt$Profile))
  
  names(longnames) <- as.character(unique(TEST.melt$Profile))
  
  ## Colors
  
  cols.list <- c(
    "salmon",
    "olivedrab4",
    "mediumaquamarine",
    "lightskyblue",
    "mediumorchid")
  
  names(cols.list) <- names(longnames)
  color.this.var <- as.character(cols.list[ as.character(unique(TEST.melt$Profile)) ])
  
  traj.CI <- TEST.melt
  
  p <- ggplot(data = traj.CI,
              aes(x = date,
                  y = mean, ymin = low.95, ymax = up.95,
                  color = Profile,
                  fill = Profile,
                  group = Profile))
  
  p <- p +  geom_ribbon(data = traj.CI,
                        aes(x = date,
                            y = mean, ymin = low.95, ymax = up.95,
                            color = Profile,
                            fill = Profile,
                            group = Profile),alpha = .2, inherit.aes=TRUE, color=FALSE)
  
  p <- p +  scale_fill_manual(values = c(color.this.var),labels = longnames) + scale_color_manual(values = c(color.this.var), labels = longnames)
  p <- p + geom_line() + geom_ribbon(alpha = 0.2, color = FALSE)
  
  p <- p + theme_bw()
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%b-%y")
  # p <- p + scale_y_continuous(limits = c(0,round((max(traj.CI$up_95)+.1),1)), breaks = seq(from = 0, to = round((max(traj.CI$up_95)+.1),1), by = round((max(traj.CI$up_95)+.1),1)/5))
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 12, face = "bold"))
  type.in <- sapply(colnames(CFR.OUT[1]), function(x) {str_split(x, "_")[[1]][1]})
  p <- p + ylab(type.in) + xlab(NULL) 
  p <- p + theme(legend.title=element_blank())
  
}


# plot.CFR.chooes = CFR.choose.profile.plot(is.CFR=1, Profile=Profile, freq.PREV.q=freq.PREV.q, CFR.OUT=CFR.OUT)
# plot.IFR.choose = CFR.choose.profile.plot(is.CFR=0, Profile=Profile, freq.PREV.q=freq.PREV.q, CFR.OUT=IFR.OUT)
# plot.IFR.choose


