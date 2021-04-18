
########################################################################################
## FUNCTION TO CALCULATE CI
########################################################################################

get.CI.SCENARIOS <- function(TEST.out, init.date.data) {
  
  library(data.table)
  init.date.data= init.date.data #"2020-02-01"
  
  ### Add CFR and IFR to the list (EXTRA STEP NOW THAT THIS IS BEING USED ALSO FOR summary_table)
  traj <- dplyr::mutate(TEST.out, I_stu=I_on+I_off, Q_stu=Q_on+Q_off, P_stu=P_on+P_off)
  
  #traj <- dplyr::mutate(TEST.out, Itot=I+A, CFRobs=(D/Idetectcum), CFRactual=(D/(Itotcum)) )
  #traj <-  dplyr::select(traj,c(1:7,CFRobs,CFRactual,vars.plus.R))
  
  # # Add variable for Test_new (daily tests)
  # for (i in 1:max(traj$run)) {
  #   traj_subset = subset(traj, run == i)
  #   traj_subset = traj_subset[order(traj_subset$time),]
  #   
  #   # add starting level for new cases and deaths
  #   traj_subset$Test_new = traj_subset$Test[1]
  # 
  #   for (j in 2:nrow(traj_subset)) {
  #     traj_subset$Test_new[j] = round(traj_subset$Test[j] - traj_subset$Test[j-1])
  #   }
  #   
  #   # include in main dataset
  #   traj$Test_new[traj$run==i] = traj_subset$Test_new
  # }
  
  
  ###
  colnames(traj)[colnames(traj) == "time"] <- "date"
  
  ## TO SAVE MEMORY
  rm(TEST.out)
  
  print("Starting CI calc")
  
  ### MELTING AND APPLYING SUMMARY STAT FUNCTIONS
  df.traj <- reshape2::melt(traj, measure.vars = c(4:ncol(traj)), variable.name = "state.name")
  df.traj_dt <- as.data.table(df.traj)
  
  traj.CI <- df.traj_dt[, list(
    N=.N,
    mean = mean(value),
    median = quantile(value, c(.5),na.rm=TRUE),
    low_95 = quantile(value, c(.025),na.rm=TRUE),
    up_95 = quantile(value, c(.975),na.rm=TRUE),
    up_50 = quantile(value,.75,na.rm=TRUE),
    low_50 = quantile(value,.25,na.rm=TRUE)),
    by = c("date", "state.name", "scenario.name")]
  traj.CI <- as.data.frame(traj.CI)
  
  ## TO ALIGN DATES: MODEL
  init.date = init.date.data #"2020-01-03"
  init.date <- as.Date(init.date)
  traj.CI[["date"]] <- traj.CI[["date"]] + init.date
  
  return(traj.CI)
  
}


########################################################################################
## FUNCTION TO PLOT SCENARIO CI
########################################################################################

plot <- plot.SCENARIOS(traj.CI=traj.CI, endDatePlot=endDatePlot, startDatePlot=startDatePlot, vars.to.plot=vars.to.plot, longnames=longnames, filter.scenarios=filter.scenarios)

## ADD SCENARIO ID

plot.SCENARIOS <- function(traj.CI, endDatePlot, startDatePlot, vars.to.plot, longnames, filter.scenarios) {
  
  #longnames <- c("Inf. high risk", "Inf. medium risk", "Inf. low risk", "Number of Tests")
  names(longnames) <- vars.to.plot #c("I_on","I_off","I_saf","Test")
  
  ## Filter only to variables of interest
  traj.CI <- traj.CI %>%  dplyr::filter(state.name %in% vars.to.plot) 
  
  if (!is.null(filter.scenarios)){
    traj.CI <- traj.CI %>% dplyr::filter(scenario.name %in% levels(traj.CI$scenario.name)[filter.scenarios])
  }
  
  ## Select only more recent dates
  init.date <- as.Date(startDatePlot)
  endDatePlot <- as.Date(endDatePlot) #startDatePlot + time.steps.4plot #- 40  # the constant 40 because the traj are not aligned to start date
  traj.CI <- traj.CI %>% dplyr::filter(date >= startDatePlot) %>% dplyr::filter(date < endDatePlot)
  
  # ### HARD CODED -- COME BACK TO GENERALIZE
  # traj.CI$scenario.id <- rep("scenario_1", times=nrow(traj.CI))
  
  ## PLOTTING
  traj.CI.line <- reshape2::melt(traj.CI[c("date","scenario.name", "state.name", "median")], id.vars = c("date", "scenario.name", "state.name"))
  traj.CI.area <- reshape2::melt(traj.CI[c("date","scenario.name", "state.name", "low_95", "low_50", "up_50", "up_95")], id.vars = c("date", "scenario.name", "state.name"))
  traj.CI.area$type <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj.CI.area$CI <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj.CI.area$variable <- NULL
  traj.CI.area <- reshape2::dcast(traj.CI.area, "date+scenario.name+state.name+CI~type")
  
  p <- ggplot(transform(traj.CI.area, state.name = factor(state.name, levels=vars.to.plot)))
  
  #####################
  ### colors and names
  #####################
  
  
  
  ## Colors
  
  # cols.list <- c(
  #   "salmon",
  #   "sandybrown",
  #   "navajowhite3",
  #   "olivedrab4",
  #   "olivedrab2",
  #   "mediumseagreen",
  #   "mediumaquamarine",
  #   "mediumturquoise",
  #   "cyan2",
  #   "lightskyblue",
  #   "steelblue2",
  #   "mediumpurple",
  #   "mediumorchid",
  #   "plum1",
  #   "violetred1",
  #   "deeppink4",
  #   "grey50",
  #   "mediumturquoise",
  #   "lightskyblue",
  #   "violetred1",
  #   "grey50",
  #   "grey50",
  #   "grey50"
  # )
  
  # names(cols.list) <- names(longnames)
  # color.this.var <- as.character(cols.list[vars.to.plot])
  
  # ## CAPACITY DATA FRAME
  # capacity.vals <- as.data.frame(matrix(NA, nrow=length(levels(traj.CI$state.name)), ncol=2))
  # capacity.vals[,1] <- levels(traj.CI$state.name)
  # rownames(capacity.vals) <- levels(traj.CI$state.name)
  # capacity.vals["Q_on",2] <- 750 
  # capacity.vals["Q_off",2] <- 750
  # capacity.vals["Q_saf",2] <-750
  # colnames(capacity.vals) <- c("state.name","capacity")
  
  ## PLOT OPTIONS
  #  p <- p + facet_grid(state.name ~ scenario.id, labeller=labeller(state.name=longnames, scenario.id=longnames.scenarios), scales='free')
  p <- p + facet_grid(state.name ~ scenario.name, labeller=labeller(state.name=longnames), scales='free')
  p <- p + geom_ribbon(data = traj.CI.area, aes_string(x = "date", ymin = "low", ymax = "up", alpha = "CI", fill = "state.name"),show.legend = c(fill=FALSE))
  p <- p + geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable", colour = "state.name"), size = 1, show.legend = c(colour=FALSE))
  
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.20, "50" = 0.50), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  
  # ## ADD CAPACITY
  # capacity.vals <- capacity.vals %>% filter(state.name %in% vars.to.plot)
  # p <- p + geom_hline(data= capacity.vals, aes(yintercept=capacity),linetype = "dashed")
  
 ## Finish plot 
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 week" , date_labels = "%d-%b-%y")
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 10, face = "bold"),
                 strip.text.y = element_text(size=10, face="bold"))
  p <- p + ylab("Numbers in Compartments") + xlab(NULL)
  p <- p + scale_y_continuous(labels = scales::comma)
  p <- p + theme(strip.background = element_rect(colour="black", fill="white", 
                                                 size=1, linetype="solid"))
  p
}


