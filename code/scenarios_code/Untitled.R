
plot.SCENARIOS <- function(traj.CI, data.in, endDatePlot, vars.to.plot, filter.scenarios) {
  
  ## Filter only to variables of interest
  traj.CI <- traj.CI %>%  dplyr::filter(state.name %in% vars.to.plot) 
  
  if (!is.null(filter.scenarios)){
    traj.CI <- traj.CI %>% dplyr::filter(scenario.id %in% levels(traj.CI$scenario.id)[filter.scenarios])
  }
  
  ## Select only more recent dates
  init.date <- as.Date("2020-03-01")
  startDatePlot <- init.date #- date.offset.4plot -1 #15
  endDatePlot <- as.Date(endDatePlot) #startDatePlot + time.steps.4plot #- 40  # the constant 40 because the traj are not aligned to start date
  traj.CI <- traj.CI %>% dplyr::filter(date >= startDatePlot) %>% dplyr::filter(date < endDatePlot)
  
  if(!is.null(data.in)){
    ## ALIGN DATES: DATA
    no_obs <- nrow(data.in)
    step <- 0:(no_obs-1)
    date <- init.date + step
    data.date <- cbind(date,data.in)
    rownames(data.date) <- step
    
    ## Select only more recent dates
    data.date <- data.date %>% dplyr::filter(date > startDatePlot) %>% dplyr::filter(date < endDatePlot)
    data <- reshape2::melt(data.date, measure.vars = c(2:ncol(data.date)), variable.name = "state.name")
  }
  
  ## PLOTTING
  traj.CI.line <- reshape2::melt(traj.CI[c("date","scenario.id",  "protect.id", "NPI.id", "state.name", "median")], id.vars = c("date", "scenario.id",  "protect.id", "NPI.id","state.name"))
  traj.CI.area <- reshape2::melt(traj.CI[c("date","scenario.id",  "protect.id", "NPI.id", "state.name", "low_95", "low_50", "up_50", "up_95")], id.vars = c("date", "scenario.id",  "protect.id", "NPI.id","state.name"))
  traj.CI.area$type <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj.CI.area$CI <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj.CI.area$variable <- NULL
  traj.CI.area <- reshape2::dcast(traj.CI.area, "date+scenario.id+state.name+CI~type")
  
  p <- ggplot(transform(traj.CI.area, state.name = factor(state.name, levels=vars.to.plot)))
  
  #####################
  ### colors and names
  #####################
  
  longnames <- c("Susceptible",
                 "New Obs. Infected",
                 "Current Obs. Infected",
                 "Cum. Obs. Infected",
                 "Current Tot. Infected",
                 "Cum. Tot. Infected",
                 "New in Hospital",
                 "Current in Hospital",
                 "Cum. in Hospital",
                 "Current in ICU",
                 "Cum. in ICU",
                 "Current Ventilation",
                 "Cum. Ventilation",
                 "New Deaths",
                 "Cum. Deaths",
                 "Recovered",
                 "R0(t)",
                 "Alpha(t)",
                 "Kappa(t)",
                 "Delta(t)",
                 "r(t)",
                 "CFR",
                 "IFR"
  )
  
  names(longnames) <-  c(
    "S",
    "I_detect_new",
    "I",
    "Idetectcum",
    "Itot",
    "Itotcum",
    "H_new",
    "Htot",
    "Htotcum",
    "Q",
    "Qcum",
    "V",
    "Vcum",
    "D_new",
    "D",
    "R",
    "Rt",
    "Alpha_t",
    "Kappa_t",
    "Delta_t",
    "r_t",
    "CFRobs",
    "CFRactual"
  )
  
  ## Colors
  
  cols.list <- c(
    "salmon",
    "sandybrown",
    "navajowhite3",
    "olivedrab4",
    "olivedrab2",
    "mediumseagreen",
    "mediumaquamarine",
    "mediumturquoise",
    "cyan2",
    "lightskyblue",
    "steelblue2",
    "mediumpurple",
    "mediumorchid",
    "plum1",
    "violetred1",
    "deeppink4",
    "grey50",
    "mediumturquoise",
    "lightskyblue",
    "violetred1",
    "grey50",
    "grey50",
    "grey50"
  )
  
  names(cols.list) <- names(longnames)
  color.this.var <- as.character(cols.list[vars.to.plot])
  
  ## CAPACITY DATA FRAME
  capacity.vals <- as.data.frame(matrix(NA, nrow=length(levels(traj.CI$state.name)), ncol=2))
  capacity.vals[,1] <- levels(traj.CI$state.name)
  rownames(capacity.vals) <- levels(traj.CI$state.name)
  capacity.vals["Htot",2] <- 4000 
  capacity.vals["Q",2] <- 2245
  capacity.vals["V",2] <-1000
  colnames(capacity.vals) <- c("state.name","capacity")
  
  ## PLOT OPTIONS
  #  p <- p + facet_grid(state.name ~ scenario.id, labeller=labeller(state.name=longnames, scenario.id=longnames.scenarios), scales='free')
  p <- p + facet_grid(state.name ~ scenario.id, labeller=labeller(state.name=longnames), scales='free')
  p <- p + geom_ribbon(data = traj.CI.area, aes_string(x = "date", ymin = "low", ymax = "up", alpha = "CI", fill = "state.name"),show.legend = c(fill=FALSE))
  p <- p + geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable", colour = "state.name"), size = 1, show.legend = c(colour=FALSE))
  
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.20, "50" = 0.50), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  
  ## ADD CAPACITY
  capacity.vals <- capacity.vals %>% filter(state.name %in% vars.to.plot)
  p <- p + geom_hline(data= capacity.vals, aes(yintercept=capacity),linetype = "dashed")
  
  ## ADD DATA
  if(!is.null(data.in)){
    p <- p + geom_point(data = data, aes_string(x = "date", y = "value"), size = .5, colour = "black")
  }
  
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%d-%b-%y")
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 8, face = "bold"))
  p <- p + ylab("Numbers in Compartments") + xlab(NULL)
  p <- p + scale_y_continuous(labels = scales::comma)
  p <- p + theme(strip.background = element_rect(colour="black", fill="white", 
                                                 size=1, linetype="solid"))
  p
}

