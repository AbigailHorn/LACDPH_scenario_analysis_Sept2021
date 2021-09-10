
data.in = la_data
colnames(data.in) = c("date","1_Infected","1_New_Infections","4_Died","4_New_Deaths","2_In_Hospital","3_In_ICU")

################################################################################
## CUMULATIVE COUNTS

#traj.CI = traj.0
traj.CI = traj.base
scenario.IDs = c("EveryoneVax.05_Rt.1.00","NooneVax_Rt.1.00")

######## Vars to plot
vars.to.plot.scenario = c("Idetectcum", "Htotcum", "Qcum", "D")
traj.CI <- traj.CI %>% dplyr::filter(state.name %in% vars.to.plot.scenario) 

######## Scenarios to plot
#scenario.IDs = c("NooneVax_Rt.0.80","NooneVax_Rt.0.90","NooneVax_Rt.1.00","NooneVax_Rt.1.10","NooneVax_Rt.1.20")
#scenario.IDs = c("NooneVax_Isolate.00","NooneVax_Rt.1.00","NooneVax_Isolate.20")
#scenario.IDs = c("EveryoneVax.05_Rt.0.80","EveryoneVax.05_Rt.0.90","EveryoneVax.05_Rt.1.00","EveryoneVax.05_Rt.1.10","EveryoneVax.05_Rt.1.20")
#scenario.IDs = "EveryoneVax.05_Rt.1.00"
#scenario.IDs = "NooneVax_Rt.1.00"

traj.CI <- traj.CI %>% dplyr::filter(scenario %in% scenario.IDs)

traj.CI$state.name <- recode_factor(traj.CI$state.name, 
                                    Idetectcum = "1_Infected",
                                    Htotcum = "2_Hospitalized",
                                    Qcum = "3_ICU",
                                    D = "4_Died")
vars.to.plot.relevel = c("1_Infected", "2_Hospitalized", "3_ICU", "4_Died")

time.steps.plot = 365
plot.SCENARIOS(traj.CI=traj.CI, data.in=data.in, time.steps.plot=365, vars.to.plot=vars.to.plot.relevel, add.capacity = FALSE)


################################################################################
## PEAK COUNTS

#traj.CI = traj.0
traj.CI = traj.base
scenario.IDs = c("EveryoneVax.05_Rt.1.00","NooneVax_Rt.1.00")

######## Vars to plot
vars.to.plot.scenario = c("I_detect_new", "Htot", "Q", "D_new")
traj.CI <- traj.CI %>% dplyr::filter(state.name %in% vars.to.plot.scenario) 

######## Scenarios to plot
#scenario.IDs = c("NooneVax_Rt.0.80","NooneVax_Rt.0.90","NooneVax_Rt.1.00","NooneVax_Rt.1.10","NooneVax_Rt.1.20")
#scenario.IDs = c("NooneVax_Isolate.00","NooneVax_Rt.1.00","NooneVax_Isolate.20")
#scenario.IDs = c("EveryoneVax.05_Rt.0.80","EveryoneVax.05_Rt.0.90","EveryoneVax.05_Rt.1.00","EveryoneVax.05_Rt.1.10","EveryoneVax.05_Rt.1.20")
#scenario.IDs = "EveryoneVax.05_Rt.1.00"
#scenario.IDs = "NooneVax_Rt.1.00"

traj.CI <- traj.CI %>% dplyr::filter(scenario %in% scenario.IDs)

traj.CI$state.name <- recode_factor(traj.CI$state.name, 
                                    I_detect_new = "1_New_Infections",
                                    Htot = "2_In_Hospital",
                                    Q = "3_In_ICU",
                                    D_new = "4_New_Deaths")
vars.to.plot.relevel = c("1_New_Infections", "2_In_Hospital", "3_In_ICU", "4_New_Deaths")

plot.SCENARIOS(traj.CI=traj.CI, data.in=data.in, time.steps.plot=100, vars.to.plot=vars.to.plot.relevel, add.capacity=TRUE)

################################################################################
## FUNCTION
plot.SCENARIOS <- function(traj.CI, data.in, time.steps.plot, vars.to.plot, add.capacity=TRUE) {
  
  ## Filter only to variables of interest
  #traj.CI <- traj.CI %>% dplyr::filter(state.name %in% vars.to.plot) 
  
  # if (!is.null(filter.scenarios)){
  #   traj.CI <- traj.CI %>% dplyr::filter(scenario.id %in% levels(traj.CI$scenario.id)[filter.scenarios])
  # }
  
  ## Select only more recent dates
  startDatePlot <- as.Date(start_date) #init.date #- date.offset.4plot -1 #15
  endDatePlot <- startDatePlot + time.steps.plot # - 40
  traj.CI <- traj.CI %>% dplyr::filter(date >= startDatePlot) %>% dplyr::filter(date < endDatePlot)
  
  if(!is.null(data.in)){
    data <- data.in %>% dplyr::filter(date > startDatePlot-1)
    data <- reshape2::melt(data, measure.vars = c(2:ncol(data)), variable.name = "state.name")
    dataCum = data %>% filter(state.name %in% c("4_Died", "1_Infected")) %>% group_by(state.name) %>% mutate(value = value - dplyr::first(value)) %>% as.data.frame()
    dataNew = data %>% filter(state.name %in% c("4_Died", "1_Infected") == FALSE) %>% as.data.frame()
    data=as.data.frame(rbind(dataCum,dataNew))
    
    data <- data %>% dplyr::filter(state.name %in% vars.to.plot) 
    
    # rebase0 <- function(data,thisVar){
    #   data_thisVar = data %>% dplyr::filter(state.name %in% thisVar) %>% mutate(value = value-min(value))
    #   return(rbind(data_thisVar,data))
    # }
    # data = rebase0(data, "D")
    # data = rebase0(data,"Idetectcum")
  }
  
  ## PLOTTING
  traj.CI.line <- reshape2::melt(traj.CI[c("date","scenario", "state.name", "median")], id.vars = c("date", "scenario", "state.name"))
  traj.CI.area <- reshape2::melt(traj.CI[c("date","scenario", "state.name", "low_95", "low_50", "up_50", "up_95")], id.vars = c("date", "scenario", "state.name"))
  traj.CI.area$type <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj.CI.area$CI <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj.CI.area$variable <- NULL
  traj.CI.area <- reshape2::dcast(traj.CI.area, "date+scenario+state.name+CI~type")
  
  #traj.CI.area$state.name = factor(traj.CI.area$state.name, levels = vars.to.plot)
  
  p <- ggplot(traj.CI.area)
  
  ## CAPACITY DATA FRAME
  capacity.vals <- as.data.frame(matrix(NA, nrow=length(levels(traj.CI$state.name)), ncol=2))
  capacity.vals[,1] <- levels(traj.CI$state.name)
  rownames(capacity.vals) <- levels(traj.CI$state.name)
  capacity.vals["2_In_Hospital",2] <- 4000 
  capacity.vals["3_In_ICU",2] <- 2245
  capacity.vals["V",2] <-1000
  colnames(capacity.vals) <- c("state.name","capacity")
  
  ## PLOT OPTIONS
  #  p <- p + facet_grid(state.name ~ scenario.id, labeller=labeller(state.name=longnames, scenario.id=longnames.scenarios), scales='free')
  p <- p + facet_grid(state.name ~ scenario, scales='free') #labeller=labeller(state.name=longnames), scales='free')
  p <- p + geom_ribbon(data = traj.CI.area, aes_string(x = "date", ymin = "low", ymax = "up", alpha = "CI", fill = "state.name"),show.legend = c(fill=FALSE))
  p <- p + geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable", colour = "state.name"), size = 1, show.legend = c(colour=FALSE))
  
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.20, "50" = 0.50), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  
  ## ADD CAPACITY
  if (add.capacity==TRUE){
    capacity.vals <- capacity.vals %>% filter(state.name %in% vars.to.plot)
    p <- p + geom_hline(data= capacity.vals, aes(yintercept=capacity),linetype = "dashed")
  }
  
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


