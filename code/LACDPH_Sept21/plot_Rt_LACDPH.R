
########################################################################################
## GET Rt(t,y)
########################################################################################

p_tab <- read.csv(path(data.dir, paste0(scenario.IDs[this.scenario],".csv")), sep=",", header=TRUE)
start_date <- p_tab$Value[which(p_tab$Var=="start_date")] %>% as.Date()
#p_tab$Value <-gsub(",","",p_tab$Value)
p_tab$Value <- as.numeric(as.character(p_tab$Value))
p_tab$Lower <- as.numeric(as.character(p_tab$Lower))
p_tab$Upper <- as.numeric(as.character(p_tab$Upper))

########## Fixed parameters
everyone_vac <- p_tab$Value[which(p_tab$Var=="everyone_vac")]
N_LA <- p_tab$Value[which(p_tab$Var=="N_LA")]
N_eligible <- p_tab$Value[which(p_tab$Var=="N_eligible")]

########## Variables to be sampled over replicates
## Define functions to generate replicates
cut01 <- function(x) pmax(pmin(x, 1), 0)
getRandomNorm <- function(rep, value) {
  replicatedVar <- 
    rnorm(rep,
          mean=p_tab$Value[which(p_tab$Var==value)],
          sd=calcStDevFromLowerUpper(mean=p_tab$Value[which(p_tab$Var==value)],
                                     lower=p_tab$Lower[which(p_tab$Var==value)],
                                     upper=p_tab$Upper[which(p_tab$Var==value)]))
  return(replicatedVar)
}
calcStDevFromLowerUpper <- function(mean=mean, lower=lower, upper=upper, CI=.95) {
  z.stat <- qnorm((CI + (1-CI)/2))
  se.u <- (upper-mean)/z.stat
  se.l <- (mean-lower)/z.stat
  mean(c(se.u, se.l))
}
## Get variables and replicates
vac_efficacy <- cut01(getRandomNorm(replicates, "vac_efficacy"))
N_vac <- getRandomNorm(replicates, "N_vac")
isolate_prop <- cut01(getRandomNorm(replicates, "isolate_prop"))
#isolate_prop[is.na(isolate_prop)] <- 0
vac_prop <- N_vac/N_LA                                        # Proportion of LAC population vaccinated

I_ini <- round(getRandomNorm(replicates, "I_ini"))
E_ini <- round(getRandomNorm(replicates, "E_ini"))
## Get S_ini based on scenario 
if(everyone_vac == 1){
  S_ini = N_LA - N_eligible*vac_efficacy    # Susceptible are all <12 yrs + those for whom vaccine was not effective  
  #S_ini = N_LA - N_vac*vac_efficacy    # Susceptible are all <12 yrs + those for whom vaccine was not effective  
}  else{
  S_ini = N_LA*(1-isolate_prop)          # Susceptible are all except those that self-isolate  
}

Alpha <- cut01(getRandomNorm(replicates, "Alpha"))
Kappa <- cut01(getRandomNorm(replicates, "Kappa"))
Delta <- cut01(getRandomNorm(replicates, "Delta"))
r <- cut01(getRandomNorm(replicates, "r"))

##########################################################################################
########## Get Rt over time

## Evaluate Rt from csv file
fn_t_readin_path <- path(data.dir, "fn_t_readin.csv")
fn_t_readin = as.data.frame(read.csv(fn_t_readin_path, sep=",",stringsAsFactors = FALSE))
Rt_dates <- as.Date(fn_t_readin$t)
Rt_t <- round(as.numeric(Rt_dates - start_date))
Rt_y <- fn_t_readin$Rt

## Get Rt in vaxed and unvaxed populations using weighting equations
Rt_uv_v <- function(Rt_eff, vac_efficacy, vac_prop, everyone_vac){
  w_v = vac_efficacy*vac_prop
  w_uv = 1-w_v
  Rt = Rt_eff #* 1/w_uv  # For LA overall, adjusted to get Rt based on behavior alone, not including size of susceptible population
  Rt_v = round(Rt/(5-4*w_v),2)   # For vaccinated
  Rt_uv = round(5*Rt_v,2)        # For unvaccinated
  #if (everyone_vac==1) return(Rt_v) else return(Rt_uv)
  #return(Rt_v)
  return(Rt_uv)
}

## Apply function to get Rt(t,y) over replicates of vac_efficacy and vac_prop
Rt_vac=sapply(Rt_y,Rt_uv_v,vac_efficacy=vac_efficacy,vac_prop=vac_prop,everyone_vac=0)               # Get Rt(t,y) over iterations - use this version if want to simulate over multiple Rt values
#Rt_vac = mapply(Rt_uv_v,Rt_eff=Rt_y,vac_efficacy=vac_efficacy, vac_prop=vac_prop)   # Get Rt(t,y) - use this version if want one Rt_y for each Rt_t

#################################
## Confidence intervals
posterior.CI <- function(posterior.var, round.by=4){
  #median = quantile(posterior.var, c(.5), na.rm=TRUE)
  low_95 = quantile(posterior.var, c(.025), na.rm=TRUE)
  low_50 = quantile(posterior.var, c(.25), na.rm=TRUE)
  mean = mean(posterior.var)
  up_50 = quantile(posterior.var, c(.75), na.rm=TRUE)
  up_95 = quantile(posterior.var, c(.975), na.rm=TRUE)
  posterior.CI <- as.data.frame(cbind(low_95,low_50,mean,up_50,up_95))
  posterior.CI <- round(posterior.CI, digits=round.by)
  return(posterior.CI)
}
ABC.par.CI <- apply(Rt_vac, MARGIN=2, FUN=posterior.CI)  # Returns a list for each time point

## Get Rt_y
Rt_y <- as.data.frame(do.call(rbind, ABC.par.CI))
colnames(Rt_y) <- names(ABC.par.CI[[1]])
Rt_dates_modify = Rt_dates
Rt_dates_modify[length(Rt_dates)] <- "2022-05-31"
Rt_y$date <- Rt_dates_modify

# Rt_y$state.name = "Rt_v"
# Rt_v = Rt_y
# 
# Rt_y$state.name = "Rt_uv"
# Rt_uv = Rt_y

Rt_y <- rbind(Rt_v, Rt_uv)

## Plot together
test = plot.together.capacity(traj.CI=Rt_y,  vars.to.plot=c("Rt_v","Rt_uv"), data.in=NULL)
test + geom_hline(yintercept =1,linetype="dashed")

###### NOTE TO USE ######
## Output first Rt_uv and second Rt_v and then bind to get Rt_y, which is what we plot here.

# 
#   data.processed <- reshape2::melt(Rty, measure.vars = c(2:ncol(data.in)), variable.name = "state.name")
#   
#   test <- test + geom_point(data = data.processed,
#                       aes(x = date, y = value,
#                           color = state.name),
#                       alpha = 0.7,
#                       inherit.aes = TRUE)

## Plot Single
traj.CI <- Rt_y
var.to.plot <- "Rt"  
time.steps.4plot <- as.numeric(as.Date(endDatePlot) - as.Date("2020-03-01"))
R_t_plot <- plot.model.single(traj.CI, data.in=NULL, init.date.data=NULL, ymax=4, var.to.plot=var.to.plot)
R_t_plot



























########################################################################################
## PLOT SINGLE
########################################################################################

plot.model.single <- function(traj.CI, data.in=NULL, init.date.data=NULL, time.steps.plot=365, ymax=NULL, var.to.plot=NULL) {
  
  ## Filter only to variable of interest
  traj.CI <- traj.CI %>%  dplyr::filter(state.name==var.to.plot)
  
  ## Select only more recent dates
  #init.date <- init.date.data
  init.date <- as.Date("2021-06-01") #as.Date(lubridate::ydm(init.date))
  startDatePlot <- init.date #- date.offset.4plot #15
  endDatePlot <- startDatePlot + 365 #- 40  # the constant 40 because the traj are not aligned to start date
  traj.CI <- traj.CI %>% dplyr::filter(date > startDatePlot-1) %>% dplyr::filter(date < endDatePlot)
  
  y.max.in <- max(traj.CI$up_95)
  
  ## Data in -- plot only for the selected variable
  if(!is.null(data.in)){
    
    if(var.to.plot %in% colnames(data.in)) {  # FIX LATER -- REMOVE TMP
      
      ## Filter only to variable of interest
      data.in<- data.in %>% dplyr::select(var.to.plot)
      
      ## ALIGN DATES: DATA
      no_obs <- nrow(data.in)
      step <- 0:(no_obs-1)
      date <- init.date + step
      data.date <- cbind(date,data.in)
      rownames(data.date) <- step
      
      ## Select only more recent dates
      data.date <- data.date %>% dplyr::filter(date > startDatePlot)
      data <- reshape2::melt(data.date, measure.vars = c(2:ncol(data.date)), variable.name = "state.name")
    }
    
    else {data.in = NULL}
  }
  
  ## PLOTTING
  #traj.CI.line <- reshape2::melt(traj.CI[c("date", "state.name", "mean", "median")], id.vars = c("date", "state.name"))
  traj.CI.line <- reshape2::melt(traj.CI[c("date", "state.name", "mean")], id.vars = c("date", "state.name"))
  traj.CI.area <- reshape2::melt(traj.CI[c("date", "state.name", "low_95", "low_50", "up_50", "up_95")], id.vars = c("date", "state.name"))
  traj.CI.area$type <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj.CI.area$CI <- sapply(traj.CI.area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj.CI.area$variable <- NULL
  traj.CI.area <- reshape2::dcast(traj.CI.area, "date+state.name+CI~type")
  
  longnames <- c("R(t) Effective",
                 "R(t)")
  
  names(longnames) <-  c("Rt_eff",
                         "Rt"
  )
  
  ## Colors
  
  cols.list <- c(
    "deepskyblue4",
    "deepskyblue1"
  )
  names(cols.list) <- names(longnames)
  color.this.var <- as.character(cols.list[var.to.plot])
  
  ### THE PLOT
  p <- ggplot(traj.CI.area)
  
  if(!is.null(ymax)){
    p <- p + facet_wrap(~state.name, labeller = labeller(state.name = longnames), scales = "fixed")
  }
  else {
    p <- p + facet_wrap(~state.name, labeller = labeller(state.name = longnames), scales = "free_y")}
  
  p <- p + geom_ribbon(data = traj.CI.area, aes_string(x = "date", ymin = "low", ymax = "up", alpha = "CI"), fill = color.this.var, show.legend = c(fill=FALSE))
  p <- p + geom_line(data = traj.CI.line, aes_string(x = "date", y = "value", linetype = "variable"), colour = color.this.var, size=1, show.legend = c(colour=FALSE))
  
  ## ADD LEGENDS
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.20, "50" = 0.50), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  
  ## ADD DATA
  if(!is.null(data.in)){
    p <- p + geom_point(data = data, aes_string(x = "date", y = "value"), size = 1, colour = "black")
  }
  
  ## FINAL THEMES AND EDITING
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%d-%b-%y")
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 12, face = "bold"))
  #  p <- p + ylab(paste0("Number  ", as.character(longnames[var.to.plot]))) + xlab(NULL)
  p <- p + ylab(NULL) + xlab(NULL)
  
  
  if(!is.null(ymax)){
    p <- p + scale_y_continuous(limits=c(0,ymax) , labels = scales::comma, breaks=seq(0,ymax,ymax/10))
  }
  else {  p <- p + scale_y_continuous(labels = scales::comma)}
  
  p
  
}



########################################################################################
## PLOT TOGETHER
########################################################################################
## USED IN:
## plot.param.t
## CFR.IFR.plots

plot.together.capacity <- function(traj.CI=traj.CI,  vars.to.plot, data.in) {
  
  ###########
  ### traj.CI
  ###########
  
  ## Select only more recent dates
  init.date <- as.Date("2021-06-01") #as.Date(lubridate::ydm(init.date))
  startDatePlot <- init.date - 15 #- date.offset.4plot #15
  endDatePlot <- startDatePlot + 365+15 #- 40  # the constant 40 because the traj are not aligned to start date
  #traj.CI <- traj.CI %>% dplyr::filter(date > startDatePlot-5) %>% dplyr::filter(date < endDatePlot)
  
  y.max.in <- max(traj.CI$up_95)
  
  ## Add title
  traj.CI$title <- "Modeled Effective R(t) in Vaccinated vs. Unvaccinated"

  
  #####################
  ### colors and names
  #####################
  
  longnames <- c("R(t) Overall",
                 "R(t) Unvax",
                 "R(t) Vax"
  )
  
  names(longnames) <-  c(
    "Rt",
    "Rt_uv",
    "Rt_v"
  )
  
  ## Colors
  
  cols.list <- c(
    ## R(t) Scenarios
    "deeppink1",
    "cornflowerblue",
    "antiquewhite4"
  )
  
  names(cols.list) <- names(longnames)
  color.this.var <- as.character(cols.list[vars.to.plot])
  
  ##################
  ### CREATE PLOT
  ##################
  
  p <- ggplot(data = traj.CI,
              aes(x = date,
                  y = mean, ymin = low_95, ymax = up_95,
                  color = state.name,
                  fill = state.name,
                  group = state.name))
  
  p <- p +  geom_ribbon(data = traj.CI,
                        aes(x = date,
                            y = mean, ymin = low_50, ymax = up_50,
                            color = state.name,
                            fill = state.name,
                            group = state.name),alpha = .5, inherit.aes=TRUE, color=FALSE)
  
  p <- p +  scale_fill_manual(values = c(color.this.var),labels = longnames) + scale_color_manual(values = c(color.this.var), labels = longnames)
  p <- p + geom_line() + geom_ribbon(alpha = 0.2, color = FALSE)
  
  if(!is.null(data.in)){
    data.processed <- reshape2::melt(data.in, measure.vars = c(2:ncol(data.in)), variable.name = "state.name")

    p <- p + geom_point(data = data.processed,
                        aes(x = date, y = value,
                            color = state.name),
                        alpha = 0.7,
                        inherit.aes = FALSE)
  }
  
  ##################
  ## FINISH PLOT
  p <- p + theme_bw() + theme(legend.title = element_blank())
  p <- p + scale_x_date(limits = as.Date(c(startDatePlot,endDatePlot)), date_breaks = "1 month" , date_labels = "%d-%b-%y")
  p <- p + scale_y_continuous(limits = c(0,y.max.in), breaks = seq(from = 0, to = y.max.in, by = y.max.in/5))
  p <- p + theme(axis.text.x = element_text(angle = 90),
                 strip.text.x = element_text(size = 12, face = "bold"))
  #  p <- p + ylab(paste0("Number  ", as.character(longnames[var.to.plot]))) + xlab(NULL)
  #p <- p + ylab("Probability") + xlab(NULL)
  p <- p + ylab("Effective R(t)") + xlab(NULL)
  #p <- p + labs(title = title.input)
  #p<-p+theme(plot.title = element_text(size = 12, hjust = 0.5, face="bold"))
  p <- p + facet_grid(. ~ title)
  
  
  p
  
  
}

