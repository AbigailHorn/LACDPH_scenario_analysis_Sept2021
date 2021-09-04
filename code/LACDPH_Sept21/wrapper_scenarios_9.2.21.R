
###################################################################################################
## Packages

library(reshape2)
library(tidyverse)
library(ggplot2)
library(plotly)
library(ggrepel)
library(bindata)
library(odin)
library(fitR)
library(knitr)
library(EasyABC)
library(gridExtra)
library(odin)
library(lubridate)
library(EasyABC)
library(gridExtra)
library(kableExtra)
library(plyr)
library(dplyr)
library(data.table)
library(scales)
library(EasyABC)
library(patchwork)

library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(network)
library(tidygraph)
library(ggraph)
library(visNetwork)
library(networkD3)
library(ggmosaic)
library(formattable)
library(DT)
library(reshape)
library(here)
library(fs)

library(MASS)
library(plotly)

lang_output <- function(x, lang) {
  cat(c(sprintf("```%s", lang), x, "```"), sep = "\n")
}
r_output <- function(x) lang_output(x, "r")

knitr::opts_chunk$set(
  fig.width = 9.5,
  fig.height = 8,
  eval=TRUE,
  echo=FALSE,
  warning=FALSE,
  cache=FALSE,
  message=FALSE,
  include=TRUE
)

code.LACDPH=here("code/LACDPH_Sept21/")
code.dir=here("code/")
data.dir=here("data/LACDPH_Sept21")
data.update.dir=here("data/risk_update")
result.dir = here("results/")
fig.dir = here("figs/")
output.dir = here("output/")
output.density.dir = here("output/epi_model_output/density")
code.paper.dir=here("code/epi_model_code")
code.risk.dir=here("code/risk_model_code/2021")
code.scenarios.dir=here("code/scenarios_code/")

###################################################################################################
## COVID INPUT DATA
## latest_data: cumulative and daily counts for "Htotcum","D","Vcum","Idetectcum","H_new","D_new"
latimes_path <- path(code.paper.dir, "LAtimes.R")
source(latimes_path)
la_data <- latimes_readin()

###################################################################################################
## LOAD EPIDEMIC MODEL
## And compile the model
path_seihqdr_model <- path(code.paper.dir, "stochastic_SEIAHQDR_Alphat_rt.R")
seihqdr_generator <- odin::odin(path_seihqdr_model)

###################################################################################################
## LOAD MODEL AND PLOTTING FUNCTIONS
model_code <- path(code.LACDPH, "model_plot_LACDPH.R")
source(model_code)

####################################################################################
## Global parameters
start_date = as.Date("2021-06-01")
start_time = 0

####################################################################################
## Get simulation data
replicates = 100
iter = 10
time.steps=500
scenario.IDs = c("NooneVax_Isolate.00", "NooneVax_Isolate.10", "NooneVax_Isolate.20", "EveryoneVax_Efficacy.850", "EveryoneVax_Efficacy.725", "EveryoneVax_Efficacy.600")
#scenario.IDs = "EveryoneVax_Efficacy.850"

traj.0 = model.output.to.plot.SIM(scenario.IDs, replicates=replicates,iter=iter,time.steps=time.steps)


####################################################################################
## PLOTTING
data.in <- la_data
time.steps.plot = as.numeric(as.Date("2022-03-31") - (start_date))
vars.to.plot <- c(#"S",
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
                   "D_new",
                   "D"
                   #"R"
)

#traj.plot = traj.0 %>% filter(scenario=="Scenario_NooneVac_Efficacy.85") #filter(scenario=="Scenario_NooneVac_Efficacy.85")

traj.CI=traj.0
plot.model.data.all(traj.CI=traj.CI, data.in=data.in, time.steps.plot=time.steps.plot, vars.to.plot=vars.to.plot)


########################################################################################
## GETTING MODEL SUMMARY STATISTICS FUNCTION (FOR TABLES)
########################################################################################

scenario.IDs <- c("EveryoneVax_Efficacy.850", "NooneVax_Isolate.10")
traj.raw <- correlated.param.SIM(scenario.IDs, replicates=replicates,iter=iter,time.steps=366)
traj.0 = model.output.to.plot.SIM(scenario.IDs, replicates=replicates,iter=iter,time.steps=time.steps)
traj.0$scenario <- as.factor(traj.0$scenario)

######## Get Peaks
traj.peak <- traj.0 %>% dplyr::filter(state.name %in% c("I","Htot","Q","D_new")) %>% dplyr::filter(date > "2021-08-01") %>% 
  dplyr::group_by(scenario,state.name) %>% 
  dplyr::summarise(
    peak_median = max(median, na.rm=TRUE),
    date_peak = match(peak_median,median),
    peak_low95 = low_95[date_peak],
    peak_up95 = up_95[date_peak]
  ) %>% as.data.frame()

data <- la_data %>% filter(date>"2021-06-01")
data <- reshape2::melt(data, measure.vars = c(2:ncol(data)), variable.name = "state.name")
dataCum = data %>% filter(state.name %in% c("D", "Idetectcum")) %>% group_by(state.name) %>% mutate(value = value - dplyr::first(value)) %>% as.data.frame()
dataNew = data %>% filter(state.name %in% c("D", "Idetectcum") == FALSE) %>% as.data.frame()
data=as.data.frame(rbind(dataCum,dataNew))
data %>% 
  dplyr::group_by(state.name) %>% 
  dplyr::summarise(
    peak = max(value, na.rm=TRUE)
  ) %>% as.data.frame()

######## Get Cumulative Counts

traj.cum <- traj.0 %>% dplyr::filter(date %in% as.Date(c("2021-09-01","2022-06-01"))) %>% 
  dplyr::filter(state.name %in% c("Idetectcum","Htotcum","Qcum","D")) %>% 
  dplyr::filter(scenario=="NooneVax_Isolate.10") %>% 
  select(c(date,state.name,median,low_95,up_95))
data %>% 
  dplyr::filter(date %in% as.Date(c("2021-09-01","2022-06-01"))) %>% 
  dplyr::filter(state.name %in% c("Idetectcum","Htotcum","Qcum","D"))


