
<br>
<br> 

# Vaccination model-based scenario analysis for LA County 

### Requested by the Los Angeles County Department of Public Health (LACDPH)

## Research questions:

- What might infection, hospitalization, and death trend lines have looked like if nobody had been vaccinated?

- What might future infection, hospitalization, and death trend lines look like if all eligible persons were now vaccinated?

- How do these scenarios compare to observed infections, hospitalizations, and deaths between June 1 – Sept 1, 2021?

<br>

### Time period of analysis:

* June 1, 2021 – June 1, 2022

<br>

## Work by: 

<br>
Abigail Horn

Research Associate, Divisions of Biostatistics and Health Behavior

<br>
Dave Conti

Professor of Biostatistics

<br>
Department of Population and Public Health Sciences

University of Southern California

<br>

## Scenarios and parameters evaluated

- Final presentation available at **[this link](https://docs.google.com/presentation/d/1qhLLlo9a5K3Huhsd3gwkPwudSSIKHuTL/edit?usp=sharing&ouid=114012102276366140518&rtpof=true&sd=true)**.

- Scenarios, parameter values, and sources found at **[this link](https://docs.google.com/spreadsheets/d/1PBXMgmchEnlYTxcWHqZR4cuUcvKFO5LuW4G3CbSLXOo/edit#gid=490110439)**.

- Differences in infection, hospitalization, ICU admission, and death rates between vaccinated and unvaccinated populations in LA County from [LACDPH study published in MMWR](https://www.cdc.gov/mmwr/volumes/70/wr/mm7034e5.htm#contribAff).

- Values of $R_{t,eff}$ for LA County overall based on **(CovidEstim)[https://covidestim.org/us/CA/06037]** estimates. 

<br>

## Guide to using code

- Use `code/LACDPH_Sept21/wrapper_scenarios_9.2.21.R` to run the modeling scenarios and extract output for summarizing tables

- Use `code/LACDPH_Sept21/plot_scenario_facetgrid.R` to plot multiple scenarios together. Output may need to be customized to particular scenarios and compartmental variables chosen for display.