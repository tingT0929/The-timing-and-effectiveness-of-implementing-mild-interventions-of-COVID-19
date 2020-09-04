# The-timing-and-effectiveness-of-implementing-mild-interventions-of-COVID-19
The timing and effectiveness of implementing mild interventions of COVID-19 in large industrial regions via a synthetic control method

## 1. Data 
### 1) Abstract
Our data source consists of the population, latitude, population density, GDP per capita, the cumulative confirmed cases and deaths of COVID-19, the confirmed cases per 100,000 for 68 counties in the United States from March 1 to March 29, 2020 and for Shenzhen from January 19 to Feburary 16, 2020. These data were collected from the offical websites. This enabled us to evaluate the treatment effects of interventions implemented in Shenzhen by a synthetic control method and delay effects of those interventions by a proposed SIHR (Susceptible, Infectious, Hospitalized, Removed) model.   

### 2) Availability
The data to reproduce our results are available.

### 3) Description
The data incorporte 1 `.csv` file and 1 `.rda` file.
- The cumulative confirmed cases (or confirmed cases per 100,000) of COVID-19, the population, latitude, population density, and GDP per capita for 68 counties of the United States and Shenzhen were collected in the `.csv` file (`szcpr-0710.csv`)
- The cumulative confirmed cases of COVID-19 in Shenzhen from January 19 to February 29 and the initial parameters from the proposed SIHR model were complied in the  `.rda` file (`results_sz.rda`)

### 4) Permissions
The data were orignially collected by the authors.

----
## 2. Code
### 1) Abstract
The codes incorported all the scripts to reproduce the analysis in the paper. 

### 2) Reporducibility
- The estimation of Susceptible (S), Infectious (I), Hospitalized(H), and Removed (R) individuals for the dynamic system (A1) in the paper by runing `Model.R`.
- The simulation of delay effects by runing `Simulation.R`.
- The evaluation of treatment effects of interventions implemented in Shenzhen by running `Synth_PCA0831.R`.

----
## 3. Paper

- https://www.medrxiv.org/content/10.1101/2020.05.18.20105

