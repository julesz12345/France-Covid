## Lab 2 - Preprocessing

#Import data
urlfile<-'https://raw.githubusercontent.com/julesz12345/France-Covid/main/countries-aggregated_csv.csv'
countries<-read.csv(urlfile)

#define the sum of active cases  (active = total - fatal - recovered) 
countries$Active_Cases=countries$Total-countries$Deaths-countries$Recovered

#add the day column 

#reshape
library(reshape2)
library(tidyverse)
names(countries)[1] <- "Date"
countries$Date=as.Date(countries$Date) #set date format
countries_m =melt(countries , id=c("Date","Country"))

#France Observations

#France
library(ggplot2)
france_m<-subset(countries_m, Country == "France")
p1=ggplot(france_m, aes(x=Date, y=value, group=variable, color=variable))
p1+ geom_line()+theme_bw()

#Data - Unmelted
france<-subset(countries, Country == "France")

#france population
N_France = 65.273*1000000
france$Population=ifelse(france$Country=="France",N_France,0)

#renaming columns
names(france)[names(france) == "ï..Date"] <- "date"
names(france)[names(france) == "Total"] <- "cases_total"
names(france)[names(france) == "Deaths"] <- "cases_fatal"
names(france)[names(france) == "Recovered"] <- "cases_recovered"
names(france)[names(france) == "Active_Cases"] <- "cases_active"
names(france)[names(france) == "Population"] <- "N"
names(france)[names(france) == "Day"] <- "day"

################################################################################

# Week 2

# Transform epidemiological data into SIR-data:
# - susceptible=N-total_cases
# - infected=cases_active 
# - removed=cases_fatal+cases_recovered

france$susceptible=france$N-france$cases_total 
france$infected=france$cases_active 
france$removed=france$cases_recovered+france$cases_fatal

################################################################################

#Just sample by data of first - second and third wave
france1 <- france[france$date >= '2020-03-01' & france$date <= '2020-06-01',]
france2 <- france[france$date >= '2020-09-26' & france$date <= '2020-12-01',]
france3 <- france[france$date >= '2021-03-18' & france$date <= '2021-05-15',]

#Peak of waves
# 1: april 6 2020
# 2: november 15 2020
# 3: april 8 2021

#add the day column
france1$day <- (1:nrow(france1))
france2$day <- (1:nrow(france2))
france3$day <- (1:nrow(france3))


#To see long term behavior of endemic, we simulate
#long-term behavior with estimated parameters

##R-LAB: STEP 1 - SET UP DIFFERENTIAL EQUATIONS

#######################################
#Pandemic Behavior 3 years
#######################################

library(deSolve)
sirmod = function(t, y, parms) {
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  alpha=parms["alpha"]
  beta = parms["beta"]
  mu = parms["mu"]
  gamma = parms["gamma"]
  N = parms["N"]
  
  ###DEFINE EQUATIONS
  dS = alpha - mu*S - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  dR = gamma * I - mu * R
  res = c(dS, dI, dR)
  
  # Return list of gradients
  list(res)
}

#Predicted behavior of the pandemic in the next 3 years
#Assuming parameters of the first wave

#Let's look at a evolution of epidemic over 1000 days
times = seq(0, 1000, by = 1/10)
parms = c(alpha=1/80, mu = 1/100, N = 1, beta =0.42 , gamma = 0.13)

#We are assuming a normalized population of N=1(we can multiply final results by France's real population: 65.273*1000000)
#birth rate of alpha=1/80
#death rate of mu=1/100

#Values of the first wave
#interaction rate: beta=0.42
#removal rate of gamma=0.13

start = c(S = 0.999, I = 0.001, R = 0)
#We assume that at beginning, 0.1% of population was infected.

#R-LAB: STEP 2 - SIMULATE BEHAVIOR FOR SIR MODEL
out = ode(y = start, times = times, func = sirmod, 
          parms = parms)
out=as.data.frame(out)
head(round(out, 3))

#Plot of long-term behavior
library(ggplot2)
ggplot(out, aes(x=time))+geom_line(aes(y=I),col="green")+geom_line(aes(y=S),col="red")+geom_line(aes(y=R),col="blue")
#Green - Infections
#Red - Susceptible
#Blue - Removed

#Green - Infections
ggplot(out, aes(x=time))+geom_line(aes(y=I),col="green", size=.5)+theme_bw()

###########
#Predicted first peak wave = 03/03/2020
#Predicted second peak wave = 05/26/2020
#Predicted third peak wave = 09/28/2020
###########

#Let's calculate long-term R0 and steady-state proportion
#of (Susceptible, Infected, Recovered}-people in a population
alpha=parms[1]
mu=parms[2]
beta=parms[4]
gamma=parms[5]
R0=as.numeric(alpha*beta/(mu*(mu+gamma)))
S_eq=as.numeric((gamma+mu)/beta)
I_eq=as.numeric((mu/beta)*(R0-1))
R_eq=as.numeric((gamma/beta)*(R0-1))

R0 #3.75
S_eq #.33
I_eq #.06
R_eq #.85

"""In the long-run, we will reach an equilibrium where
0.06 of population is currently infected
0.85 had the infection (and either recovered or died from it)
0.33 never got infection (are susceptible)"""


###R-LAB: SIR-DEMOGRAPHY MODEL

#Let's estimate effect of social distancing and hygiene by
#changing previous estimates, with beta= 0.15 where this coefficient is the beta of the second wave (instead of 0.42 - first wave)

times = seq(0, 1000, by = 1/10)
parms = c(alpha=1/80, mu = 1/100, N = 1, beta =0.15 , gamma = 0.13)
start = c(S = 0.999, I = 0.001, R = 0)
out2 = ode(y = start, times = times, func = sirmod, 
            parms = parms)
out2=as.data.frame(out2)
head(round(out2, 3))
out$I2=out2$I

ggplot(out, aes(x=time))+geom_line(aes(y=I), group = 1, col="red")+geom_line(aes(y=I2), group = 1,col="blue")+theme_bw()

#red no measures
#blue with measures

###########
#Predicted peak first wave with measures = 27/07/2020
###########

###With social distancing, we delayed second wave, and made each wave much smaller. 
###Also, in long run, endemic equilibrium of infected is reduced 
###If we remove social distancing midway, blue curve would converge towards green curve


#######################################
#Immigration
#######################################

#Let's estimate model assuming that 1% of immigrants
#come infected (e.g.,1% of travelers arriving to France
#               bring covid-19)
#In this case, k=0.99

sirmod = function(t, y, parms) {
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  alpha=parms["alpha"]
  beta = parms["beta"]
  mu = parms["mu"]
  gamma = parms["gamma"]
  N = parms["N"]
  ###DEFINE EQUATIONS
  dS = alpha*0.99 - mu*S - beta * S * I/N
  dI = alpha*0.01+beta * S * I/N - (mu + gamma) * I
  dR = gamma * I - mu * R
  res = c(dS, dI, dR)
  # Return list of gradients
  list(res)
}

times = seq(0, 1000, by = 1/10)
parms = c(alpha=1/80, mu = 1/100, N = 1, beta =0.42 , gamma = 0.13)
start = c(S = 0.999, I = 0.001, R = 0)
out3 = ode(y = start, times = times, func = sirmod, 
           parms = parms)
out3=as.data.frame(out3)
head(round(out3, 3))
out$I3=out3$I
ggplot(out, aes(x=time))+geom_line(aes(y=I),col="green")+geom_line(aes(y=I3),col="blue")+theme_bw()

#In this modified model, first wave is larger, but subsequent waves are attenuated
#Conclusion, adding infected traveller make initial spike larger but subsequent spikes smaller

#####
#SIRS models capture dynamics of diseases where reinfection is possible
#Xi is the rate at which people lose their immunity.
#1/Xi=time immunity lasts

library(deSolve)
sirmod = function(t, y, parms) {
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  # Pull parameter values from parms vector
  beta = parms["beta"]
  gamma = parms["gamma"]
  xi = parms["xi"]
  N = parms["N"]
  
  # Define equations
  dS = - beta * S * I/N + xi*R
  dI = beta * S * I/N - ( gamma) * I
  dR = gamma * I - xi*R
  res = c(dS, dI, dR)
  # Return list of gradients
  list(res)
}

#Assuming 200 days of immunity ~ 6 months

times = seq(0, 1000, by = 1/10)
parms = c( N = 1, beta =0.42 , gamma = 0.13, xi=0.005)
start = c(S = 0.999, I = 0.001, R = 0)
out = ode(y = start, times = times, func = sirmod,  parms = parms)
out=as.data.frame(out)
head(round(out, 3))

ggplot(out, aes(x=time))+geom_line(aes(y=I),col="green")+theme_bw()

immunity=200 #days
times = seq(0, 1000, by = 1/10)
parms = c( N = 1, beta =0.42 , gamma = 0.13, xi=1/immunity)
start = c(S = 0.999, I = 0.001, R = 0)
out = ode(y = start, times = times, func = sirmod,  parms = parms)
out=as.data.frame(out)
head(round(out, 3))

ggplot(out, aes(x=time))+geom_line(aes(y=I),col="green")+theme_bw()

#Endemic equilibrium
#Let's estimate steady-state equilibrium under both cases

N=1
####xi=1/0.01111 = 90 days immunity ~ 3 months
xi=0.01111
I_eq=as.numeric((N-gamma/beta)/(1+gamma/xi))
I_eq #output = 0.054

#With 90 days of immunity, 5.4% of the
#people would be infected in long run
#endemic equilibrium

####xi=1/0.005 = 200 days immunity
xi=0.005
I_eq=as.numeric((N-gamma/beta)/(1+gamma/xi))
I_eq #output = 0.255

#With 200 days of immunity, 2.5% would
#be infected in the endemic equilibrium

##################Vaccines

library(deSolve)
sirmod = function(t, y, parms) {
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  # Pull parameter values from params vector
  alpha=parms["alpha"]
  beta = parms["beta"]
  mu = parms["mu"]
  gamma = parms["gamma"]
  rho = parms["rho"]
  N = parms["N"]
# Define equations
dS = alpha - mu*S - beta * S * I/N-rho*S
dI = beta * S * I/N - (mu + gamma) * I
dR = gamma * I - mu * R+rho*S
res = c(dS, dI, dR)
# Return list of gradients
list(res)
}

times = seq(0, 1000, by = 1/10)
parms = c(alpha=1/80, mu = 1/100, N = 1, beta =0.17 , gamma = 0.13, rho=0.01*0.95)
start = c(S = 0.999, I = 0.001, R = 0)
out = ode(y = start, times = times, func = sirmod,
           parms = parms)
out=as.data.frame(out)
head(round(out, 3))
ggplot(out, aes(x=time))+geom_line(aes(y=I),col="green")+theme_bw()

