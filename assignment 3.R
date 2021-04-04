
#CS112 ASSIGNMENT 3
#Question 1

#install.packages("cobalt")
#install.packages("MatchIt")
#install.packages("Matching")
#install.packages("ggplot2")
#install.packages("gridExtra")
library(cobalt)
library(MatchIt)
library(foreign)
library(gridExtra)
library(Matching)
library(ggplot2)
lead_data
lead_data <- read.dta("https://wps.pearsoned.com/wps/media/objects/11422/11696965/data3eu/lead_mortality.dta")

#Step 2
m.out1 <- matchit(lead ~ ph + temperature + infrate , data = lead_data, method = NULL)
bal.tab(m.out1, un = TRUE, m.threshold = .1, v.threshold = 2)
balanceplot_ph <- bal.plot(m.out1, "ph", which = "unadjusted")
balanceplot_temperature <- bal.plot(m.out1, "temperature", which = "unadjusted")
grid.arrange(balanceplot_ph, balanceplot_temperature) 


#Step 3 Prima Facie Treatment Effect
treatment <- subset(lead_data, lead == 1)
control <- subset(lead_data, lead == 0)

ph_treatment_mean <- mean(treatment$ph)
ph_control_mean <- mean(control$ph)

ph_treatment_mean - ph_control_mean # equal to 0.249

temp_treatment_mean <- mean(treatment$temperature)
temp_control_mean <- mean(control$temperature)

temp_treatment_mean - temp_control_mean # equal to 0.93


#Step 5 Propensity Score Matching

#Define variables
X=cbind(lead_data$ph,lead_data$temperature)
La = cbind(lead_data$ph)
Ka = cbind(lead_data$temperature)
Y=cbind(lead_data$infrate)
Tr=cbind(lead_data$lead)

#Propensity score model
glm1 <- glm(Tr ~ X, data = lead_data, family="binomial")
summary(glm1)

#Average treatment on the treated effect
mout1 <- Match(Y=Y,Tr=Tr, X=glm1$fitted)
summary(mout1)

#Check balance
mb <- MatchBalance(Tr~X, match.out = mout1, nboots=1000, data =lead_data)
MatchBalance(lead~infrate+ph+temperature, match.out = mout1, nboots=1000, data =lead_data)

#Visualizations
m.out2 <- matchit(lead ~ ph + temperature, data = lead_data, method = "nearest")
balanceplot_ph2 <- bal.plot(m.out2, "ph", which = "both")
balanceplot_temperature2 <- bal.plot(m.out2, "temperature", which = "both")
grid.arrange(balanceplot_ph2, balanceplot_temperature2) 


#Step 7 Rosenbaum's Sensitivity Method
library(rbounds)
psens(mout1, Gamma=5, GammaInc=.05)


#Step 8 Genetic Matching

#Define covariates and treatment variables, then run the GenMatch
gen1 <- GenMatch(Tr=Tr, X=X, M=1, estimand="ATT", pop.size=1000, max.generation=1000, wait.generation=50, caliper=0.1)

#Use weights from GenMatch to match dataset using Match
mgen1 <- Match(Tr=Tr, X=X, M=1, Y=Y, estimand="ATT", Weight.matrix=gen1, caliper=0.1)
summary(mgen1)

#Get results of the Match using MatchBalance
mgenbalance <- MatchBalance(lead ~ temperature + ph, data = lead_data, match.out = mgen1, nboots=10)
mgenbalance

#Data visualizations
m.out3 <- matchit(lead ~ ph + temperature, data = lead_data, method = "genetic")
balanceplot_ph3 <- bal.plot(m.out3, "ph", which = "both")
balanceplot_temperature3 <- bal.plot(m.out3, "temperature", which = "both")
grid.arrange(balanceplot_ph3, balanceplot_temperature3) 

#Step 9 Genetic Matching -> Changes in caliper
#Narrow caliper on temperature and wide caliper on  pH
genetic_match_caliper <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=30, max.generation=10, wait.generation=5, caliper=c(0.01, 100000))
matched_genetic_match_narrow <- Match(Tr=Tr, X=X, estimand="ATT", M=1, Weight.matrix=genetic_match_caliper, caliper=c(0.01, 100000))
summary(matched_genetic_match_caliper)
mb_matched_genetic_match_caliper <- MatchBalance(lead ~ temperature + ph, data = lead_data, match.out = matched_genetic_match_narrow, nboots=100)

#Step 10 Exact option vs caliper 
genetic_match_exact <- GenMatch(Tr=Tr, X=X, estimand="ATT", M=1, pop.size=30, max.generation=10, wait.generation=5)
matched_genetic_match_exact <- Match(Tr=Tr, X=X, estimand="ATT", M=1, Weight.matrix=genetic_match, exact = 1)
summary(matched_genetic_match_exact)
mb_matched_genetic_match_exact <- MatchBalance(lead ~ temperature + ph, data = lead_data, match.out = matched_genetic_match_exact, nboots=100)

#Step 12 Rosenbaum's Sensitivity Method
psens(mgen1, Gamma=2.5, GammaInc=0.05)


#########################


#Question 2 

library(Matching)
demo(GerberGreenImai)




