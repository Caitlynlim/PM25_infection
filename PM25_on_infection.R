#Preparation
rm(list = ls())
setwd("/Users/CaitlynLim/Documents/EC4383/Assignment 2")
library(AER)
library(stargazer)
library(data.table)
library(dplyr)

# load data set (written in Stata's format)
library(rio)                                         # load required package
library(haven)
endata <- data.frame(as_factor(read_dta("EnvironmentalHealthData.dta"))) # read data set
endata <- spread_attrs(droplevels(gather_attrs(endata))) # cleaning up
endata2 <- read.csv("EnvironmentalHealthData.csv")
data("Endata", package = "AER")

sum(is.na(endata)) #No NAs

#Calculating baseline mean
as.numeric(endata$pm25)
mean(endata$pm25[159:167])#Mean over 2015 Haze
r1 <- mean(endata$resp_infect[3:10])  #Mean respiratory over 2012
r2 <- mean(endata$resp_infect[55:62]) #Mean respiratory over 2013
r3 <- mean(endata$resp_infect[107:114]) #Mean respiratory over 2014
mean((r1+r2+r3)/3) #Mean respiratory over non-2015 haze constituting as baseline

p1 <- mean(endata$pm25[3:10])  #Mean PM25 over 2012
p2 <- mean(endata$pm25[55:62]) #Mean PM25 over 2013
p3 <- mean(endata$pm25[107:114]) #Mean PM25 over 2014
mean((p1+p2+p3)/3) #Mean PM25 over non-2015 haze constituting as baseline

mean(endata$pm25) #Mean over total period 

#Number of respiratory infection vs all 4 illnesses
sum(endata$resp_infect)
mean(endata$resp_infect)
mean(endata$four_categ)
diffFourRes <- mean(endata$four_categ) - mean(endata$resp_infect)
(mean(endata$resp_infect)/mean(endata$four_categ))*100
endata$other_dies <- (endata$four_categ - endata$resp_infect)
cor(endata$other_dies, endata$resp_infect)


#Listing my variables
y1 <- cbind(endata$resp_infect) #dependent variable
y2 <- cbind(endata$pm25) #endogeneous variable 
x1 <- cbind(endata$schoolHol, endata$month, endata$year, endata$tp, endata$ws, endata$pp,
            endata$dewPointDep, endata$hm) #exogeneous variables 
z <- cbind(endata$pos_tpGrad_0to78m) #instrument 
zalt <- cbind(endata$idwddw_frp, endata$pos_tpGrad_0to78m) #alternative instruments

#OLS Regression
ols <- lm(y1 ~ y2 + x1 + zalt + endata$month:endata$tp + endata$month:endata$pp + endata$schoolHol:endata$month)
summary(ols, vcov = sandwich, df = Inf, diagnostics = TRUE) 
instrFtest <- waldtest(ols,.~.-zalt)
print(instrFtest)

#IV F-stats test
first_stage <- lm(endata$pm25~endata$pm25 + endata$schoolHol + endata$month + endata$year + endata$tp +
                    endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
                    endata$month:endata$pp + endata$schoolHol:endata$month 
                  + endata$idwddw_frp + endata$pos_tpGrad_0to78m,
                  data=endata)
instrFtest <- waldtest(first_stage,.~.-endata$idwddw_frp-endata$pos_tpGrad_0to78m)
print(instrFtest)

#(MAIN REGRESSION) IV regression with interaction
iv2 <- ivreg(endata$resp_infect ~ endata$pm25  + endata$schoolHol + endata$month + endata$year + endata$tp + 
               endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
               endata$month:endata$pp + endata$schoolHol:endata$month | .-endata$pm25 + endata$idwddw_frp + 
               endata$pos_tpGrad_0to78m)
summary(iv2, vcov = sandwich, df = Inf, diagnostics = TRUE)

#IV regression quadratic - found insignificant 
IVPM25 <- endata$pm25
IVPM25sq <- endata$pm25^2

iv5 <- ivreg(endata$resp_infect ~ IVPM25  + IVPM25sq + endata$schoolHol + endata$month + endata$year + endata$tp + 
               endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
               endata$month:endata$pp + endata$schoolHol:endata$month | .-IVPM25 + endata$idwddw_frp + 
               endata$pos_tpGrad_0to78m)
summary(iv5, vcov = sandwich, df = Inf, diagnostics = TRUE)

#IV regression without interaction
iv1 <- ivreg(y1 ~ y2 + x1 | x1 + zalt)
summary(iv1, vcov = sandwich, df = Inf, diagnostics = TRUE)

#IV regression with trend as covariate 
iv4 <- ivreg(endata$resp_infect ~ endata$pm25  + endata$schoolHol + endata$month + endata$year + endata$tp + 
               endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
               endata$month:endata$pp + endata$schoolHol:endata$month | .-endata$pm25 + endata$idwddw_frp + 
               endata$pos_tpGrad_0to78m)
summary(iv5, vcov = sandwich, df = Inf, diagnostics = TRUE)

#IV regression dropping adjusted radiative power
iv3 <- ivreg(endata$resp_infect ~ endata$pm25 + endata$schoolHol + endata$month + endata$year + endata$tp + 
               endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
               endata$month:endata$pp + endata$schoolHol:endata$month | .-endata$pm25 + 
               endata$pos_tpGrad_0to78m)
summary(iv3, vcov = sandwich, df = Inf, diagnostics = TRUE)

#demonstrate data in table form
stargazer(ols, 
          type = "text", 
          dep.var.labels = "Number of visits to the clinic for respiratory related illness", 
          covariate.labels = c("PM2.5", "School Holidays", "Month", "Year", "Temperature",
                               "Wind speed", "Precipitation", "Dew Point", "Humidity",
                               "Direction-Bearing weighted FRP", "Atmospheric Inversion", 
                               "Month:Temperature", "Month:Precipitation", "Month:SchoolHolidays"), 
          digits = 2
)

stargazer(iv2, 
          type = "text", 
          dep.var.labels = "Number of visits to the clinic for respiratory related illness", 
          covariate.labels = c("PM2.5", "School Holidays", "Month", "Year", "Temperature",
                               "Wind speed", "Precipitation", "Dew Point", "Humidity",
                               "Month:Temperature", "Month:Precipitation", "Month:SchoolHolidays"), 
          digits = 2
)

#Robustness Check

#test with lagged respiratory infection found insignificant
iv4 <- ivreg(endata2$lagResp ~ endata$pm25 + endata$schoolHol + endata$month + endata$year + endata$tp +
               endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
               endata$month:endata$pp + endata$schoolHol:endata$month | .-endata$pm25 + endata$idwddw_frp + 
               endata$pos_tpGrad_0to78m, is.na = TRUE)
summary(iv4, vcov = sandwich, df = Inf, diagnostics = TRUE)

#test with placebo Instruments
first_stage <- lm(endata$pm25~endata$pm25 + endata$schoolHol + endata$month + endata$year + endata$tp +
                    endata$ws + endata$pp + endata$dewPointDep + endata$hm + endata$month:endata$tp + 
                    endata$month:endata$pp + endata$schoolHol:endata$month 
                  + endata2$frpPlacebo + endata2$posPlacebo,
                  data=endata)
instrFtest <- waldtest(first_stage,.~.-endata2$frpPlacebo-endata2$posPlacebo)
print(instrFtest)







#summary(lm(endata$resp_infect ~ endata$pm25 + endata$month +endata$schoolHol +endata$ws
#+endata$tpmax +endata$hm+ endata$pp))

#hazeJO <-subset(endata, endata$trend<=167 & endata$trend>=150)
