
#--packages required
library(doBy)
library(sqldf)
library(xtable)
library(plyr)
library(dplyr)
library(ggplot2)
library(lattice)
library(MASS)
library(car)
library(gridExtra)
library(lme4)
library(vioplot)
library(broman) 
library(gamlss)
library(lubridate)

# -- Master tables

#############
#data sets. Pat=description of patients, PSA: events, Arzt: arzt description, Chronic diseases
data=read.csv("f0008_01_pat.csv", stringsAsFactors=FALSE) 
data2=read.csv("f0008_02_psa.csv", stringsAsFactors=FALSE) 
data3=read.csv("f0008_03_arzt.csv", stringsAsFactors=FALSE) 
cronic_d=read.csv("f0008_04_chronic_disease.csv", stringsAsFactors=FALSE)
  
data<-orderBy(~pat_id+study_start_dt, data)
data2<-orderBy(~pat_id, data2)
data_all<-merge(data, data2, all=TRUE)

#db eligible but it is not enough we should check age and >= 90 days of follow up
data_el<-subset(data_all, elig_vect_1=="1-1-1-1-1-1-1")

#data cleaning for missing
data_el$arzt_region[data_el$arzt_region=="\\N"] <-""
data_el$arzt_region_code[data_el$arzt_region_code=="\\N"] <-""
data_el$arzt_doby[data_el$arzt_doby=="\\N"] <-""
data_el$arzt_doby<-as.numeric(as.character(data_el$arzt_doby))
data_el$employ_status[data_el$employ_status=="Selbständig "]<-"Selbständig"
data_el$event_prostcancer_dt[data_el$event_prostcancer_dt=="\\N"] <-""
data_el$event_urinarytract_dt[data_el$event_urinarytract_dt=="\\N"] <-""
data_el$event_prosthyper_dt[data_el$event_prosthyper_dt=="\\N"] <-""
data_el$event_atc_g04_dt[data_el$event_atc_g04_dt=="\\N"] <-""

data_el$arzt_region_code[data_el$arzt_region_code=="AGR"|data_el$arzt_region_code=="RE"|data_el$arzt_region_code=="MIX"]<-"OTHER"
data_el$arzt_id<-as.character(data_el$arzt_id)

#first disease event
data_el$event_first_disease_dt<-pmin(as.Date(data_el$event_prostcancer_dt, format="%Y-%m-%d"), as.Date(data_el$event_urinarytract_dt, format="%Y-%m-%d"), as.Date(data_el$event_prosthyper_dt, format="%Y-%m-%d"), as.Date(data_el$event_atc_g04_dt, format="%Y-%m-%d"), na.rm=TRUE)

#mortality rate (for simulation of death cases)
#we start with ch male demography (years 2009 and 2010)
mortality=read.csv("f0008_ch_demography.csv") 
mortality1<-subset(mortality, survey_year==2009 & sex2=="m")
mortality2<-subset(mortality, survey_year==2010 & sex2=="m")
mort_2009_m <- reshape(mortality1, direction = "wide", idvar="age_y", timevar="fact2", drop=c("canton", "survey_year", "sex2"))
mort_2010_m <- reshape(mortality2, direction = "wide", idvar="age_y", timevar="fact2", drop=c("canton", "survey_year", "sex2"))
mort_2009_m$year<-2009
mort_2010_m$year<-2010
#we put together year 2009 and 2010
mort_0910<-join(mort_2009_m, mort_2010_m, type="full")
mort_0910$rate<-mort_0910$n.death_cases/mort_0910$n.pop_middle
#we extract age between [50,80] and we use year 2010  
mort_0910<-subset(mort_0910, age_y>=50 & age_y <=80 & year==2010)

#Parametric estimation of hazard rates (Gompertz distribution y=b*exp(cx) where x is age). 
#The hazard function is increasing from b at time zero to infinite at infinite time.
#In the Gompertz distribution the log of the hazard is linear in x.
#This distribution  provides  a  remarkably  close  fit  to  adult  mortality  in contemporary developed countries.

#We use NLS function (algorithm port, starting values for b=0.2 and c=0.1)
fm4DNase1 <- nls(rate ~ b* exp(c*age_y),
                 data=mort_0910,
                 start = list(b = 0.2, c = 0.1),
                 algorithm = "port")

#to see estimated coefficients, std and p values 
summary(fm4DNase1)

#verifying fitting graphically
plot(mort_0910$age_y, mort_0910$rate)
lines(mort_0910$age_y, predict(fm4DNase1), col = 2)

#we generate time of death t on the basis of Gompertz distribution, using the relationship with uniform distribution. 
#In fact if we generate U from Uniform distribution we have that t=log(1-log(U)/b)/c is Gompertz distributed
set.seed(50)
y2 <- log(1-log(runif(12722))/summary(fm4DNase1)[10]$coefficients[1,1])/summary(fm4DNase1)[10]$coefficients[2,1]

#simulating death for cases where we don't have a consultation in the last trimester and where cases are not diseased
data_sim<-subset(data_el, as.Date(last_kons_dt)<"2016-10-01" & is.na(event_first_disease_dt))
data_sim<-orderBy(~pat_id, data_sim)
dead_sim<-data.frame(pat_id=unique(data_sim$pat_id), y2, age_dc=tapply(data_sim$pat_age_end_2, data_sim$pat_id, mean))
dead_sim$dead<-as.numeric(dead_sim$y2<=75 & dead_sim$y2>=55 & dead_sim$age_dc<=y2)                     

data_el<-merge(dead_sim, data_el, all=TRUE)

#we adjust follow up on the basis of simulated results
data_el$fup_end_dt_1<-ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), data_el$fup_end_dt_1)
data_el$fup_end_dt_2<-ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), data_el$fup_end_dt_2)
data_el$fup_end_dt_3<-ifelse(data_el$dead==1, data_el$last_kons_dt, data_el$study_end_dt)
data_el$fup_end_dt_3<-ifelse(is.na(data_el$dead), ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), data_el$study_end_dt), data_el$fup_end_dt_3)

data_el$fup_y_1<-(as.Date(data_el$fup_end_dt_1)-as.Date(data_el$fup_start_dt_1))/365
data_el$fup_y_2<-(as.Date(data_el$fup_end_dt_2)-as.Date(data_el$fup_start_dt_2))/365
data_el$fup_y_3<-(as.Date(data_el$fup_end_dt_3)-as.Date(data_el$fup_start_dt_3))/365

data_el$pat_age_end_1<-data_el$pat_age_start_1+floor(as.numeric(data_el$fup_y_1))
data_el$pat_age_end_2<-data_el$pat_age_start_2+floor(as.numeric(data_el$fup_y_2))
data_el$pat_age_end_3<-data_el$pat_age_start_3+floor(as.numeric(data_el$fup_y_3))

#eligibility: age (start of end of study) between [55,75] and days of follow up > 90
data_el$elig_2_1<-as.numeric((data_el$pat_age_start_1>=55 & data_el$pat_age_start_1<=75 & as.numeric(data_el$fup_y_1)>=90/365) | (data_el$pat_age_end_1>=55 & data_el$pat_age_end_1<=75 & as.numeric(data_el$fup_y_1)>=90/365)) 
data_el$elig_2_2<-as.numeric((data_el$pat_age_start_2>=55 & data_el$pat_age_start_2<=75 & as.numeric(data_el$fup_y_2)>=90/365) | (data_el$pat_age_end_2>=55 & data_el$pat_age_end_2<=75 & as.numeric(data_el$fup_y_2)>=90/365)) 
data_el$elig_2_3<-as.numeric((data_el$pat_age_start_3>=55 & data_el$pat_age_start_3<=75 & as.numeric(data_el$fup_y_3)>=90/365) | (data_el$pat_age_end_3>=55 & data_el$pat_age_end_3<=75 & as.numeric(data_el$fup_y_3)>=90/365)) 

#chronic disease dataset
chr<-subset(cronic_d, select=c(pat_id, cd_count, effective_dt, expiry_dt, arzt_id))
chr$arzt_id<-as.character(chr$arzt_id)

chr1<-subset(chr, pat_id %in% el2$pat_id)
chr2<-subset(chr, pat_id %in% n_psa_id$pat_id)
chr2$dt<-cut(as.Date(chr2$effective_dt), as.Date(quarter_starts))
chr2<-subset(chr2, select=c("pat_id", "dt", "cd_count"))

#export (file to be updated)
write.csv(data_el, file="f0008_data.csv")


