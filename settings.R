
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

# -- Master tables

#############
#data sets. Pat=description of patients, PSA: events, Arzt: arzt description
data=read.csv("f0008_01_pat.csv", stringsAsFactors=FALSE) 
data2=read.csv("f0008_02_psa.csv", stringsAsFactors=FALSE) 
data3=read.csv("f0008_03_arzt.csv", stringsAsFactors=FALSE) 

data<-orderBy(~pat_id+study_start_dt, data)
data2<-orderBy(~pat_id, data2)
data_all<-merge(data, data2, all=TRUE)
#eligible
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
mortality=read.csv("f0008_ch_demography.csv") 
mortality1<-subset(mortality, survey_year==2009 & sex2=="m")
mortality2<-subset(mortality, survey_year==2010 & sex2=="m")
mort_2009_m <- reshape(mortality1, direction = "wide", idvar="age_y", timevar="fact2", drop=c("canton", "survey_year", "sex2"))
mort_2010_m <- reshape(mortality2, direction = "wide", idvar="age_y", timevar="fact2", drop=c("canton", "survey_year", "sex2"))
mort_2009_m$year<-2009
mort_2010_m$year<-2010
mort_0910<-join(mort_2009_m, mort_2010_m, type="full")
mort_0910$rate<-mort_0910$n.death_cases/mort_0910$n.pop_middle
mort_0910<-subset(mort_0910, age_y>=50 & age_y <=80 & year==2010)

#Parametric estimation of hazard rates (Gompertz distribution)
fm4DNase1 <- nls(rate ~ b* exp(c*age_y),
                 data=mort_0910,
                 start = list(b = 0.2, c = 0.1),
                 algorithm = "port")
summary(fm4DNase1)

#verifying fitting
plot(mort_0910$age_y, mort_0910$rate)
lines(mort_0910$age_y, predict(fm4DNase1), col = 2)

#generate time of death
set.seed(50)
y2 <- log(1-log(runif(12722))/summary(fm4DNase1)[10]$coefficients[1,1])/summary(fm4DNase1)[10]$coefficients[2,1]

#simulating death for cases where we don't have a consultation in the last trimester and where cases are not diseased
data_sim<-subset(data_el, as.Date(last_kons_dt)<"2016-10-01" & is.na(event_first_disease_dt))
data_sim<-orderBy(~pat_id, data_sim)
dead_sim<-data.frame(pat_id=unique(data_sim$pat_id), y2, age_dc=tapply(data_sim$pat_age_end_2, data_sim$pat_id, mean))
dead_sim$dead<-as.numeric(dead_sim$y2<=75 & dead_sim$y2>=55 & dead_sim$age_dc<=y2)                     

data_el<-merge(dead_sim, data_el, all=TRUE)


data_el$fup_end_dt_1<-ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), data_el$fup_end_dt_1)
data_el$fup_end_dt_2<-ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), data_el$fup_end_dt_2)
#data_el$fup_end_dt_3<-ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), ifelse(is.na(data_el$dead), data_el$study_end_dt, ))
data_el$fup_end_dt_3<-ifelse(data_el$dead==1, data_el$last_kons_dt, data_el$study_end_dt)
data_el$fup_end_dt_3<-ifelse(is.na(data_el$dead), ifelse(!is.na(data_el$event_first_disease_dt), as.character.Date(data_el$event_first_disease_dt), data_el$study_end_dt), data_el$fup_end_dt_3)

data_el$fup_y_1<-(as.Date(data_el$fup_end_dt_1)-as.Date(data_el$fup_start_dt_1))/365
data_el$fup_y_2<-(as.Date(data_el$fup_end_dt_2)-as.Date(data_el$fup_start_dt_2))/365
data_el$fup_y_3<-(as.Date(data_el$fup_end_dt_3)-as.Date(data_el$fup_start_dt_3))/365

data_el$pat_age_end_1<-data_el$pat_age_start_1+floor(data_el$fup_y_1)
data_el$pat_age_end_2<-data_el$pat_age_start_2+floor(data_el$fup_y_2)
data_el$pat_age_end_3<-data_el$pat_age_start_3+floor(data_el$fup_y_3)

#export
write.csv(data_el, file="f0008_data.csv")


