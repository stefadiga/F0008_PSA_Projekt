
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

# -- Master tables

#############
#Study End
#data=read.csv("f0008_psa_2_studyend.csv") 
data=read.csv("f0008_01_pat.csv") 
data2=read.csv("f0008_02_psa.csv") 
data3=read.csv("f0008_03_arzt.csv") 

data<-orderBy(~pat_id+study_start_dt, data)
data2<-orderBy(~pat_id, data2)
data_all<-merge(data, data2, all=TRUE)


#data cleaning for missing
data$arzt_region<-as.character(data$arzt_region)
data$arzt_region[data$arzt_region=="\\N"] <-""

data$arzt_region_code<-as.character(data$arzt_region_code)
data$arzt_region_code[data$arzt_region_code=="\\N"] <-""

data$practice_type[data$practice_type==''] <- ""

data$arzt_doby<-as.character(data$arzt_doby)
data$arzt_doby[data$arzt_doby=="\\N"] <-""
data$arzt_doby<-as.numeric(as.character(data$arzt_doby))

data$psa_value<-as.character(data$psa_value)
data$psa_value[data$psa_value=="\\N"] <-""
data$psa_value<-as.numeric(as.character(data$psa_value))

data$employ_status<-as.character(data$employ_status)
data$employ_status[data$employ_status=="Selbst?ndig "]<-"Selbst?ndig"

#mortality rate
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

weib<-fitdistr(mort_0910$rate, densfun=dweibull,start=list(scale=1,shape=5))

curve(dweibull(x, scale=weib$estimate[1], shape=weib$estimate[2]),
 from=50, to=80, add=TRUE)

fm4DNase1 <- nls(rate ~ b* exp(c*age_y),
                 data=mort_0910,
                 start = list(b = 0.2, c = 0.1),
                 algorithm = "port")
summary(fm4DNase1)

plot(mort_0910$age_y, mort_0910$rate)
lines(mort_0910$age_y, predict(fm4DNase1), col = 2)

rnos<-runif(100) 
which(rnos<= 0.1) 
y2 <- (-log(summary(fm4DNase1)[10]$coefficients[1,1]) + log(runif(1000))) /summary(fm4DNase1)[10]$coefficients[2,1]

#arzt consultation - tutti i pazienti del dottore
arzt_cons=read.csv("f0008_arzt.csv") 


#adding sum of measurement by patient
n_psa_id<-data.frame(pat_id=unique(data$pat_id), n_psa=tapply(data$psa_n, data$pat_id, sum), n_arzt=tapply(data$arzt_id, data$pat_id, n_distinct), n_py=tapply(data$py, data$pat_id, sum), 
                     age_m=tapply(data$age, data$pat_id, mean), psa_mean=tapply(data$psa_value, data$pat_id, mean, na.rm=TRUE), fr=tapply(data$psa_n, data$pat_id, sum)/tapply(data$py, data$pat_id, sum))
#db by patient
n_psa_id<-join(n_psa_id, data, match ="first")

#Data Frame with sum of measurement by start_dt
n_psa_start_dt_s<-data.frame(start_dt=unique(data$start_dt), n_psa_start_dt=tapply(data$psa_n, data$start_dt, sum), ir=tapply(data$psa_n, data$start_dt, sum)/tapply(data$py, data$start_dt, sum))

############
#Observed
data1=read.csv("f0008_psa_1_observed.csv") 
data1<-orderBy(~pat_id+start_dt, data1)

#data cleaning for missing
data1$arzt_region<-as.character(data1$arzt_region)
data1$arzt_region[data1$arzt_region=="\\N"] <-""

data1$arzt_region_code<-as.character(data1$arzt_region_code)
data1$arzt_region_code[data1$arzt_region_code=="\\N"] <-""

data1$practice_type[data1$practice_type==''] <- ""

data1$arzt_doby<-as.character(data1$arzt_doby)
data1$arzt_doby[data1$arzt_doby=="\\N"] <-""
data1$arzt_doby<-as.numeric(as.character(data1$arzt_doby))

data1$psa_value<-as.character(data1$psa_value)
data1$psa_value[data1$psa_value=="\\N"] <-""
data1$psa_value<-as.numeric(as.character(data1$psa_value))

data1$employ_status<-as.character(data1$employ_status)
data1$employ_status[data1$employ_status=="Selbst?ndig "]<-"Selbst?ndig"

#adding sum of measurement by patient
n_psa1_id<-data.frame(pat_id=unique(data1$pat_id), n_psa=tapply(data1$psa_n, data1$pat_id, sum), n_arzt=tapply(data1$arzt_id, data1$pat_id, n_distinct), n_py=tapply(data1$py, data1$pat_id, sum), 
                     age_m=tapply(data1$age, data1$pat_id, mean), psa_mean=tapply(data1$psa_value, data1$pat_id, mean, na.rm=TRUE), fr=tapply(data1$psa_n, data1$pat_id, sum)/tapply(data1$py, data1$pat_id, sum))
#db by patient
n_psa1_id<-join(n_psa_id, data1, match ="first")

#Data Frame with sum of measurement by start_dt
n_psa_start_dt_s1<-data.frame(start_dt=unique(data1$start_dt), n_psa_start_dt=tapply(data1$psa_n, data1$start_dt, sum), ir=tapply(data1$psa_n, data1$start_dt, sum)/tapply(data1$py, data1$start_dt, sum))

