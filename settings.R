
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
#Study End
#data=read.csv("f0008_psa_2_studyend.csv") 
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

data_el$event_first_disease_dt<-pmin(as.Date(data_el$event_prostcancer_dt, format="%Y-%m-%d"), as.Date(data_el$event_urinarytract_dt, format="%Y-%m-%d"), as.Date(data_el$event_prosthyper_dt, format="%Y-%m-%d"), as.Date(data_el$event_atc_g04_dt, format="%Y-%m-%d"), na.rm=TRUE)
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

#weib<-fitdistr(mort_0910$rate, densfun=dweibull,start=list(scale=1,shape=5))

#curve(dweibull(x, scale=weib$estimate[1], shape=weib$estimate[2]),
# from=50, to=80, add=TRUE)

fm4DNase1 <- nls(rate ~ b* exp(c*age_y),
                 data=mort_0910,
                 start = list(b = 0.2, c = 0.1),
                 algorithm = "port")
summary(fm4DNase1)

plot(mort_0910$age_y, mort_0910$rate)
lines(mort_0910$age_y, predict(fm4DNase1), col = 2)

set.seed(50)
y2 <- log(1-log(runif(12722))/summary(fm4DNase1)[10]$coefficients[1,1])/summary(fm4DNase1)[10]$coefficients[2,1]

#simulating death
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

write.csv(data_el, file="f0008_data.csv")


# v is the indicator vector of first patient's id appearance
v <- c(1, rep(0, length(data_el$pat_id)-1))
for (i in 2:length(data_el$pat_id))
  if ((data_el$pat_id[i] %in% 1:data_el$pat_id[i-1]) == FALSE)
    v[i] <- 1;


quarter_starts <- seq(from=as.Date("2009-01-01"), to=as.Date("2017-01-01"), by="quarter")

cev <- function(end, start) {
  
  py <- rep(0, (length(quarter_starts)-1));
  events <- rep(0, (length(quarter_starts)-1));
  
  for (i in 1:(length(quarter_starts)-1)) {
    
    Qstart <- rep(as.character(quarter_starts[i]), length(data_el$pat_id));
    Qend <- rep(as.character(quarter_starts[i+1]), length(data_el$pat_id));
    
    py[i] <- sum(v * (pmax(0,
                              as.Date(pmin(end, Qend)) -
                              as.Date(pmax(start, Qstart)))) / 365);
    
    events[i] <- sum((data_el$valid_from_dt >= Qstart) &
                     (data_el$valid_from_dt <= Qend) &
                     (data_el$valid_from_dt >= start) &
                     (data_el$valid_from_dt <= end), na.rm = T);
    
  }
return(cbind(py, events))}

inc<-data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1), cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3))  
#arzt consultation - tutti i pazienti del dottore
arzt_cons=read.csv("f0008_arzt.csv") 


#adding sum of measurement by patient
n_psa_id<-data.frame(pat_id=unique(data_el$pat_id), n_psa=tapply(data_el$psa, data_el$pat_id, length)) 
#                     age_m=tapply(data$age, data$pat_id, mean), psa_mean=tapply(data$psa_value, data$pat_id, mean, na.rm=TRUE), fr=tapply(data$psa_n, data$pat_id, sum)/tapply(data$py, data$pat_id, sum))
#db by patient
n_psa_id<-join(n_psa_id, data_el, match ="first")

#Data Frame with sum of measurement by start_dt
n_psa_id$trim1<-cut(as.Date(n_psa_id$fup_end_dt_1), seq(from=as.Date('2009-1-1'),to=as.Date('2017-1-1'),by='quarter'), right=TRUE)
n_psa_id$trim2<-cut(as.Date(n_psa_id$fup_end_dt_2), seq(from=as.Date('2009-1-1'),to=as.Date('2017-1-1'),by='quarter'), right=TRUE)
n_psa_id$trim3<-cut(as.Date(n_psa_id$fup_end_dt_3), seq(from=as.Date('2009-1-1'),to=as.Date('2017-1-1'),by='quarter'), right=TRUE)

n_psa_start_dt<-data.frame(start_dt_1=unique(n_psa_id$trim1), n_psa_start_dt_1=tapply(n_psa_id$n_psa, n_psa_id$trim1, sum), py1=tapply(n_psa_id$fup_y_1, n_psa_id$trim1, sum), ir1=tapply(n_psa_id$n_psa, n_psa_id$trim1, sum)/tapply(n_psa_id$fup_y_1, n_psa_id$trim1, sum),
                            start_dt_2=unique(n_psa_id$trim2), n_psa_start_dt_2=tapply(n_psa_id$n_psa, n_psa_id$trim2, sum), py2=tapply(n_psa_id$fup_y_2, n_psa_id$trim2, sum), ir2=tapply(n_psa_id$n_psa, n_psa_id$trim2, sum)/tapply(n_psa_id$fup_y_2, n_psa_id$trim2, sum),
                            start_dt_3=unique(n_psa_id$trim3), n_psa_start_dt_3=tapply(n_psa_id$n_psa, n_psa_id$trim3, sum), py3=tapply(n_psa_id$fup_y_3, n_psa_id$trim3, sum), ir3=tapply(n_psa_id$n_psa, n_psa_id$trim3, sum)/tapply(n_psa_id$fup_y_3, n_psa_id$trim3, sum))
n_psa_start_dt<-orderBy(~start_dt_1, n_psa_start_dt)

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

