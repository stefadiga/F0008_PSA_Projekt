
#--packages required
library(doBy)
library(sqldf)
library(xtable)
library(plyr)
library(dplyr)
library(ggplot2)
library(lattice)


# -- Master tables

#############
#Study End
data=read.csv("f0008_psa_2_studyend.csv") 

data<-orderBy(~pat_id+start_dt, data)


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
data$employ_status[data$employ_status=="Selbständig "]<-"Selbständig"

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


data1$psa_value<-as.numeric(as.character(data1$psa_value))
data1$arzt_doby<-as.numeric(as.character(data1$arzt_doby))

data1$arzt_region[data1$arzt_region=="\\N"] <-""
data1$arzt_region_code[data1$arzt_region_code=="\\N"] <-""
data1$practice_type[data1$practice_type==''] <- ""

#adding sum of measurement by patient
n_psa_id1<-data.frame(pat_id=unique(data1$pat_id), n_psa=tapply(data1$psa_n, data1$pat_id, sum), n_arzt=tapply(data1$arzt_id, data1$pat_id, n_distinct), n_py=tapply(data1$py, data1$pat_id, sum), fr=tapply(data1$psa_n, data1$pat_id, sum)/tapply(data1$py, data1$pat_id, sum))
data1<-merge(data1, n_psa_id)

#Data Frame with sum of measurement by start_dt
n_psa_start_dt_s1<-data.frame(start_dt=unique(data1$start_dt), n_psa_start_dt=tapply(data1$psa_n, data1$start_dt, sum), ir=tapply(data1$psa_n, data1$start_dt, sum)/tapply(data1$py, data1$start_dt, sum))

