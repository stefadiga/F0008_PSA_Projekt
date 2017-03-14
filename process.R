
# -----TO DO-----------------------------------------------
# Mixed model PSA Value
# ----------------------------------------------------------------

#to review
a00<-subset(n_psa_id, n_psa>=6)

a001a0_3 <- subset(data, pat_id %in% a00$pat_id & !is.na(psa_value))


xyplot(psa_value~start_dt|as.factor(pat_id), data=a001a0_3,
       type=c("p","g","r"),col="dark blue",col.line="black",
       xlab="start_dt",
       ylab="PSA value")


# -------------------------------
# Incidence rates
# ----------------------------------

# v is the indicator vector of first patient's id appearance
v <- c(1, rep(0, length(data_el$pat_id)-1))
for (i in 2:length(data_el$pat_id))
  if ((data_el$pat_id[i] %in% 1:data_el$pat_id[i-1]) == FALSE)
    v[i] <- 1;


quarter_starts <- seq(from=as.Date("2009-01-01"), to=as.Date("2017-01-01"), by="quarter")

cev <- function(end, start, subgroup) {
  
  py <- rep(0, (length(quarter_starts)-1));
  events <- rep(0, (length(quarter_starts)-1));
  
  for (i in 1:(length(quarter_starts)-1)) {
    
    
    Qstart <- rep(as.character(quarter_starts[i]), length(data_el$pat_id));
    Qend <- rep(as.character(quarter_starts[i+1]), length(data_el$pat_id));
    
    py[i] <- sum(v * subgroup * (pmax(0,
                           as.Date(pmin(end, Qend)) -
                             as.Date(pmax(start, Qstart)))) / 365);
    
    events[i] <- sum(subgroup * (data_el$valid_from_dt >= Qstart) &
                       (data_el$valid_from_dt <= Qend) &
                       (data_el$valid_from_dt >= start) &
                       (data_el$valid_from_dt <= end), na.rm = T);
    
  }
  return(cbind(py, events))}

inc<-data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, 1), cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, 1), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, 1))  

inc_arzt_sex<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_sex=="f"),sex="f", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_sex=="f"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_sex=="f")), 
          data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_sex=="m"),sex="m", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_sex=="m"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_sex=="m")))

inc_arzt_prac<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Doppelpraxis"),type="Doppelpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Doppelpraxis"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Doppelpraxis")), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Einzelpraxis"),type="Einzelpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Einzelpraxis"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Einzelpraxis")),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Gruppenpraxis"),type="Gruppenpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Gruppenpraxis"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Gruppenpraxis")))

inc_arzt_reg<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="CEN"),region="CEN", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="CEN"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="CEN")), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="IND"),region="IND-Ter", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="IND"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="IND")),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="PERI"),region="PERI", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="PERI"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="PERI")),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="SUB"),region="SUB", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="SUB"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="SUB")))

inc_arzt_emp<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="Angestellt"), employ="Angestellt", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="Angestellt"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Angestellt")), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="Selbständig"), employ="Selbständig", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="Selbständig"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Selbständig")))

#Plot Incidence
i1<-qplot(as.Date(inc$dt), inc$events/inc$py,  geom=c("point", "smooth"),
            ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed")
i2<-qplot(as.Date(inc$dt), inc$events.1/inc$py.1,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End")
i3<-qplot(as.Date(inc$dt), inc$events.2/inc$py.2,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed")
grid.arrange(i1, i2, i3, ncol=2)

#verifying sum py per dt and per patients
#Observed
sum(inc$py)
sum(data_el$fup_y_1*v)
#Study End
sum(inc$py.1)
sum(data_el$fup_y_2*v)
#Mixed
sum(inc$py.2)
sum(data_el$fup_y_3*v)

#by sex
sum(inc_arzt_sex$py *(inc_arzt_sex$sex=="f"))
sum(data_el$fup_y_1*v*(data_el$arzt_sex=="f"))
sum(inc_arzt_sex$py*(inc_arzt_sex$sex=="m"))
sum(data_el$fup_y_1*v*(data_el$arzt_sex=="m"))
#Study End
sum(inc_arzt_sex$py.1 *(inc_arzt_sex$sex=="f"))
sum(data_el$fup_y_2*v*(data_el$arzt_sex=="f"))
sum(inc_arzt_sex$py.1*(inc_arzt_sex$sex=="m"))
sum(data_el$fup_y_2*v*(data_el$arzt_sex=="m"))
#Mixed
sum(inc_arzt_sex$py.2 *(inc_arzt_sex$sex=="f"))
sum(data_el$fup_y_3*v*(data_el$arzt_sex=="f"))
sum(inc_arzt_sex$py.2*(inc_arzt_sex$sex=="m"))
sum(data_el$fup_y_3*v*(data_el$arzt_sex=="m"))

#by practice
sum(inc_arzt_prac$py *(inc_arzt_prac$type=="Doppelpraxis"))
sum(data_el$fup_y_1*v*(data_el$practice_type=="Doppelpraxis"))
sum(inc_arzt_prac$py*(inc_arzt_prac$type=="Einzelpraxis"))
sum(data_el$fup_y_1*v*(data_el$practice_type=="Einzelpraxis"))
sum(inc_arzt_prac$py*(inc_arzt_prac$type=="Gruppenpraxis"))
sum(data_el$fup_y_1*v*(data_el$practice_type=="Gruppenpraxis"))

#Study End
sum(inc_arzt_prac$py.1 *(inc_arzt_prac$type=="Doppelpraxis"))
sum(data_el$fup_y_2*v*(data_el$practice_type=="Doppelpraxis"))
sum(inc_arzt_prac$py.1*(inc_arzt_prac$type=="Einzelpraxis"))
sum(data_el$fup_y_2*v*(data_el$practice_type=="Einzelpraxis"))
sum(inc_arzt_prac$py.1*(inc_arzt_prac$type=="Gruppenpraxis"))
sum(data_el$fup_y_2*v*(data_el$practice_type=="Gruppenpraxis"))

#Mixed
sum(inc_arzt_prac$py.2 *(inc_arzt_prac$type=="Doppelpraxis"))
sum(data_el$fup_y_3*v*(data_el$practice_type=="Doppelpraxis"))
sum(inc_arzt_prac$py.2*(inc_arzt_prac$type=="Einzelpraxis"))
sum(data_el$fup_y_3*v*(data_el$practice_type=="Einzelpraxis"))
sum(inc_arzt_prac$py.2*(inc_arzt_prac$type=="Gruppenpraxis"))
sum(data_el$fup_y_3*v*(data_el$practice_type=="Gruppenpraxis"))


#by sex
i1<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events/inc_arzt_sex$py), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + theme(legend.title=element_blank())

i2<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events.1/inc_arzt_sex$py.1), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + theme(legend.title=element_blank())

i3<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events.2/inc_arzt_sex$py.2), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + theme(legend.title=element_blank())

grid.arrange(i1, i2, i3, ncol=2)

#by practice
i1<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events/inc_arzt_prac$py), color=inc_arzt_prac$type, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Practice type")

i2<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events.1/inc_arzt_prac$py.1), color=inc_arzt_prac$type, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + theme(legend.title=element_blank())

i3<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events.2/inc_arzt_prac$py.2), color=inc_arzt_prac$type, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + theme(legend.title=element_blank())

grid.arrange(i1, i2, i3, ncol=2)

#by region
#to do add AGR?
i1<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events/inc_arzt_reg$py), color=inc_arzt_reg$region, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + theme(legend.title=element_blank())

i2<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events.1/inc_arzt_reg$py.1), color=inc_arzt_reg$region, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + theme(legend.title=element_blank())

i3<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events.2/inc_arzt_reg$py.2), color=inc_arzt_reg$region, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + theme(legend.title=element_blank())

grid.arrange(i1, i2, i3, ncol=2)

#by employ status
i1<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events/inc_arzt_emp$py), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + theme(legend.title=element_blank())

i2<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events.1/inc_arzt_emp$py.1), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + theme(legend.title=element_blank())

i3<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events.2/inc_arzt_emp$py.2), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + theme(legend.title=element_blank())

grid.arrange(i1, i2, i3, ncol=2)

#Plot Frequency
qplot(as.Date(inc$dt),inc$events, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests")

n_psa_id<-data.frame(pat_id=unique(data_el$pat_id), n_psa=tapply(!is.na(data_el$psa), data_el$pat_id, sum)) 
n_psa_id<-join(n_psa_id, data_el, match ="first")
n_psa_id$n_inc1<-n_psa_id$n_psa/as.numeric(n_psa_id$fup_y_1)
n_psa_id$n_inc2<-n_psa_id$n_psa/as.numeric(n_psa_id$fup_y_2)
n_psa_id$n_inc3<-n_psa_id$n_psa/as.numeric(n_psa_id$fup_y_3)

######To DO
####
#tabular ohne missing
table(practice_type, exclude= "")
table(employ_status, exclude= "")


#study end tables (categories of frequency rates per patient)
n_psa_id$fr_c<-cut(n_psa_id$fr, c(0,0.1, 0.15,0.25,0.3,0.5,1, 5), right=FALSE)
n_psa_id$c_age_m <- cut(n_psa_id$age_m, c(55,57,60,65,70,76), right=FALSE)
n_psa_id$c_age_dum <- as.numeric(n_psa_id$age_m>=62)
n_psa_id$c_doby <- cut(n_psa_id$arzt_doby, c(1945,1955,1965,1990), right=FALSE, dig.lab=4)

table(n_psa_id$fr_c)
table(n_psa_id$c_age_m)
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$c_age_dum))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$c_age_dum))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_region_code, exclude=""))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$practice_type, exclude=""))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$practice_type, exclude=""))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_sex))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_sex))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$employ_status, exclude=""))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$c_doby))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$c_doby))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$employ_status))
table(n_psa_id$employ_status)
#mean value of PSA

#observed tables (categories of frequency rates per patient)
n_psa_id$fr_c<-cut(n_psa_id$fr, c(0,0.1,0.5,1,5), right=FALSE)
n_psa_id$c_age_m <- cut(n_psa_id$age_m, c(55,57,60,65,70,76), right=FALSE)
n_psa_id$c_age_dum <- as.numeric(n_psa_id$age_m>=62)
n_psa_id$c_doby <- cut(n_psa_id$arzt_doby, c(1945,1955,1965,1990), right=FALSE, dig.lab=4)

table(n_psa_id$fr_c)
table(n_psa_id$c_age_m)
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$c_age_dum))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$c_age_dum))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_region_code, exclude=""))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$practice_type, exclude=""))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$practice_type, exclude=""))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_sex))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_sex))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$employ_status, exclude=""))
ftable(xtabs(~n_psa_id$fr_c+ n_psa_id$c_doby))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$c_doby))
chisq.test(xtabs(~n_psa_id$fr_c+ n_psa_id$employ_status))
table(n_psa_id$employ_status)


#graphs
barplot(xtabs(~n_psa_id$fr_c+ n_psa_id$c_doby))
barplot(xtabs(~n_psa_id$fr_c+ n_psa_id$employ_status, exclude=""))
barplot(xtabs(~n_psa_id$fr_c+ n_psa_id$practice_type, exclude=""))
barplot(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_region_code, exclude=""))
barplot(xtabs(~n_psa_id$fr_c+ n_psa_id$arzt_sex))
