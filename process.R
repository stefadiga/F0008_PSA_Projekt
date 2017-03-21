
# -----TO DO-----------------------------------------------
# Mixed model PSA Value
# ----------------------------------------------------------------

#adding sum of measurement by patient
n_psa_id<-data.frame(pat_id=unique(data_el$pat_id), n_psa=tapply(data_el$psa, data_el$pat_id, length)) 
#                     age_m=tapply(data$age, data$pat_id, mean), psa_mean=tapply(data$psa_value, data$pat_id, mean, na.rm=TRUE), fr=tapply(data$psa_n, data$pat_id, sum)/tapply(data$py, data$pat_id, sum))
#db by patient
n_psa_id<-join(n_psa_id, data_el, match ="first")
a00<-subset(n_psa_id, n_psa>=6)

a001a0_3 <- subset(data_el, pat_id %in% a00$pat_id & !is.na(psa))


xyplot(psa~as.Date(valid_from_dt)|as.factor(pat_id), data=a001a0_3,
       type=c("p","g","r"),col="dark blue",col.line="black", ylim=c(0,8),
       xlab="start_dt",
       ylab="PSA value")

fm_mix<-lmer(psa~as.Date(valid_from_dt)+(as.Date(valid_from_dt)|pat_id), a001a0_3)

plot(fm_mix, sqrt(abs(resid(.)))~fitted(.), type=c("p", "smooth"))
------------------------------
# Incidence rates
# ----------------------------------

# v is the indicator vector of first patient's id appearance

data_el<-orderBy(~pat_id, data_el)

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
                             as.Date(pmax(start, Qstart)))) / 365, na.rm = T);
    
    events[i] <- sum(subgroup * (data_el$valid_from_dt >= Qstart) &
                       (data_el$valid_from_dt < Qend) &
                       (data_el$valid_from_dt >= start) &
                       (data_el$valid_from_dt < end), na.rm = T);
    
  }
  return(cbind(py, events))}


#Data Set for Graphs
inc<-data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$elig_2_1), cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$elig_2_3))  

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

data_el$c_doby <- cut(data_el$arzt_doby, c(1940,1950,1960,1970, 1980), right=FALSE, labels=FALSE)

inc_arzt_doby<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==2),doby="[1950,1960)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==2)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==3),doby="[1960,1970)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==3), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==4),doby="[1970,1980)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==4), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==4)))

data_el$c_age <- cut(data_el$pat_age_start_3, c(55,60,65,70,75), right=FALSE)

inc_pat_age<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[55,60)"),age="[55,60)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[55,60)"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[55,60)")), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[60,65)"),age="[60,65)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[60,65)"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[60,65)")),
                   data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[65,70)"),age="[65,70)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[65,70)"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[65,70)")),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[70,75)"),age="[70,75)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[70,75)"), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[70,75)")))

cev_age <- function(end, start, class_age) {
  # end & start are dates for patients' stay in the study period
  # for instance class_age <- c(50, 60, 70) should partition the possible ages in
  # [0, 49], [50, 59], [60, 69], [70, 120]
  
  res <- NULL;
  
  n_classes <- length(class_age);
  
  for (j in 1:(n_classes+1)) {
    
    # computes start and end instants for class j
    if (j == 1) {
      class_start <- 0; class_end <- class_age[1]-1;
    } else {
      if (j == (n_classes+1)) {
        class_start <- class_age[n_classes]; class_end <- 120;
      } else {
        class_start <- class_age[j-1]; class_end <- class_age[j]-1;
      }
    }
    
    py <- rep(0, (length(quarter_starts)-1));
    events <- rep(0, (length(quarter_starts)-1));
    
    # loops for all quarters
    for (i in 1:(length(quarter_starts)-1)) {
      
      # quarter i start and end instants
      Qstart <- rep(as.character(quarter_starts[i]), length(data_el$pat_id));
      Qend <- rep(as.character(quarter_starts[i+1]), length(data_el$pat_id));
      
      # current ages for patients in the database
      current_age <- floor(2009 - data_el$doby + # age at start of study period
                             i / 4);               # years age at current quarter
      
      # indicator function of patients belonging to class j
      grouping <- (current_age >= class_start) & (current_age <= class_end);
      
      py[i] <- sum(v * grouping * (pmax(0,
                                        as.Date(pmin(end, Qend)) -
                                          as.Date(pmax(start, Qstart)))) / 365, na.rm = T);
      
      events[i] <- sum(grouping * (data_el$valid_from_dt >= Qstart) &
                         (data_el$valid_from_dt < Qend) &
                         (data_el$valid_from_dt >= start) &
                         (data_el$valid_from_dt < end), na.rm = T);
      
    }
    
    res <- cbind(res, py, events);
   
  }
  
  return(res);
}

class_age <-  c(50, 55, 60, 65, 70, 75);
date()
prova_2 <- cev_age(data_el$fup_end_dt_3, data_el$fup_start_dt_3, class_age);
date()

inc_pat_age_c<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,1], events=prova_2[,2], class_age=rep(50, 32)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,3], events=prova_2[,4], class_age=rep(55, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,5], events=prova_2[,6], class_age=rep(60, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,7], events=prova_2[,8], class_age=rep(65, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,9], events=prova_2[,10], class_age=rep(70, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,11], events=prova_2[,12], class_age=rep(75, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,13], events=prova_2[,14], class_age=rep(76, 32)))
             
# total number of events:
sum(prova[, 2*(1:(length(class_age)+1))])

#CHECK

#verifying sum py per dt and per patients
#Observed
sum(inc$py)
sum(data_el$fup_y_1*data_el$elig_2_1*v)
#Study End
sum(inc$py.1)
sum(data_el$fup_y_2*data_el$elig_2_2*v)
#Mixed
sum(inc$py.2)
sum(data_el$fup_y_3*data_el$elig_2_3*v)

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

#Employ status
#Observed
sum(inc_arzt_emp$py *(inc_arzt_emp$employ=="Angestellt"))
sum(data_el$fup_y_1*v*(data_el$employ_status=="Angestellt"))
sum(inc_arzt_emp$py*(inc_arzt_emp$employ=="Selbständig"))
sum(data_el$fup_y_1*v*(data_el$employ_status=="Selbständig"))
#Study End
sum(inc_arzt_emp$py.1 *(inc_arzt_emp$employ=="Angestellt"))
sum(data_el$fup_y_2*v*(data_el$employ_status=="Angestellt"))
sum(inc_arzt_emp$py.1*(inc_arzt_emp$employ=="Selbständig"))
sum(data_el$fup_y_2*v*(data_el$employ_status=="Selbständig"))
#Mixed
sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="Angestellt"))
sum(data_el$fup_y_3*v*(data_el$employ_status=="Angestellt"))
sum(inc_arzt_emp$py.2*(inc_arzt_emp$employ=="Selbständig"))
sum(data_el$fup_y_3*v*(data_el$employ_status=="Selbständig"))

#Arzt DOBY
#Observed
sum(inc_arzt_doby$py *(inc_arzt_doby$doby=="[1950,1960)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==2), na.rm=T)
sum(inc_arzt_doby$py*(inc_arzt_doby$doby=="[1960,1970)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==3), na.rm=T)
sum(inc_arzt_doby$py*(inc_arzt_doby$doby=="[1970,1980)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==4), na.rm=T)
#Study End
sum(inc_arzt_doby$py.1 *(inc_arzt_doby$doby=="[1950,1960)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==2), na.rm=T)
sum(inc_arzt_doby$py.1*(inc_arzt_doby$doby=="[1960,1970)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==3), na.rm=T)
sum(inc_arzt_doby$py.1*(inc_arzt_doby$doby=="[1970,1980)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==4), na.rm=T)
#Mixed
sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1950,1960)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==2), na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby=="[1960,1970)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==3), na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby=="[1970,1980)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==4), na.rm=T)

#Pat Age
#Observed
sum(inc_pat_age$py *(inc_pat_age$age=="[55,60)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[55,60)"), na.rm=T)
sum(inc_pat_age$py*(inc_pat_age$age=="[60,65)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[60,65)"), na.rm=T)
sum(inc_pat_age$py*(inc_pat_age$age=="[65,70)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[65,70)"), na.rm=T)
sum(inc_pat_age$py*(inc_pat_age$age=="[70,75)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[70,75)"), na.rm=T)
#Study End
sum(inc_pat_age$py.1 *(inc_pat_age$age=="[55,60)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[55,60)"), na.rm=T)
sum(inc_pat_age$py.1*(inc_pat_age$age=="[60,65)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[60,65)"), na.rm=T)
sum(inc_pat_age$py.1*(inc_pat_age$age=="[65,70)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[65,70)"), na.rm=T)
sum(inc_pat_age$py.1*(inc_pat_age$age=="[70,75)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[70,75)"), na.rm=T)
#Mixed
sum(inc_pat_age$py.2 *(inc_pat_age$age=="[55,60)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[55,60)"), na.rm=T)
sum(inc_pat_age$py.2*(inc_pat_age$age=="[60,65)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[60,65)"), na.rm=T)
sum(inc_pat_age$py.2*(inc_pat_age$age=="[65,70)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[65,70)"), na.rm=T)
sum(inc_pat_age$py.2*(inc_pat_age$age=="[70,75)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[70,75)"), na.rm=T)

#events
#Observed
sum(inc$events)
sum(as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1), na.rm=T)
#Study End
sum(inc$events.1)
sum(as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2), na.rm=T)
#Mixed
sum(inc$events.2)
sum(as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3), na.rm=T)

#Plot Incidence
i1<-qplot(as.Date(inc$dt), inc$events/inc$py,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed")
i2<-qplot(as.Date(inc$dt), inc$events.1/inc$py.1,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End")
i3<-qplot(as.Date(inc$dt), inc$events.2/inc$py.2,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed")
grid.arrange(i1, i2, i3, ncol=2)


#by sex
i1<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events/inc_arzt_sex$py), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Sex")

i2<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events.1/inc_arzt_sex$py.1), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Sex")

i3<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events.2/inc_arzt_sex$py.2), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Sex")

grid.arrange(i1, i2, i3, ncol=2)

#by practice
i1<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events/inc_arzt_prac$py), color=inc_arzt_prac$type, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Practice type")

i2<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events.1/inc_arzt_prac$py.1), color=inc_arzt_prac$type, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Practice type")

i3<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events.2/inc_arzt_prac$py.2), color=inc_arzt_prac$type, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Practice type")

grid.arrange(i1, i2, i3, ncol=2)

#by region
#to do add AGR?
i1<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events/inc_arzt_reg$py), color=inc_arzt_reg$region, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Region")

i2<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events.1/inc_arzt_reg$py.1), color=inc_arzt_reg$region, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Region")

i3<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events.2/inc_arzt_reg$py.2), color=inc_arzt_reg$region, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Region")

grid.arrange(i1, i2, i3, ncol=2)

#by employ status
i1<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events/inc_arzt_emp$py), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Employ status")

i2<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events.1/inc_arzt_emp$py.1), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Employ status")

i3<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events.2/inc_arzt_emp$py.2), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Employ status")

grid.arrange(i1, i2, i3, ncol=2)

#by arzt doby
i1<-qplot(as.Date(inc_arzt_doby$dt), as.numeric(inc_arzt_doby$events/inc_arzt_doby$py), color=inc_arzt_doby$doby, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Arzt DOBY")

i2<-qplot(as.Date(inc_arzt_doby$dt), as.numeric(inc_arzt_doby$events.1/inc_arzt_doby$py.1), color=inc_arzt_doby$doby, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Arzt DOBY")

i3<-qplot(as.Date(inc_arzt_doby$dt), as.numeric(inc_arzt_doby$events.2/inc_arzt_doby$py.2), color=inc_arzt_doby$doby, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Arzt DOBY")

grid.arrange(i1, i2, i3, ncol=2)

#by age
i1<-qplot(as.Date(inc_pat_age_c$dt[inc_pat_age_c$class_age>=55]), as.numeric(inc_pat_age_c$events[inc_pat_age_c$class_age>=55]/inc_pat_age_c$py[inc_pat_age_c$class_age>=55]), color=inc_pat_age_c$class_age[inc_pat_age_c$class_age>=55], geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Pat Age")

#i2<-qplot(as.Date(inc_pat_age$dt), as.numeric(inc_pat_age$events.1/inc_pat_age$py.1), color=inc_pat_age$age, geom=c("point", "smooth"),
      #    ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Pat Age")

#i3<-qplot(as.Date(inc_pat_age$dt), as.numeric(inc_pat_age$events.2/inc_pat_age$py.2), color=inc_pat_age$age, geom=c("point", "smooth"),
       #   ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Pat Age")

#grid.arrange(i1, i2, i3, ncol=2)

inc_pat_age_c$inc<-as.numeric(inc_pat_age_c$events/inc_pat_age_c$py)
inc_pat_age_c <- subset(inc_pat_age_c, inc_pat_age_c$class_age >=60 & inc_pat_age_c$class_age < 76)

i1<-qplot(data=inc_pat_age_c, x=as.Date(dt), y=inc, geom=c("point","smooth"),ylim=c(0,0.3),
      xlab="dt", main="Mixed follow up") + 
   facet_grid(.~as.character(class_age))
  
inc_pat_age$inc<-as.numeric(inc_pat_age$events.2/inc_pat_age$py.2)

i2<-qplot(data=inc_pat_age, x=as.Date(dt), y=inc, geom=c("point","smooth"),ylim=c(0,0.3),
      xlab="dt", main="Mixed follow up") + 
    facet_grid(.~age)
grid.arrange(i1, i2, ncol=2)

#Plot Frequency
qplot(as.Date(inc$dt),inc$events, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests")

#PSA

med<-tapply(data_el$psa, cut(as.Date(data_el$valid_from_dt), as.Date(quarter_starts)), median)
setq<-tapply(data_el$psa, cut(as.Date(data_el$valid_from_dt), as.Date(quarter_starts)), quantile, 0.75)
novq<-tapply(data_el$psa, cut(as.Date(data_el$valid_from_dt), as.Date(quarter_starts)), quantile, 0.90)

quant<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq)))          

qplot(data=quant, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
   ylim(c(0,8)), xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
ylab(expression(paste('PSA (' , mu,g/l,')')))
  
#by sex
med_f<-tapply(data_el$psa[data_el$arzt_sex=="f"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_sex=="f"]), as.Date(quarter_starts)), median) 
med_m<-tapply(data_el$psa[data_el$arzt_sex=="m"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_sex=="m"]), as.Date(quarter_starts)), median) 
setq_f<-tapply(data_el$psa[data_el$arzt_sex=="f"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_sex=="f"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_m<-tapply(data_el$psa[data_el$arzt_sex=="m"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_sex=="m"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_f<-tapply(data_el$psa[data_el$arzt_sex=="f"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_sex=="f"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_m<-tapply(data_el$psa[data_el$arzt_sex=="m"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_sex=="m"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_sex<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_f), sex="f"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_m), sex="m"),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_f), sex="f"),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_f), sex="m"),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_f), sex="f"),          
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_m), sex="m"))

qplot(data=quant_sex, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ sex)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

#region
med_C<-tapply(data_el$psa[data_el$arzt_region_code=="CEN"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="CEN"]), as.Date(quarter_starts)), median) 
med_I<-tapply(data_el$psa[data_el$arzt_region_code=="IND"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="IND"]), as.Date(quarter_starts)), median) 
med_O<-tapply(data_el$psa[data_el$arzt_region_code=="OTHER"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), median) 
med_P<-tapply(data_el$psa[data_el$arzt_region_code=="PERI"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="PERI"]), as.Date(quarter_starts)), median) 
med_S<-tapply(data_el$psa[data_el$arzt_region_code=="SUB"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="SUB"]), as.Date(quarter_starts)), median) 
setq_C<-tapply(data_el$psa[data_el$arzt_region_code=="CEN"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="CEN"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_I<-tapply(data_el$psa[data_el$arzt_region_code=="IND"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="IND"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_O<-tapply(data_el$psa[data_el$arzt_region_code=="OTHER"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_P<-tapply(data_el$psa[data_el$arzt_region_code=="PERI"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="PERI"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_S<-tapply(data_el$psa[data_el$arzt_region_code=="SUB"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="SUB"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_C<-tapply(data_el$psa[data_el$arzt_region_code=="CEN"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="CEN"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_I<-tapply(data_el$psa[data_el$arzt_region_code=="IND"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="IND"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_O<-tapply(data_el$psa[data_el$arzt_region_code=="OTHER"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_P<-tapply(data_el$psa[data_el$arzt_region_code=="PERI"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="PERI"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_S<-tapply(data_el$psa[data_el$arzt_region_code=="SUB"], cut(as.Date(data_el$valid_from_dt[data_el$arzt_region_code=="SUB"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_reg<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_C), region="CEN"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_I), region="IND"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_O), region="OTHER"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_P), region="PERI"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_S), region="SUB"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_C), region="CEN"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_I), region="IND"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_O), region="OTHER"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_P), region="PERI"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_S), region="SUB"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_C), region="CEN"), 
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_I), region="IND"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_O), region="OTHER"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_P), region="PERI"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_S), region="SUB"))         
                 
qplot(data=quant_reg, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ region)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

#by employ status
med_a<-tapply(data_el$psa[data_el$employ_status=="Angestellt"], cut(as.Date(data_el$valid_from_dt[data_el$employ_status=="Angestellt"]), as.Date(quarter_starts)), median) 
med_s<-tapply(data_el$psa[data_el$employ_status=="Selbständig"], cut(as.Date(data_el$valid_from_dt[data_el$employ_status=="Selbständig"]), as.Date(quarter_starts)), median) 
setq_a<-tapply(data_el$psa[data_el$employ_status=="Angestellt"], cut(as.Date(data_el$valid_from_dt[data_el$employ_status=="Angestellt"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_s<-tapply(data_el$psa[data_el$employ_status=="Selbständig"], cut(as.Date(data_el$valid_from_dt[data_el$employ_status=="Selbständig"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_a<-tapply(data_el$psa[data_el$employ_status=="Angestellt"], cut(as.Date(data_el$valid_from_dt[data_el$employ_status=="Angestellt"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_s<-tapply(data_el$psa[data_el$employ_status=="Selbständig"], cut(as.Date(data_el$valid_from_dt[data_el$employ_status=="Selbständig"]), as.Date(quarter_starts)),quantile, 0.90) 

quant_emp<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a), employ="Angestellt"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_s), employ="Selbständig"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a), employ="Angestellt"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_s), employ="Selbständig"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a), employ="Angestellt"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_s), employ="Selbständig"))

qplot(data=quant_emp, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ employ)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

#by praxis
med_d<-tapply(data_el$psa[data_el$practice_type=="Doppelpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), median) 
med_e<-tapply(data_el$psa[data_el$practice_type=="Einzelpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), median) 
med_g<-tapply(data_el$psa[data_el$practice_type=="Gruppenpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), median) 
setq_d<-tapply(data_el$psa[data_el$practice_type=="Doppelpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_e<-tapply(data_el$psa[data_el$practice_type=="Einzelpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_g<-tapply(data_el$psa[data_el$practice_type=="Gruppenpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_d<-tapply(data_el$psa[data_el$practice_type=="Doppelpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_e<-tapply(data_el$psa[data_el$practice_type=="Einzelpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_g<-tapply(data_el$psa[data_el$practice_type=="Gruppenpraxis"], cut(as.Date(data_el$valid_from_dt[data_el$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_pra<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_g), practice="Gruppenpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_g), practice="Gruppenpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_g), practice="Gruppenpraxis"))

qplot(data=quant_pra, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ practice)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

######Dummy tables
####
#tabular ohne missing
table(data_el$practice_type, exclude= "")
table(employ_status, exclude= "")

#to remark patients with age <= 55 & > = 76 with elig_vec 1-1-1-1-1-1-1
summary(data_el$pat_age_end_1)
table(data_el$dead*v)


######### Report tables
########


#N PSA
count_psa<-as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_1 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1)]), data_el$pat_id[data_el$elig_2_1 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName = "n_pat")

count_psa_2<-as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_2 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2)]), data_el$pat_id[data_el$elig_2_2 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName = "n_pat")

count_psa_3<-as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)]), data_el$pat_id[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName ="n_pat")


#check with length(unique(data_el$pat_id*data_el$elig_2_1) 
el<-subset(data_el, elig_2_1==1)
el1<-subset(data_el, elig_2_2==1)
el2<-subset(data_el, elig_2_3==1)

n_distinct(el$arzt_id[unique(el$pat_id)])

tab1<-data.frame(method=c("Observed", "Study Period", "Mixed"), n_pat=c(n_distinct(el$pat_id), n_distinct(el1$pat_id), n_distinct(el2$pat_id)), psa=c(sum(as.numeric(count_psa$n_psa)*count_psa$n_pat),sum(as.numeric(count_psa_2$n_psa)*count_psa_2$n_pat), sum(as.numeric(count_psa_3$n_psa)*count_psa_3$n_pat)), py=c(sum(data_el$fup_y_1*data_el$elig_2_1*v), sum(data_el$fup_y_2*data_el$elig_2_2*v), sum(data_el$fup_y_3*data_el$elig_2_3*v)), n_arzt=c(n_distinct(el$arzt_id[unique(el$pat_id)]), n_distinct(el1$arzt_id[unique(el1$pat_id)]),n_distinct(el2$arzt_id[unique(el2$pat_id)])))

#study end tables (categories of frequency rates per patient)
n_psa_id$fr_c<-cut(n_psa_id$fr, c(0,0.1, 0.15,0.25,0.3,0.5,1, 5), right=FALSE)
n_psa_id$c_age_dum <- as.numeric(n_psa_id$age_m>=62)

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




