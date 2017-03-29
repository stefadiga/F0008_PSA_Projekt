
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

inc_arzt_sex<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_sex=="f" & data_el$elig_2_1),sex="f", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_sex=="f" & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_sex=="f" & data_el$elig_2_3)), 
          data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_sex=="m" & data_el$elig_2_1),sex="m", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_sex=="m" & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_sex=="m" & data_el$elig_2_3)))

inc_arzt_prac<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Doppelpraxis" & data_el$elig_2_1),type="Doppelpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Doppelpraxis"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Doppelpraxis" & data_el$elig_2_3)), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Einzelpraxis" & data_el$elig_2_1),type="Einzelpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Einzelpraxis"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Einzelpraxis" & data_el$elig_2_3)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Gruppenpraxis" & data_el$elig_2_1),type="Gruppenpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Gruppenpraxis"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Gruppenpraxis" & data_el$elig_2_3)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="" & data_el$elig_2_1),type="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="" & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="" & data_el$elig_2_3)))

inc_arzt_reg<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="CEN"& data_el$elig_2_1),region="CEN", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="CEN"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="CEN"& data_el$elig_2_3)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="IND"& data_el$elig_2_1),region="IND-Ter", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="IND"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="IND"& data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="PERI"& data_el$elig_2_1),region="PERI", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="PERI"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="PERI"& data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="SUB"& data_el$elig_2_1),region="SUB", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="SUB"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="SUB"& data_el$elig_2_3)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="OTHER"& data_el$elig_2_1),region="OTHER", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="OTHER"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="OTHER"& data_el$elig_2_3)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code==""& data_el$elig_2_1),region="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code==""& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code==""& data_el$elig_2_3)))

inc_arzt_emp<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="Angestellt"& data_el$elig_2_1), employ="Angestellt", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="Angestellt"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Angestellt"& data_el$elig_2_3)), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="Selbständig"& data_el$elig_2_1), employ="Selbständig", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="Selbständig"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Selbständig"& data_el$elig_2_3)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="AssistentIn"& data_el$elig_2_1), employ="AssistentIn", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="AssistentIn"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="AssistentIn"& data_el$elig_2_3)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status==""& data_el$elig_2_1), employ="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status==""& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status==""& data_el$elig_2_3)))

data_el$c_doby <- cut(data_el$arzt_doby, c(1940,1950,1960,1970,1980, 1990), right=FALSE, labels=FALSE)

inc_arzt_doby<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==1 & data_el$elig_2_1),doby="[1940,1950)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==1 & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==1 & data_el$elig_2_3)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==2 & data_el$elig_2_1),doby="[1950,1960)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==2 & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==2 & data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==3 & data_el$elig_2_1),doby="[1960,1970)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==3 & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==3 & data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==4 & data_el$elig_2_1),doby="[1970,1980)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==4 & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==4 & data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==5 & data_el$elig_2_1),doby="[1980,1990)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==5 & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==5 & data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, is.na(data_el$c_doby) & data_el$elig_2_1),doby="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, is.na(data_el$c_doby) & data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, is.na(data_el$c_doby) & data_el$elig_2_3)))

data_el$c_age <- cut(data_el$pat_age_start_3, c(55,60,65,70,75), right=FALSE)

inc_pat_age<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[55,60)"& data_el$elig_2_1),age="[55,60)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[55,60)"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[55,60)"& data_el$elig_2_3)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[60,65)"& data_el$elig_2_1),age="[60,65)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[60,65)"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[60,65)"& data_el$elig_2_3)),
                   data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[65,70)"& data_el$elig_2_1),age="[65,70)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[65,70)"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[65,70)"& data_el$elig_2_3)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[70,75)"& data_el$elig_2_1),age="[70,75)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[70,75)"& data_el$elig_2_2), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[70,75)"& data_el$elig_2_3)))

#DB EL2!!!!
el2<-subset(data_el, elig_2_3==1)
el2<-orderBy(~pat_id, el2)

v2 <- c(1, rep(0, length(el2$pat_id)-1))
for (i in 2:length(el2$pat_id))
  if ((el2$pat_id[i] %in% 1:el2$pat_id[i-1]) == FALSE)
    v2[i] <- 1;

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
    ca <- rep(0, (length(quarter_starts)-1));
    # loops for all quarters
    for (i in 1:(length(quarter_starts)-1)) {
      
      # quarter i start and end instants
      Qstart <- rep(as.character(quarter_starts[i]), length(el2$pat_id));
      Qend <- rep(as.character(quarter_starts[i+1]), length(el2$pat_id));
      
      # current ages for patients in the database
      current_age <- floor(2009 - el2$doby + # age at start of study period
                             i / 4);               # years age at current quarter
      
      # indicator function of patients belonging to class j
      grouping <- (current_age >= class_start) & (current_age <= class_end);
      ca[i]<-sum(grouping)
      py[i] <- sum(v2 * grouping * (pmax(0,
                                        as.Date(pmin(end, Qend)) -
                                          as.Date(pmax(start, Qstart)))) / 365, na.rm = T);
      
      events[i] <- sum(grouping * (el2$valid_from_dt >= Qstart) &
                         (el2$valid_from_dt < Qend) &
                         (el2$valid_from_dt >= start) &
                         (el2$valid_from_dt < end), na.rm = T);
      
    }
    
    res <- cbind(res, py, events, ca);
   
  }
  
  return(res);
}



class_age <-  c(55, 60, 65, 70, 75);
date()
prova_2 <- cev_age(el2$fup_end_dt_3, el2$fup_start_dt_3, class_age);
date()

inc_pat_age_c<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,1], events=prova_2[,2], cases=prova_2[,3], class_age=rep(55, 32)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,4], events=prova_2[,5], cases=prova_2[,6], class_age=rep(60, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,7], events=prova_2[,8], cases=prova_2[,9], class_age=rep(65, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,10], events=prova_2[,11], cases=prova_2[,12], class_age=rep(70, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,13], events=prova_2[,14], cases=prova_2[,15], class_age=rep(75, 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,16], events=prova_2[,17], cases=prova_2[,18], class_age=rep(76, 32))) 

             
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
sum(data_el$fup_y_1*v*(data_el$arzt_sex=="f")*data_el$elig_2_1)
sum(inc_arzt_sex$py*(inc_arzt_sex$sex=="m"))
sum(data_el$fup_y_1*v*(data_el$arzt_sex=="m")*data_el$elig_2_1)
#Study End
sum(inc_arzt_sex$py.1 *(inc_arzt_sex$sex=="f"))
sum(data_el$fup_y_2*v*(data_el$arzt_sex=="f")*data_el$elig_2_2)
sum(inc_arzt_sex$py.1*(inc_arzt_sex$sex=="m"))
sum(data_el$fup_y_2*v*(data_el$arzt_sex=="m")*data_el$elig_2_2)
#Mixed
sum(inc_arzt_sex$py.2 *(inc_arzt_sex$sex=="f"))
sum(data_el$fup_y_3*v*(data_el$arzt_sex=="f")*data_el$elig_2_3)
sum(inc_arzt_sex$py.2*(inc_arzt_sex$sex=="m"))
sum(data_el$fup_y_3*v*(data_el$arzt_sex=="m")*data_el$elig_2_3)

#by practice
sum(inc_arzt_prac$py *(inc_arzt_prac$type=="Doppelpraxis"))
sum(data_el$fup_y_1*v*(data_el$practice_type=="Doppelpraxis")*data_el$elig_2_1)
sum(inc_arzt_prac$py*(inc_arzt_prac$type=="Einzelpraxis"))
sum(data_el$fup_y_1*v*(data_el$practice_type=="Einzelpraxis")*data_el$elig_2_1)
sum(inc_arzt_prac$py*(inc_arzt_prac$type=="Gruppenpraxis"))
sum(data_el$fup_y_1*v*(data_el$practice_type=="Gruppenpraxis")*data_el$elig_2_1)
sum(inc_arzt_prac$py*(inc_arzt_prac$type==""))
sum(data_el$fup_y_1*v*(data_el$practice_type=="")*data_el$elig_2_1)

#Study End
sum(inc_arzt_prac$py.1 *(inc_arzt_prac$type=="Doppelpraxis"))
sum(data_el$fup_y_2*v*(data_el$practice_type=="Doppelpraxis")*data_el$elig_2_2)
sum(inc_arzt_prac$py.1*(inc_arzt_prac$type=="Einzelpraxis"))
sum(data_el$fup_y_2*v*(data_el$practice_type=="Einzelpraxis")*data_el$elig_2_2)
sum(inc_arzt_prac$py.1*(inc_arzt_prac$type=="Gruppenpraxis"))
sum(data_el$fup_y_2*v*(data_el$practice_type=="Gruppenpraxis")*data_el$elig_2_2)
sum(inc_arzt_prac$py.1*(inc_arzt_prac$type==""))
sum(data_el$fup_y_2*v*(data_el$practice_type=="")*data_el$elig_2_2)

#Mixed
sum(inc_arzt_prac$py.2 *(inc_arzt_prac$type=="Doppelpraxis"))
sum(data_el$fup_y_3*v*(data_el$practice_type=="Doppelpraxis")*data_el$elig_2_3)
sum(inc_arzt_prac$py.2*(inc_arzt_prac$type=="Einzelpraxis"))
sum(data_el$fup_y_3*v*(data_el$practice_type=="Einzelpraxis")*data_el$elig_2_3)
sum(inc_arzt_prac$py.2*(inc_arzt_prac$type=="Gruppenpraxis"))
sum(data_el$fup_y_3*v*(data_el$practice_type=="Gruppenpraxis")*data_el$elig_2_3)
sum(inc_arzt_prac$py.2*(inc_arzt_prac$type==""))
sum(data_el$fup_y_3*v*(data_el$practice_type=="")*data_el$elig_2_3)

#Region Type
sum(inc_arzt_reg$py.2 *(inc_arzt_reg$region=="CEN"))
sum(data_el$fup_y_3*v*(data_el$arzt_region_code=="CEN")*data_el$elig_2_3)
sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="IND-Ter"))
sum(data_el$fup_y_3*v*(data_el$arzt_region_code=="IND")*data_el$elig_2_3)
sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="PERI"))
sum(data_el$fup_y_3*v*(data_el$arzt_region_code=="PERI")*data_el$elig_2_3)
sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="SUB"))
sum(data_el$fup_y_3*v*(data_el$arzt_region_code=="SUB")*data_el$elig_2_3)
sum(inc_arzt_reg$py.2*(inc_arzt_reg$region==""))
sum(data_el$fup_y_3*v*(data_el$arzt_region_code=="")*data_el$elig_2_3)

#Employ status
#Observed
sum(inc_arzt_emp$py *(inc_arzt_emp$employ=="Angestellt"))
sum(data_el$fup_y_1*v*(data_el$employ_status=="Angestellt")*data_el$elig_2_1)
sum(inc_arzt_emp$py*(inc_arzt_emp$employ=="Selbständig"))
sum(data_el$fup_y_1*v*(data_el$employ_status=="Selbständig")*data_el$elig_2_1)
#Study End
sum(inc_arzt_emp$py.1 *(inc_arzt_emp$employ=="Angestellt"))
sum(data_el$fup_y_2*v*(data_el$employ_status=="Angestellt")*data_el$elig_2_2)
sum(inc_arzt_emp$py.1*(inc_arzt_emp$employ=="Selbständig"))
sum(data_el$fup_y_2*v*(data_el$employ_status=="Selbständig")*data_el$elig_2_2)
#Mixed
sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="Angestellt"))
sum(data_el$fup_y_3*v*(data_el$employ_status=="Angestellt")*data_el$elig_2_3)
sum(inc_arzt_emp$py.2*(inc_arzt_emp$employ=="Selbständig"))
sum(data_el$fup_y_3*v*(data_el$employ_status=="Selbständig")*data_el$elig_2_3)
sum(inc_arzt_emp$py.2*(inc_arzt_emp$employ=="AssistentIn"))
sum(data_el$fup_y_3*v*(data_el$employ_status=="AssistentIn")*data_el$elig_2_3)

#Arzt DOBY
#Observed
sum(inc_arzt_doby$py *(inc_arzt_doby$doby=="[1940,1950)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==1)*data_el$elig_2_1, na.rm=T)
sum(inc_arzt_doby$py*(inc_arzt_doby$doby=="[1950,1960)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==2)*data_el$elig_2_1, na.rm=T)
sum(inc_arzt_doby$py*(inc_arzt_doby$doby=="[1960,1970)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==3)*data_el$elig_2_1, na.rm=T)
sum(inc_arzt_doby$py*(inc_arzt_doby$doby=="[1970,1980)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==4)*data_el$elig_2_1, na.rm=T)
sum(inc_arzt_doby$py*(inc_arzt_doby$doby=="[1980,1990)"))
sum(data_el$fup_y_1*v*(data_el$c_doby==4)*data_el$elig_2_1, na.rm=T)

#Study End
sum(inc_arzt_doby$py.1 *(inc_arzt_doby$doby=="[1940,1950)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==1)*data_el$elig_2_2, na.rm=T)
sum(inc_arzt_doby$py.1*(inc_arzt_doby$doby=="[1950,1960)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==2)*data_el$elig_2_2, na.rm=T)
sum(inc_arzt_doby$py.1*(inc_arzt_doby$doby=="[1960,1970)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==3)*data_el$elig_2_2, na.rm=T)
sum(inc_arzt_doby$py.1*(inc_arzt_doby$doby=="[1970,1980)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==4)*data_el$elig_2_2, na.rm=T)
sum(inc_arzt_doby$py.1*(inc_arzt_doby$doby=="[1980,1990)"))
sum(data_el$fup_y_2*v*(data_el$c_doby==5)*data_el$elig_2_2, na.rm=T)

#Mixed
sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1940,1950)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==1)*data_el$elig_2_3, na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby=="[1950,1960)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==2)*data_el$elig_2_3, na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby=="[1960,1970)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==3)*data_el$elig_2_3, na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby=="[1970,1980)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==4)*data_el$elig_2_3, na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby=="[1980,1990)"))
sum(data_el$fup_y_3*v*(data_el$c_doby==5)*data_el$elig_2_3, na.rm=T)
sum(inc_arzt_doby$py.2*(inc_arzt_doby$doby==""))
sum(data_el$fup_y_3*v*(is.na(data_el$c_doby))*data_el$elig_2_3, na.rm=T)

#Pat Age
#Observed
sum(inc_pat_age$py *(inc_pat_age$age=="[55,60)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[55,60)")*data_el$elig_2_1, na.rm=T)
sum(inc_pat_age$py*(inc_pat_age$age=="[60,65)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[60,65)")*data_el$elig_2_1, na.rm=T)
sum(inc_pat_age$py*(inc_pat_age$age=="[65,70)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[65,70)")*data_el$elig_2_1, na.rm=T)
sum(inc_pat_age$py*(inc_pat_age$age=="[70,75)"))
sum(data_el$fup_y_1*v*(data_el$c_age=="[70,75)")*data_el$elig_2_1, na.rm=T)
#Study End
sum(inc_pat_age$py.1 *(inc_pat_age$age=="[55,60)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[55,60)")*data_el$elig_2_2, na.rm=T)
sum(inc_pat_age$py.1*(inc_pat_age$age=="[60,65)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[60,65)")*data_el$elig_2_2, na.rm=T)
sum(inc_pat_age$py.1*(inc_pat_age$age=="[65,70)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[65,70)")*data_el$elig_2_2, na.rm=T)
sum(inc_pat_age$py.1*(inc_pat_age$age=="[70,75)"))
sum(data_el$fup_y_2*v*(data_el$c_age=="[70,75)")*data_el$elig_2_2, na.rm=T)
#Mixed
sum(inc_pat_age$py.2 *(inc_pat_age$age=="[55,60)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[55,60)")*data_el$elig_2_3, na.rm=T)
sum(inc_pat_age$py.2*(inc_pat_age$age=="[60,65)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[60,65)")*data_el$elig_2_3, na.rm=T)
sum(inc_pat_age$py.2*(inc_pat_age$age=="[65,70)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[65,70)")*data_el$elig_2_3, na.rm=T)
sum(inc_pat_age$py.2*(inc_pat_age$age=="[70,75)"))
sum(data_el$fup_y_3*v*(data_el$c_age=="[70,75)")*data_el$elig_2_3, na.rm=T)

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
graphics.off()
windows(width=15, height=10)
i1<-qplot(as.Date(inc$dt), inc$events/inc$py,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") +
     geom_errorbar(aes(ymin=exp(log(inc$events/inc$py)- qnorm(0.975)*sqrt(inc$events/inc$py)), ymax=exp(log(inc$events/inc$py)+ qnorm(0.975)*sqrt(inc$events/inc$py))))

i2<-qplot(as.Date(inc$dt), inc$events.1/inc$py.1,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End")+
     geom_errorbar(aes(ymin=exp(log(inc$events.1/inc$py.1)- qnorm(0.975)*sqrt(inc$events.1/inc$py.1)), ymax=exp(log(inc$events.1/inc$py.1)+ qnorm(0.975)*sqrt(inc$events.1/inc$py.1))))

i3<-qplot(as.Date(inc$dt), inc$events.2/inc$py.2,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed")+
    geom_errorbar(aes(ymin=exp(log(inc$events.2/inc$py.2)- qnorm(0.975)*sqrt(inc$events.2/inc$py.2)), ymax=exp(log(inc$events.2/inc$py.2)+ qnorm(0.975)*sqrt(inc$events.2/inc$py.2))))

grid.arrange(i1, i2, i3, ncol=2)

savePlot("plot_a00xa1_ir.jpg",type="jpg")
save.image()


#by sex
#i1<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events/inc_arzt_sex$py), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Sex")

#i2<-qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events.1/inc_arzt_sex$py.1), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Sex")

graphics.off()
windows(width=15, height=10)

qplot(as.Date(inc_arzt_sex$dt), as.numeric(inc_arzt_sex$events.2/inc_arzt_sex$py.2), color=inc_arzt_sex$sex, geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Sex")
savePlot("plot_a00xa2_ir.jpg", type="jpg")
save.image()

#grid.arrange(i1, i2, i3, ncol=2)

#by practice

#i1<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events/inc_arzt_prac$py), color=inc_arzt_prac$type, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Practice type")

#i2<-qplot(as.Date(inc_arzt_prac$dt), as.numeric(inc_arzt_prac$events.1/inc_arzt_prac$py.1), color=inc_arzt_prac$type, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Practice type")
graphics.off()
windows(width=15, height=10)
qplot(as.Date(inc_arzt_prac$dt[inc_arzt_prac$type!=""]), as.numeric(inc_arzt_prac$events.2[inc_arzt_prac$type!=""]/inc_arzt_prac$py.2[inc_arzt_prac$type!=""]), color=inc_arzt_prac$type[inc_arzt_prac$type!=""], geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Practice type")
savePlot("plot_a00xa4_ir.jpg",type="jpg")
save.image()

#grid.arrange(i1, i2, i3, ncol=2)

#by region
#to do add AGR?
#i1<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events/inc_arzt_reg$py), color=inc_arzt_reg$region, geom=c("point", "smooth"),
 #         ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Region")

#i2<-qplot(as.Date(inc_arzt_reg$dt), as.numeric(inc_arzt_reg$events.1/inc_arzt_reg$py.1), color=inc_arzt_reg$region, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Region")
graphics.off()
windows(width=15, height=10)

qplot(as.Date(inc_arzt_reg$dt[inc_arzt_reg$region!=""]), as.numeric(inc_arzt_reg$events.2[inc_arzt_reg$region!=""]/inc_arzt_reg$py.2[inc_arzt_reg$region!=""]), color=inc_arzt_reg$region[inc_arzt_reg$region!=""], geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Region")
savePlot("plot_a00xa3_ir.jpg",type="jpg")
save.image()

#grid.arrange(i1, i2, i3, ncol=2)

#by employ status
#i1<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events/inc_arzt_emp$py), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Employ status")

#i2<-qplot(as.Date(inc_arzt_emp$dt), as.numeric(inc_arzt_emp$events.1/inc_arzt_emp$py.1), color=inc_arzt_emp$employ, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Employ status")
graphics.off()
windows(width=15, height=10)

qplot(as.Date(inc_arzt_emp$dt[inc_arzt_emp$employ!=""]), as.numeric(inc_arzt_emp$events.2[inc_arzt_emp$employ!=""]/inc_arzt_emp$py.2[inc_arzt_emp$employ!=""]), color=inc_arzt_emp$employ[inc_arzt_emp$employ!=""], geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Employ status")
savePlot("plot_a00xa5_ir.jpg",type="jpg")
save.image()

#grid.arrange(i1, i2, i3, ncol=2)

#by arzt doby
#i1<-qplot(as.Date(inc_arzt_doby$dt), as.numeric(inc_arzt_doby$events/inc_arzt_doby$py), color=inc_arzt_doby$doby, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Arzt DOBY")

#i2<-qplot(as.Date(inc_arzt_doby$dt), as.numeric(inc_arzt_doby$events.1/inc_arzt_doby$py.1), color=inc_arzt_doby$doby, geom=c("point", "smooth"),
#          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Arzt DOBY")
graphics.off()
windows(width=15, height=10)

qplot(as.Date(inc_arzt_doby$dt[inc_arzt_doby$doby!=""]), as.numeric(inc_arzt_doby$events.2[inc_arzt_doby$doby!=""]/inc_arzt_doby$py.2[inc_arzt_doby$doby!=""]), color=inc_arzt_doby$doby[inc_arzt_doby$doby!=""], geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Arzt DOBY")
savePlot("plot_a00xa6_ir.jpg",type="jpg")
save.image()

#grid.arrange(i1, i2, i3, ncol=2)

#by age

#i1<-qplot(as.Date(inc_pat_age_c$dt[inc_pat_age_c$class_age>=55]), as.numeric(inc_pat_age_c$events[inc_pat_age_c$class_age>=55]/inc_pat_age_c$py[inc_pat_age_c$class_age>=55]), color=inc_pat_age_c$class_age[inc_pat_age_c$class_age>=55], geom=c("point", "smooth"),
 #         ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") + labs(color = "Pat Age")

#i2<-qplot(as.Date(inc_pat_age$dt), as.numeric(inc_pat_age$events.1/inc_pat_age$py.1), color=inc_pat_age$age, geom=c("point", "smooth"),
      #    ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End") + labs(color = "Pat Age")
graphics.off()
windows(width=15, height=10)

qplot(as.Date(inc_pat_age$dt), as.numeric(inc_pat_age$events.2/inc_pat_age$py.2), color=inc_pat_age$age, geom=c("point", "smooth"),
         ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed") + labs(color = "Pat Age")
savePlot("plot_a00xa7.0_ir.jpg",type="jpg")
save.image()

#grid.arrange(i1, i2, i3, ncol=2)

inc_pat_age_c$inc<-as.numeric(inc_pat_age_c$events/inc_pat_age_c$py)
inc_pat_age_c <- subset(inc_pat_age_c, inc_pat_age_c$class_age >=60 & inc_pat_age_c$class_age < 76)

graphics.off()
windows(width=15, height=10)

qplot(data=inc_pat_age_c, x=as.Date(dt), y=inc, geom=c("point","smooth"),ylim=c(0,0.3),
      xlab="dt", main="Mixed follow up") + 
   facet_grid(.~as.character(class_age))

savePlot("plot_a00xa7_ir.jpg",type="jpg")
save.image()

inc_pat_age$inc<-as.numeric(inc_pat_age$events.2/inc_pat_age$py.2)

i2<-qplot(data=inc_pat_age, x=as.Date(dt), y=inc, geom=c("point","smooth"),ylim=c(0,0.3),
      xlab="dt", main="Mixed follow up") + 
    facet_grid(.~age)
grid.arrange(i1, i2, ncol=2)

#Plot Frequency
qplot(as.Date(inc$dt),inc$events, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests")

#PSA

med<-tapply(el2$psa, cut(as.Date(el2$valid_from_dt), as.Date(quarter_starts)), median)
setq<-tapply(el2$psa, cut(as.Date(el2$valid_from_dt), as.Date(quarter_starts)), quantile, 0.75)
novq<-tapply(el2$psa, cut(as.Date(el2$valid_from_dt), as.Date(quarter_starts)), quantile, 0.90)


quant<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq)))          

graphics.off()
windows(width=15, height=10)

qplot(data=quant, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
   ylim(c(0,8)), xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("a00xa0_3_1_psa_quantile.jpg",type="jpg")
save.image()

#by sex
med_f<-tapply(el2$psa[el2$arzt_sex=="f"], cut(as.Date(el2$valid_from_dt[el2$arzt_sex=="f"]), as.Date(quarter_starts)), median) 
med_m<-tapply(el2$psa[el2$arzt_sex=="m"], cut(as.Date(el2$valid_from_dt[el2$arzt_sex=="m"]), as.Date(quarter_starts)), median) 
setq_f<-tapply(el2$psa[el2$arzt_sex=="f"], cut(as.Date(el2$valid_from_dt[el2$arzt_sex=="f"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_m<-tapply(el2$psa[el2$arzt_sex=="m"], cut(as.Date(el2$valid_from_dt[el2$arzt_sex=="m"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_f<-tapply(el2$psa[el2$arzt_sex=="f"], cut(as.Date(el2$valid_from_dt[el2$arzt_sex=="f"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_m<-tapply(el2$psa[el2$arzt_sex=="m"], cut(as.Date(el2$valid_from_dt[el2$arzt_sex=="m"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_sex<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_f), sex="f"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_m), sex="m"),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_f), sex="f"),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_f), sex="m"),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_f), sex="f"),          
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_m), sex="m"))

graphics.off()
windows(width=15, height=10)

qplot(data=quant_sex, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ sex)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("plot_a002a2_3_psa_quantile",type="jpg")
save.image()

#region
med_C<-tapply(el2$psa[el2$arzt_region_code=="CEN"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="CEN"]), as.Date(quarter_starts)), median) 
med_I<-tapply(el2$psa[el2$arzt_region_code=="IND"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="IND"]), as.Date(quarter_starts)), median) 
med_O<-tapply(el2$psa[el2$arzt_region_code=="OTHER"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), median) 
med_P<-tapply(el2$psa[el2$arzt_region_code=="PERI"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="PERI"]), as.Date(quarter_starts)), median) 
med_S<-tapply(el2$psa[el2$arzt_region_code=="SUB"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="SUB"]), as.Date(quarter_starts)), median) 
setq_C<-tapply(el2$psa[el2$arzt_region_code=="CEN"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="CEN"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_I<-tapply(el2$psa[el2$arzt_region_code=="IND"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="IND"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_O<-tapply(el2$psa[el2$arzt_region_code=="OTHER"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_P<-tapply(el2$psa[el2$arzt_region_code=="PERI"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="PERI"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_S<-tapply(el2$psa[el2$arzt_region_code=="SUB"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="SUB"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_C<-tapply(el2$psa[el2$arzt_region_code=="CEN"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="CEN"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_I<-tapply(el2$psa[el2$arzt_region_code=="IND"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="IND"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_O<-tapply(el2$psa[el2$arzt_region_code=="OTHER"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_P<-tapply(el2$psa[el2$arzt_region_code=="PERI"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="PERI"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_S<-tapply(el2$psa[el2$arzt_region_code=="SUB"], cut(as.Date(el2$valid_from_dt[el2$arzt_region_code=="SUB"]), as.Date(quarter_starts)), quantile, 0.90) 

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

graphics.off()
windows(width=15, height=10)

qplot(data=quant_reg, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up", ylim=c(0,15)) + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ region)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("plot_a002a3_3_psa_quantile",type="jpg")
save.image()

#by employ status
med_a<-tapply(el2$psa[el2$employ_status=="Angestellt"], cut(as.Date(el2$valid_from_dt[el2$employ_status=="Angestellt"]), as.Date(quarter_starts)), median) 
med_s<-tapply(el2$psa[el2$employ_status=="Selbständig"], cut(as.Date(el2$valid_from_dt[el2$employ_status=="Selbständig"]), as.Date(quarter_starts)), median) 
setq_a<-tapply(el2$psa[el2$employ_status=="Angestellt"], cut(as.Date(el2$valid_from_dt[el2$employ_status=="Angestellt"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_s<-tapply(el2$psa[el2$employ_status=="Selbständig"], cut(as.Date(el2$valid_from_dt[el2$employ_status=="Selbständig"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_a<-tapply(el2$psa[el2$employ_status=="Angestellt"], cut(as.Date(el2$valid_from_dt[el2$employ_status=="Angestellt"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_s<-tapply(el2$psa[el2$employ_status=="Selbständig"], cut(as.Date(el2$valid_from_dt[el2$employ_status=="Selbständig"]), as.Date(quarter_starts)),quantile, 0.90) 

quant_emp<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a), employ="Angestellt"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_s), employ="Selbständig"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a), employ="Angestellt"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_s), employ="Selbständig"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a), employ="Angestellt"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_s), employ="Selbständig"))

graphics.off()
windows(width=15, height=10)

qplot(data=quant_emp, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ employ)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("plot_a002a5_3_psa_quantile.jpg",type="jpg")
save.image()

#by praxis
med_d<-tapply(el2$psa[el2$practice_type=="Doppelpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), median) 
med_e<-tapply(el2$psa[el2$practice_type=="Einzelpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), median) 
med_g<-tapply(el2$psa[el2$practice_type=="Gruppenpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), median) 
setq_d<-tapply(el2$psa[el2$practice_type=="Doppelpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_e<-tapply(el2$psa[el2$practice_type=="Einzelpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_g<-tapply(el2$psa[el2$practice_type=="Gruppenpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_d<-tapply(el2$psa[el2$practice_type=="Doppelpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_e<-tapply(el2$psa[el2$practice_type=="Einzelpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_g<-tapply(el2$psa[el2$practice_type=="Gruppenpraxis"], cut(as.Date(el2$valid_from_dt[el2$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_pra<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_g), practice="Gruppenpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_g), practice="Gruppenpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_g), practice="Gruppenpraxis"))

graphics.off()
windows(width=15, height=10)


qplot(data=quant_pra, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ practice)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("plot_a002a4_3_psa_quantile.jpg",type="jpg")
save.image()

#by arzt age
med_a1<-tapply(el2$psa[el2$c_doby==1], cut(as.Date(el2$valid_from_dt[el2$c_doby==1]), as.Date(quarter_starts)), median) 
med_a2<-tapply(el2$psa[el2$c_doby==2], cut(as.Date(el2$valid_from_dt[el2$c_doby==2]), as.Date(quarter_starts)), median) 
med_a3<-tapply(el2$psa[el2$c_doby==3], cut(as.Date(el2$valid_from_dt[el2$c_doby==3]), as.Date(quarter_starts)), median) 
med_a4<-tapply(el2$psa[el2$c_doby==4], cut(as.Date(el2$valid_from_dt[el2$c_doby==4]), as.Date(quarter_starts)), median) 
med_a5<-tapply(el2$psa[el2$c_doby==5], cut(as.Date(el2$valid_from_dt[el2$c_doby==5]), as.Date(quarter_starts)), median) 

setq_a1<-tapply(el2$psa[el2$c_doby==1], cut(as.Date(el2$valid_from_dt[el2$c_doby==1]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a2<-tapply(el2$psa[el2$c_doby==2], cut(as.Date(el2$valid_from_dt[el2$c_doby==2]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a3<-tapply(el2$psa[el2$c_doby==3], cut(as.Date(el2$valid_from_dt[el2$c_doby==3]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a4<-tapply(el2$psa[el2$c_doby==4], cut(as.Date(el2$valid_from_dt[el2$c_doby==4]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a5<-tapply(el2$psa[el2$c_doby==5], cut(as.Date(el2$valid_from_dt[el2$c_doby==5]), as.Date(quarter_starts)), quantile, 0.75) 

novq_a1<-tapply(el2$psa[el2$c_doby==1], cut(as.Date(el2$valid_from_dt[el2$c_doby==1]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a2<-tapply(el2$psa[el2$c_doby==2], cut(as.Date(el2$valid_from_dt[el2$c_doby==2]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a3<-tapply(el2$psa[el2$c_doby==3], cut(as.Date(el2$valid_from_dt[el2$c_doby==3]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a4<-tapply(el2$psa[el2$c_doby==4], cut(as.Date(el2$valid_from_dt[el2$c_doby==4]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a5<-tapply(el2$psa[el2$c_doby==5], cut(as.Date(el2$valid_from_dt[el2$c_doby==5]), as.Date(quarter_starts)), quantile, 0.90) 

quant_adob<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a1), arzt_doby="1940"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a2), arzt_doby="1950"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a3), arzt_doby="1960"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a4), arzt_doby="1970"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a5), arzt_doby="1980"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a1), arzt_doby="1940"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a2), arzt_doby="1950"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a3), arzt_doby="1960"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a4), arzt_doby="1970"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a5), arzt_doby="1980"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a1), arzt_doby="1940"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a2), arzt_doby="1950"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a3), arzt_doby="1960"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a4), arzt_doby="1970"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a5), arzt_doby="1980"))
graphics.off()
windows(width=15, height=10)


qplot(data=quant_adob, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ arzt_doby)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("plot_a002a6_3_psa_quantile.jpg",type="jpg")
save.image()

#by pat age
med_p1<-tapply(el2$psa[el2$c_age=="[55,60)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[55,60)"]), as.Date(quarter_starts)), median) 
med_p2<-tapply(el2$psa[el2$c_age=="[60,65)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[60,65)"]), as.Date(quarter_starts)), median) 
med_p3<-tapply(el2$psa[el2$c_age=="[65,70)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[65,70)"]), as.Date(quarter_starts)), median) 
med_p4<-tapply(el2$psa[el2$c_age=="[70,75)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[70,75)"]), as.Date(quarter_starts)), median) 
setq_p1<-tapply(el2$psa[el2$c_age=="[55,60)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[55,60)"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_p2<-tapply(el2$psa[el2$c_age=="[60,65)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[60,65)"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_p3<-tapply(el2$psa[el2$c_age=="[65,70)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[65,70)"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_p4<-tapply(el2$psa[el2$c_age=="[70,75)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[70,75)"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_p1<-tapply(el2$psa[el2$c_age=="[55,60)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[55,60)"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_p2<-tapply(el2$psa[el2$c_age=="[60,65)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[60,65)"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_p3<-tapply(el2$psa[el2$c_age=="[65,70)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[65,70)"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_p4<-tapply(el2$psa[el2$c_age=="[70,75)"], cut(as.Date(el2$valid_from_dt[el2$c_age=="[70,75)"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_pdob<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_p1), p_age="[55,60)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_p2), p_age="[60,65)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_p3), p_age="[65,70)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_p4), p_age="[70,75)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_p1), p_age="[55,60)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_p2), p_age="[60,65)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_p3), p_age="[65,70)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_p4), p_age="[70,75)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_p1), p_age="[55,60)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_p2), p_age="[60,65)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_p3), p_age="[65,70)"),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_p4), p_age="[70,75)"))

graphics.off()
windows(width=15, height=10)


qplot(data=quant_pdob, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ p_age)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))

savePlot("plot_a002a7_3_psa_quantile.jpg",type="jpg")
save.image()

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

el<-subset(data_el, elig_2_1==1)
el1<-subset(data_el, elig_2_2==1)
el2<-subset(data_el, elig_2_3==1)


count_psa<-rbind(data.frame(n_psa=0, n_pat=n_distinct(el$pat_id)-sum(inc$events)), as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_1 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1)]), data_el$pat_id[data_el$elig_2_1 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName = "n_pat"))

count_psa_2<-rbind(data.frame(n_psa=0, n_pat=n_distinct(el1$pat_id)-sum(inc$events.1)), as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_2 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2)]), data_el$pat_id[data_el$elig_2_2 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName = "n_pat"))

count_psa_3<-rbind(data.frame(n_psa=0, n_pat=n_distinct(el2$pat_id)-sum(inc$events.2)), as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)]), data_el$pat_id[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName ="n_pat"))


tab1<-data.frame(method=c("Observed", "Study Period", "Mixed"), n_pat=c(n_distinct(el$pat_id), n_distinct(el1$pat_id), n_distinct(el2$pat_id)), psa=c(sum(as.numeric(count_psa$n_psa)*count_psa$n_pat),sum(as.numeric(count_psa_2$n_psa)*count_psa_2$n_pat), sum(as.numeric(count_psa_3$n_psa)*count_psa_3$n_pat)), py=c(sum(data_el$fup_y_1*data_el$elig_2_1*v), sum(data_el$fup_y_2*data_el$elig_2_2*v), sum(data_el$fup_y_3*data_el$elig_2_3*v)), n_arzt=c(n_distinct(el$arzt_id), n_distinct(el1$arzt_id),n_distinct(el2$arzt_id)))
tab1$ir<-as.numeric(tab1$psa)/as.numeric(tab1$py)

#tab2<-data.frame(method=c("Observed", "Study Period", "Mixed"), n_pat=c(n_distinct(el$pat_id[el$arzt_sex=="m"]),n_distinct(el1$pat_id[el1$arzt_sex=="m"]), n_distinct(el2$pat_id[el2$arzt_sex=="m"])))   
tab2<-data.frame(method=c("Mixed", "Mixed"), arzt_sex=c("f", "m"), n_pat=c(n_distinct(el2$pat_id[el2$arzt_sex=="f"]),n_distinct(el2$pat_id[el2$arzt_sex=="m"])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$arzt_sex=="f"]),n_distinct(el2$arzt_id[el2$arzt_sex=="m"])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$arzt_sex=="f" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_sex=="m" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_sex$py.2 *(inc_arzt_sex$sex=="f")), sum(inc_arzt_sex$py.2 *(inc_arzt_sex$sex=="m"))))   
tab2$ir<-as.numeric(tab2$n_psa)/as.numeric(tab2$py)

tab3<-data.frame(method=c("Mixed", "Mixed", "Mixed", "Mixed", "Mixed", "Mixed"), region=c("CEN", "IND-Ter", "SUB", "PERI", "OTHER", "NA"), n_pat=c(n_distinct(el2$pat_id[el2$arzt_region_code=="CEN"]), n_distinct(el2$pat_id[el2$arzt_region_code=="IND"]), n_distinct(el2$pat_id[el2$arzt_region_code=="SUB"]), n_distinct(el2$pat_id[el2$arzt_region_code=="PERI"]), n_distinct(el2$pat_id[el2$arzt_region_code=="OTHER"]), n_distinct(el2$pat_id[el2$arzt_region_code==""])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$arzt_region_code=="CEN"]), n_distinct(el2$arzt_id[el2$arzt_region_code=="IND"]), n_distinct(el2$arzt_id[el2$arzt_region_code=="SUB"]),  n_distinct(el2$arzt_id[el2$arzt_region_code=="PERI"]), n_distinct(el2$arzt_id[el2$arzt_region_code=="OTHER"]), n_distinct(el2$arzt_id[el2$arzt_region_code==""])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$arzt_region_code=="CEN" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="IND" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="SUB" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="PERI" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="OTHER" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_reg$py.2 *(inc_arzt_reg$region=="CEN")), sum(inc_arzt_reg$py.2 *(inc_arzt_reg$region=="IND-Ter")), sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="SUB")), sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="PERI")),sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="OTHER")), sum(inc_arzt_reg$py.2*(inc_arzt_reg$region==""))))   
tab3$ir<-as.numeric(tab3$n_psa)/as.numeric(tab3$py)

tab4<-data.frame(method=c("Mixed", "Mixed", "Mixed", "Mixed"), practice_type=c("Einzelpraxis", "Gruppenpraxis", "Doppelpraxis", "NA"), n_pat=c(n_distinct(el2$pat_id[el2$practice_type=="Einzelpraxis"]), n_distinct(el2$pat_id[el2$practice_type=="Gruppenpraxis"]), n_distinct(el2$pat_id[el2$practice_type=="Doppelpraxis"]), n_distinct(el2$pat_id[el2$practice_type==""])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$practice_type=="Einzelpraxis"]), n_distinct(el2$arzt_id[el2$practice_type=="Gruppenpraxis"]), n_distinct(el2$arzt_id[el2$practice_type=="Doppelpraxis"]),  n_distinct(el2$arzt_id[el2$practice_type==""])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$practice_type=="Einzelpraxis" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$practice_type=="Gruppenpraxis" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$practice_type=="Doppelpraxis" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$practice_type=="" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_prac$py.2 *(inc_arzt_prac$type=="Einzelpraxis")), sum(inc_arzt_prac$py.2 *(inc_arzt_prac$type=="Gruppenpraxis")), sum(inc_arzt_prac$py.2*(inc_arzt_prac$type=="Doppelpraxis")), sum(inc_arzt_prac$py.2*(inc_arzt_prac$type==""))))   
tab4$ir<-as.numeric(tab4$n_psa)/as.numeric(tab4$py)

tab5<-data.frame(method=c("Mixed", "Mixed", "Mixed", "Mixed"), employ_status=c("Angestellt", "Selbständig", "Assistentin", "NA"), n_pat=c(n_distinct(el2$pat_id[el2$employ_status=="Angestellt"]),n_distinct(el2$pat_id[el2$employ_status=="Selbständig"]), n_distinct(el2$pat_id[el2$employ_status=="AssistentIn"]), n_distinct(el2$pat_id[el2$employ_status==""])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$employ_status=="Angestellt"]), n_distinct(el2$arzt_id[el2$employ_status=="Selbständig"]), n_distinct(el2$arzt_id[el2$employ_status=="AssistentIn"]), n_distinct(el2$arzt_id[el2$employ_status==""])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$employ_status=="Angestellt" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$employ_status=="Selbständig" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$employ_status=="AssistentIn" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$employ_status=="" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="Angestellt")), sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="Selbständig")), sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="AssistentIn")), sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ==""))))   
tab5$ir<-as.numeric(tab5$n_psa)/as.numeric(tab5$py)

el2$c_doby <- cut(el2$arzt_doby, c(1940,1950,1960,1970,1980, 1990), right=FALSE, labels=FALSE)

tab6<-data.frame(method=c("Mixed", "Mixed", "Mixed", "Mixed", "Mixed", "Mixed"), arzt_doby=c("[1940,1950)", "[1950,1960)", "[1960,1970)", "[1970,1980)", "[1980,1990)", "NA"), n_pat=c(as.vector(table(el2$c_doby*v2))[2:6], n_distinct(el2$pat_id[is.na(el2$c_doby)])),
                 n_arzt=c(as.vector(tapply(el2$arzt_id, el2$c_doby, n_distinct)), n_distinct(el2$arzt_id[is.na(el2$c_doby)])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$c_doby==1 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==2 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==3 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==4 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==5 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[is.na(el2$c_doby) & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1940,1950)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1950,1960)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1960,1970)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1970,1980)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1980,1990)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby==""))))   
tab6$ir<-as.numeric(tab6$n_psa)/as.numeric(tab6$py)

#patient age
el2$c_age <- cut(el2$pat_age_start_3, c(55,60,65,70,75), right=FALSE)


tab7<-data.frame(method=c("Mixed", "Mixed", "Mixed", "Mixed"), pat_age=c("[55,60)", "[60,65)", "[65,70)","[70,75)"), n_pat=c(as.vector(tapply(el2$pat_id, el2$c_age, n_distinct))),
                 n_arzt=c(as.vector(tapply(el2$arzt_id, el2$c_age, n_distinct))), 
                 n_psa=c(sum(!is.na(el2$psa[el2$c_age=="[55,60)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[60,65)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[65,70)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[70,75)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_pat_age$py.2 *(inc_pat_age$age=="[55,60)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[60,65)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[65,70)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[70,75)"))))   
tab7$ir<-as.numeric(tab7$n_psa)/as.numeric(tab7$py)

age_m<-(el2$pat_age_start_3*as.numeric(el2$fup_y_3)*v2)
#weighted age (age with weight the contribution to the study)
sum(age_m)/sum(as.numeric(el2$fup_y_3)*v2)


sum(inc_pat_age_c$py[inc_pat_age_c$class_age==60])
sum(inc_pat_age_c$py[inc_pat_age_c$class_age==65])
sum(inc_pat_age_c$py[inc_pat_age_c$class_age==70])
sum(inc_pat_age_c$py[inc_pat_age_c$class_age==75])

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


class_age_1 <-  c(55, 56, 60, 61, 65, 66, 70, 71, 75, 76)
prova_3 <- cev_age(el2$fup_end_dt_3, el2$fup_start_dt_3, class_age_1)

inc_pat_age_c1<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,4], events=prova_3[,5], cases=prova_3[,6], class_age=rep(55, 32)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,10], events=prova_3[,11], cases=prova_3[,12], class_age=rep(60, 32)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,16], events=prova_3[,17], cases=prova_3[,18], class_age=rep(65, 32)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,22], events=prova_3[,23], cases=prova_3[,24], class_age=rep(70, 32)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,28], events=prova_3[,29], cases=prova_3[,30], class_age=rep(75, 32))) 

c55<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c60<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c65<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c70<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c75<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
