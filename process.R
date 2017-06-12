
# -----TO DO-----------------------------------------------
# Mixed model PSA Value
# ----------------------------------------------------------------

#adding sum of measurement by patient

n_psa_id<-subset(el2, as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3))
n_psa_id$dt<-cut(as.Date(n_psa_id$valid_from_dt), as.Date(quarter_starts))

n_psa_id2<-data.frame(pat_id=unique(n_psa_id$pat_id), n_psa=tapply(n_psa_id$psa, n_psa_id$pat_id, length))


#chronic diseases

#group_by(chr2, pat_id, dt)
#summarise(group_by(chr2, pat_id, dt), min(cd_count), max(cd_count))

chr2<-distinct(chr2, pat_id, dt, .keep_all = TRUE)

n_psa_id_c<-merge(n_psa_id, chr2, all.x=TRUE)

#db chronic disease for PSA
#n_psa_id_c<-filter(group_by(n_psa_id_c, pat_id, valid_from_dt,psa), psa == max(psa))

a00<-subset(n_psa_id2, n_psa>=6)

a001a0_3 <- subset(el2, pat_id %in% a00$pat_id & !is.na(psa))

a001a0_4<-merge(el2, n_psa_id2, all=TRUE)
#n_kons / length fup
a001a0_4$rate_p<-a001a0_4$n_kons/as.numeric(a001a0_4$fup_y_3)
# n_psa / n_kons
a001a0_4$rate_p1<-a001a0_4$n_psa/a001a0_4$n_kons

n_psa_arzt<- n_psa_id %>% group_by(arzt_id) %>% 
       summarize(min_dt = min(as.Date(fup_start_dt_3))
                , max_dt = max(as.Date(fup_end_dt_3))
                , diff0 = as.numeric(max_dt-min_dt)/365.25
                , n_psa_a = sum(!is.na(psa), na.rm=TRUE)
                , rate = n_psa_a/diff0
                )

# to obtain the same using tapply (for example for n_psa) use 
#n_psa_arzt<-data.frame(arzt_id=unique(n_psa_id$arzt_id), n_psa_a=tapply(n_psa_id$psa, n_psa_id$arzt_id, length))


#hist(n_psa_arzt$rate)
#hist(subset(n_psa_arzt, rate<=100)$rate)
n_psa_arzt %>% group_by(a20=floor(rate/10)*10) %>% summarize(n=n())
n_psa_arzt$c_arzt_f <- cut(n_psa_arzt$rate, c(0,10,65), right=FALSE, labels=FALSE)

#DB with Arzt and patient rate screening
a001a0_4<-merge(a001a0_4, n_psa_arzt, all=TRUE)
#a001a0_4$n_psa_a[is.na(a001a0_4$n_psa_a)] <- 0
a001a0_4$c_pat_k <- cut(a001a0_4$rate_p, c(0,10,800), right=FALSE, labels=FALSE)
a001a0_4$c_pat_k1 <- cut(a001a0_4$rate_p1, c(0,0.1,1.1), right=FALSE, labels=FALSE)

------------------------------
# Incidence rates
# ----------------------------------
# CEV - function for incidence rates calculation (general or subgroups)
# note that we need to recalculate vector v (one observation per patient) that is computed 
# outside the function. So, check the value of it before doing computation (normal is computed for db data_el)
# CEV_age - function for incidence rates for age classes of patient during study
# CEV_cron - function for incidences rates for chronic diseases
data_el<-orderBy(~pat_id, data_el)

# v is the indicator vector of first patient's id appearance (necessary for incidence calculation and function CEV)

v <- c(1, rep(0, length(data_el$pat_id)-1))
for (i in 2:length(data_el$pat_id))
  if ((data_el$pat_id[i] %in% 1:data_el$pat_id[i-1]) == FALSE)
    v[i] <- 1;

#definition of trimesters between 01/01/2009 and 31/12/2016
quarter_starts <- seq(from=as.Date("2009-01-01"), to=as.Date("2017-01-01"), by="quarter")

#CEV(end of follow up, start of follow up, subgroup category, dB)
#check the value of v since is not calculated within the function

cev <- function(end, start, subgroup, data) {
  
  py <- rep(0, (length(quarter_starts)-1));
  events <- rep(0, (length(quarter_starts)-1));
  
  for (i in 1:(length(quarter_starts)-1)) {
    
    
    Qstart <- rep(as.character(quarter_starts[i]), length(data$pat_id));
    Qend <- rep(as.character(quarter_starts[i+1]), length(data$pat_id));
    
    py[i] <- sum(v * subgroup * (pmax(0,
                           as.Date(pmin(end, Qend)) -
                             as.Date(pmax(start, Qstart)))) / 365, na.rm = T);
    
    events[i] <- sum(subgroup * (data$valid_from_dt >= Qstart) &
                       (data$valid_from_dt < Qend) &
                       (data$valid_from_dt >= start) &
                       (data$valid_from_dt < end), na.rm = T);
    
  }
  return(cbind(py, events))}



#Data Set for Graphs
inc<-data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$elig_2_1, data_el), cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$elig_2_3, data_el))  
#cut_off periods
inc$c_dt<-findInterval(as.Date(inc$dt), c(as.numeric(as.Date("2011-11-01")), as.numeric(as.Date("2012-07-01")), as.numeric(as.Date("2014-05-01"))))

inc_arzt_sex<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_sex=="f" & data_el$elig_2_1, data_el), sex="f", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_sex=="f" & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_sex=="f" & data_el$elig_2_3, data_el)), 
          data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_sex=="m" & data_el$elig_2_1, data_el),sex="m", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_sex=="m" & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_sex=="m" & data_el$elig_2_3, data_el)))

inc_arzt_sex$c_dt<-findInterval(as.Date(inc_arzt_sex$dt), as.numeric(as.Date("2014-05-01")))

inc_arzt_prac<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Doppelpraxis" & data_el$elig_2_1, data_el),type="Doppelpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Doppelpraxis"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Doppelpraxis" & data_el$elig_2_3, data_el)), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Einzelpraxis" & data_el$elig_2_1, data_el),type="Einzelpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Einzelpraxis"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Einzelpraxis" & data_el$elig_2_3, data_el)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="Gruppenpraxis" & data_el$elig_2_1, data_el),type="Gruppenpraxis", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="Gruppenpraxis"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Gruppenpraxis" & data_el$elig_2_3, data_el)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$practice_type=="" & data_el$elig_2_1, data_el),type="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$practice_type=="" & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="" & data_el$elig_2_3, data_el)))

inc_arzt_prac_g<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$practice_type=="Einzelpraxis" & data_el$elig_2_3, data_el),type="Einzelpraxis"),
                       data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, (data_el$practice_type=="Gruppenpraxis" | data_el$practice_type=="Doppelpraxis") & data_el$elig_2_3, data_el),type="Gruppen/Doppelpraxis"))
inc_arzt_prac_g$c_dt<-findInterval(as.Date(inc_arzt_prac_g$dt), as.numeric(as.Date("2014-05-01")))


inc_arzt_reg<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="CEN"& data_el$elig_2_1, data_el),region="CEN", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="CEN"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="CEN"& data_el$elig_2_3, data_el)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="IND"& data_el$elig_2_1, data_el),region="IND-Ter", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="IND"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="IND"& data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="PERI"& data_el$elig_2_1, data_el),region="PERI", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="PERI"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="PERI"& data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="SUB"& data_el$elig_2_1, data_el),region="SUB", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="SUB"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="SUB"& data_el$elig_2_3, data_el)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code=="OTHER"& data_el$elig_2_1, data_el),region="OTHER", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code=="OTHER"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="OTHER"& data_el$elig_2_3, data_el)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$arzt_region_code==""& data_el$elig_2_1, data_el),region="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$arzt_region_code==""& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code==""& data_el$elig_2_3, data_el)))

inc_arzt_reg_g<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="CEN"& data_el$elig_2_3, data_el),region="CEN"), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="SUB"& data_el$elig_2_3, data_el), region="SUB"),
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$arzt_region_code=="PERI"& data_el$elig_2_3, data_el), region="PERI"))
inc_arzt_reg_g$c_dt<-findInterval(as.Date(inc_arzt_reg_g$dt), as.numeric(as.Date("2014-05-01")))

                    
inc_arzt_emp<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="Angestellt"& data_el$elig_2_1, data_el), employ="Angestellt", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="Angestellt"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Angestellt"& data_el$elig_2_3, data_el)), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="Selbständig"& data_el$elig_2_1, data_el), employ="Selbständig", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="Selbständig"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Selbständig"& data_el$elig_2_3, data_el)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status=="AssistentIn"& data_el$elig_2_1, data_el), employ="AssistentIn", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status=="AssistentIn"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="AssistentIn"& data_el$elig_2_3, data_el)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$employ_status==""& data_el$elig_2_1, data_el), employ="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$employ_status==""& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status==""& data_el$elig_2_3, data_el)))

inc_arzt_emp_g<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Angestellt"& data_el$elig_2_3, data_el), employ="Angestellt"),
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$employ_status=="Selbständig"& data_el$elig_2_3, data_el), employ="Selbständig"))
inc_arzt_emp_g$c_dt<-findInterval(as.Date(inc_arzt_emp_g$dt), as.numeric(as.Date("2014-05-01")))

data_el$c_doby <- cut(data_el$arzt_doby, c(1940,1950,1960,1970,1980, 1990), right=FALSE, labels=FALSE)

inc_arzt_doby<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==1 & data_el$elig_2_1, data_el),doby="[1940,1950)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==1 & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==1 & data_el$elig_2_3, data_el)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==2 & data_el$elig_2_1, data_el),doby="[1950,1960)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==2 & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==2 & data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==3 & data_el$elig_2_1, data_el),doby="[1960,1970)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==3 & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==3 & data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==4 & data_el$elig_2_1, data_el),doby="[1970,1980)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==4 & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==4 & data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_doby==5 & data_el$elig_2_1, data_el),doby="[1980,1990)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_doby==5 & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==5 & data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, is.na(data_el$c_doby) & data_el$elig_2_1, data_el),doby="", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, is.na(data_el$c_doby) & data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, is.na(data_el$c_doby) & data_el$elig_2_3, data_el)))

inc_arzt_doby_g<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, (data_el$c_doby==1 | data_el$c_doby==2) & data_el$elig_2_3, data_el),doby="[1940,1960)"),
                       data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_doby==3 & data_el$elig_2_3, data_el),doby="[1960,1970)"),
                       data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, (data_el$c_doby==4 | data_el$c_doby==5) & data_el$elig_2_3, data_el),doby="[1970,1990)"))
inc_arzt_doby_g$c_dt<-findInterval(as.Date(inc_arzt_doby_g$dt), as.numeric(as.Date("2014-05-01")))

data_el$c_age <- cut(data_el$pat_age_start_3, c(55,60,65,70,75,76), right=FALSE)

inc_pat_age<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[55,60)"& data_el$elig_2_1, data_el),age="[55,60)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[55,60)"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[55,60)"& data_el$elig_2_3, data_el)), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[60,65)"& data_el$elig_2_1, data_el),age="[60,65)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[60,65)"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[60,65)"& data_el$elig_2_3, data_el)),
                   data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[65,70)"& data_el$elig_2_1, data_el),age="[65,70)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[65,70)"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[65,70)"& data_el$elig_2_3, data_el)),
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[70,75)"& data_el$elig_2_1, data_el),age="[70,75)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[70,75)"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[70,75)"& data_el$elig_2_3, data_el)),
                  data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(data_el$fup_end_dt_1, data_el$fup_start_dt_1, data_el$c_age=="[75,76)"& data_el$elig_2_1, data_el),age="[75+)", cev(data_el$fup_end_dt_2, data_el$fup_start_dt_2, data_el$c_age=="[75,76)"& data_el$elig_2_2, data_el), cev(data_el$fup_end_dt_3, data_el$fup_start_dt_3, data_el$c_age=="[75,76)"& data_el$elig_2_3, data_el)))
inc_pat_age$c_dt<-findInterval(as.Date(inc_pat_age$dt), as.numeric(as.Date("2014-05-01")))


#DB patient Eligible
el2<-subset(data_el, elig_2_3==1)
el2<-orderBy(~pat_id, el2)

v2 <- c(1, rep(0, length(el2$pat_id)-1))
for (i in 2:length(el2$pat_id))
  if ((el2$pat_id[i] %in% 1:el2$pat_id[i-1]) == FALSE)
    v2[i] <- 1;

#CEV_age used only with el2 DB (end of fu, start of fu, age classes)
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
      ca[i]<-sum(grouping*v2)
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



class_age <-  c(55, 60, 65, 70, 75, 76);
date()
prova_2 <- cev_age(el2$fup_end_dt_3, el2$fup_start_dt_3, class_age);
date()

inc_pat_age_c<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,1], events=prova_2[,2], cases=prova_2[,3], class_age=rep("[-55)", 32)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,4], events=prova_2[,5], cases=prova_2[,6], class_age=rep("[55,60)", 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,7], events=prova_2[,8], cases=prova_2[,9], class_age=rep("[60,65)", 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,10], events=prova_2[,11], cases=prova_2[,12], class_age=rep("[65,70)", 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,13], events=prova_2[,14], cases=prova_2[,15], class_age=rep("[70,75)", 32)), 
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2[,16], events=prova_2[,17], cases=prova_2[,18], class_age=rep("[75+)", 32))) 
inc_pat_age_c$c_dt<-findInterval(as.Date(inc_pat_age_c$dt), as.numeric(as.Date("2014-05-01")))

#DB eligible
el<-subset(data_el, elig_2_1==1)
el1<-subset(data_el, elig_2_2==1)
el2<-subset(data_el, elig_2_3==1)

el2$c_doby <- cut(el2$arzt_doby, c(1940,1950,1960,1970,1980, 1990), right=FALSE, labels=FALSE)
el2$c_age <- cut(el2$pat_age_start_3, c(55,60,65,70,75,76), right=FALSE)

#number of people in classe age (general population)
mort_2010_m$c_age <- cut(mort_2010_m$age_y, c(55,60,65,70,75,76), right=FALSE)

#weight for adjusted
wei_ad<-tapply(mort_2010_m$n.pop_middle, mort_2010_m$c_age, sum, na.rm=TRUE)/sum(tapply(mort_2010_m$n.pop_middle, mort_2010_m$c_age, sum, na.rm=TRUE))


class_age_1 <-  c(55, 56, 60, 61, 65, 66, 70, 71, 75, 76)
prova_3 <- cev_age(el2$fup_end_dt_3, el2$fup_start_dt_3, class_age_1)

inc_pat_age_c1<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,4], events=prova_3[,5], cases=prova_3[,6], class_age=rep("[55, 60)", 32)),
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,10], events=prova_3[,11], cases=prova_3[,12], class_age=rep("[60, 65)", 32)), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,16], events=prova_3[,17], cases=prova_3[,18], class_age=rep("[65, 70)", 32)), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,22], events=prova_3[,23], cases=prova_3[,24], class_age=rep("[70, 75)", 32)), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,28], events=prova_3[,29], cases=prova_3[,30], class_age=rep("[75+)", 32))) 

# total number of events:
sum(prova[, 2*(1:(length(class_age)+1))])

####### Incidence for high/low rates
a001a0_4<-orderBy(~pat_id, a001a0_4)
# we recompute v because we use db a001a04
v <- c(1, rep(0, length(a001a0_4$pat_id)-1))
for (i in 2:length(a001a0_4$pat_id))
  if ((a001a0_4$pat_id[i] %in% 1:a001a0_4$pat_id[i-1]) == FALSE)
    v[i] <- 1;


inc_arzt_high<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(a001a0_4$fup_end_dt_3, a001a0_4$fup_start_dt_3, a001a0_4$c_arzt_f==1, a001a0_4), freq="l"), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(a001a0_4$fup_end_dt_3, a001a0_4$fup_start_dt_3, a001a0_4$c_arzt_f==2, a001a0_4), freq="h"))
inc_arzt_high$c_dt<-findInterval(as.Date(inc_arzt_high$dt), as.numeric(as.Date("2014-05-01")))

inc_paz_high<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(a001a0_4$fup_end_dt_3, a001a0_4$fup_start_dt_3, a001a0_4$c_pat_k==1, a001a0_4), freq="l"), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(a001a0_4$fup_end_dt_3, a001a0_4$fup_start_dt_3, a001a0_4$c_pat_k==2, a001a0_4), freq="h"))
inc_paz_high$c_dt<-findInterval(as.Date(inc_paz_high$dt), as.numeric(as.Date("2014-05-01")))

inc_paz1_high<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(a001a0_4$fup_end_dt_3, a001a0_4$fup_start_dt_3, a001a0_4$c_pat_k1==1, a001a0_4), freq="l"), 
                     data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]),cev(a001a0_4$fup_end_dt_3, a001a0_4$fup_start_dt_3, a001a0_4$c_pat_k1==2, a001a0_4), freq="h"))
inc_paz1_high$c_dt<-findInterval(as.Date(inc_paz1_high$dt), as.numeric(as.Date("2014-05-01")))


#chronic disease

# counter of cronic diseases, one row per patient, one column per quarter
n_quarters <- length(quarter_starts);
cronic_disease_counter <- matrix(0, length(el2$pat_id), n_quarters);
quarter_limits <- c(quarter_starts[2:n_quarters], quarter_starts[n_quarters]);
for (i in 1:length(el2$pat_id)) {
  pat_line <- rep(0, n_quarters);
  onset_dates <- chr$effective_dt[chr$pat_id == el2$pat_id[i]];
  jump <- chr$cd_count[chr$pat_id == el2$pat_id[i]];
  if (length(onset_dates) > 0)
    for (j in 1:length(onset_dates)) {
      pat_line <- pat_line + (jump[j]-pat_line)*(onset_dates[j] < quarter_limits);
    }
  cronic_disease_counter[i, ] <- pat_line;
}

v2 <- c(1, rep(0, length(el2$pat_id)-1))
for (i in 2:length(el2$pat_id))
  if ((el2$pat_id[i] %in% 1:el2$pat_id[i-1]) == FALSE)
    v2[i] <- 1;

cev_cronic <- function(end, start, class_cronic) {
  # end & start are dates for patients' stay in the study period
  # for instance class_cronic <- c(4, 11) should partition the possible number of
  # cronic diseases in [0, 3], [4, 10], [11, 20]
  
  res <- NULL;
  
  n_classes <- length(class_cronic);
  
  for (j in 1:(n_classes+1)) {
    
    # computes start and end cutoffs for class j
    if (j == 1) {
      class_start <- 0; class_end <- class_cronic[1]-1;
    } else {
      if (j == (n_classes+1)) {
        class_start <- class_cronic[n_classes]; class_end <- 120;
      } else {
        class_start <- class_cronic[j-1]; class_end <- class_cronic[j]-1;
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
      
      #      # current ages for patients in the database
      #      current_age <- floor(2009 - el2$doby + # age at start of study period
      #                             i / 4);               # years age at current quarter
      
      # indicator function of patients belonging to class j
      #      grouping <- (current_age >= class_start) & (current_age <= class_end);
      grouping <- (cronic_disease_counter[, i] >= class_start) &
        (cronic_disease_counter[, i] <= class_end);
      
      ca[i]<-sum(grouping*v2)
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

class_cronic <- c(1, 4, 11);
date()
prova_2b <- cev_cronic(el2$fup_end_dt_3, el2$fup_start_dt_3, class_cronic);
date()

inc_pat_cron<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2b[,1], events=prova_2b[,2], cases=prova_2b[,3], cronic=rep("0", 32)),
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2b[,4], events=prova_2b[,5], cases=prova_2b[,6], cronic=rep("[1,3]", 32)), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2b[,7], events=prova_2b[,8], cases=prova_2b[,9], cronic=rep("[4,10]", 32)), 
                    data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_2b[,10], events=prova_2b[,11], cases=prova_2b[,12], cronic=rep("[11+)", 32))) 
inc_pat_cron$c_dt<-findInterval(as.Date(inc_pat_cron$dt), as.numeric(as.Date("2014-05-01")))

count_dis<-data.frame(pat_id=el2$pat_id, arzt_id=el2$arzt_id, psa=el2$psa, data_psa=el2$valid_from_dt, fup_start=el2$fup_start_dt_3, fup_end=el2$fup_end_dt_3, count_cron=cronic_disease_counter[,33])
count_dis$c_cron <- cut(count_dis$count_cron, c(0,1,4,11,21), right=FALSE)
count_dis$count_cron_start<-cronic_disease_counter[,5]
count_dis_psa<-subset(count_dis, !is.na(psa), select=c("pat_id", "count_cron"))
count_dis_psa<-distinct(count_dis_psa)
n_psa_id_c<-join(n_psa_id_c, count_dis_psa)
n_psa_id_c$cd_count[is.na(n_psa_id_c$cd_count)]<-n_psa_id_c$count_cron[is.na(n_psa_id_c$cd_count)]
n_psa_id_c$c_cron <- cut(n_psa_id_c$cd_count, c(0,1,4,11,21), right=FALSE)

#CHECKING
##########################
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
sum(as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3), na.rm=T)

#N PSA

el<-subset(data_el, elig_2_1==1)
el1<-subset(data_el, elig_2_2==1)
el2<-subset(data_el, elig_2_3==1)

el2$c_doby <- cut(el2$arzt_doby, c(1940,1950,1960,1970,1980, 1990), right=FALSE, labels=FALSE)
el2$c_age <- cut(el2$pat_age_start_3, c(55,60,65,70,75,76), right=FALSE)



class_age_1 <-  c(55, 56, 60, 61, 65, 66, 70, 71, 75, 76)
prova_3 <- cev_age(el2$fup_end_dt_3, el2$fup_start_dt_3, class_age_1)

inc_pat_age_c1<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,4], events=prova_3[,5], cases=prova_3[,6], class_age=rep("[55, 60)", 32)),
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,10], events=prova_3[,11], cases=prova_3[,12], class_age=rep("[60, 65)", 32)), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,16], events=prova_3[,17], cases=prova_3[,18], class_age=rep("[65, 70)", 32)), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,22], events=prova_3[,23], cases=prova_3[,24], class_age=rep("[70, 75)", 32)), 
                      data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), py=prova_3[,28], events=prova_3[,29], cases=prova_3[,30], class_age=rep("[75+)", 32))) 

#patient at studyage (check)
c55<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==55 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c60<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==60 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c65<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==65 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c70<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==70 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )
c75<-sum(inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2009-01-01"],inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2010-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2011-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2012-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2013-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2014-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2015-01-01"], inc_pat_age_c1$cases[inc_pat_age_c1$class_age==75 & as.character(inc_pat_age_c1$dt)=="2016-01-01"] )

el2$ind55<-as.numeric(floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))- el2$doby)>=55 & floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))-el2$doby) <=59 | (floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))- el2$doby)>=55 & floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))-el2$doby) <=59))
el2$ind60<-as.numeric(floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))- el2$doby)>=60 & floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))-el2$doby) <=64 | (floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))- el2$doby)>=60 & floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))-el2$doby) <=64))
el2$ind65<-as.numeric(floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))- el2$doby)>=65 & floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))-el2$doby) <=69 | (floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))- el2$doby)>=65 & floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))-el2$doby) <=69))
el2$ind70<-as.numeric(floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))- el2$doby)>=70 & floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))-el2$doby) <=74 | (floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))- el2$doby)>=70 & floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))-el2$doby) <=74))
el2$ind75<-as.numeric(floor(as.numeric(format(as.Date(el2$fup_end_dt_3), '%Y'))- el2$doby)>=75  | (floor(as.numeric(format(as.Date(el2$fup_start_dt_3), '%Y'))- el2$doby)>=75))



age_m<-(el2$pat_age_start_3*as.numeric(el2$fup_y_3)*v2)
#weighted age (age with weight the contribution to the study)
sum(age_m)/sum(as.numeric(el2$fup_y_3)*v2)


######################Not reported (old PSA quantiles)
#PSA

med<-tapply(n_psa_id$psa, cut(as.Date(n_psa_id$valid_from_dt), as.Date(quarter_starts)), median)
setq<-tapply(n_psa_id$psa, cut(as.Date(n_psa_id$valid_from_dt), as.Date(quarter_starts)), quantile, 0.75)
novq<-tapply(n_psa_id$psa, cut(as.Date(n_psa_id$valid_from_dt), as.Date(quarter_starts)), quantile, 0.90)


quant<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq)),
             data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq)))          


#by sex
med_f<-tapply(n_psa_id$psa[n_psa_id$arzt_sex=="f"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_sex=="f"]), as.Date(quarter_starts)), median) 
med_m<-tapply(n_psa_id$psa[n_psa_id$arzt_sex=="m"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_sex=="m"]), as.Date(quarter_starts)), median) 
setq_f<-tapply(n_psa_id$psa[n_psa_id$arzt_sex=="f"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_sex=="f"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_m<-tapply(n_psa_id$psa[n_psa_id$arzt_sex=="m"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_sex=="m"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_f<-tapply(n_psa_id$psa[n_psa_id$arzt_sex=="f"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_sex=="f"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_m<-tapply(n_psa_id$psa[n_psa_id$arzt_sex=="m"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_sex=="m"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_sex<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_f), sex="f"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_m), sex="m"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_f), sex="f"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_f), sex="m"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_f), sex="f"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_m), sex="m"))

#region
med_C<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="CEN"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="CEN"]), as.Date(quarter_starts)), median) 
med_I<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="IND"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="IND"]), as.Date(quarter_starts)), median) 
med_O<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="OTHER"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), median) 
med_P<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="PERI"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="PERI"]), as.Date(quarter_starts)), median) 
med_S<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="SUB"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="SUB"]), as.Date(quarter_starts)), median) 
setq_C<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="CEN"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="CEN"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_I<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="IND"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="IND"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_O<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="OTHER"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_P<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="PERI"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="PERI"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_S<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="SUB"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="SUB"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_C<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="CEN"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="CEN"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_I<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="IND"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="IND"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_O<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="OTHER"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="OTHER"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_P<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="PERI"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="PERI"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_S<-tapply(n_psa_id$psa[n_psa_id$arzt_region_code=="SUB"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$arzt_region_code=="SUB"]), as.Date(quarter_starts)), quantile, 0.90) 

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


#by employ status
med_a<-tapply(n_psa_id$psa[n_psa_id$employ_status=="Angestellt"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$employ_status=="Angestellt"]), as.Date(quarter_starts)), median) 
med_s<-tapply(n_psa_id$psa[n_psa_id$employ_status=="Selbständig"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$employ_status=="Selbständig"]), as.Date(quarter_starts)), median) 
setq_a<-tapply(n_psa_id$psa[n_psa_id$employ_status=="Angestellt"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$employ_status=="Angestellt"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_s<-tapply(n_psa_id$psa[n_psa_id$employ_status=="Selbständig"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$employ_status=="Selbständig"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_a<-tapply(n_psa_id$psa[n_psa_id$employ_status=="Angestellt"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$employ_status=="Angestellt"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_s<-tapply(n_psa_id$psa[n_psa_id$employ_status=="Selbständig"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$employ_status=="Selbständig"]), as.Date(quarter_starts)),quantile, 0.90) 

quant_emp<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_a), employ="Angestellt"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_s), employ="Selbständig"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_a), employ="Angestellt"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_s), employ="Selbständig"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_a), employ="Angestellt"),          
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_s), employ="Selbständig"))


#by praxis
med_d<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Doppelpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), median) 
med_e<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Einzelpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), median) 
med_g<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Gruppenpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), median) 
setq_d<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Doppelpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_e<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Einzelpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_g<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Gruppenpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_d<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Doppelpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Doppelpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_e<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Einzelpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Einzelpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_g<-tapply(n_psa_id$psa[n_psa_id$practice_type=="Gruppenpraxis"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$practice_type=="Gruppenpraxis"]), as.Date(quarter_starts)), quantile, 0.90) 

quant_pra<-rbind(data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.5, value=as.vector(med_g), practice="Gruppenpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.75, value=as.vector(setq_g), practice="Gruppenpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_d), practice="Doppelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_e), practice="Einzelpraxis"),
                 data.frame(dt=as.character(quarter_starts[1:(length(quarter_starts)-1)]), psa_quant=0.9, value=as.vector(novq_g), practice="Gruppenpraxis"))


#by arzt age
med_a1<-tapply(n_psa_id$psa[n_psa_id$c_doby==1], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==1]), as.Date(quarter_starts)), median) 
med_a2<-tapply(n_psa_id$psa[n_psa_id$c_doby==2], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==2]), as.Date(quarter_starts)), median) 
med_a3<-tapply(n_psa_id$psa[n_psa_id$c_doby==3], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==3]), as.Date(quarter_starts)), median) 
med_a4<-tapply(n_psa_id$psa[n_psa_id$c_doby==4], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==4]), as.Date(quarter_starts)), median) 
med_a5<-tapply(n_psa_id$psa[n_psa_id$c_doby==5], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==5]), as.Date(quarter_starts)), median) 

setq_a1<-tapply(n_psa_id$psa[n_psa_id$c_doby==1], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==1]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a2<-tapply(n_psa_id$psa[n_psa_id$c_doby==2], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==2]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a3<-tapply(n_psa_id$psa[n_psa_id$c_doby==3], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==3]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a4<-tapply(n_psa_id$psa[n_psa_id$c_doby==4], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==4]), as.Date(quarter_starts)), quantile, 0.75) 
setq_a5<-tapply(n_psa_id$psa[n_psa_id$c_doby==5], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==5]), as.Date(quarter_starts)), quantile, 0.75) 

novq_a1<-tapply(n_psa_id$psa[n_psa_id$c_doby==1], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==1]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a2<-tapply(n_psa_id$psa[n_psa_id$c_doby==2], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==2]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a3<-tapply(n_psa_id$psa[n_psa_id$c_doby==3], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==3]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a4<-tapply(n_psa_id$psa[n_psa_id$c_doby==4], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==4]), as.Date(quarter_starts)), quantile, 0.90) 
novq_a5<-tapply(n_psa_id$psa[n_psa_id$c_doby==5], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_doby==5]), as.Date(quarter_starts)), quantile, 0.90) 

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

#by pat age
med_p1<-tapply(n_psa_id$psa[n_psa_id$c_age=="[55,60)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[55,60)"]), as.Date(quarter_starts)), median) 
med_p2<-tapply(n_psa_id$psa[n_psa_id$c_age=="[60,65)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[60,65)"]), as.Date(quarter_starts)), median) 
med_p3<-tapply(n_psa_id$psa[n_psa_id$c_age=="[65,70)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[65,70)"]), as.Date(quarter_starts)), median) 
med_p4<-tapply(n_psa_id$psa[n_psa_id$c_age=="[70,75)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[70,75)"]), as.Date(quarter_starts)), median) 
setq_p1<-tapply(n_psa_id$psa[n_psa_id$c_age=="[55,60)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[55,60)"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_p2<-tapply(n_psa_id$psa[n_psa_id$c_age=="[60,65)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[60,65)"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_p3<-tapply(n_psa_id$psa[n_psa_id$c_age=="[65,70)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[65,70)"]), as.Date(quarter_starts)), quantile, 0.75) 
setq_p4<-tapply(n_psa_id$psa[n_psa_id$c_age=="[70,75)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[70,75)"]), as.Date(quarter_starts)), quantile, 0.75) 
novq_p1<-tapply(n_psa_id$psa[n_psa_id$c_age=="[55,60)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[55,60)"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_p2<-tapply(n_psa_id$psa[n_psa_id$c_age=="[60,65)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[60,65)"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_p3<-tapply(n_psa_id$psa[n_psa_id$c_age=="[65,70)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[65,70)"]), as.Date(quarter_starts)), quantile, 0.90) 
novq_p4<-tapply(n_psa_id$psa[n_psa_id$c_age=="[70,75)"], cut(as.Date(n_psa_id$valid_from_dt[n_psa_id$c_age=="[70,75)"]), as.Date(quarter_starts)), quantile, 0.90) 

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

