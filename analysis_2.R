# 1. Baseline characteristics of patients study start 01/01/2010 - Mixed follow-up
#db of start follow_up >= 01/2010
quarter_starts_e <- seq(from=as.Date("2010-01-01"), to=as.Date("2017-01-01"), by="quarter")

el2_2010<-subset(el2, fup_end_dt_3>="2010-01-01", select=c(pat_id, dead, arzt_id, doby, n_kons, fup_start_dt_3, fup_end_dt_3, fup_y_3, arzt_doby, arzt_sex, practice_type, employ_status, arzt_region, arzt_region_code, psa, valid_from_dt, valid_to_dt))

n_base<-distinct(el2_2010, pat_id, .keep_all=TRUE) 

count_dis_e<-distinct(subset(count_dis, select=c("pat_id", "count_cron", "count_cron_start")))

n_base<-join(n_base, count_dis_e)

n_base$cat_c_start<-cut(n_base$count_cron_start, c(0,1,4,11), right=FALSE)

n_base$age<-floor(2010-n_base$doby)
n_base$arzt_age<-floor(2010-n_base$arzt_doby)

library(nlme)
par.plot()

el2_2010<-orderBy(~pat_id, el2_2010)

v_e2 <- c(1, rep(0, length(el2_2010$pat_id)-1))
for (i in 2:length(el2_2010$pat_id))
  if ((el2_2010$pat_id[i] %in% 1:el2_2010$pat_id[i-1]) == FALSE)
    v_e2[i] <- 1;

el2_2010$py_2012<-v_e2*(pmax(0,
              as.Date(pmin(el2_2010$fup_end_dt_3, "2012-07-01")) -
                as.Date(pmax(el2_2010$fup_start_dt_3, "2010-01-01")))/365)
el2_2010$py_2014<-v_e2*(pmax(0,
                      as.Date(pmin(el2_2010$fup_end_dt_3, "2014-05-01")) -
                        as.Date(pmax(el2_2010$fup_start_dt_3, "2012-07-01")))/365)
el2_2010$py_2016<-v_e2*(pmax(0,
                      as.Date(pmin(el2_2010$fup_end_dt_3, "2017-01-01")) -
                        as.Date(pmax(el2_2010$fup_start_dt_3, "2014-05-01")))/365)

el2_2010$ev_2012<-(el2_2010$valid_from_dt >= "2010-01-01") &
  (el2_2010$valid_from_dt < "2012-07-01") &
  (el2_2010$valid_from_dt >= el2_2010$fup_start_dt_3) &
  (el2_2010$valid_from_dt < el2_2010$fup_end_dt_3)

el2_2010$ev_2014<-(el2_2010$valid_from_dt >= "2012-07-01") &
  (el2_2010$valid_from_dt < "2014-05-01") &
  (el2_2010$valid_from_dt >= el2_2010$fup_start_dt_3) &
  (el2_2010$valid_from_dt < el2_2010$fup_end_dt_3)

el2_2010$ev_2016<-(el2_2010$valid_from_dt >= "2014-05-01") &
  (el2_2010$valid_from_dt < "2017-01-01") &
  (el2_2010$valid_from_dt >= el2_2010$fup_start_dt_3) &
  (el2_2010$valid_from_dt < el2_2010$fup_end_dt_3)

n_psa_id_e<-subset(el2_2010, pmax(ev_2012, ev_2014, ev_2016)==TRUE)
n_psa_id_e$dt<-cut(as.Date(n_psa_id_e$valid_from_dt), as.Date(quarter_starts_e))

n_psa_id2_e<-data.frame(pat_id=unique(n_psa_id_e$pat_id), n_psa=tapply(n_psa_id_e$psa, n_psa_id_e$pat_id, length))

inc_arzt<-data.frame(cbind(arzt_id=as.vector(rownames(tapply(el2_2010$py_2012, el2_2010$arzt_id, sum))), py_1=as.numeric(tapply(el2_2010$py_2012, el2_2010$arzt_id, sum)), py_2=as.numeric(tapply(el2_2010$py_2014, el2_2010$arzt_id, sum)), py_3=as.numeric(tapply(el2_2010$py_2016, el2_2010$arzt_id, sum)), ev_1=as.numeric(tapply(el2_2010$ev_2012, el2_2010$arzt_id, sum, na.rm=TRUE)), ev_2=as.numeric(tapply(el2_2010$ev_2014, el2_2010$arzt_id, sum, na.rm=TRUE)), ev_3=as.numeric(tapply(el2_2010$ev_2016, el2_2010$arzt_id, sum, na.rm=TRUE))), stringsAsFactors=FALSE)
inc_arzt$py<-as.numeric(inc_arzt$py_1)+as.numeric(inc_arzt$py_2)+as.numeric(inc_arzt$py_3)
inc_arzt$ev<-as.numeric(inc_arzt$ev_1)+as.numeric(inc_arzt$ev_2)+as.numeric(inc_arzt$ev_3)
inc_arzt$inc<-inc_arzt$ev/inc_arzt$py
inc_arzt$inc1<-as.numeric(inc_arzt$ev_1)/as.numeric(inc_arzt$py_1)
inc_arzt$inc2<-as.numeric(inc_arzt$ev_2)/as.numeric(inc_arzt$py_2)
inc_arzt$inc3<-as.numeric(inc_arzt$ev_3)/as.numeric(inc_arzt$py_3)

dt<-c(as.Date("2010-01-01"), as.Date("2012-07-01"), as.Date("2014-05-01"), as.Date("2017-01-01"))

tab1_c<-data.frame(Characteristic=c("Age(yr)", "Mean (sd)", "Median", "Inter-Range", "PSA", "Yes", "N_PSA", "N_PSA/patient", "N_PSA/pat Range", "Chronic diseases", "[0,1)", "[1,4)", "[4,11)", "Length FUP (yr)", "Mean (sd)", "Median", "Inter_range", "Konsultations", "Mean (sd)", "Median", "Inter_range"), 
                   Value=c("", paste(summary(n_base$age)[4]," (", format(sd(n_base$age), digits=2), ")", sep=""), summary(n_base$age)[3], paste("(", summary(n_base$age)[2], ",", summary(n_base$age)[5], ")", sep=""), "n (%)", paste(length(n_psa_id2_e$pat_id), " (", format(length(n_psa_id2_e$pat_id)/length(n_base$pat_id)*100, digits=2), "%)", sep=""), 
                           sum(as.numeric(inc_arzt$ev_1), as.numeric(inc_arzt$ev_2), as.numeric(inc_arzt$ev_3)), format(sum(as.numeric(inc_arzt$ev_1), as.numeric(inc_arzt$ev_2), as.numeric(inc_arzt$ev_3))/length(n_psa_id2_e$pat_id), digits=2), paste("(", summary(n_psa_id2_e$n_psa)[2], ",", summary(n_psa_id2_e$n_psa)[6], ")", sep=""),
                           "n (%)", paste(as.vector(table(n_base$cat_c_start)), " (", format(as.vector(table(n_base$cat_c_start))/sum((table(n_base$cat_c_start)))*100, digits=2), ")", sep=""),
                           "", paste(summary(as.numeric(n_base$fup_y_3))[4], " (", format(sd(n_base$fup_y_3), digits=2), ")", sep=""), summary(as.numeric(n_base$fup_y_3))[3], paste("(", summary(as.numeric(n_base$fup_y_3))[2], ",", summary(as.numeric(n_base$fup_y_3))[5], ")", sep=""),
                           "", paste(summary(n_base$n_kons)[4], " (", format(sd(n_base$n_kons), digits=2), ")", sep=""), summary(n_base$n_kons)[3], paste("(", summary(n_base$n_kons)[2], ",", summary(n_base$n_kons)[5], ")", sep="")))
#Arzt characteristics
# age at 2010: we have 75!
tab2_c<-data.frame(Characteristic=c("Age(yr)", "Mean (sd)", "Median", "Inter range", "N_pat", "N_pat/arzt", "Median", "Inter-Range","Sex","f", "m", "PSA", "Yes", "N_PSA", "N_PSA/arzt", "Median", "Inter Range", "Practice","Doppelpraxis", "Einzelpraxis", "Gruppenpraxis", "missing", "Employ status", "Angestellt", "AssistentIn", "Selbständig", "missing", "Region", "CEN", "SUB", "PERI", "IND", "OTHER", "missing"), 
                   Value=c("", paste(summary(n_base$arzt_age)[4], " (", format(sd(n_base$arzt_age, na.rm=TRUE), digits = 2), ")", sep="") , summary(n_base$arzt_age)[3], paste("(", summary(n_base$arzt_age)[2], ",", summary(n_base$arzt_age)[5], ")", sep=""), 
                           length(n_base$pat_id), summary(tapply(n_base$pat_id, n_base$arzt_id, length))[4], summary(tapply(n_base$pat_id, n_base$arzt_id, length))[3], paste("(", summary(tapply(n_base$pat_id, n_base$arzt_id, length))[2], ",", summary(tapply(n_base$pat_id, n_base$arzt_id, length))[5], ")", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$arzt_sex=="f"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_sex=="f"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_sex=="m"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_sex=="m"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
                           "n (%)", paste(length(inc_arzt$arzt_id[inc_arzt$ev>0]), " (", format(length(inc_arzt$arzt_id[inc_arzt$ev>0])/length(inc_arzt$arzt_id)*100, digits=2), "%)", sep=""), sum(inc_arzt$ev), format(sum(inc_arzt$ev)/length(inc_arzt$arzt_id[inc_arzt$ev>0]), digits=2), summary(inc_arzt$ev)[3], paste("(", summary(inc_arzt$ev)[2], ", ", summary(inc_arzt$ev)[5], ")", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$practice_type=="Doppelpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Doppelpraxis"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type=="Einzelpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Einzelpraxis"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type=="Gruppenpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Gruppenpraxis"]))/length(unique(n_base$arzt_id))*100, digits=3), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type==""])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$employ_status=="Angestellt"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="Angestellt"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status=="AssistentIn"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="AssistentIn"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status=="Selbständig"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="Selbständig"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status==""])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="CEN"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="CEN"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="SUB"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="SUB"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="PERI"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="PERI"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="IND"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="IND"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="OTHER"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region=="OTHER"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code==""])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep="")))

#PSA characteristics
# to add: distribution of chronic diseases
tab3_c<-data.frame(Period=c("[01/01/2010-01/07/2012)", "[01/07/2012-01/05/2014)", "[01/05/2014-01/01/2017)", "Total"),
           N_PSA=c(sum(as.numeric(inc_arzt$ev_1)), sum(as.numeric(inc_arzt$ev_2)), sum(as.numeric(inc_arzt$ev_3)), sum(inc_arzt$ev)),
           PY=c(sum(as.numeric(inc_arzt$py_1)), sum(as.numeric(inc_arzt$py_2)), sum(as.numeric(inc_arzt$py_3)), sum(inc_arzt$py)),
           ir=c(sum(as.numeric(inc_arzt$ev_1))/sum(as.numeric(inc_arzt$py_1)), sum(as.numeric(inc_arzt$ev_2))/sum(as.numeric(inc_arzt$py_2)), sum(as.numeric(inc_arzt$ev_3))/sum(as.numeric(inc_arzt$py_3)), sum(as.numeric(inc_arzt$ev))/sum(as.numeric(inc_arzt$py))),
           PSA_mean_value=c(summary(n_psa_id_e$psa[n_psa_id_e$ev_2012==TRUE])[4], summary(n_psa_id_e$psa[n_psa_id_e$ev_2014==TRUE])[4], summary(n_psa_id_e$psa[n_psa_id_e$ev_2016==TRUE])[4], summary(n_psa_id_e$psa)[4]))

cc<-reshape(inc_arzt, direction="long", varying=list(names(inc_arzt)[11:13]), v.names="Inc", 
            idvar=c("arzt_id"), timevar="Period", times=1:3, drop=c("py_1", "py", "py_2", "py_3", "ev_1", "ev_2", "ev_3", "ev", "inc"))
cc<-orderBy(~arzt_id, cc)

par.plot(Inc~Period,data=cc,subject=arzt_id)

#Mixed model
ff<-join(n_psa_id_e, inc_arzt)
ff<-join(ff, n_psa_id2_e)

ff$c_dt<-findInterval(as.Date(ff$dt), c(as.numeric(as.Date("2010-01-01")), as.numeric(as.Date("2012-07-01")), as.numeric(as.Date("2014-05-01"))))
ff$ev_per<-as.numeric(ff$ev_1)*(ff$c_dt==1)+as.numeric(ff$ev_2)*(ff$c_dt==2)+as.numeric(ff$ev_3)*(ff$c_dt==3)

ff$age_psa<-floor(year(ff$valid_from_dt)-ff$doby)
tapply(ff$psa, ff$c_dt, mean)
tapply(ff$psa, ff$age_psa, mean)
tapply(ff$psa, ff$ev_per, mean)
tapply(ff$psa, ff$n_psa, mean)
tapply(ff$psa, ff$ch_per, mean)
plot(tapply(ff$psa, ff$n_psa, mean))
summary(lm(psa~n_psa, ff))

ggplot(ff, aes(arzt_id, log(psa))) + 
  geom_boxplot(aes(fill = factor(c_dt)))+
  theme_bw()

kruskal.test(psa~factor(arzt_id), data=ff)
kruskal.test(log(psa)~factor(c_dt), data=ff)
kruskal.test(log(psa)~factor(pat_id), data=ff)


fm_mix<-lmer(log(psa)~ev+(1|arzt_id), ff)
plot(fm_mix, sqrt(abs(resid(.)))~fitted(.), type=c("p", "smooth"))
fm_mix<-lmer(log(psa)~ev_per+age_psa+factor(c_dt)+(1|arzt_id/pat_id), ff)
coef(fm_mix)

n_q <- length(dt);
cronic_disease_counter_e <- matrix(0, length(ff$pat_id), n_q);
quarter_li <- c(dt[2:n_q], dt[n_q]);
for (i in 1:length(ff$pat_id)) {
  pat_line_e <- rep(0, n_q);
  onset_dates_e <- chr$effective_dt[chr$pat_id == ff$pat_id[i]];
  jump_e <- chr$cd_count[chr$pat_id == ff$pat_id[i]];
  if (length(onset_dates_e) > 0)
    for (j in 1:length(onset_dates_e)) {
      pat_line_e <- pat_line_e + (jump_e[j]-pat_line_e)*(onset_dates_e[j] < quarter_li);
    }
  cronic_disease_counter_e[i, ] <- pat_line_e;
}
ff$ch_per<-as.numeric(cronic_disease_counter_e[,2])*(ff$c_dt==1)+as.numeric(cronic_disease_counter_e[,3])*(ff$c_dt==2)+as.numeric(cronic_disease_counter_e[,4])*(ff$c_dt==3)
