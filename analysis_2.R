# 1. Baseline characteristics of patients study start 01/01/2010 - Mixed follow-up
#db of start follow_up >= 01/2010

el2_2010<-subset(el2, fup_end_dt_3>="2010-01-01", select=c(pat_id, dead, arzt_id, doby, n_kons, fup_start_dt_3, fup_end_dt_3, fup_y_3, arzt_doby, arzt_sex, practice_type, employ_status, arzt_region, arzt_region_code, psa, valid_from_dt, valid_to_dt))

n_base<-distinct(el2_2010, pat_id, .keep_all=TRUE) 

el2_2010<-merge(el2_2010, inc_pat, all=TRUE)
n_base<-join(n_base, inc_pat)
n_base$n_psa<-n_base$ev_1+n_base$ev_2+n_base$ev_3

count_dis_e<-distinct(subset(count_dis, select=c("pat_id", "count_cron", "count_cron_start", "count_cron_2", "count_cron_3")))
n_base<-join(n_base, count_dis_e)

n_base$cat_c_start<-cut(n_base$count_cron_start, c(0,1,4,11,20), right=FALSE)
count_dis_e$cat_c_2<-cut(count_dis_e$count_cron_2, c(0,1,4,11,20), right=FALSE)
count_dis_e$cat_c_3<-cut(count_dis_e$count_cron_3, c(0,1,4,11, 20), right=FALSE)
count_dis_e$cat_c_4<-cut(count_dis_e$count_cron, c(0,1,4,11, 20), right=FALSE)

n_base$age_start<-floor(year(as.Date(n_base$fup_start_dt_3))-n_base$doby)
n_base$arzt_age<-floor(2010-n_base$arzt_doby)

#n_psa_id_e<-subset(el2_2010, (el2_2010$age_psa>=55 & el2_2010$age_psa<=75) & as.Date(el2_2010$valid_from_dt) < as.Date(el2_2010$fup_end_dt_3))
#n_psa_id_e$dt<-cut(as.Date(n_psa_id_e$valid_from_dt), as.Date(quarter_starts_e))

#n_psa_id2_e<-data.frame(pat_id=unique(n_psa_id_e$pat_id), n_psa=tapply(n_psa_id_e$psa, n_psa_id_e$pat_id, length))

#dt<-c(as.Date("2010-01-01"), as.Date("2012-07-01"), as.Date("2014-05-01"), as.Date("2017-01-01"))
#Patient characteristics
tab1_c<-data.frame(Characteristic=c("N_Patients", "Baseline", "Age(yr)", "Mean (sd)", "Median", "Inter-Range", "PSA", "Yes", "N_PSA", "N_PSA/patient", "N_PSA/pat Range", "PSA value",  "Mean (sd)", "Median", "Inter-range", "Chronic diseases", "[0,1)", "[1,4)", "[4,11)", ">= 11", "Length FUP (yr)", "Mean (sd)", "Median", "Inter_range", "Konsultations", "Mean (sd)", "Median", "Inter_range", "Konsult / length FUP", "Mean (sd)", "Median", "Inter_range", "N_PSA / Konsult", "Mean (sd)", "Median", "Inter_range" ), 
                   Total=c(length(n_base$pat_id), "", "", paste(summary(n_base$age_start)[4]," (", format(sd(n_base$age_start), digits=2), ")", sep=""), summary(n_base$age_start)[3], paste("(", summary(n_base$age_start)[2], ",", summary(n_base$age_start)[5], ")", sep=""), "n (%)", paste(length(n_base$n_psa[n_base$n_psa>0]), " (", format(length(n_base$n_psa[n_base$n_psa>0])/length(n_base$pat_id)*100, digits=2), "%)", sep=""), 
                           sum(inc$events), format(sum(inc$events)/length(n_base$n_psa[n_base$n_psa>0]), digits=2), paste("(", summary(n_base$n_psa)[2], ",", summary(n_base$n_psa)[6], ")", sep=""),
                           "", paste(summary(el2_2010$psa[el2_2010$c_dt>=1])[4], " (", format(sd(el2_2010$psa[el2_2010$c_dt>=1], na.rm=TRUE), digits=2), ")", sep=""), summary(el2_2010$psa[el2_2010$c_dt>=1])[3], paste("(", summary(el2_2010$psa[el2_2010$c_dt>=1])[2], ",", summary(el2_2010$psa[el2_2010$c_dt>=1])[5], ")", sep=""),
                          "n (%)", paste(as.vector(table(n_base$cat_c_start)), " (", format(as.vector(table(n_base$cat_c_start))/sum((table(n_base$cat_c_start)))*100, digits=2), ")", sep=""),
                           "", paste(summary(as.numeric(n_base$fup_y_3))[4], " (", format(sd(n_base$fup_y_3), digits=2), ")", sep=""), summary(as.numeric(n_base$fup_y_3))[3], paste("(", summary(as.numeric(n_base$fup_y_3))[2], ",", summary(as.numeric(n_base$fup_y_3))[5], ")", sep=""),
                           "", paste(summary(n_base$n_kons)[4], " (", format(sd(n_base$n_kons), digits=2), ")", sep=""), summary(n_base$n_kons)[3], paste("(", summary(n_base$n_kons)[2], ",", summary(n_base$n_kons)[5], ")", sep=""),
                           "", paste(format(summary(n_base$n_kons/as.numeric(n_base$fup_y_3))[4], digits=2), " (", format(sd(n_base$n_kons/as.numeric(n_base$fup_y_3)), digits=3), ")", sep=""), format(summary(n_base$n_kons/as.numeric(n_base$fup_y_3))[3], digits=2), paste("(", format(summary(n_base$n_kons/as.numeric(n_base$fup_y_3))[2], digits=2), ",", format(summary(n_base$n_kons/as.numeric(n_base$fup_y_3))[5], digits=2), ")", sep=""),
                           "", paste(format(summary(n_base$n_psa/n_base$n_kons)[4], digits=2), " (", format(sd(n_base$n_psa/n_base$n_kons), digits=2), ")", sep=""), format(summary(n_base$n_psa/n_base$n_kons)[3], digits=2), paste("(", format(summary(n_base$n_psa/n_base$n_kons)[2], digits=2), ",", format(summary(n_base$n_psa/n_base$n_kons)[5], digits=2), ")", sep="")),
                   Time1=c(length(n_base$pat_id[year(n_base$fup_start_dt_3)<=2012 & year(n_base$fup_end_dt_3)>2010]), "", "", paste(summary(n_base$age[year(n_base$fup_start_dt_3)<=2012 & year(n_base$fup_end_dt_3)>2010])[4]," (", format(sd(n_base$age[year(n_base$fup_start_dt_3)<=2012 & year(n_base$fup_end_dt_3)>2010]), digits=2), ")", sep=""), summary(n_base$age[year(n_base$fup_start_dt_3)<=2012 & year(n_base$fup_end_dt_3)>2010])[3], paste("(", summary(n_base$age[year(n_base$fup_start_dt_3)<=2012 & year(n_base$fup_end_dt_3)>2010])[2], ",", summary(n_base$age[year(n_base$fup_start_dt_3)<=2012 & year(n_base$fup_end_dt_3)>2010])[5], ")", sep=""), "n (%)", paste(length(n_base$ev_1[n_base$ev_1>0]), " (", format(length(n_base$ev_1[n_base$ev_1>0])/length(n_base$pat_id)*100, digits=2), "%)", sep=""), sum(n_base$ev_1), format(sum(n_base$ev_1)/length(n_base$n_psa[n_base$ev_1>0]), digits=2), paste("(", summary(n_base$ev_1)[2], ",", summary(n_base$ev_1)[6], ")", sep=""),
                           "", paste(tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`1`[4], " (", format(sd(el2_2010$psa[el2_2010$c_dt==1], na.rm=TRUE), digits=2), ")", sep=""), tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`1`[3], paste("(", tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`1`[2], ",", tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`1`[5], ")", sep=""),
                           "n (%)", paste(as.vector(table(count_dis_e$cat_c_2)), " (", format(as.vector(table(count_dis_e$cat_c_2))/sum((table(count_dis_e$cat_c_2)))*100, digits=2), ")", sep=""),
                           "", "", "", "", "", "", "", "", "", "", "", "", "",
                            paste(format(summary(n_base$ev_1/n_base$n_kons)[4], digits=2), " (", format(sd(n_base$ev_1/n_base$n_kons), digits=2), ")", sep=""), format(summary(n_base$ev_1/n_base$n_kons)[3], digits=2), paste("(", format(summary(n_base$ev_1/n_base$n_kons)[2], digits=2), ",", format(summary(n_base$ev_1/n_base$n_kons)[5], digits=2), ")", sep="")
                           ),
                   Time2=c(length(n_base$pat_id[year(n_base$fup_start_dt_3)<=2014 & year(n_base$fup_end_dt_3)>2012]), "", "", "", "", "", paste(length(n_base$ev_2[n_base$ev_2>0]), " (", format(length(n_base$ev_2[n_base$ev_2>0])/length(n_base$pat_id)*100, digits=2), "%)", sep=""), sum(n_base$ev_2), format(sum(n_base$ev_2)/length(n_base$n_psa[n_base$ev_2>0]), digits=2), paste("(", summary(n_base$ev_2)[2], ",", summary(n_base$ev_2)[6], ")", sep=""),
                           "", paste(tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`2`[4], " (", format(sd(el2_2010$psa[el2_2010$c_dt==2], na.rm=TRUE), digits=2), ")", sep=""), tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`2`[3], paste("(", tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`2`[2], ",", tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`2`[5], ")", sep=""),
                           "n (%)", paste(as.vector(table(count_dis_e$cat_c_3)), " (", format(as.vector(table(count_dis_e$cat_c_3))/sum((table(count_dis_e$cat_c_3)))*100, digits=2), ")", sep=""),
                           "", "", "", "", "", "", "", "", "", "", "", "", "",
                           paste(format(summary(n_base$ev_2/n_base$n_kons)[4], digits=2), " (", format(sd(n_base$ev_2/n_base$n_kons), digits=2), ")", sep=""), format(summary(n_base$ev_2/n_base$n_kons)[3], digits=2), paste("(", format(summary(n_base$ev_2/n_base$n_kons)[2], digits=2), ",", format(summary(n_base$ev_2/n_base$n_kons)[5], digits=2), ")", sep="")
                   ),
                   Time3=c(length(n_base$pat_id[year(n_base$fup_start_dt_3)<=2017 & year(n_base$fup_end_dt_3)>2014]), "", "", "", "", "", paste(length(n_base$ev_3[n_base$ev_3>0]), " (", format(length(n_base$ev_3[n_base$ev_3>0])/length(n_base$pat_id)*100, digits=2), "%)", sep=""), sum(n_base$ev_3), format(sum(n_base$ev_3)/length(n_base$n_psa[n_base$ev_3>0]), digits=2), paste("(", summary(n_base$ev_3)[2], ",", summary(n_base$ev_3)[6], ")", sep=""),
                           "", paste(tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`3`[4], " (", format(sd(el2_2010$psa[el2_2010$c_dt==3], na.rm=TRUE), digits=2), ")", sep=""), tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`3`[3], paste("(", tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`3`[2], ",", tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], summary)$`3`[5], ")", sep=""),
                           "n (%)", paste(as.vector(table(count_dis_e$cat_c_4)), " (", format(as.vector(table(count_dis_e$cat_c_4))/sum((table(count_dis_e$cat_c_4)))*100, digits=2), ")", sep=""),
                           "", "", "", "", "", "", "", "", "", "", "", "", "",
                           paste(format(summary(n_base$ev_3/n_base$n_kons)[4], digits=2), " (", format(sd(n_base$ev_3/n_base$n_kons), digits=2), ")", sep=""), format(summary(n_base$ev_3/n_base$n_kons)[3], digits=2), paste("(", format(summary(n_base$ev_3/n_base$n_kons)[2], digits=2), ",", format(summary(n_base$ev_3/n_base$n_kons)[5], digits=2), ")", sep="")
                   ))
#Mean SD
#Arzt characteristics
# age at 2010: we have 75!
tab2_c<-data.frame(Characteristic=c("Age(yr)", "Mean (sd)", "Median", "Inter range", "N_pat", "N_pat/arzt", "Mean (sd)", "Median", "Inter-Range","Sex","f", "m", "PSA", "Yes", "N_PSA", "N_PSA/arzt", "Mean (sd)", "Median", "Inter Range", "Konsultations", "Mean", "Median", "Inter-range","N_PSA / length FUP", "Mean (sd)", "Median", "Inter_range","Practice","Doppelpraxis", "Einzelpraxis", "Gruppenpraxis", "missing", "Employ status", "Angestellt", "AssistentIn", "Selbständig", "missing", "Region", "CEN", "SUB", "PERI", "IND", "OTHER", "missing"), 
                   Value=c("Baseline", paste(summary(n_base$arzt_age)[4], " (", format(sd(n_base$arzt_age, na.rm=TRUE), digits = 2), ")", sep="") , summary(n_base$arzt_age)[3], paste("(", summary(n_base$arzt_age)[2], ",", summary(n_base$arzt_age)[5], ")", sep=""), 
                           length(n_base$pat_id), "", paste(summary(tapply(n_base$pat_id, n_base$arzt_id, length))[4], " (", format(sd(tapply(n_base$pat_id, n_base$arzt_id, length)), digits=2), ")", sep=""), summary(tapply(n_base$pat_id, n_base$arzt_id, length))[3], paste("(", summary(tapply(n_base$pat_id, n_base$arzt_id, length))[2], ",", summary(tapply(n_base$pat_id, n_base$arzt_id, length))[5], ")", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$arzt_sex=="f"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_sex=="f"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_sex=="m"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_sex=="m"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
                           "n (%)", paste(n_distinct(n_base$arzt_id[n_base$n_psa>0]), " (", format(n_distinct(n_base$arzt_id[n_base$n_psa>0])/n_distinct(n_base$arzt_id)*100, digits=2), "%)", sep=""), sum(n_base$n_psa), "", paste(format(summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[4], digits=2), " (", format(sd(tapply(n_base$n_psa, n_base$arzt_id, sum)), digits=2), ")", sep=""), summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[3], paste("(", summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[2], ", ", summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[5], ")", sep=""),
                           "", paste(summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[4], " (", format(sd(tapply(n_base$n_kons, n_base$arzt_id, sum)), digits=2), ")", sep=""), summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[3], paste("(", summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[2], ",", summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[5], ")", sep=""),
                           "", paste(summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[4], " (", format(sd(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum)), digits=2), ")", sep=""), summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[3], paste("(", summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[2], ",", summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[5], ")", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$practice_type=="Doppelpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Doppelpraxis"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type=="Einzelpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Einzelpraxis"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type=="Gruppenpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Gruppenpraxis"]))/length(unique(n_base$arzt_id))*100, digits=3), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type==""])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$employ_status=="Angestellt"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="Angestellt"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status=="AssistentIn"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="AssistentIn"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status=="Selbständig"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="Selbständig"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status==""])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
                           "n (%)", paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="CEN"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="CEN"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="SUB"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="SUB"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="PERI"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="PERI"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="IND"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="IND"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="OTHER"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region=="OTHER"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code==""])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep="")), 
Time1=c("", "", "", "",  
        length(n_base$pat_id), "", paste(summary(tapply(n_base$pat_id, n_base$arzt_id, length))[4], " (", format(sd(tapply(n_base$pat_id, n_base$arzt_id, length)), digits=2), ")", sep=""), summary(tapply(n_base$pat_id, n_base$arzt_id, length))[3], paste("(", summary(tapply(n_base$pat_id, n_base$arzt_id, length))[2], ",", summary(tapply(n_base$pat_id, n_base$arzt_id, length))[5], ")", sep=""),
        "n (%)", paste(length(unique(n_base$arzt_id[n_base$arzt_sex=="f"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_sex=="f"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_sex=="m"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_sex=="m"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
        "n (%)", paste(n_distinct(n_base$arzt_id[n_base$n_psa>0]), " (", format(n_distinct(n_base$arzt_id[n_base$n_psa>0])/n_distinct(n_base$arzt_id)*100, digits=2), "%)", sep=""), sum(n_base$n_psa), "", paste(format(summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[4], digits=2), " (", format(sd(tapply(n_base$n_psa, n_base$arzt_id, sum)), digits=2), ")", sep=""), summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[3], paste("(", summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[2], ", ", summary(tapply(n_base$n_psa, n_base$arzt_id, sum))[5], ")", sep=""),
        "", paste(summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[4], " (", format(sd(tapply(n_base$n_kons, n_base$arzt_id, sum)), digits=2), ")", sep=""), summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[3], paste("(", summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[2], ",", summary(tapply(n_base$n_kons, n_base$arzt_id, sum))[5], ")", sep=""),
        "", paste(summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[4], " (", format(sd(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum)), digits=2), ")", sep=""), summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[3], paste("(", summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[2], ",", summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(as.numeric(n_base$fup_y_3), n_base$arzt_id, sum))[5], ")", sep=""),
        "n (%)", paste(length(unique(n_base$arzt_id[n_base$practice_type=="Doppelpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Doppelpraxis"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type=="Einzelpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Einzelpraxis"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type=="Gruppenpraxis"])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type=="Gruppenpraxis"]))/length(unique(n_base$arzt_id))*100, digits=3), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$practice_type==""])), "(", format(length(unique(n_base$arzt_id[n_base$practice_type==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
        "n (%)", paste(length(unique(n_base$arzt_id[n_base$employ_status=="Angestellt"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="Angestellt"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status=="AssistentIn"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="AssistentIn"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status=="Selbständig"])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status=="Selbständig"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$employ_status==""])), "(", format(length(unique(n_base$arzt_id[n_base$employ_status==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""),
        "n (%)", paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="CEN"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="CEN"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="SUB"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="SUB"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="PERI"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="PERI"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="IND"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code=="IND"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code=="OTHER"])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region=="OTHER"]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep=""), paste(length(unique(n_base$arzt_id[n_base$arzt_region_code==""])), "(", format(length(unique(n_base$arzt_id[n_base$arzt_region_code==""]))/length(unique(n_base$arzt_id))*100, digits=2), "%)", sep="")))

#PSA characteristics
# to add: distribution of chronic diseases
tab3_c<-data.frame(Period=c("[01/01/2010-01/07/2012)", "[01/07/2012-01/05/2014)", "[01/05/2014-01/01/2017)", "Total"),
                   N_PSA=c(sum(inc$events[inc$c_dt==1]), sum(inc$events[inc$c_dt==2]), sum(inc$events[inc$c_dt==3]), sum(inc$events)),
                   PY=c(sum(inc$py[inc$c_dt==1]), sum(inc$py[inc$c_dt==2]), sum(inc$py[inc$c_dt==3]), sum(inc$py)),
                   ir=c(sum(inc$events[inc$c_dt==1])/sum(inc$py[inc$c_dt==1]), sum(inc$events[inc$c_dt==2])/sum(inc$py[inc$c_dt==2]), sum(inc$events[inc$c_dt==3])/sum(inc$py[inc$c_dt==3]), sum(inc$events)/sum(inc$py)),
                   ir_mean_pat=c(summary(n_base$ev_1/n_base$py_1)[4], summary(n_base$ev_2/n_base$py_2)[4], summary(n_base$ev_3/n_base$py_3)[4], summary(n_base$n_psa/(n_base$py_1+n_base$py_2+n_base$py_3))[4]),
                   ir_mean_arzt=c(summary(tapply(n_base$ev_1, n_base$arzt_id, sum)/tapply(n_base$py_1, n_base$arzt_id, sum))[4], summary(tapply(n_base$ev_2, n_base$arzt_id, sum)/tapply(n_base$py_2, n_base$arzt_id, sum))[4], summary(tapply(n_base$ev_3, n_base$arzt_id, sum)/tapply(n_base$py_3, n_base$arzt_id, sum))[4], summary(tapply(n_base$n_psa, n_base$arzt_id, sum)/tapply(n_base$py_1+n_base$py_2+n_base$py_3, n_base$arzt_id, sum))[4]))

#cc<-reshape(inc_arzt, direction="long", varying=list(names(inc_arzt)[11:13]), v.names="Inc", 
 #           idvar=c("arzt_id"), timevar="Period", times=1:3, drop=c("py_1", "py", "py_2", "py_3", "ev_1", "ev_2", "ev_3", "ev", "inc"))
#cc<-orderBy(~arzt_id, cc)

#par.plot(Inc~Period,data=cc,subject=arzt_id)

################
#Plot Incidence
##############

#Overall 
graphics.off()
windows(width=15, height=10)
qplot(as.Date(inc$dt), inc$events/inc$py,  geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests")+
  geom_errorbar(aes(ymin=exp(log(inc$events/inc$py)- qnorm(0.975)/sqrt(inc$events)), ymax=exp(log(inc$events/inc$py)+ qnorm(0.975)/sqrt(inc$events))), width=.1) +
  geom_line(size=0.1) +
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()
savePlot("plot_overall_ir.jpg",type="jpg")
save.image()


#Mixed model
#ff<-join(n_psa_id_e, inc_arzt)
#ff<-join(ff, n_psa_id2_e)

el2_2010$c_dt<-findInterval(as.Date(el2_2010$valid_from_dt), c(as.numeric(as.Date("2010-01-01")), as.numeric(as.Date("2012-07-01")), as.numeric(as.Date("2014-05-01"))))
el2_2010$ev_per<-as.numeric(el2_2010$ev_1)*(el2_2010$c_dt==1)+as.numeric(el2_2010$ev_2)*(el2_2010$c_dt==2)+as.numeric(el2_2010$ev_3)*(el2_2010$c_dt==3)

el2_2010$age_psa<-floor(year(el2_2010$valid_from_dt)-el2_2010$doby)
el2_2010<-join(el2_2010, count_dis_e)
el2_2010$ch_per<-el2_2010$count_cron_2*(el2_2010$c_dt==1)+el2_2010$count_cron_3*(el2_2010$c_dt==2)+el2_2010$count_cron*(el2_2010$c_dt==3)

tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$c_dt[el2_2010$c_dt>=1], mean)
tapply(el2_2010$psa, el2_2010$age_psa, mean)
tapply(el2_2010$psa, el2_2010$ev_per, mean)
tapply(el2_2010$psa, cut(el2_2010$ch_per, c(0,1,4,11, 20), right=FALSE), mean)
#interesting
tapply(el2_2010$psa[el2_2010$c_dt>=1], el2_2010$ev_1[el2_2010$c_dt>=1]+el2_2010$ev_2[el2_2010$c_dt>=1]+el2_2010$ev_3[el2_2010$c_dt>=1], mean)
#summary(lm(psa~n_psa, el2_2010))

ggplot(el2_2010[el2_2010$c_dt>=1], aes(arzt_id[el2_2010$c_dt>=1], log(psa)[el2_2010$c_dt>=1])) + 
  geom_boxplot(aes(fill = factor(c_dt[el2_2010$c_dt>=1])))+
  theme_bw()


kruskal.test(psa~factor(arzt_id), data=el2_2010)
kruskal.test(log(psa)~factor(c_dt), data=el2_2010)
kruskal.test(log(psa)~factor(pat_id), data=el2_2010)


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




















#####
el2_2010<-orderBy(~pat_id, el2_2010)

v_e2 <- c(1, rep(0, length(el2_2010$pat_id)-1))
for (i in 2:length(el2_2010$pat_id))
  if ((el2_2010$pat_id[i] %in% 1:el2_2010$pat_id[i-1]) == FALSE)
    v_e2[i] <- 1;

dt<-as.character(c(as.Date("2010-01-01"), as.Date("2012-07-01"), as.Date("2014-05-01"), as.Date("2017-01-01")))

age_1012<-floor(year(as.Date(pmax(el2_2010$fup_start_dt_3, "2010-01-01"))) - el2_2010$doby)>=55 & floor(year(as.Date(pmax(el2_2010$fup_start_dt_3, "2010-01-01"))) - el2_2010$doby)<=75 
age_1214<-floor(year(as.Date(pmax(el2_2010$fup_start_dt_3, "2012-07-01"))) - el2_2010$doby)>=55 & floor(year(as.Date(pmax(el2_2010$fup_start_dt_3, "2010-01-01"))) - el2_2010$doby)<=75 
age_1416<-floor(year(as.Date(pmax(el2_2010$fup_start_dt_3, "2014-05-01"))) - el2_2010$doby)>=55 & floor(year(as.Date(pmax(el2_2010$fup_start_dt_3, "2010-01-01"))) - el2_2010$doby)<=75 

el2_2010$py_2012<-v_e2*age_1012*(pmax(0,
              as.Date(pmin(el2_2010$fup_end_dt_3, "2012-07-01")) -
                as.Date(pmax(el2_2010$fup_start_dt_3, "2010-01-01")))/365)
el2_2010$py_2014<-v_e2*age_1214*(pmax(0,
                      as.Date(pmin(el2_2010$fup_end_dt_3, "2014-05-01")) -
                        as.Date(pmax(el2_2010$fup_start_dt_3, "2012-07-01")))/365)
el2_2010$py_2016<-v_e2*age_1416*(pmax(0,
                      as.Date(pmin(el2_2010$fup_end_dt_3, "2017-01-01")) -
                        as.Date(pmax(el2_2010$fup_start_dt_3, "2014-05-01")))/365)

el2_2010$ev_2012<-(el2_2010$valid_from_dt >= "2010-01-01") & (el2_2010$age_psa >= 55) & (el2_2010$age_psa <= 75) &
  (el2_2010$valid_from_dt < "2012-07-01") &
  (el2_2010$valid_from_dt >= el2_2010$fup_start_dt_3) &
  (el2_2010$valid_from_dt < el2_2010$fup_end_dt_3)

el2_2010$ev_2014<-(el2_2010$valid_from_dt >= "2012-07-01") & (el2_2010$age_psa >= 55) & (el2_2010$age_psa <= 75) &
  (el2_2010$valid_from_dt < "2014-05-01") &
  (el2_2010$valid_from_dt >= el2_2010$fup_start_dt_3) &
  (el2_2010$valid_from_dt < el2_2010$fup_end_dt_3)

el2_2010$ev_2016<-(el2_2010$valid_from_dt >= "2014-05-01") & (el2_2010$age_psa >= 55) & (el2_2010$age_psa <= 75) &
  (el2_2010$valid_from_dt < "2017-01-01") &
  (el2_2010$valid_from_dt >= el2_2010$fup_start_dt_3) &
  (el2_2010$valid_from_dt < el2_2010$fup_end_dt_3)


inc_arzt<-data.frame(cbind(arzt_id=as.vector(rownames(tapply(el2_2010$py_2012, el2_2010$arzt_id, sum))), py_1=as.numeric(tapply(el2_2010$py_2012, el2_2010$arzt_id, sum)), py_2=as.numeric(tapply(el2_2010$py_2014, el2_2010$arzt_id, sum)), py_3=as.numeric(tapply(el2_2010$py_2016, el2_2010$arzt_id, sum)), ev_1=as.numeric(tapply(el2_2010$ev_2012, el2_2010$arzt_id, sum, na.rm=TRUE)), ev_2=as.numeric(tapply(el2_2010$ev_2014, el2_2010$arzt_id, sum, na.rm=TRUE)), ev_3=as.numeric(tapply(el2_2010$ev_2016, el2_2010$arzt_id, sum, na.rm=TRUE))), stringsAsFactors=FALSE)
inc_arzt$py<-as.numeric(inc_arzt$py_1)+as.numeric(inc_arzt$py_2)+as.numeric(inc_arzt$py_3)
inc_arzt$ev<-as.numeric(inc_arzt$ev_1)+as.numeric(inc_arzt$ev_2)+as.numeric(inc_arzt$ev_3)
inc_arzt$inc<-inc_arzt$ev/inc_arzt$py
inc_arzt$inc1<-as.numeric(inc_arzt$ev_1)/as.numeric(inc_arzt$py_1)
inc_arzt$inc2<-as.numeric(inc_arzt$ev_2)/as.numeric(inc_arzt$py_2)
inc_arzt$inc3<-as.numeric(inc_arzt$ev_3)/as.numeric(inc_arzt$py_3)

