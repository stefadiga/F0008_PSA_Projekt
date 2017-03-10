
# ----------------------------------------------------------------
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

#plot(n_psa_start_dt_s$start_dt,n_psa_start_dt_s$n_psa_start_dt)
#lines(lowess(n_psa_start_dt_s$start_dt,n_psa_start_dt_s$n_psa_start_dt), col="blue")

# Incidence
i1<-qplot(as.Date(inc$dt), inc$events/inc$py,  geom=c("point", "smooth"),
            xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed")
i2<-qplot(as.Date(inc$dt), inc$events.1/inc$py.1,  geom=c("point", "smooth"),
          xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End")
i3<-qplot(as.Date(inc$dt), inc$events.2/inc$py.2,  geom=c("point", "smooth"),
          xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed")
grid.arrange(i1, i2, i3, ncol=2)
#ggplot(n_psa_start_dt_s, aes(as.Date(n_psa_start_dt_s$start_dt),n_psa_start_dt_s$n_psa_start_dt)) + geom_point() +
#  geom_smooth(method = 'loess', aes(x=as.Date(n_psa_start_dt_s$start_dt), y=n_psa_start_dt_s$n_psa_start_dt), se = FALSE) 

qplot(as.Date(n_psa_start_dt$start_dt_1),n_psa_start_dt$ir1, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study end")

#Plot Observed
qplot(as.Date(n_psa_start_dt$start_dt_3),n_psa_start_dt$n_psa_start_dt_3, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests: Observed")

qplot(as.Date(n_psa_start_dt$start_dt_3),n_psa_start_dt$ir3, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed")

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
