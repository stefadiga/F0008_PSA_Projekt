
# ----------------------------------------------------------------
# Mixed model PSA Value
# ----------------------------------------------------------------

a001a0_3 <- subset(data, !is.na(data$psa_value) & data$n_psa>=6)


xyplot(psa_value~start_dt|as.factor(pat_id), data=a001a0_3,
       type=c("p","g","r"),col="dark blue",col.line="black",
       xlab="start_dt",
       ylab="PSA value")


# -------------------------------

#plot(n_psa_start_dt_s$start_dt,n_psa_start_dt_s$n_psa_start_dt)
#lines(lowess(n_psa_start_dt_s$start_dt,n_psa_start_dt_s$n_psa_start_dt), col="blue")

# Plot Study End
qplot(as.Date(n_psa_start_dt_s$start_dt),n_psa_start_dt_s$n_psa_start_dt, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests")

#ggplot(n_psa_start_dt_s, aes(as.Date(n_psa_start_dt_s$start_dt),n_psa_start_dt_s$n_psa_start_dt)) + geom_point() +
#  geom_smooth(method = 'loess', aes(x=as.Date(n_psa_start_dt_s$start_dt), y=n_psa_start_dt_s$n_psa_start_dt), se = FALSE) 

qplot(as.Date(n_psa_start_dt_s$start_dt),n_psa_start_dt_s$ir, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests")

#Plot Observed
qplot(as.Date(n_psa_start_dt_s1$start_dt),n_psa_start_dt_s1$n_psa_start_dt, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests")

qplot(as.Date(n_psa_start_dt_s1$start_dt),n_psa_start_dt_s1$ir, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests")

#tabular ohne missing
table(practice_type, exclude= "")
table(employ_status, exclude= "")



