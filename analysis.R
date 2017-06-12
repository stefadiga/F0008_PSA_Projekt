# DB to use
######
# data_el: complete db
# el2: eligible for follow - up mixed method
# inc: overall incidence rates for the three methods
# inc_arzt_"" e inc_pat_"" db: for incidence rates graphs by categories _"" and computation
# a001a0_3: db for PSA value by patient graph
# a001a0_4: db with Arzt screening and patient rate calculation 
# n_psa_id for log(psa) graphs
# n_psa_id_c for log(psa) chronic diseases
# tab# for report summary
# quant_"" for quantile calculation (old version)


#######
# PSA by patient ID

graphics.off()
windows(width=15, height=10)
xyplot(log(psa)~as.Date(valid_from_dt)|as.factor(pat_id), data=a001a0_3,
       type=c("p","g","r"),col="dark blue",col.line="black", 
       xlab="start_dt",
       ylab="log PSA value",
       main="PSA value (log) by patient id")
savePlot("plot_psa_id.jpg",type="jpg")
save.image()

#Mixed model
fm_mix<-lmer(log(psa)~(as.Date(valid_from_dt))+as.Date(valid_from_dt)|pat_id, a001a0_3)
plot(fm_mix, sqrt(abs(resid(.)))~fitted(.), type=c("p", "smooth"))

################
#Plot Incidence
##############

#Overall (three methods - Fig 3 report)
graphics.off()
windows(width=15, height=10)
i1<-qplot(as.Date(inc$dt), inc$events/inc$py,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Observed") +
  geom_errorbar(aes(ymin=exp(log(inc$events/inc$py)- qnorm(0.975)/sqrt(inc$events)), ymax=exp(log(inc$events/inc$py)+ qnorm(0.975)/sqrt(inc$events))), width=.1) +
  geom_line(size=0.1) +
  theme_bw()
i2<-qplot(as.Date(inc$dt), inc$events.1/inc$py.1,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Study End")+
  geom_errorbar(aes(ymin=exp(log(inc$events.1/inc$py.1)- qnorm(0.975)/sqrt(inc$events.1)), ymax=exp(log(inc$events.1/inc$py.1)+ qnorm(0.975)/sqrt(inc$events.1))), width=.1) +
  geom_line(size=0.1) +
  theme_bw()
i3<-qplot(as.Date(inc$dt), inc$events.2/inc$py.2,  geom=c("point", "smooth"),
          ylim=c(0,0.3), xlab="Start_dt", ylab="Incidence of PSA Test", main="PSA tests: Mixed")+
  geom_errorbar(aes(ymin=exp(log(inc$events.2/inc$py.2)- qnorm(0.975)/sqrt(inc$events.2)), ymax=exp(log(inc$events.2/inc$py.2)+ qnorm(0.975)/sqrt(inc$events.2))), width=.1) +
  geom_line(size=0.1) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()
grid.arrange(i1, i2, i3, ncol=2)

savePlot("plot_a00xa1_ir.jpg",type="jpg")
save.image()

#number of patients that enter in the study
sump<-cumsum(tapply(el2$pat_id, cut(as.Date(el2$fup_start_dt_3), as.Date(quarter_starts)), n_distinct))

graphics.off()
windows(width=15, height=10)

qplot(x=as.Date(quarter_starts[1:31]), y=as.vector(sump[1:31])) + ylim(c(0,22000))+
  xlab("Start_dt")+ ylab("Number of patients") + 
  labs(title ="Number of patients entering in the study")+
  geom_point() +
  geom_line()+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()
savePlot("plot_a00xa1cum_ir.jpg", type="jpg")
save.image()
####### Mixed Model
#1 # by sex
graphics.off()
windows(width=15, height=10)

ggplot(inc_arzt_sex, aes(as.Date(dt), as.numeric(events.2/py.2), color=sex,linetype=sex))+
      ylim(c(0,0.3))+ xlab("Start_dt")+ ylab("Incidence of PSA Test") + 
  labs(color = "Sex", title = paste("Incidence rate PSA tests by sex",  # Add a multi-line title
                                  "III Model (Mixed)",
                                  "Error bars represent 95% Confidence Intervals",
                                  sep = "\n")) + 
  geom_errorbar(aes(ymin=exp(log(inc_arzt_sex$events.2/inc_arzt_sex$py.2)- qnorm(0.975)/sqrt(inc_arzt_sex$events.2)), ymax=exp(log(inc_arzt_sex$events.2/inc_arzt_sex$py.2)+ qnorm(0.975)/sqrt(inc_arzt_sex$events.2)))) +
  geom_ribbon(data=inc_arzt_sex, aes(ymin=exp(log(events.2/py.2)- qnorm(0.975)/sqrt(events.2)), ymax=exp(log(events.2/py.2)+ qnorm(0.975)/sqrt(events.2)), fill=sex), alpha=0.3) +
  geom_point() +
  geom_line(show.legend = FALSE)+
  #geom_line(lwd=.8, aes(linetype=sex))+
  #geom_line(position=position_dodge(0.5)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Sex",
                       breaks=c("f", "m"),
                       labels=c("f", "m"))+
  scale_fill_manual(values = alpha(c("pink", "cadetblue1"), 0.3), name  ="Sex",
                 breaks=c("f", "m"),
                 labels=c("f", "m"))+
  theme_bw()
savePlot("plot_a00xa2_ir.jpg", type="jpg")
save.image()

summary(pois.b<-glm(events.2~sex,offset=log(py.2),data=inc_arzt_sex,family="poisson"))


# 2 # by practice
graphics.off()
windows(width=15, height=10)
ggplot(inc_arzt_prac_g, aes(as.Date(dt), as.numeric(events/py), color=type, linetype=type))+ 
      ylim(c(0,0.3)) + 
  xlab("Start_dt") + 
  ylab("Incidence of PSA Test") + 
  labs(color = "Practice type", title = paste("Incidence rate PSA tests by practice",  # Add a multi-line title
                                            "III Model (Mixed)",
                                            "Error bars represent 95% Confidence Intervals",
                                            sep = "\n")) +
  geom_point()+
  #geom_point(position = position_dodge(width = 70))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_prac_g$events/inc_arzt_prac_g$py)- qnorm(0.975)/sqrt(inc_arzt_prac_g$events)), ymax=exp(log(inc_arzt_prac_g$events/inc_arzt_prac_g$py)+ qnorm(0.975)/sqrt(inc_arzt_prac_g$events))), width=.1, position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_arzt_prac_g$events/inc_arzt_prac_g$py)- qnorm(0.975)/sqrt(inc_arzt_prac_g$events)), ymax=exp(log(inc_arzt_prac_g$events/inc_arzt_prac_g$py)+ qnorm(0.975)/sqrt(inc_arzt_prac_g$events)))) +
  geom_ribbon(data=inc_arzt_prac_g, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=type), alpha=0.3) +
  geom_line(show.legend = FALSE)+
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Practice type",
                breaks=c("Einzelpraxis", "Gruppen/Doppelpraxis"),
                labels=c("Einzelpraxis", "Gruppen/Doppelpraxis"))+
  scale_fill_manual(values = c("yellow", "green"), name  ="Practice type",
                    breaks=c("Einzelpraxis", "Gruppen/Doppelpraxis"),
                    labels=c("Einzelpraxis", "Gruppen/Doppelpraxis"))+
  theme_bw()
savePlot("plot_a00xa4_ir.jpg",type="jpg")
save.image()

summary(pois.b<-glm(events~type,offset=log(py),data=inc_arzt_prac_g,family="poisson"))

# 3 # by region
graphics.off()
windows(width=15, height=10)

ggplot(inc_arzt_reg_g, aes(as.Date(dt), as.numeric(events/py),color=region, linetype=region)) +
      ylim(c(0,0.3)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Region", title = paste("Incidence rate PSA tests by region",  # Add a multi-line title
                                                                    "III Model (Mixed)",
                                                                    "Error bars represent 95% Confidence Intervals",
                                                                    sep = "\n"))+
 #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
                        #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))) +
  geom_ribbon(data=inc_arzt_reg_g, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=region), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Region",
                 breaks=c("CEN", "SUB", "PERI"),
                 labels=c("CEN", "SUB", "PERI"))+
  scale_fill_manual(values = c("pink", "green", "lightblue"), name  ="Region",
                    breaks=c("CEN", "SUB", "PERI"),
                    labels=c("CEN", "SUB", "PERI"))+
  theme_bw()
savePlot("plot_a00xa3_ir.jpg",type="jpg")
save.image()

summary(pois.b<-glm(events~region,offset=log(py),data=inc_arzt_reg_g,family="poisson"))

# 4 # by employ status
graphics.off()
windows(width=15, height=10)

ggplot(inc_arzt_emp_g, aes(as.Date(dt), as.numeric(events/py),color=employ, linetype=employ)) +
  ylim(c(0,0.3)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Employ", title = paste("Incidence rate PSA tests by employ status",  # Add a multi-line title
                                                                      "III Model (Mixed)",
                                                                      "Error bars represent 95% Confidence Intervals",
                                                                      sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_arzt_emp_g$events/inc_arzt_emp_g$py)- qnorm(0.975)/sqrt(inc_arzt_emp_g$events)), ymax=exp(log(inc_arzt_emp_g$events/inc_arzt_emp_g$py)+ qnorm(0.975)/sqrt(inc_arzt_emp_g$events)))) +
  geom_ribbon(data=inc_arzt_emp_g, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=employ), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Employ",
                 breaks=c("Angestellt", "Selbständig"),
                 labels=c("Angestellt", "Selbständig"))+
  scale_fill_manual(values = c("orange", "lightblue"), name  ="Employ",
                    breaks=c("Angestellt", "Selbständig"),
                    labels=c("Angestellt", "Selbständig"))+
  theme_bw()
savePlot("plot_a00xa5_ir.jpg",type="jpg")
save.image()
summary(pois.b<-glm(events~employ,offset=log(py),data=inc_arzt_emp_g,family="poisson"))

# 5 # by arzt doby
graphics.off()
windows(width=15, height=10)

ggplot(inc_arzt_doby_g, aes(as.Date(dt), as.numeric(events/py),color=doby, linetype=doby)) +
  ylim(c(0,0.3)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Arzt Jahrgang", title = paste("Incidence rate PSA tests by Doctor's year of birth",  # Add a multi-line title
                                                                      "III Model (Mixed)",
                                                                      "Error bars represent 95% Confidence Intervals",
                                                                      sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_arzt_doby_g$events/inc_arzt_doby_g$py)- qnorm(0.975)/sqrt(inc_arzt_doby_g$events)), ymax=exp(log(inc_arzt_doby_g$events/inc_arzt_doby_g$py)+ qnorm(0.975)/sqrt(inc_arzt_doby_g$events)))) +
  geom_ribbon(data=inc_arzt_doby_g, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=doby), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Arzt Jahrgang",
                 breaks=c("[1940,1960)", "[1960,1970)", "[1970,1990)"),
                 labels=c("[1940,1960)", "[1960,1970)", "[1970,1990)"))+
  scale_fill_manual(values = c("pink", "green", "lightblue"), name  ="Arzt Jahrgang",
                    breaks=c("[1940,1960)", "[1960,1970)", "[1970,1990)"),
                    labels=c("[1940,1960)", "[1960,1970)", "[1970,1990)"))+
  theme_bw()
savePlot("plot_a00xa6_ir.jpg",type="jpg")
save.image()
summary(pois.b<-glm(events~doby,offset=log(py),data=inc_arzt_doby_g,family="poisson"))

# 6 # by age (at study entry)
graphics.off()
windows(width=15, height=10)

ggplot(inc_pat_age, aes(as.Date(dt), as.numeric(events.2/py.2),color=age, linetype=age)) +
  ylim(c(0,0.3)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Alter Patient", title = paste("Incidence rate PSA tests by patient's age at study entry",  # Add a multi-line title
                                                                             "III Model (Mixed)",
                                                                             "Error bars represent 95% Confidence Intervals",
                                                                             sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_pat_age$events.2/inc_pat_age$py.2)- qnorm(0.975)/sqrt(inc_pat_age$events.2)), ymax=exp(log(inc_pat_age$events.2/inc_pat_age$py.2)+ qnorm(0.975)/sqrt(inc_pat_age$events.2)))) +
  geom_ribbon(data=inc_pat_age, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=age), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Alter Patient",
                 breaks=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"),
                 labels=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"))+
  scale_fill_manual(values = c("pink", "yellow", "green", "lightblue", "violet"), name  ="Alter Patient",
                    breaks=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"),
                    labels=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"))+
  theme_bw()

savePlot("plot_a00xa7.0_ir.jpg",type="jpg")
save.image()
summary(pois.b<-glm(events.2~age,offset=log(py.2),data=inc_pat_age,family="poisson"))

# 7 # age during study
inc_pat_age_c<-subset(inc_pat_age_c, class_age!="[-55)")
inc_pat_age_c$inc<-as.numeric(inc_pat_age_c$events/inc_pat_age_c$py)

graphics.off()
windows(width=15, height=10)


ggplot(inc_pat_age_c, aes(as.Date(dt), inc, linetype=class_age, color=class_age)) +
  ylim(c(0,0.3)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Patient Alter", title = paste("Incidence rate PSA tests by patient's age during study",  # Add a multi-line title
                                                                             "III Model (Mixed)",
                                                                             "Error bars represent 95% Confidence Intervals",
                                                                             sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_pat_age_c$inc)- qnorm(0.975)/sqrt(inc_pat_age_c$events)), ymax=exp(log(inc_pat_age_c$inc)+ qnorm(0.975)/sqrt(inc_pat_age_c$events)))) +
  geom_ribbon(data=inc_pat_age_c, aes(ymin=exp(log(inc)- qnorm(0.975)/sqrt(events)), ymax=exp(log(inc)+ qnorm(0.975)/sqrt(events)), fill=class_age), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Patient Alter",
                 breaks=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"),
                 labels=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"))+
  scale_fill_manual(values = c("pink", "yellow", "green", "lightblue", "violet"), name  ="Patient Alter",
                    breaks=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"),
                    labels=c("[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75+)"))+
  
 facet_grid(.~class_age)+
theme_bw()
savePlot("plot_a00xa7_ir.jpg",type="jpg")
save.image()
summary(pois.b<-glm(events~class_age,offset=log(py),data=inc_pat_age_c,family="poisson"))

# 8 # by arzt freq
graphics.off()
windows(width=15, height=10)

ggplot(inc_arzt_high, aes(as.Date(dt), y=events/py, color=freq, linetype=freq)) +
  ylim(c(0,0.45)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Screening doctor rate", title = paste("Incidence rate PSA tests by frequency of screening",  # Add a multi-line title
                                                                             "III Model (Mixed) - High rate doctor= number of PSA >= 10 per year of follow up",
                                                                             "Error bars represent 95% Confidence Intervals",
                                                                             sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_arzt_high$events/inc_arzt_high$py)- qnorm(0.975)/sqrt(inc_arzt_high$events)), ymax=exp(log(inc_arzt_high$events/inc_arzt_high$py)+ qnorm(0.975)/sqrt(inc_arzt_high$events)))) +
  geom_ribbon(data=inc_arzt_high, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=freq), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Screening doctor rate",
                 breaks=c("l", "h"),
                 labels=c("l", "h"))+
  scale_fill_manual(values = c("blue", "red"), name  ="Screening doctor rate",
                    breaks=c("l", "h"),
                    labels=c("l", "h"))+
  
   theme_bw()
savePlot("plot_a00xa8_ir.jpg",type="jpg")
save.image()

summary(pois.b<-glm(events~freq,offset=log(py),data=inc_arzt_high,family="poisson"))


# 9 # by pat freq

graphics.off()
windows(width=15, height=10)

ggplot(inc_paz_high, aes(as.Date(dt), y=events/py, color=freq, linetype=freq)) +
  ylim(c(0,0.45)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Konsultation patient rate", title = paste("Incidence rate PSA tests by frequency of patient consultation",  # Add a multi-line title
                                                                              "III Model (Mixed) - High rate patient = number of konsultation >= 10 per year of follow up",
                                                                              "Error bars represent 95% Confidence Intervals",
                                                                              sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_paz_high$events/inc_paz_high$py)- qnorm(0.975)/sqrt(inc_paz_high$events)), ymax=exp(log(inc_paz_high$events/inc_paz_high$py)+ qnorm(0.975)/sqrt(inc_paz_high$events)))) +
  geom_ribbon(data=inc_paz_high, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=freq), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Konsultation patient rate",
                 breaks=c("l", "h"),
                 labels=c("l", "h"))+
  scale_fill_manual(values = c("blue", "red"), name  ="Konsultation patient rate",
                    breaks=c("l", "h"),
                    labels=c("l", "h"))+
  
  theme_bw()
savePlot("plot_a00xa9_ir.jpg",type="jpg")
save.image()

summary(pois.b<-glm(events~freq,offset=log(py),data=inc_paz_high,family="poisson"))

graphics.off()
windows(width=15, height=10)

ggplot(inc_paz1_high, aes(as.Date(dt), y=events/py, color=freq, linetype=freq)) +
  ylim(c(0,0.45)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Screening patient rate", title = paste("Incidence rate PSA tests by frequency of patient screening",  # Add a multi-line title
                                                                                         "III Model (Mixed) - High rate patient = number of PSA / n.Konsultation >=0.1",
                                                                                         "Error bars represent 95% Confidence Intervals",
                                                                                         sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_paz1_high$events/inc_paz1_high$py)- qnorm(0.975)/sqrt(inc_paz1_high$events)), ymax=exp(log(inc_paz1_high$events/inc_paz1_high$py)+ qnorm(0.975)/sqrt(inc_paz1_high$events)))) +
  geom_ribbon(data=inc_paz1_high, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=freq), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Screening patient rate",
                 breaks=c("l", "h"),
                 labels=c("l", "h"))+
  scale_fill_manual(values = c("blue", "red"), name  ="Screening patient rate",
                    breaks=c("l", "h"),
                    labels=c("l", "h"))+
  
  theme_bw()
savePlot("plot_a00xa9.2_ir.jpg",type="jpg")
save.image()

# 10 # by chronic disease
graphics.off()
windows(width=15, height=10)

ggplot(inc_pat_cron, aes(as.Date(dt), as.numeric(events/py),color=cronic, linetype=cronic)) +
  ylim(c(0,0.3)) +
  #geom_point(position = position_dodge(width = 70))+  
  geom_point()+
  geom_line(show.legend = FALSE)+
  xlab("Start_dt") +
  ylab("Incidence of PSA Test") +labs(color = "Chronic diseases", title = paste("Incidence rate PSA tests by patient's number of chronic diseases",  # Add a multi-line title
                                                                             "III Model (Mixed)",
                                                                             "Error bars represent 95% Confidence Intervals",
                                                                             sep = "\n"))+
  #geom_errorbar(aes(ymin=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)- qnorm(0.975)/sqrt(inc_arzt_reg_g$events)), ymax=exp(log(inc_arzt_reg_g$events/inc_arzt_reg_g$py)+ qnorm(0.975)/sqrt(inc_arzt_reg_g$events)))
  #, width=.1,position = position_dodge(width = 70)) +
  geom_errorbar(aes(ymin=exp(log(inc_pat_cron$events/inc_pat_cron$py)- qnorm(0.975)/sqrt(inc_pat_cron$events)), ymax=exp(log(inc_pat_cron$events/inc_pat_cron$py)+ qnorm(0.975)/sqrt(inc_pat_cron$events)))) +
  geom_ribbon(data=inc_pat_cron, aes(ymin=exp(log(events/py)- qnorm(0.975)/sqrt(events)), ymax=exp(log(events/py)+ qnorm(0.975)/sqrt(events)), fill=cronic), alpha=0.3) +
  #geom_line(size=0.1,position = position_dodge(width = 70)) +
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  scale_linetype(name  ="Chronic diseases",
                 breaks=c("0", "[1,3]", "[4,10]", "[11+)"),
                 labels=c("0", "[1,3]", "[4,10]", "[11+)"))+
  scale_fill_manual(values = c("pink", "yellow", "green", "lightblue"), name  ="Chronic diseases",
                    breaks=c("0", "[1,3]", "[4,10]", "[11+)"),
                    labels=c("0", "[1,3]", "[4,10]", "[11+)"))+
  theme_bw()

savePlot("plot_a00xa10.0_ir.jpg",type="jpg")
save.image()

#glm comparison
summary(pois.b<-glm(events~factor(cronic),offset=log(py),data=subset(inc_pat_cron, cronic!="[11+)"),family="poisson"))


#Other options

#Plot Frequency
qplot(as.Date(inc$dt),inc$events, geom=c("point", "smooth"),
      xlab="Start_dt", ylab="Frequency of PSA Test", main="PSA tests")


######## PSA
graphics.off()
windows(width=15, height=10)

qplot(data=quant, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","line"),
      ylim(c(0,8)), xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  ylab(expression(paste('PSA (' , mu,g/l,')'))) +
  theme_bw()

savePlot("a00xa0_3_1_psa_quantile.jpg",type="jpg")
save.image()

#Boxplot
graphics.off()
windows(width=15, height=10)

boxplot(n_psa_id$psa,
        notch = TRUE, col=c("red"), ylim=c(0,15), main="Boxplot of PSA")

savePlot("a00xa0_3_1_2_psa_quantile.jpg",type="jpg")
save.image()

graphics.off()
windows(width=15, height=10)

boxplot(log(n_psa_id$psa),
        notch = TRUE, col=c("red"), main="Boxplot of log PSA")
savePlot("a00xa0_3_1_3_psa_quantile.jpg",type="jpg")
save.image()

# PSA hist plot

graphics.off()
windows(width=15, height=10)

hist(log(n_psa_id$psa), density=20, breaks=20, prob=TRUE, xlab="log(PSA)", col="green", main="Histogram of log PSA (Normal distribution is plotted over histogram)")
curve(dnorm(x, mean=mean(log(n_psa_id$psa)), sd=sqrt(var(log(n_psa_id$psa)))), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#vioplot(log(n_psa_id$psa), col="gold")
#title("Violin Plots of log PSA")

savePlot("a00xa0_log_psa_hist.jpg",type="jpg")
save.image()

#boxplot by trimester
graphics.off()
windows(width=15, height=10)
boxplot(log(n_psa_id$psa)~cut(as.Date(n_psa_id$valid_from_dt), as.Date(quarter_starts)), notch=T)

savePlot("a00xa0_3_1_5_psa_quantile.jpg",type="jpg")
save.image()

# by arzt sex
graphics.off()
windows(width=15, height=10)
ggplot(n_psa_id, aes(dt, log(psa))) + 
geom_boxplot(aes(fill = arzt_sex)) +
  #facet_grid(.~arzt_sex)+
#geom_boxplot(aes(fill = arzt_sex))+
#geom_boxplot(aes(fill = arzt_sex), notch = TRUE, notchwidth = 1)+
xlab("Start_dt") +
ylab("log(PSA)") +
theme_bw()
savePlot("a00xa0_3_1_5_b_psa_quantile.jpg",type="jpg")
save.image()

# 1 arzt sex plot + smooth
graphics.off()
windows(width=15, height=10)

ggplot(n_psa_id, aes(y=log(psa), x=as.Date(valid_from_dt), color=arzt_sex))+
         ylim(c(-3,+2))+
         ylab("log PSA")+ 
         xlab("dt")+
        #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()
# facet_grid(.~arzt_sex)
savePlot("a00xa0_3_1_6_psa_quantile.jpg",type="jpg")
save.image()


#by sex
graphics.off()
windows(width=15, height=10)

qplot(data=quant_sex, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ sex)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))+
theme_bw()
savePlot("plot_a002a2_3_psa_quantile",type="jpg")
save.image()

# by region

n_psa_id_s<- subset(n_psa_id, arzt_region_code=="CEN"|arzt_region_code=="SUB"|arzt_region_code=="PERI")

ggplot(n_psa_id_s, aes(cut(as.Date(valid_from_dt), as.Date(quarter_starts)), log(psa))) + 
  geom_boxplot(aes(fill = arzt_region_code))+
  ylab("log PSA")+ 
  xlab("dt")+
  theme_bw()

graphics.off()
windows(width=15, height=10)

# 2 # by region

ggplot(n_psa_id_s, aes(y=log(psa), x=as.Date(valid_from_dt), color=arzt_region_code))+
  ylab("log PSA")+ 
  xlab("dt")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a3_3_psa_quantile",type="jpg")
save.image()

qplot(data=quant_reg, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", main="Study period follow up", ylim=c(0,15)) + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ region)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))+
  theme_bw()


#by employ status

ggplot(n_psa_id, aes(cut(as.Date(valid_from_dt), as.Date(quarter_starts)), log(psa))) + 
  geom_boxplot(aes(fill = employ_status))+
  theme_bw()


qplot(data=quant_emp, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ employ)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))+
  theme_bw()

n_psa_id_e<-subset(n_psa_id, employ_status!="")

graphics.off()
windows(width=15, height=10)

# 3 # by employ status

ggplot(n_psa_id_e, aes(y=log(psa), x=as.Date(valid_from_dt), color=employ_status))+
  ylab("log PSA")+ 
  xlab("dt")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a5_3_psa_quantile.jpg",type="jpg")
save.image()

#by praxis
graphics.off()
windows(width=15, height=10)

ggplot(n_psa_id, aes(cut(as.Date(valid_from_dt), as.Date(quarter_starts)), log(psa))) + 
  geom_boxplot(aes(fill = practice_type))+
  theme_bw()

qplot(data=quant_pra, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ practice)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))+
theme_bw()

n_psa_id$p<-n_psa_id$practice_type
n_psa_id$p[n_psa_id$p=="Einzelpraxis"| n_psa_id$p=="Gruppenpraxis"] <-"Einzel/Gruppen Praxis"
n_psa_id_p<-subset(n_psa_id, practice_type!="")

# 4 by praxis
graphics.off()
windows(width=15, height=10)

ggplot(n_psa_id_p, aes(y=log(psa), x=as.Date(valid_from_dt), color=p))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Praxis type")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a4_3_psa_quantile.jpg",type="jpg")
save.image()

#by arzt age
graphics.off()
windows(width=15, height=10)

ggplot(n_psa_id, aes(cut(as.Date(valid_from_dt), as.Date(quarter_starts)), log(psa))) + 
  geom_boxplot(aes(fill = factor(c_doby)))+
  theme_bw()

qplot(data=quant_adob, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ arzt_doby)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))+
  theme_bw()

n_psa_id$art_a[n_psa_id$c_doby==1 | n_psa_id$c_doby==2] <-"[1940, 1960)"
n_psa_id$art_a[n_psa_id$c_doby==3] <-"[1960, 1970)"
n_psa_id$art_a[n_psa_id$c_doby==4 | n_psa_id$c_doby==5] <-"[1970, 1990)"


n_psa_id_sa<-subset(n_psa_id, !is.na(art_a))
graphics.off()
windows(width=15, height=10)

# 5 by arzt age
ggplot(n_psa_id_sa, aes(y=log(psa), x=as.Date(valid_from_dt), color=art_a))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Jahrgang Arzt")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a6_3_psa_quantile.jpg",type="jpg")
save.image()

#by patient age at study entry

ggplot(n_psa_id, aes(cut(as.Date(valid_from_dt), as.Date(quarter_starts)), log(psa))) + 
  geom_boxplot(aes(fill = c_age))+
  theme_bw()

qplot(data=quant_pdob, x=as.Date(dt), y=value, colour=as.character(psa_quant), geom=c("point","smooth"),
      xlab="dt", ylim=c(0,10), main="Study period follow up") + labs(color = "Quantile") + 
  scale_colour_discrete(labels=c("Median", "75%-Quant", "90%-Quant")) +
  facet_grid(.~ p_age)+
  ylab(expression(paste('PSA (' , mu,g/l,')')))+
theme_bw()

n_psa_id_pt_a<-subset(n_psa_id, !is.na(c_age))

graphics.off()
windows(width=15, height=10)

# 6 by patienta age (at beginning)
ggplot(n_psa_id_pt_a, aes(y=log(psa), x=as.Date(valid_from_dt), color=c_age))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Alter Patient")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a7_3_psa_quantile.jpg",type="jpg")
save.image()

# by arzt screening rate
a001a0_4_b<-subset(a001a0_4, !is.na(c_arzt_f))
a001a0_4_b$c_arzt_f[a001a0_4_b$c_arzt_f==1]<-"l"
a001a0_4_b$c_arzt_f[a001a0_4_b$c_arzt_f==2]<-"h"

graphics.off()
windows(width=15, height=10)

# 7 arzt screening
ggplot(a001a0_4_b, aes(y=log(psa), x=as.Date(valid_from_dt), color=c_arzt_f))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Arzt screening rate")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a8_3_psa_log.jpg",type="jpg")
save.image()

# by pat Kons rate
a001a0_4_bp<-subset(a001a0_4, !is.na(c_pat_k))
a001a0_4_bp$c_pat_k[a001a0_4_bp$c_pat_k==1]<-"l"
a001a0_4_bp$c_pat_k[a001a0_4_bp$c_pat_k==2]<-"h"

graphics.off()
windows(width=15, height=10)

# 8 by pat Kons rate

ggplot(a001a0_4_bp, aes(y=log(psa), x=as.Date(valid_from_dt), color=c_pat_k))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Pat kons rate")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a9_3_psa_log.jpg",type="jpg")
save.image()

#by pat screening rate
a001a0_4_bps<-subset(a001a0_4, !is.na(c_pat_k1))
a001a0_4_bps$c_pat_k1[a001a0_4_bps$c_pat_k1==1]<-"l"
a001a0_4_bps$c_pat_k1[a001a0_4_bps$c_pat_k1==2]<-"h"

graphics.off()
windows(width=15, height=10)

# 9 by pat screening rate

ggplot(a001a0_4_bps, aes(y=log(psa), x=as.Date(valid_from_dt), color=c_pat_k1))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Pat screening rate")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a10_3_psa_log.jpg",type="jpg")
save.image()


graphics.off()
windows(width=15, height=10)

# 10 by chronic diseases

ggplot(n_psa_id_c, aes(y=log(psa), x=as.Date(valid_from_dt), color=c_cron))+
  ylab("log PSA")+ 
  xlab("dt")+
  labs(color="Chronic disease")+
  #geom_point()+
  #geom_smooth(method = "lm", formula = y ~ splines::bs(x, 5))+
  geom_smooth(method=loess)+
  #stat_smooth(method="gam", formula=y~s(x,k=10, bs="cs"), se=FALSE, size=1)+
  geom_vline(xintercept=as.numeric(as.Date("2011-11-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2012-07-01")), linetype=4)+
  geom_vline(xintercept=as.numeric(as.Date("2014-05-01")), linetype=4)+
  theme_bw()

savePlot("plot_a002a10_4_psa_log.jpg",type="jpg")
save.image()

#######
#Report tables
######

#n psa
count_psa<-rbind(data.frame(n_psa=0, n_pat=n_distinct(el$pat_id)-sum(inc$events)), as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_1 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1)]), data_el$pat_id[data_el$elig_2_1 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_1)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName = "n_pat"))
count_psa_2<-rbind(data.frame(n_psa=0, n_pat=n_distinct(el1$pat_id)-sum(inc$events.1)), as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_2 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2)]), data_el$pat_id[data_el$elig_2_2 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_2)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName = "n_pat"))
count_psa_3<-rbind(data.frame(n_psa=0, n_pat=n_distinct(el2$pat_id)-sum(inc$events.2)), as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)]), data_el$pat_id[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName ="n_pat"))
#by arzt
count_psa_3a<-rbind(as.data.frame(table(tapply(!is.na(data_el$psa[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)]), data_el$arzt_id[data_el$elig_2_3 & as.Date(data_el$valid_from_dt) < as.Date(data_el$fup_end_dt_3)], sum, na.rm=T), dnn=c("n_psa")), stringsAsFactors = FALSE, responseName ="n_arzt"))
count_psa_3a<-rbind(data.frame(n_psa=0, n_arzt=n_distinct(el2$arzt_id)-sum(count_psa_3a$n_arzt)), count_psa_3a)

#function for reporting p-values
format_pval <- function(x){
  if (x < .001) return(paste('<', '.001'))
  if (x > .250) return(paste('>', '.250'))
  paste('=', myround(x, 3))   # 3 = no. of digits to round p value to if .001 < p < .250.
}

#incidence overall
tab1<-data.frame(Follow_up=c("Observed", "Study Period", "Mixed"), n_pat=c(n_distinct(el$pat_id), n_distinct(el1$pat_id), n_distinct(el2$pat_id)), psa=c(sum(as.numeric(count_psa$n_psa)*count_psa$n_pat),sum(as.numeric(count_psa_2$n_psa)*count_psa_2$n_pat), sum(as.numeric(count_psa_3$n_psa)*count_psa_3$n_pat)), py=c(sum(as.numeric(data_el$fup_y_1)*data_el$elig_2_1*v), sum(as.numeric(data_el$fup_y_2)*data_el$elig_2_2*v), sum(as.numeric(data_el$fup_y_3)*data_el$elig_2_3*v)), n_arzt=c(n_distinct(el$arzt_id), n_distinct(el1$arzt_id),n_distinct(el2$arzt_id)))
tab1$ir<-as.numeric(tab1$psa)/as.numeric(tab1$py)
#pairwisw poisson rate ratio (reference Observed)
tab1$rate_ratio<-c("-", poisson.test(tab1$psa[1:2], tab1$py[1:2])[5], poisson.test(c(tab1$psa[1], tab1$psa[3]), c(tab1$py[1], tab1$py[3]))[5])
tab1$p_value<-c("-", format_pval(poisson.test(tab1$psa[1:2], tab1$py[1:2])[3]), format_pval(poisson.test(c(tab1$psa[1], tab1$psa[3]), c(tab1$py[1], tab1$py[3]))[3]))

py_date<-c(sum(inc$py.2 *(inc$c_dt==0)), sum(inc$py.2 *(inc$c_dt==1)), sum(inc$py.2 *(inc$c_dt==2)), sum(inc$py.2 *(inc$c_dt==3)))   
psa_date<-c(sum(inc$events.2 *(inc$c_dt==0)), sum(inc$events.2 *(inc$c_dt==1)), sum(inc$events.2 *(inc$c_dt==2)), sum(inc$events.2 *(inc$c_dt==3)))   
tab1$ir_adj<-c("-", "-", format(ageadjust.direct(tab1$psa[3], tab1$py[3], rate = tab8$ir, tapply(mort_2010_m$n.pop_middle, mort_2010_m$c_age, sum, na.rm=TRUE), conf.level = 0.95)[2], digits=4))
tab1$cint_adj<-c("-", "-", paste("(", format(ageadjust.direct(tab1$psa[3], tab1$py[3], rate = tab8$ir, tapply(mort_2010_m$n.pop_middle, mort_2010_m$c_age, sum, na.rm=TRUE), conf.level = 0.95)[3], digits=4), ",", format(ageadjust.direct(tab1$psa[3], tab1$py[3], rate = tab8$ir, tapply(mort_2010_m$n.pop_middle, mort_2010_m$c_age, sum, na.rm=TRUE), conf.level = 0.95)[4], digits=4), ")", sep=""))

#incidence overall by time
tab1.b<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed"), period=c("[01/01/09, 01/11/11)", "[01/11/11, 01/07/12)", "[01/07/12, 01/05/14)", "[01/05/14, 31/12/16]"), 
                   py=py_date, psa=psa_date)
tab1.b$ir<-as.numeric(tab1.b$psa)/as.numeric(tab1.b$py)
resp1<-summary(pois.b<-glm(events.2~factor(c_dt),offset=log(py.2),data=inc,family="poisson"))
tab1.b$rate_ratio<-c("-", format(exp(-as.vector(resp1$coefficients[,1][2:4])), digits=4))
tab1.b$p_value<-c("-", as.character(lapply(as.vector(resp1$coefficients[,4][2:4]), format_pval)))
tab1.b$cint<-c("-", paste("(", cbind(format(exp(-as.vector(confint(glm(events.2~factor(c_dt),offset=log(py.2),data=inc,family="poisson"))[,2][2:4])), digits=4)), ",", cbind(format(exp(-as.vector(confint(glm(events.2~factor(c_dt),offset=log(py.2),data=inc,family="poisson"))[,1][2:4])), digits=4)), ")", sep=""))

#incidence sex
tab2<-data.frame(Follow_up=c("Mixed", "Mixed"), arzt_sex=c("f", "m"), n_pat=c(n_distinct(el2$pat_id[el2$arzt_sex=="f"]),n_distinct(el2$pat_id[el2$arzt_sex=="m"])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$arzt_sex=="f"]),n_distinct(el2$arzt_id[el2$arzt_sex=="m"])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$arzt_sex=="f" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_sex=="m" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_sex$py.2 *(inc_arzt_sex$sex=="f")), sum(inc_arzt_sex$py.2 *(inc_arzt_sex$sex=="m"))))   
tab2$ir<-as.numeric(tab2$n_psa)/as.numeric(tab2$py)
tab2$rate_ratio<-c("-", poisson.test(tab2$n_psa,tab2$py)[5])
tab2$p_value<-c("-", format_pval(poisson.test(tab2$n_psa,tab2$py)[3]))

py_date_f<-c(sum(inc_arzt_sex$py.2*(inc_arzt_sex$c_dt==0)*(inc_arzt_sex$sex=="f")), sum(inc_arzt_sex$py.2*(inc_arzt_sex$c_dt==1)*(inc_arzt_sex$sex=="f")))
psa_date_f<-c(sum(inc_arzt_sex$events.2 *(inc_arzt_sex$c_dt==0)*(inc_arzt_sex$sex=="f")), sum(inc_arzt_sex$events.2 *(inc_arzt_sex$c_dt==1)*(inc_arzt_sex$sex=="f")))   

tab2.b<-data.frame(Follow_up=c("Mixed", "Mixed"), period=c("[01/01/09, 01/05/14)", "[01/05/14, 31/12/16]"), 
                   py=py_date_f, psa=psa_date_f)
tab2.b$ir<-as.numeric(tab2.b$psa)/as.numeric(tab2.b$py)
tab2.b$rate_ratio<-c("-", poisson.test(tab2.b$psa,tab2.b$py)[5])
tab2.b$p_value<-c("-", format_pval(poisson.test(tab2.b$psa,tab2.b$py)[3]))

cbind(summary(glm(events~sex+sex:factor(c_dt),offset=log(py),data=inc_arzt_sex,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~sex+sex:factor(c_dt),offset=log(py),data=inc_arzt_sex,family="poisson"))$coefficients[, 4], format_pval))

#Region
tab3<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed", "Mixed", "Mixed"), region=c("CEN", "IND-Ter", "SUB", "PERI", "OTHER", "NA"), n_pat=c(n_distinct(el2$pat_id[el2$arzt_region_code=="CEN"]), n_distinct(el2$pat_id[el2$arzt_region_code=="IND"]), n_distinct(el2$pat_id[el2$arzt_region_code=="SUB"]), n_distinct(el2$pat_id[el2$arzt_region_code=="PERI"]), n_distinct(el2$pat_id[el2$arzt_region_code=="OTHER"]), n_distinct(el2$pat_id[el2$arzt_region_code==""])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$arzt_region_code=="CEN"]), n_distinct(el2$arzt_id[el2$arzt_region_code=="IND"]), n_distinct(el2$arzt_id[el2$arzt_region_code=="SUB"]),  n_distinct(el2$arzt_id[el2$arzt_region_code=="PERI"]), n_distinct(el2$arzt_id[el2$arzt_region_code=="OTHER"]), n_distinct(el2$arzt_id[el2$arzt_region_code==""])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$arzt_region_code=="CEN" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="IND" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="SUB" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="PERI" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="OTHER" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$arzt_region_code=="" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_reg$py.2 *(inc_arzt_reg$region=="CEN")), sum(inc_arzt_reg$py.2 *(inc_arzt_reg$region=="IND-Ter")), sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="SUB")), sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="PERI")),sum(inc_arzt_reg$py.2*(inc_arzt_reg$region=="OTHER")), sum(inc_arzt_reg$py.2*(inc_arzt_reg$region==""))))   
tab3$ir<-as.numeric(tab3$n_psa)/as.numeric(tab3$py)

resp2<-summary(pois.b<-glm(events~region,offset=log(py),data=inc_arzt_reg_g,family="poisson"))
tab3$rate_ratio<-c("-", " ", format(exp(-as.vector(resp2$coefficients[,1][2:3])), digits=4), " ", " ")
tab3$p_value<-c("-", " ", as.character(lapply(as.vector(resp2$coefficients[,4][2:3]), format_pval)), " ", " ")

cbind(summary(glm(events~region+region:factor(c_dt),offset=log(py),data=inc_arzt_reg_g,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~region+region:factor(c_dt),offset=log(py),data=inc_arzt_reg_g,family="poisson"))$coefficients[, 4], format_pval))


#poisson.test(c(tab3$n_psa[1], tab3$n_psa[3]),c(tab3$py[1], tab3$py[3]))
#poisson.test(c(tab3$n_psa[1], tab3$n_psa[4]),c(tab3$py[1], tab3$py[4]))
#poisson.test(c(tab3$n_psa[4], tab3$n_psa[3]),c(tab3$py[4], tab3$py[3]))

#pairwise.prop.test(x=c(tab3$n_psa[1], tab3$n_psa[3:4]),n=c(tab3$py[1], tab3$py[3:4]),p.adjust.method="bonferroni")

#Praxis
tab4<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed"), practice=c("Einzelpraxis", "Gruppenpraxis", "Doppelpraxis", "NA"), n_pat=c(n_distinct(el2$pat_id[el2$practice_type=="Einzelpraxis"]), n_distinct(el2$pat_id[el2$practice_type=="Gruppenpraxis"]), n_distinct(el2$pat_id[el2$practice_type=="Doppelpraxis"]), n_distinct(el2$pat_id[el2$practice_type==""])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$practice_type=="Einzelpraxis"]), n_distinct(el2$arzt_id[el2$practice_type=="Gruppenpraxis"]), n_distinct(el2$arzt_id[el2$practice_type=="Doppelpraxis"]),  n_distinct(el2$arzt_id[el2$practice_type==""])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$practice_type=="Einzelpraxis" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$practice_type=="Gruppenpraxis" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$practice_type=="Doppelpraxis" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$practice_type=="" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_prac$py.2 *(inc_arzt_prac$type=="Einzelpraxis")), sum(inc_arzt_prac$py.2 *(inc_arzt_prac$type=="Gruppenpraxis")), sum(inc_arzt_prac$py.2*(inc_arzt_prac$type=="Doppelpraxis")), sum(inc_arzt_prac$py.2*(inc_arzt_prac$type==""))))   
tab4$ir<-as.numeric(tab4$n_psa)/as.numeric(tab4$py)

tab4$rate_ratio<-c("-", poisson.test(c(tab4$n_psa[1], tab4$n_psa[2]+tab4$n_psa[3]),c(tab4$py[1], tab4$py[2]+tab4$py[3]))[5], " * ", " ")
tab4$p_value<-c("-", format_pval(poisson.test(c(tab4$n_psa[1], tab4$n_psa[2]+tab4$n_psa[3]),c(tab4$py[1], tab4$py[2]+tab4$py[3]))[3]), " * ", " ")

cbind(summary(glm(events~type+type:factor(c_dt),offset=log(py),data=inc_arzt_prac_g,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~type+type:factor(c_dt),offset=log(py),data=inc_arzt_prac_g,family="poisson"))$coefficients[, 4], format_pval))

#Employ status
tab5<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed"), employ=c("Angestellt", "Selbständig", "Assistentin", "NA"), n_pat=c(n_distinct(el2$pat_id[el2$employ_status=="Angestellt"]),n_distinct(el2$pat_id[el2$employ_status=="Selbständig"]), n_distinct(el2$pat_id[el2$employ_status=="AssistentIn"]), n_distinct(el2$pat_id[el2$employ_status==""])),
                 n_arzt=c(n_distinct(el2$arzt_id[el2$employ_status=="Angestellt"]), n_distinct(el2$arzt_id[el2$employ_status=="Selbständig"]), n_distinct(el2$arzt_id[el2$employ_status=="AssistentIn"]), n_distinct(el2$arzt_id[el2$employ_status==""])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$employ_status=="Angestellt" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$employ_status=="Selbständig" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$employ_status=="AssistentIn" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$employ_status=="" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="Angestellt")), sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="Selbständig")), sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ=="AssistentIn")), sum(inc_arzt_emp$py.2 *(inc_arzt_emp$employ==""))))   
tab5$ir<-as.numeric(tab5$n_psa)/as.numeric(tab5$py)
tab5$rate_ratio<-c("-", poisson.test(c(tab5$n_psa[1], tab5$n_psa[2]),c(tab5$py[1], tab5$py[2]))[5], "", "")
tab5$p_value<-c("-", format_pval(poisson.test(c(tab5$n_psa[1], tab5$n_psa[2]),c(tab5$py[1], tab5$py[2]))[3]), "", "")

cbind(summary(glm(events~employ+employ:factor(c_dt),offset=log(py),data=inc_arzt_emp_g,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~employ+employ:factor(c_dt),offset=log(py),data=inc_arzt_emp_g,family="poisson"))$coefficients[, 4], format_pval))

#Arzt doby

tab6<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed", "Mixed", "Mixed"), arzt_doby=c("[1940,1950)", "[1950,1960)", "[1960,1970)", "[1970,1980)", "[1980,1990)", "NA"), n_pat=c(as.vector(table(el2$c_doby*v2))[2:6], n_distinct(el2$pat_id[is.na(el2$c_doby)])),
                 n_arzt=c(as.vector(tapply(el2$arzt_id, el2$c_doby, n_distinct)), n_distinct(el2$arzt_id[is.na(el2$c_doby)])), 
                 n_psa=c(sum(!is.na(el2$psa[el2$c_doby==1 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==2 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==3 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==4 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_doby==5 & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[is.na(el2$c_doby) & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1940,1950)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1950,1960)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1960,1970)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1970,1980)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby=="[1980,1990)")), sum(inc_arzt_doby$py.2 *(inc_arzt_doby$doby==""))))   
tab6$ir<-as.numeric(tab6$n_psa)/as.numeric(tab6$py)
resp3<-summary(pois.b<-glm(events~doby,offset=log(py),data=inc_arzt_doby_g,family="poisson"))
tab6$rate_ratio<-c("-", "-", format(exp(-as.vector(resp3$coefficients[,1][2:3])), digits=4), " * ", " ")
tab6$p_value<-c("-", "-", as.character(lapply(as.vector(resp3$coefficients[,4][2:3]), format_pval)), " * ", " ")
cbind(summary(glm(events~doby+doby:factor(c_dt),offset=log(py),data=inc_arzt_doby_g,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~doby+doby:factor(c_dt),offset=log(py),data=inc_arzt_doby_g,family="poisson"))$coefficients[, 4], format_pval))

#pairwise.prop.test(x=c(tab6$n_psa[1]+tab6$n_psa[2], tab6$n_psa[3], tab6$n_psa[4]+tab6$n_psa[5]),n=c(tab6$py[1]+tab6$py[2], tab6$py[3], tab6$py[4]+tab6$py[5]),p.adjust.method="bonferroni")
#poisson.test(c(tab6$n_psa[1]+tab6$n_psa[2], tab6$n_psa[3]),c(tab6$py[1]+tab6$py[2], tab6$py[3]))
#poisson.test(c(tab6$n_psa[3], tab6$n_psa[4]+tab6$n_psa[5]),c(tab6$py[3], tab6$py[4]+tab6$py[5]))
#poisson.test(c(tab6$n_psa[1]+tab6$n_psa[2], tab6$n_psa[4]+tab6$n_psa[5]),c(tab6$py[1]+tab6$py[2], tab6$py[4]+tab6$py[5]))

#patient age at study entry
tab7<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed", "Mixed"), pat_age=c("[55,60)", "[60,65)", "[65,70)","[70,75)", "[75+)"), n_pat=c(as.vector(tapply(el2$pat_id, el2$c_age, n_distinct))),
                 n_arzt=c(as.vector(tapply(el2$arzt_id, el2$c_age, n_distinct))), 
                 n_psa=c(sum(!is.na(el2$psa[el2$c_age=="[55,60)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[60,65)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[65,70)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[70,75)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)])), sum(!is.na(el2$psa[el2$c_age=="[75,76)" & as.Date(el2$valid_from_dt) < as.Date(el2$fup_end_dt_3)]))),
                 py=c(sum(inc_pat_age$py.2 *(inc_pat_age$age=="[55,60)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[60,65)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[65,70)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[70,75)")), sum(inc_pat_age$py.2 *(inc_pat_age$age=="[75+)"))))   
tab7$ir<-as.numeric(tab7$n_psa)/as.numeric(tab7$py)
resp4<-summary(pois.b<-glm(events.2~age,offset=log(py.2),data=inc_pat_age,family="poisson"))
tab7$rate_ratio<-c("-", format(exp(-as.vector(resp4$coefficients[,1][2:5])), digits=4))
tab7$p_value<-c("-", as.character(lapply(as.vector(resp4$coefficients[,4][2:5]), format_pval)))
tab7$cint<-c("-", paste("(", cbind(format(exp(-as.vector(confint(glm(events.2~age,offset=log(py.2),data=inc_pat_age,family="poisson"))[,2][2:5])), digits=4)), ",", cbind(format(exp(-as.vector(confint(glm(events.2~age,offset=log(py.2),data=inc_pat_age,family="poisson"))[,1][2:5])), digits=4)), ")", sep=""))

cbind(summary(glm(events.2~age+age:factor(c_dt),offset=log(py.2),data=inc_pat_age,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events.2~age+age:factor(c_dt),offset=log(py.2),data=inc_pat_age,family="poisson"))$coefficients[, 4], format_pval))

#patient age during the study
tab8<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed", "Mixed"), pat_age=c("[55,60)", "[60,65)", "[65,70)","[70,75)", "[75+)"), n_pat=c(as.vector(tapply(el2$pat_id, el2$ind55, n_distinct))[2], as.vector(tapply(el2$pat_id, el2$ind60, n_distinct))[2], as.vector(tapply(el2$pat_id, el2$ind65, n_distinct))[2], as.vector(tapply(el2$pat_id, el2$ind70, n_distinct))[2], as.vector(tapply(el2$pat_id, el2$ind75, n_distinct))[2]),
                 n_arzt=c(as.vector(tapply(el2$arzt_id, el2$ind55, n_distinct))[2], as.vector(tapply(el2$arzt_id, el2$ind60, n_distinct))[2], as.vector(tapply(el2$arzt_id, el2$ind65, n_distinct))[2], as.vector(tapply(el2$arzt_id, el2$ind70, n_distinct))[2], as.vector(tapply(el2$arzt_id, el2$ind75, n_distinct))[2]), 
                 n_psa=c(sum(inc_pat_age_c$events[inc_pat_age_c$class_age=="[55,60)"]), sum(inc_pat_age_c$events[inc_pat_age_c$class_age=="[60,65)"]), sum(inc_pat_age_c$events[inc_pat_age_c$class_age=="[65,70)"]), sum(inc_pat_age_c$events[inc_pat_age_c$class_age=="[70,75)"]), sum(inc_pat_age_c$events[inc_pat_age_c$class_age=="[75+)"])),
                 py=c(sum(inc_pat_age_c$py[inc_pat_age_c$class_age=="[55,60)"]), sum(inc_pat_age_c$py[inc_pat_age_c$class_age=="[60,65)"]), sum(inc_pat_age_c$py[inc_pat_age_c$class_age=="[65,70)"]), sum(inc_pat_age_c$py[inc_pat_age_c$class_age=="[70,75)"]), sum(inc_pat_age_c$py[inc_pat_age_c$class_age=="[75+)"])))   
tab8$ir<-as.numeric(tab8$n_psa)/as.numeric(tab8$py)
tab8$ir_adj<-tab8$ir*wei_ad
resp5<-summary(pois.b<-glm(events~class_age,offset=log(py),data=inc_pat_age_c,family="poisson"))
tab8$rate_ratio<-c("-", format(exp(-as.vector(resp5$coefficients[,1][2:5])), digits=4))
tab8$p_value<-c("-", as.character(lapply(as.vector(resp5$coefficients[,4][2:5]), format_pval)))
tab8$cint<-c("-", paste("(", cbind(format(exp(-as.vector(confint(glm(events~class_age,offset=log(py),data=inc_pat_age_c,family="poisson"))[,2][2:5])), digits=4)), ",", cbind(format(exp(-as.vector(confint(glm(events~class_age,offset=log(py),data=inc_pat_age_c,family="poisson"))[,1][2:5])), digits=4)), ")", sep=""))

sum(tab8$ir*wei_ad)

cbind(summary(glm(events~class_age+class_age:factor(c_dt),offset=log(py),data=inc_pat_age_c,family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~class_age+class_age:factor(c_dt),offset=log(py),data=inc_pat_age_c,family="poisson"))$coefficients[, 4], format_pval))

#Screening doctor rate
tab9<-data.frame(Follow_up=c("Mixed", "Mixed"), arzt_freq=c("l", "h"), n_pat=as.vector(tapply(a001a0_4$pat_id, a001a0_4$c_arzt_f, n_distinct)),
                 n_arzt=as.vector(table(n_psa_arzt$c_arzt_f)), 
                 n_psa=c(sum(inc_arzt_high$events[inc_arzt_high$freq=="l"]), sum(inc_arzt_high$events[inc_arzt_high$freq=="h"])), 
                 py=c(sum(inc_arzt_high$py[inc_arzt_high$freq=="l"]), sum(inc_arzt_high$py [inc_arzt_high$freq=="h"])))   
tab9$ir<-as.numeric(tab9$n_psa)/as.numeric(tab9$py)
tab9$rate_ratio<-c("-", poisson.test(tab9$n_psa, tab9$py)[5])
tab9$p_value<-c("-", format_pval(poisson.test(tab9$n_psa, tab9$py)[3]))
cbind(summary(glm(events~freq+freq:factor(c_dt),offset=log(py),data=inc_arzt_high, family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~freq+freq:factor(c_dt),offset=log(py),data=inc_arzt_high,family="poisson"))$coefficients[, 4], format_pval))


#Screening patient rate
tab10<-data.frame(Follow_up=c("Mixed", "Mixed"), pat_rate=c("l", "h"), n_pat=as.vector(tapply(a001a0_4$pat_id, a001a0_4$c_pat_k, n_distinct)),
                 n_arzt=as.vector(tapply(a001a0_4$arzt_id, a001a0_4$c_pat_k, n_distinct)), 
                 n_psa=c(sum(inc_paz_high$events[inc_paz_high$freq=="l"]), sum(inc_paz_high$events[inc_paz_high$freq=="h"])), 
                 py=c(sum(inc_paz_high$py[inc_paz_high$freq=="l"]), sum(inc_paz_high$py [inc_paz_high$freq=="h"])))   
tab10$ir<-as.numeric(tab10$n_psa)/as.numeric(tab10$py)
tab10$rate_ratio<-c("-", poisson.test(tab10$n_psa,tab10$py)[5])
tab10$p_value<-c("-", format_pval(poisson.test(tab10$n_psa, tab10$py)[3]))
cbind(summary(glm(events~freq+freq:factor(c_dt),offset=log(py),data=inc_paz_high, family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~freq+freq:factor(c_dt),offset=log(py),data=inc_paz_high,family="poisson"))$coefficients[, 4], format_pval))

tab11<-data.frame(Follow_up=c("Mixed", "Mixed"), pat_rate=c("l", "h"), n_pat=as.vector(tapply(a001a0_4$pat_id, a001a0_4$c_pat_k1, n_distinct)),
                  n_arzt=as.vector(tapply(a001a0_4$arzt_id, a001a0_4$c_pat_k1, n_distinct)), 
                  n_psa=c(sum(inc_paz1_high$events[inc_paz1_high$freq=="l"]), sum(inc_paz1_high$events[inc_paz1_high$freq=="h"])), 
                  py=c(sum(inc_paz1_high$py[inc_paz1_high$freq=="l"]), sum(inc_paz1_high$py [inc_paz1_high$freq=="h"])))   

tab11$ir<-as.numeric(tab11$n_psa)/as.numeric(tab11$py)
tab11$rate_ratio<-c("-", poisson.test(tab11$n_psa,tab11$py)[5])
tab11$p_value<-c("-", format_pval(poisson.test(tab11$n_psa, tab11$py)[3]))
cbind(summary(glm(events~freq+freq:factor(c_dt),offset=log(py),data=inc_paz1_high, family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~freq+freq:factor(c_dt),offset=log(py),data=inc_paz1_high,family="poisson"))$coefficients[, 4], format_pval))


#Chronic disease
tab12<-data.frame(Follow_up=c("Mixed", "Mixed", "Mixed", "Mixed"), chronic_dis=c("0", "[1,3]", "[4,10]","[11+)"), n_pat=as.vector(tapply(count_dis$pat_id, count_dis$c_cron, n_distinct)),
                 n_arzt=as.vector(tapply(count_dis$arzt_id, count_dis$c_cron, n_distinct)), 
                 n_psa=as.vector(tapply(inc_pat_cron$events, inc_pat_cron$cronic, sum)), 
                 py=as.vector(tapply(inc_pat_cron$py, inc_pat_cron$cronic, sum)))
tab12$ir<-as.numeric(tab12$n_psa)/as.numeric(tab12$py)
resp6<-summary(pois.b<-glm(events~factor(cronic),offset=log(py),data=subset(inc_pat_cron, py!=0),family="poisson", na.action=na.exclude))
tab12$rate_ratio<-c("-", format(exp(-as.vector(resp6$coefficients[,1][2:4])), digits=4))
tab12$p_value<-c("-", as.character(lapply(as.vector(resp6$coefficients[,4][2:4]), format_pval)))
cbind(summary(glm(events~cronic+cronic:factor(c_dt),offset=log(py),data=subset(inc_pat_cron, py!=0), family="poisson"))$coefficients[, 1:3], p_value=lapply(summary(glm(events~cronic+cronic:factor(c_dt),offset=log(py),data=subset(inc_pat_cron, py!=0),family="poisson"))$coefficients[, 4], format_pval))

#Chronic disease
#number of patient by categories
#inc_pat_cron$cases[32]+inc_pat_cron$cases[64]+inc_pat_cron$cases[96]+inc_pat_cron$cases[128]

#########
#Not reported

#autoregressive
#not autoregressive
ts.plot(ts(inc$events.2/inc$py.2, frequency = 4, start = c(2009,1)))
fit<-arima(ts(inc$events.2/inc$py.2, frequency = 4, start = c(2009,1)),order=c(1,0,0))
tsdiag(fit)
pred<-predict(fit, n.ahead=10)
plot(ts(inc$events.2/inc$py.2, frequency = 4, start = c(2009,1)), xlim=c(2009,2020))
lines(pred$pred,col="blue", lwd=5)
#linear model
lm(inc$events.2[2:32]/inc$py.2[2:32]~as.Date(inc$dt[2:32]))
plot(y=inc$events.2[2:32]/inc$py.2[2:32],x=as.Date(inc$dt[2:32]))
lines(y=fitted(lm(inc$events.2[2:32]/inc$py.2[2:32]~as.Date(inc$dt[2:32]))), x=as.Date(inc$dt[2:32]))
summary(lm(inc$events.2[2:32]/inc$py.2[2:32]~as.Date(inc$dt[2:32])))

#Call:
#lm(formula = inc$events.2[2:32]/inc$py.2[2:32] ~ as.Date(inc$dt[2:32]))

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.033152 -0.007598 -0.002525  0.006373  0.046112 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            3.598e-01  5.616e-02   6.406 5.26e-07 ***
#  as.Date(inc$dt[2:32]) -1.856e-05  3.571e-06  -5.198 1.46e-05 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.01624 on 29 degrees of freedom
#Multiple R-squared:  0.4823,	Adjusted R-squared:  0.4645 
#F-statistic: 27.02 on 1 and 29 DF,  p-value: 1.462e-05

#change of 2 per mille over 10 years - not relevant
