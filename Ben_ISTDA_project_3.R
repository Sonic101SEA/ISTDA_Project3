metabric <- read.csv("metabric-analytical.csv")

data(metabric)
summary(metabric)

metabric$surv_status <- factor(metabric$surv_status, levels=c("deceased", "living"), c("1","0"))
metabric$er_status <- factor(metabric$er_status, levels=c("positive", "negative"), c("1", "0"))
metabric$her2_status <- factor(metabric$her2_status, levels=c("positive", "negative"), c("1", "0"))
metabric$pr_status <- factor(metabric$pr_status, levels=c("positive", "negative"), labels=c("1", "0"))

km_fit_1 <- survfit(Surv(surv_months, surv_status) ~ breast_surgery, data = metabric)
plot(km_fit_1)

plot(survfit(Surv(surv_months, surv_status) ~ er_status, data = metabric))
plot(survfit(Surv(surv_months, surv_status) ~ her2_status, data = metabric))
plot(survfit(Surv(surv_months, surv_status_1) ~ pr_status, data = metabric), col = 2:3, main = "Survival curve for Progesterone receptor postive individuals")
legend("topleft", c("Progesterone receptor positive", "Progesterone recptor negative"), col = 2:3, lty=1)

metabric$age_yr <- metabric$age_diagnosis + (metabric$surv_months/12)

plot(survfit(Surv(age_yr, surv_status) ~ er_status, data = metabric), col = 2:3, main = "Survival curve for Estrogene receptor postive individuals from birth")
legend("topleft", c("ER positive", "ER negative"), col = 2:3, lty=1)

plot(survfit(Surv(surv_months, surv_status_1) ~ er_status, data = metabric), col = 2:3, fun = "cloglog", main = "Complentary log-log ER")
legend("topleft", c("ER positive", "ER negative"), col = 2:3, lty=1)
plot(survfit(Surv(age_yr, surv_status) ~ er_status, data = metabric), col = 2:3, fun= "cloglog", main = "Complentary log-log ER from birth")
legend("topleft", c("ER positive", "ER negative"), col = 2:3, lty=1)

