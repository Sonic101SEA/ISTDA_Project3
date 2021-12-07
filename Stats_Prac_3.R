path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
metabric <- readRDS(here::here("data", "metabric-analytical.rds"))
head(metabric)
View(metabric)


## recode vital_status
metabric <- metabric %>% 
  mutate(death = ifelse(vital_status == 'died of disease', 1, 0))

View(metabric)

##############################################
## Fitting base model - Aim 1

km_fit_0 <- survfit(Surv(surv_months, death) ~ 1, data = metabric)
# plot(km_fit_0)

## survival prob at 60 months
summary(km_fit_0, times = c(60))
## get median
km_fit_0

## plot KM
km_fit_0_plot <- 
  ggsurvplot(km_fit_0, 
           title = "Kaplan-Meier plot for overall survivability", 
           xlab = "Time in months",
           ylab = "Overall survival probability",
           size = 0.3, 
           censor.size = 0.1,
           surv.median.line = "hv",
           )
km_fit_0_plot$plot + geom_vline(xintercept = 60)

## plot KM by molecular status

km_fit_her <- survfit(Surv(surv_months, death) ~ her2_status, data = metabric)
km_fit_er <- survfit(Surv(surv_months, death) ~ er_status, data = metabric)
km_fit_pr <- survfit(Surv(surv_months, death) ~ pr_status, data = metabric)

### TO DO: Convert to ggsurvplots with title and 
## different colours for each plot and line:
plot(km_fit_her)
plot(km_fit_er)
plot(km_fit_pr)



############################################
##############  ------- Aim 3 -------------

###### ignore this for now
# metabric %>%filter(er_status == 'positive') %>%  ggplot(aes(npi)) + geom_density()
# metabric %>%filter(her2_status == 'positive') %>%  ggplot(aes(npi)) + geom_density()
# metabric %>%filter(pr_status == 'positive') %>%  ggplot(aes(npi)) + geom_density()
# 
# metabric %>%filter(er_status == 'positive') %>%  ggplot(aes(tumor_size)) + geom_density()
# metabric %>%filter(her2_status == 'positive') %>%  ggplot(aes(tumor_size)) + geom_density()
# metabric %>%filter(pr_status == 'positive') %>%  ggplot(aes(tumor_size)) + geom_density()
# 
# table(metabric$er_status, metabric$her2_status)


######
#### Recode exposure variables into one single representative variable:

metabric <- metabric %>% mutate(er_status_binary = ifelse(er_status == "positive",1,0))
metabric <- metabric %>% mutate(her2_status_binary = ifelse(her2_status == "positive",1,0))
metabric <- metabric %>% mutate(pr_status_binary = ifelse(pr_status == "positive",1,0))
View(metabric)

metabric <- metabric %>% mutate(composite_status = paste0(er_status_binary,
                                                          her2_status_binary,
                                                          pr_status_binary))
table(metabric$composite_status)

## recode the small group 011 into the same group as 001 - 
## we will just refer to this group as the PR group

metabric <- metabric %>% mutate(composite_status = ifelse(composite_status == '011',
                                                          '001',
                                                          composite_status ))
table(metabric$composite_status)


### --------------- Fit a COX regression

## simple cox regression
fit_cox_base <- coxph(Surv(surv_months, death) ~ composite_status, data = metabric)
summary(fit_cox_base)

## diagnostics
surv_comp <- survfit(Surv(surv_months, death) ~ composite_status, data = metabric)
plot(surv_comp, col = c(1:7), main = "Kaplan Meier")
plot(surv_comp, fun = 'cloglog', col = c(1:7), main = 'cloglog')
legend("topleft", c("stand", "test"), col = c(4,6), lty =  1 ) 


cox.zph(fit_cox_base)
plot(cox.zph(fit_cox_base))
ggcoxzph(cox.zph(fit_cox_base))

# Need to justify why not using treatment -- correlated to exposure 
# ## multiple cox regression - then redo diagnostics

fit_cox_full <- coxph(Surv(surv_months, death) ~ composite_status + 
                        age_diagnosis + 
                        cellularity + 
                        npi, data = metabric)
summary(fit_cox_full)


cox.zph(fit_cox_full)
plot(cox.zph(fit_cox_full))
ggcoxzph(cox.zph(fit_cox_full))

## stratify by <60, -120 and 120+? 