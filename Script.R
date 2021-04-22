# Analyses for Final Stats Paper

# Import libraries
library(tidyverse)
library(lme4)
library(ggeffects)
library(lavaan) 
library(semPlot) 
library(simr)

# Import relevant data
setwd("/Users/eheffernan/OneDrive - University of Toronto/PSY2002/FinalPaper/Dems")

# Load cleaned data

rep_summary <- read_csv('ProlificStudySummary.csv') %>%
  select(-X1)

rep_trials <- read_csv('CleanTrialData.csv') %>%
  select(-X1)

rep_trials <- left_join(rep_trials, rep_summary %>% select(Participant, Age, Sex), by = 'Participant')

orig_dem <- read_csv('OriginalStudyDems.csv') %>%
  select(-X1) 

orig_trials <- read_csv('PSY100Trials.csv') %>%
  rename(Sex = Gender)


# Welch's t-test to confirm that ages are different
stats::t.test(orig_dem$Age, (rep_summary %>% filter(Type == 'Exception'))$Age)

# Original GLM that assesses impact of manipulation on categorization accuracy
learn_glm <- function(df)
{
  model.fit <- glmer(Learn_Corr~Count*Type*Condition+(1|Participant), 
                     data=filter(df,Count<=36 & Block == "Learn"),
                     family=binomial,
                     control=glmerControl(optimizer='optimx',calc.derivs=FALSE,
                                          optCtrl=list(method='bobyqa',starttests=FALSE,kkt=FALSE)))  
}

# Results for original and replication studies
learn.fit.orig <- learn_glm(orig_trials)
summary(learn.fit.orig)

learn.fit.rep <- learn_glm(rep_trials)
summary(learn.fit.rep)

# Plot results for visualization
# Original
learn_plot_orig <- ggpredict(learn.fit.orig, terms = c("Count", "Condition", "Type"), condition = "Late", se = TRUE)
learn_plot_orig <- plot(learn_plot_orig) + theme_classic() + 
  ggtitle('Categorization Performance in Heffernan and Mack (2021)') +
  theme(text = element_text(size=15)) + 
  ylab("Likelihood of Correct Response") + xlab("Repetition") + 
  scale_color_manual(values=c("#6096FD", "#FB7B8E")) +
  scale_fill_manual(values=c("#6096FD", "#FB7B8E")) +
  theme(strip.background = element_blank()) +
  ggsave("learn_acc_orig.png", units="in", width=8, height=4, dpi=300)
plot(learn_plot_orig)

# Replication
learn_plot_rep <- ggpredict(learn.fit.rep, terms = c("Count", "Condition", "Type"), condition = "Late", se = TRUE)
learn_plot_rep <- plot(learn_plot_rep) + theme_classic() + 
  ggtitle('Categorization Performance in Replication Study') +
  theme(text = element_text(size=15)) + 
  ylab("Likelihood of Correct Response") + xlab("Repetition") + 
  scale_color_manual(values=c("#6096FD", "#FB7B8E")) +
  scale_fill_manual(values=c("#6096FD", "#FB7B8E")) +
  theme(strip.background = element_blank()) +
  ggsave("learn_acc_rep.png", units="in", width=8, height=4, dpi=300)

# Standardize data to get model convergence
# Z-score age
rep_trials <- rep_trials %>%
  filter(Block == 'Learn') %>%
  mutate(Age_z = (Age - mean(Age))/sd(Age)) 

# Create variable with 'buckets' for repetition (0 = initial, 1 = mid, 3 = end)
rep_trials <- rep_trials %>%
  group_by(Participant) %>%
    mutate(Count_bucket = case_when( 
      Count <= 12 ~ 0,
      Count <= 24 ~ 1,
      Count <= 36 ~ 2,
      Count <= 48 ~ 3,
      Count <= 60 ~ 4,
      TRUE ~ 5)) %>%
    ungroup() 

# Add age to learn glm
learn.fit.age <- glmer(Learn_Corr~Count*Type*Condition*Age_z+(1|Participant), 
                   data=filter(rep_trials,Count<=36 & Block == "Learn"),
                   family=binomial,
                   control=glmerControl(optimizer='optimx',calc.derivs=FALSE,
                                        optCtrl=list(method='bobyqa',starttests=FALSE,kkt=FALSE)))  
summary(learn.fit.age)


# Look at impact of condition on RT

rt_fit_rep <- glmer(Learn_RT ~ Condition * Count_bucket + (1|Participant), 
                data = rep_trials %>% filter(Block == "Learn" & Type == 'Exception'),
                family = inverse.gaussian(link = "identity"),
                glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 4e5))) #similar results w elmer model log(rt)

summary(rt_fit_rep)

# Calculate effect sizes 
fixef(rt_fit_rep)
effect_count <- powerSim(rt_fit_rep, fixed("Count_bucket"), nsim=50)
effect_cond <- powerSim(rt_fit_rep, fixed("ConditionLate"), nsim=50)
effect_condInterCount <- powerSim(rt_fit_rep, fixed("ConditionLate:Count_bucket"), nsim=50)

# Visualize 
rt_rep_data <- ggpredict(rt_fit_rep, terms = c("Count_bucket", "Condition"), se = TRUE)
rt_rep_plot <- plot(rt_rep_data) + theme_classic()


# Include age in RT analysis
rt_fit_age <- glmer(Learn_RT ~ Condition * Count_bucket * Age_z + (1|Participant), 
                                  data = rep_trials %>% filter(Block == "Learn" & Type == 'Exception'),
                                  family = inverse.gaussian(link = "identity"),
                                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 4e5))) #similar results w elmer model log(rt)

summary(rt_fit_age)

# Visualize Results
rt_age_data <- ggpredict(rt_fit_age, terms = c("Count_bucket", "Condition", "Age_z"), se = TRUE)
rt_age_plot <- plot(rt_age_data) + theme_classic() +
  ggtitle('Reaction Time Across Age and Condition') +
  theme(text = element_text(size=15)) + 
  ylab("Reaction Time") + xlab("Repetition") + 
  scale_color_manual(values=c("#6096FD", "#FB7B8E")) +
  scale_fill_manual(values=c("#6096FD", "#FB7B8E")) +
  theme(strip.background = element_blank()) +
  ggsave("learn_acc_rep.png", units="in", width=8, height=4, dpi=300)
  

ageCountEffect <- powerSim(rt_fit_age, fixed("Count_bucket"),nsim=50)

ageInterCondEffect <- powerSim(rt_fit_age, fixed("ConditionLate:Age_z"), nsim=50)

AgeInterCondInterCountEffect <- powerSim(rt_fit_age, fixed("ConditionLate:Count_bucket:Age_z"),nsim=50)

# Look at accuracy faceted by sex (supplementary analysis; not included in final paper)
learn_glm_sex <- function(df)
{
  model.fit <- glmer(Learn_Corr~Count*Condition*Sex*Type +(1|Participant), 
                     data=filter(df,Count<=36 & Block == "Learn"),
                     family=binomial,
                     control=glmerControl(optimizer='optimx',calc.derivs=FALSE,
                                          optCtrl=list(method='bobyqa',starttests=FALSE,kkt=FALSE)))  
}

fit.sex.orig <- learn_glm_sex(orig_trials)
summary(fit.sex.orig)
sex_plot_orig <- ggpredict(fit.sex.orig, terms = c("Count", "Condition", "Type", "Sex"), se = TRUE)
sex_plot_orig <- plot(sex_plot_orig) + theme_classic()
plot(sex_plot)

fit.sex.rep <- learn_glm_sex(rep_trials)
summary(fit.sex.rep)
sex_plot <- ggpredict(fit.sex.rep, terms = c("Count", "Condition", "Type", "Sex"), se = TRUE)
sex_plot <- plot(sex_plot) + theme_classic()
plot(sex_plot)

# PATH ANALYSIS----
basic_model = '
  Learn_RT ~ a*Condition 
  Learn_Corr ~ b*Condition + c*Learn_RT
  
  Condition.through.Learn_RT:=a*c
'
basic_model_rep <- basic_model %>%
  sem(data = rep_trials %>% filter(Type == 'Exception' & Block == 'Learn'))

# Print results
basic_model_rep %>%
  summary(fit.measures = TRUE, 
          standardized=TRUE
  )

age_model = '
  Learn_RT ~ a*Condition + d*Age_z
  Learn_Corr ~ b*Condition + e*Age_z + c*Learn_RT
  
  Condition.through.Learn_RT:=a*c
  Age.through.Learn_RT:=d*c
  
  Age_z ~~ 0*Condition
'
 
 
 age_model_rep <- age_model %>%
   sem(data = rep_trials %>% filter(Type == 'Exception' & Block == "Learn"))
 
 # Print results
 age_model_rep %>%
   summary(fit.measures = TRUE, 
           standardized=TRUE
   )

# Plot path models
locations_base = matrix(c(0,.5, 
                          0,-.5, 
                          -.5,0), 
                   ncol=2, 
                   byrow=2
)

labels_base = c("Reaction\nTime",
           "Accuracy",
           "Condition"
)

p1 <- basic_model_rep %>%
  semPaths(whatLabels="est", 
           nodeLabels = labels_base, 
           layout = locations_base, 
           edge.label.cex = 1.25,
           edge.label.position=.5,
           label.cex = 1.25,
           sizeMan = 12
  )

# Plot path model w age
locations_age = matrix(c(.5,.5, 
                          .5,-.5, 
                          -.5,.5,
                          -.5, -.5
                          ), 
                        ncol=2, 
                        byrow=2
)

labels_age = c("Reaction\nTime",
                "Accuracy",
                "Condition",
                "Age_z"
)

p2 <- age_model_rep %>%
  semPaths(whatLabels="est", 
           nodeLabels = labels_age, 
           layout = locations_age, 
           edge.label.cex = 1.25,
           edge.label.position=.4,
           label.cex = 1.25,
           sizeMan = 12
  )


