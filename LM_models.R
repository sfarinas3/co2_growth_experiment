# install.packages(c('lme4','lmer','lmerTest','emmeans','Matrix','AICcmodavg','dplyr',
#                    'performance','broom.mixed','pbkrtest','multcomp','ggeffects','MuMIn','data.tables'))

setwd("C:/Users/sfari/OneDrive/Documents/GitHub/co2_growth_experiment")
library(lmerTest)
library(AICcmodavg)
library(broom.mixed)
library(performance)
library(emmeans)
library(multcomp)
library(ggeffects)
library(ggplot2)
library(MuMIn)
library(data.table)
library(dplyr)

## NOTE: AIC models must be fit using ML method


#Import data for herbs
ht  = read.csv('herb_transplants_model_rdy.csv')


#Run model for each species
species = c('MIGU','NAOF','POMO','VEAN')
ht[,'treatment'] <- as.factor(ht[,'treatment'])

#Lists to catch all the model results dfs
df_ls_ht = list()
aic_ls_ht = list()
ph_ls_ht = list()

for (sp in species) {
  
print('')
print('********************')
print('Running:')
print(sp)
df = subset(ht, species == sp)

#Set up null model with only an intercept
null_model = lm('total_biomass ~ 1', data = df)

#Set up model with treatment as factor
model = lm('total_biomass ~ treatment', data = df)

#Print model results
aic_r = AIC(null_model, model)
print(aic_r)
print(summary(model))
print(anova(model, null_model))
print(MuMIn::r.squaredGLMM(model))

#Capture AIC results
aic = as.data.frame(aic_r)
aic = cbind(model_type = rownames(aic), aic)
aic$species = sp
aic_ls_ht = c(aic_ls_ht, list(aic))

#Do post-hoc test
ph = summary(glht(model, linfct = mcp(treatment='Tukey')))
print(ph)
g = ggpredict(model, terms = 'treatment')

#Capture posthoc test
phdf = as.data.frame(ph$test$coefficients)
phdf = cbind(coefficients = ph$test$coefficients, phdf)
phdf = cbind(comparison = rownames(phdf), phdf)
phdf = cbind(tstat = ph$test$tstat, phdf)
phdf = cbind(stderr = ph$test$sigma, phdf)
phdf = cbind(pvalues = ph$test$pvalues, phdf)
phdf$species = sp
ph_ls_ht = c(ph_ls_ht, list(phdf))

#Save figure of predicted biomass by treatment
g = ggpredict(model, terms = 'treatment')
plot(g) + labs(title = paste0(sp,' Predicted Biomass By Treatment'))
ggsave(paste0('figures/',sp,'_predicted_bmass.jpg'), width=12, height=10)

#Create dataframe from the results and append to list
cfs = as.data.frame(summary(model)$coefficients)
cfs = cbind(parameter = rownames(cfs), cfs)
cfs$species = sp
df_ls_ht = c(df_ls_ht, list(cfs))
}


#----------------------------


#Import data for woody
wt  = read.csv('woody_transplants_model_rdy.csv')

#Run model for each species
species = c('CEOC','POFR','SAGO','TACH')

#List to catch all the model results dfs
df_ls_wt = list()
aic_ls_wt = list()
ph_ls_wt = list()

for (sp in species) {

print('')
print('********************')
print('Running:')
print(sp)
df = subset(wt, species == sp)

#Set up null model with random effects only
null_model = lmer('total_biomass ~ (1|chamber) + (1|chamber:shelf)', data = df, REML=FALSE)

#Set up model with no random effects and treatment as factor
model_no_randomeff = lm('total_biomass ~ treatment', data = df)

#Set up model with treatment as factor and random effects
model = lmer('total_biomass ~ treatment + (1|chamber) + (1|chamber:shelf)', data = df, REML=FALSE)

#Print model results
aic_r = AIC(null_model, model_no_randomeff, model)
print(aic_r)
print(tidy(model))
print(anova(model, null_model))
print(anova(model, model_no_randomeff))
print(MuMIn::r.squaredGLMM(model))

#Capture AIC results
aic = as.data.frame(aic_r)
aic = cbind(model_type = rownames(aic), aic)
aic$species = sp
aic_ls_wt = c(aic_ls_wt, list(aic))

#Do post-hoc test
ph = summary(glht(model, linfct = mcp(treatment='Tukey')))
print(ph)

#Capture posthoc test
phdf = as.data.frame(ph$test$coefficients)
phdf = cbind(coefficients = ph$test$coefficients, phdf)
phdf = cbind(comparison = rownames(phdf), phdf)
phdf = cbind(tstat = ph$test$tstat, phdf)
phdf = cbind(stderr = ph$test$sigma, phdf)
phdf = cbind(pvalues = ph$test$pvalues, phdf)
phdf$species = sp
ph_ls_wt = c(ph_ls_wt, list(phdf))

#Save figure of predicted biomass by treatment
g = ggpredict(model, terms = 'treatment')
plot(g) + labs(title = paste0(sp,' Predicted Biomass By Treatment'))
ggsave(paste0('figures/',sp,'_predicted_bmass.jpg'), width=12, height=10)

#Create dataframe from the model results and append to list
cfs = as.data.frame(coef(summary(model)))
cfs = cbind(parameter = rownames(cfs), cfs)
cfs$species = sp
df_ls_wt = c(df_ls_wt, list(cfs))
}

#Concatenate all the model result dfs in the list together into one df
fin_ls = c(df_ls_ht, df_ls_wt)
fin_df = rbindlist(fin_ls, fill=TRUE)
write.csv(fin_df, 'results/model_coeff_results.csv')

fin_ph_ls = c(ph_ls_ht, ph_ls_wt)
fin_ph_df = rbindlist(fin_ph_ls, fill=TRUE)
write.csv(fin_ph_df, 'results/posthoc_test_results.csv')

fin_aic_ls = c(aic_ls_ht, aic_ls_wt)
fin_aic_df = rbindlist(fin_aic_ls, fill=TRUE)
write.csv(fin_aic_df, 'results/aic_results.csv')