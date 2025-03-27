install.packages(c('lme4','lmer','lmerTest','emmeans','Matrix','AICcmodavg','dplyr',
                   'performance','broom.mixed','pbkrtest','multcomp','ggeffects','MuMIn','data.tables'))

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
library(performance)

## NOTE: AIC models must be fit using ML method


#Import data for herbs
ht  = read.csv('./data/herb_transplants_model_rdy.csv')


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
#plot(model)

#Print model results
aic = as.data.frame(compare_performance(null_model, model))
print(aic)
print(summary(model))

#Capture AIC results
aic$species = sp
aic_ls_ht = rbind(aic_ls_ht, aic)

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
ggsave(paste0('./figures/',sp,'_predicted_bmass.jpg'), width=12, height=10)

#Create dataframe from the results and append to list
cfs = as.data.frame(summary(model)$coefficients)
cfs = cbind(parameter = rownames(cfs), cfs)
cfs$species = sp
df_ls_ht = c(df_ls_ht, list(cfs))
}

 
#----------------------------


#Import data for woody
wt  = read.csv('./data/woody_transplants_model_rdy.csv')

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
plot(model)
ggsave(paste0('./figures/',sp,'_qq_plot.jpg'), width=12, height=10)


#Print model results
aic = as.data.frame(compare_performance(null_model, model_no_randomeff, model))
print(aic)
print(tidy(model))

#Capture AIC results
aic$ICC = NULL
aic$species = sp
aic_ls_wt = rbind(aic_ls_wt, aic)

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
ggsave(paste0('./figures/',sp,'_predicted_bmass.jpg'), width=12, height=10)

#Create dataframe from the model results and append to list
cfs = as.data.frame(coef(summary(model)))
cfs = cbind(parameter = rownames(cfs), cfs)
cfs$species = sp
df_ls_wt = c(df_ls_wt, list(cfs))
}

#Concatenate all the model result dfs in the list together into one df
fin_ls = c(df_ls_ht, df_ls_wt)
fin_df = rbindlist(fin_ls, fill=TRUE)
write.csv(fin_df, './results/model_coeff_results.csv')

fin_ph_ls = c(ph_ls_ht, ph_ls_wt)
fin_ph_df = rbindlist(fin_ph_ls, fill=TRUE)
write.csv(fin_ph_df, './results/posthoc_test_results.csv')

write.csv(aic_ls_ht, './results/aic_herb_results.csv')
write.csv(aic_ls_wt, './results/aic_woody_results.csv')

#AIC Comparison 
#https://stats.stackexchange.com/questions/232465/how-to-compare-models-on-the-basis-of-aic
