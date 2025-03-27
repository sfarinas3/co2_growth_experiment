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
ht['total_biomass_log'] = log10(ht['total_biomass'])

#Run model for each species
species = c('MIGU','NAOF','POMO','VEAN')
ht[,'treatment'] <- as.factor(ht[,'treatment'])

#Lists to catch all the model results dfs
df_ls_ht = list()
aic_ls_ht = list()
ph_ls_ht = list()
sw_ls_ht = list()

for (sp in species) {
  
  print('')
  print('********************')
  print('Running:')
  print(sp)
  df = subset(ht, species == sp)
  
  #Set up null model with only an intercept
  null_model = lm('total_biomass_log ~ 1', data = df)
  
  #Set up model with treatment as factor
  model = lm('total_biomass_log ~ treatment', data = df)
  #plot(model)
  
  #Print model results
  aic = as.data.frame(compare_performance(null_model, model))
  print(aic)
  print(summary(model))
  
  #Capture AIC results
  aic$species = sp
  aic_ls_ht = rbind(aic_ls_ht, aic)
  
  #Shapiro-Wilk test for testing if residuals are normally distributed
  sw_nm = shapiro.test(model[['residuals']])
  sw = shapiro.test(model[['residuals']])
  sw_df = data.frame(
    model = c('null_model','model'),
    w_statistic = c(sw_nm$statistic, sw$statistic),
    p_value = c(sw_nm$p.value, sw$p.value)
  )
  sw_df$species = sp
  sw_ls_ht = c(sw_ls_ht, list(sw_df))
  
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
  
  #Plot test on residuals and save
  png(filename=sprintf('./figures/%s_residuals_check.jpg',sp), width=900, height=800)
  cm = check_model(model)
  plot(cm)
  dev.off()
}


#----------------------------


#Import data for woody
wt  = read.csv('./data/woody_transplants_model_rdy.csv')
wt['total_biomass_log'] = log10(wt['total_biomass'])

#Run model for each species
species = c('CEOC','POFR','SAGO','TACH')

#List to catch all the model results dfs
df_ls_wt = list()
aic_ls_wt = list()
ph_ls_wt = list()
sw_ls_wt = list()

for (sp in species) {
  
  print('')
  print('********************')
  print('Running:')
  print(sp)
  df = subset(wt, species == sp)
  
  #Set up null model with random effects only
  null_model = lmer('total_biomass_log ~ (1|chamber) + (1|chamber:shelf)', data = df, REML=FALSE)
  
  #Set up model with no random effects and treatment as factor
  model_no_randomeff = lm('total_biomass_log ~ treatment', data = df)
  
  #Set up model with treatment as factor and random effects
  model = lmer('total_biomass_log ~ treatment + (1|chamber) + (1|chamber:shelf)', data = df, REML=FALSE)
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
  
  #Shapiro-Wilk test for testing if residuals are normally distributed
  sw_nm = shapiro.test(residuals(null_model))
  sw_nre = shapiro.test(residuals(model_no_randomeff))
  sw = shapiro.test(residuals(model))
  sw_df = data.frame(
    model = c('null_model','model_no_randomeff','model'),
    w_statistic = c(sw_nm$statistic, sw_nre$statistic, sw$statistic),
    p_value = c(sw_nm$p.value, sw_nre$p.value, sw$p.value)
  )
  sw_df$species = sp
  sw_ls_wt = c(sw_ls_wt, list(sw_df))
  
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
  
  #Plot test on residuals and save
  png(filename=sprintf('./figures/%s_residuals_check.jpg',sp), width=900, height=1200)
  cm = check_model(model)
  plot(cm)
  dev.off()
}

#Concatenate all the model result dfs in the list together into one df
fin_ls = c(df_ls_ht, df_ls_wt)
fin_df = rbindlist(fin_ls, fill=TRUE)
write.csv(fin_df, './results/model_coeff_results.csv')

fin_ph_ls = c(ph_ls_ht, ph_ls_wt)
fin_ph_df = rbindlist(fin_ph_ls, fill=TRUE)
write.csv(fin_ph_df, './results/posthoc_test_results.csv')

sw_ls = c(sw_ls_ht, sw_ls_wt)
fin_sw_df = rbindlist(sw_ls, fill=TRUE)
write.csv(fin_sw_df, './results/shapiro_wilks_test_results.csv')

write.csv(aic_ls_ht, './results/aic_herb_results.csv')
write.csv(aic_ls_wt, './results/aic_woody_results.csv')