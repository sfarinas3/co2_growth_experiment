install.packages(c('lme4','lmer','lmerTest','emmeans','Matrix'))
setwd("~/Documents/GitHub/co2_growth")
library(lmerTest)

#Import data for herbs
ht  = read.csv('herb_transplants_model_rdy.csv')

naof = subset(ht, species == 'NAOF')
model = lm('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl', data = naof)
summary(model)
 
migu = subset(ht, species == 'MIGU')
model = lm('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl', data = migu)
summary(model)

naof = subset(ht, species == 'NAOF')
model = lm('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl', data = naof)
summary(model)

pomo = subset(ht, species == 'POMO')
model = lm('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl', data = pomo)
summary(model)


#Import data for woody
wt  = read.csv('woody_transplants_model_rdy.csv')

ceoc = subset(wt, species == 'CEOC')
model = lmer('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl + (1|chamber) + (1|chamber:shelf)', data = ceoc)
summary(model)

pofr = subset(wt, species == 'POFR')
model = lmer('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl + (1|chamber) + (1|chamber:shelf)', data = pofr)
summary(model)

sago = subset(wt, species == 'SAGO')
model = lmer('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl + (1|chamber) + (1|chamber:shelf)', data = sago)
summary(model)

tach = subset(wt, species == 'TACH')
model = lmer('total_biomass ~ co2_lvl + water_lvl + co2_lvl * water_lvl + (1|chamber) + (1|chamber:shelf)', data = tach)
summary(model)
