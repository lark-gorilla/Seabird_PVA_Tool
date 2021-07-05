# Round 4 PVA

# load libraries and functions

rm(list=objects())

# install these packages if you don't already have
library(popbio)
library(readxl)
library(writexl)

setwd("C:/niras_rproj/Seabird_PVA_Tool/R/Rpackage/") # change to where you saved

#loads functions to run PCA
ff <- list.files(".", pattern="functions") ; for(k in 1:length(ff)){ source(ff[k])}
modeoptions <- read.csv("ModeOptions.csv") ## Added Version 2.1


## ##################################################

lookup_dat<-read_xlsx('C:/R4/RIAA/PVA/PVA_lookup.xlsx')


# choose species and PVA_site to run
table(lookup_dat$Species)
table(lookup_dat$PVA_site)

# example with Kittiwakes at Flamborough
sp_site<-lookup_dat[lookup_dat$Species=='Kittiwake' & lookup_dat$PVA_site=='Flamborough and Filey Coast',]
# have a look
sp_site

# need to convert absolute impact (n bird/yr mortality) to relative impact (% change in pop)
# here we take mortalities apportioned to the SPA and / by total SPA meta population -
# would be better to use total SPA population that went into apportioning analyses 

# CRM
mn_rel_imp_CRM<-(sp_site$mn_ann_mort_CRM[1]/2)/sum(sp_site$Count_bp) # what units are absolute mortality, assumed n birds so divide by 2 to get bp?
sd_rel_imp_CRM<-(sp_site$sd_ann_mort_CRM[1]/2)/sum(sp_site$Count_bp) 
# displacement
mn_rel_imp_disp<-(sp_site$mn_ann_mort_disp[1]/2)/sum(sp_site$Count_bp)
sd_rel_imp_disp<-(sp_site$sd_ann_mort_disp[1]/2)/sum(sp_site$Count_bp) 
# total sum() CRM + displacement
mn_rel_imp_both<-((sp_site$mn_ann_mort_CRM[1]+sp_site$mn_ann_mort_disp[1])/2)/sum(sp_site$Count_bp)
sd_rel_imp_both<-((sp_site$sd_ann_mort_CRM[1]+sp_site$sd_ann_mort_disp[1])/2)/sum(sp_site$Count_bp) 

# select what impact(s) running for, here we choose collision risk
my_impact_mean<-mn_rel_imp_CRM
my_impact_sd<-sd_rel_imp_CRM
## ##################################################################
# # code to run PVA with shiny app specifications
## ##################################################################



run1 <- nepva.simplescenarios(model.envstoch = "betagamma", # survival and productivity (model.prodmax=T) env stochasticity from beta distribution (yes - NIRAS) 
                                  model.demostoch = TRUE, # model demographic stochasticity (yes - NIRAS)
                                  model.dd = "nodd", # ' nodd' = no density dependence (density independence - NIRAS)
                                  model.prodmax = TRUE, # productivity rates constrained to be <= maximum brood size
                                  mbs =  sp_site$mbs[1], # maximum brood size
                                  afb = sp_site$afb[1], # age at first breeding
                                  npop = nrow(sp_site), # number of populations (useful if >1 col per SPA) 
                                  nscen = 1, # number of impact scenarios
                                  sim.n = 5000, sim.seed = 1898, nburn = 0, # n simulations, seed number and n burn - NIRAS spec
                                  demobase.specify.as.params = FALSE, # enter empirical values for prod and surv rather than estimate
                                  demobase.splitpops = TRUE, # different demographic rates for each subpopulation, ignored if npop=1
                                  demobase.splitimmat = FALSE, # different demographic rates specified for immatures? (no - NIRAS)
                                  demobase.prod = data.frame(Mean=sp_site$mn_base_prod, SD = sp_site$sd_base_prod), # baseline productivity mean(s) and SD(s)
                                  demobase.survadult = data.frame(Mean =  sp_site$mn_base_adsurv, SD = sp_site$sd_base_adsurv), # baseline survival mean(s) and SD(s)
                                  inipop.years = as.numeric(sp_site$yr), # year(s) when inital count was made
                                  inipop.inputformat = "breeding.pairs", # initial population size entered as breeding pair (SMP data) 
                                  inipop.vals = sp_site$Count_bp, # initial population value(s)
                                  impacts.relative = TRUE, # specify relative impacts (% change in population from impact)
                                  impacts.splitimmat = FALSE, # different impacts specified for immatures? (no - NIRAS)
                                  impacts.year.start = max(as.numeric(sp_site$yr))+1, # when to start impact (inipop.years+1 - NIRAS)
                                  impacts.year.end = max(as.numeric(sp_site$yr))+61,# when does impact end? (impacts.year.start+60 - NIRAS) 
                                  impacts.splitpops = FALSE,  # different impact rates for each subpopulation, ignored if npop=1
                                  impacts.scennames = c("imp1"), # naming impact for reference
                                  impacts.prod.mean = 0, # impacts on productivity, val per population (use 0, no prod impacts - NIRAS)
                                  impacts.survadult.mean = my_impact_mean,# impacts on adult survival # choose CRM/disp/CRM+disp
                                  impacts.provideses = TRUE, # whether to provide SEs for impacts (yes, I think -NIRAS)
                                  impacts.prod.se = 0, # se for productivity impact, 0 or NA?
                                  impacts.survadult.se = my_impact_sd, # se for survival impact # choose CRM/disp/CRM+disp
                                  output.agetype = "breeding.adults", #  output population as breeding adults
                                  output.year.end = max(as.numeric(sp_site$yr))+61, #when to end model (impacts.year.start+60 - NIRAS)
                                  output.year.start = max(as.numeric(sp_site$yr))+1, # when to start model report from, (inipop.years?)
                                  output.popsize.target = NULL, # don't set output pop target
                                  output.popsize.qe = NULL,# don't set output extinction risk target
                                  silent = FALSE, # report progress, errors
                                  changetablenames = TRUE) # match shiny table names


# view output run1
run1
# write it out
write_xlsx(run1, 'myloc/pva_run.xlsx')

# to interpret see pdfs in 'Seabird_PVA_Tool/Documentation'

## ##################################################################
# # create master lookup sheet - OLD
## ##################################################################

popz<-read_xlsx('C:/R4/RIAA/apportioning/SMPdata_last_count_yr_max.xlsx')

p1<-filter(popz, Master.site %in% c('Farne Islands SPA', 'Forth Islands SPA' ) | Site %in% c('Flamborough Head and Bempton Cliffs', 
                                                                                             'Filey 1', 'Filey 2', 'Filey 3', 'Orford Ness 1', 'Hollesley Marsh', 'Havergate Island', 'Butley River RSPB'))                                                                                  ))
write_xlsx(p1, 'C:/R4/RIAA/PVA/PVA_lookup.xlsx')


example7 <- nepva.calcdefaults(Species = "Herring Gull", 
                               poolregtype.BS = "Global", poolregion.BS = "Global", 
                               sourcepop.surv = "National", lookup.dir = 'C:/niras_rproj/Seabird_PVA_Tool/R/Rpackage/')

# backup code

example1 <- nepva.simplescenarios(model.envstoch = "betagamma", # survival and productivity (model.prodmax=T) env stochasticity from beta distribution (yes - NIRAS) 
                                  model.demostoch = TRUE, # model demographic stochasticity (yes - NIRAS)
                                  model.dd = "nodd", # ' nodd' = no density dependence (density independence - NIRAS)
                                  model.prodmax = TRUE, # productivity rates constrained to be <= maximum brood size
                                  mbs = 4, # maximum brood size
                                  afb = 5, # age at first breeding
                                  npop = 1, # number of populations (useful if >1 col per SPA) 
                                  nscen = 1, # number of impact scenarios
                                  sim.n = 5000, sim.seed = 1898, nburn = 0, # n simulations, seed number and n burn - NIRAS spec
                                  demobase.specify.as.params = FALSE, # enter empirical values for prod and surv rather than estimate
                                  demobase.splitpops = FALSE, # different demographic rates for each subpopulation, ignored if npop=1
                                  demobase.splitimmat = FALSE, # different demographic rates specified for immatures? (no - NIRAS)
                                  demobase.prod = data.frame(Mean=0.5, SD=0.05), # baseline productivity mean(s) and SD(s)
                                  demobase.survadult = data.frame(Mean =  0.88, SD = 0.03), # baseline survival mean(s) and SD(s)
                                  inipop.years = c(2012), # year(s) when inital count was made
                                  inipop.inputformat = "breeding.pairs", # initial population size entered as breeding adults 
                                  inipop.vals = c(791), # initial population value(s)
                                  impacts.relative = TRUE, # specify relative impacts (% change in population from impact)
                                  impacts.splitimmat = FALSE, # different impacts specified for immatures? (no - NIRAS)
                                  impacts.year.start = 2013, # when to start impact (inipop.years+1 - NIRAS)
                                  impacts.year.end = 2020,# when does impact end? (impacts.year.start+60 - NIRAS) 
                                  impacts.splitpops = FALSE,  # different impact rates for each subpopulation, ignored if npop=1
                                  impacts.scennames = c("imp1"), # naming impact for reference
                                  impacts.prod.mean = c(0), # impacts on productivity, val per population (use 0, no prod impacts - NIRAS)
                                  impacts.survadult.mean = c(+0.05),# impacts on adult survival, val per population
                                  impacts.provideses = TRUE, # whether to provide SEs for impacts (yes, I think -NIRAS)
                                  impacts.prod.se = c(0), # se for productivity impact, 0 or NA?
                                  impacts.survadult.se = c(0.01), # se for survival impact
                                  output.agetype = "breeding.adults", #  output population as breeding adults
                                  output.year.end = 2020, #when to end model (impacts.year.start+60 - NIRAS)
                                  output.year.start = 2013, # when to start model report from, (inipop.years?)
                                  output.popsize.target = NULL, # don't set output pop target
                                  output.popsize.qe = NULL,# don't set output extinction risk target
                                  silent = FALSE, # report progress, errors
                                  changetablenames = TRUE) # match shiny table names
