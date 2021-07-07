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

lookup_dat<-read_xlsx('C:/R4/RIAA/PVA/PVA_lookup_collapsed.xlsx')

# have a look
lookup_dat

PVA_run<-NULL
for(k in unique(paste(lookup_dat$Species, lookup_dat$PVA_site, sep='@')))

{  
sp_site_both<-lookup_dat[lookup_dat$Species==unlist(strsplit(k, "@"))[1] & lookup_dat$PVA_site==unlist(strsplit(k, "@"))[2],]

for(m in c('site', 'incombo')) # run for in_combo and site separately
  {
  sp_site<-sp_site_both[sp_site_both$mort_type==m,]
  # need to convert absolute impact (n bird/yr mortality) to relative impact (% change in pop)
  # here we take mortalities apportioned to the SPA and / by total SPA meta population -
  # would be better to use total SPA population that went into apportioning analyses 
  
  # CRM
  mn_rel_imp_CRM<-(sp_site$mn_ann_mort_CRM/2)/sum(sp_site$Count_bp) # what units are absolute mortality, assumed n birds so divide by 2 to get bp?
  #sd_rel_imp_CRM<-(sp_site$sd_ann_mort_CRM/2)/sum(sp_site$Count_bp) 
  # displacement
  mn_rel_imp_disp<-(sp_site$mn_ann_mort_disp/2)/sum(sp_site$Count_bp)
  #sd_rel_imp_disp<-(sp_site$sd_ann_mort_disp/2)/sum(sp_site$Count_bp) 
  # total sum() CRM + displacement
  mn_rel_imp_both<-((sp_site$mn_ann_mort_CRM+sp_site$mn_ann_mort_disp)/2)/sum(sp_site$Count_bp)
  #sd_rel_imp_both<-((sp_site$sd_ann_mort_CRM+sp_site$sd_ann_mort_disp)/2)/sum(sp_site$Count_bp) 
  
  # make impact list, not running combined (both)
  impact_mean_list<-NULL
  impact_name_list<-NULL
  if(sum(mn_rel_imp_CRM)>0){impact_mean_list<-c(impact_mean_list,mn_rel_imp_CRM)}
  if(sum(mn_rel_imp_disp)>0){impact_mean_list<-c(impact_mean_list,mn_rel_imp_disp)}
  
  if(sum(mn_rel_imp_CRM)>0){impact_name_list<-c(impact_name_list,paste('collision', 1:length(mn_rel_imp_CRM), sep='_'))}
  if(sum(mn_rel_imp_disp)>0){impact_name_list<-c(impact_name_list,paste('displacement', 1:length(mn_rel_imp_disp), sep='_'))}
  
  print(sp_site)
  print(impact_mean_list)
  print(impact_name_list)
  
  p_row<-sp_site[1,]# use first row of sp_site to lookup duplicate params in model
  
  # lookup species' immature survival rates from tool at NATIONAL level
  imm_surv<-nepva.calcdefaults(Species = p_row$Species, 
                     poolregtype.BS = "Global", poolregion.BS = "Global", 
                     sourcepop.surv = "National", lookup.dir = 'C:/niras_rproj/Seabird_PVA_Tool/R/Rpackage/')$demobase.survimmat

## ##################################################################
# # code to run PVA with shiny app specifications
## ##################################################################

run1 <- nepva.simplescenarios(model.envstoch = "betagamma", # survival and productivity (model.prodmax=T) env stochasticity from beta distribution (yes - NIRAS) 
                                  model.demostoch = TRUE, # model demographic stochasticity (yes - NIRAS)
                                  model.dd = "nodd", # ' nodd' = no density dependence (density independence - NIRAS)
                                  model.prodmax = TRUE, # productivity rates constrained to be <= maximum brood size
                                  mbs =  p_row$mbs, # maximum brood size
                                  afb = p_row$afb, # age at first breeding
                                  npop = 1, # number of populations (useful if >1 col per SPA) 
                                  nscen = length(impact_name_list), # number of impact scenarios, could do more if upp/low and if site then in-combo?
                                  sim.n = 5000, sim.seed = 1898, nburn = 0, # n simulations, seed number and n burn - NIRAS spec
                                  demobase.specify.as.params = FALSE, # enter empirical values for prod and surv rather than estimate
                                  demobase.splitpops = FALSE, # different demographic rates for each subpopulation, ignored if npop=1
                                  demobase.splitimmat = TRUE, # different demographic rates specified for immatures? (yes, National - NIRAS)
                                  demobase.prod = data.frame(Mean=p_row$mn_base_prod, SD = p_row$sd_base_prod), # baseline productivity mean(s) and SD(s) # select first row as we batch run site and incombo assessments together
                                  demobase.survadult = data.frame(Mean =  p_row$mn_base_adsurv, SD = p_row$sd_base_adsurv), # baseline survival mean(s) and SD(s)
                                  demobase.survimmat = imm_surv, # baseline survival mean(s) and SD(s) for different immature year groups 
                                  inipop.years = as.numeric(p_row$yr), # year(s) when inital count was made
                                  inipop.inputformat = "breeding.pairs", # initial population size entered as breeding pair (SMP data) 
                                  inipop.vals = p_row$Count_bp, # initial population value(s)
                                  impacts.relative = TRUE, # specify relative impacts (% change in population from impact)
                                  impacts.splitimmat = FALSE, # different impacts specified for immatures? (no - NIRAS)
                                  impacts.year.start = as.numeric(p_row$yr)+1, # when to start impact (inipop.years+1 - NIRAS)
                                  impacts.year.end = as.numeric(p_row$yr)+61,# when does impact end? (impacts.year.start+60 - NIRAS) 
                                  impacts.splitpops = FALSE,  # different impact rates for each subpopulation, ignored if npop=1
                                  impacts.scennames = impact_name_list, # naming impact for reference
                                  impacts.prod.mean = 0, # impacts on productivity, val per population (use 0, no prod impacts - NIRAS)
                                  impacts.survadult.mean = impact_mean_list,# impacts on adult survival # choose CRM/disp/CRM+disp
                                  impacts.provideses = FALSE, # whether to provide SEs for impacts (yes, I think -NIRAS) - not provided by disp so set to false
                                  impacts.prod.se = 0, # se for productivity impact, 0 or NA?, impacts.provideses = FALSE
                                  impacts.survadult.se = 0, # se for survival impact # choose CRM/disp/CRM+disp, set 0 for moment, impacts.provideses = FALSE
                                  output.agetype = "breeding.adults", #  output population as breeding adults
                                  output.year.end = as.numeric(p_row$yr)+61, #when to end model (impacts.year.start+60 - NIRAS)
                                  output.year.start = as.numeric(p_row$yr), # when to start model report from, (inipop.years?)
                                  output.popsize.target = NULL, # don't set output pop target
                                  output.popsize.qe = NULL,# don't set output extinction risk target
                                  silent = TRUE, # suppress progress text, plots
                                  changetablenames = TRUE) # match shiny table names


# view output run1
 d1<-data.frame(Species=unlist(strsplit(k, "@"))[1], PVA_site=unlist(strsplit(k, "@"))[2], mort_type=m, run1)
 print(d1)
# bind runs together
 PVA_run<-rbind(PVA_run, d1)
  } #close m loop
}# close k loop
write.csv(PVA_run, 'C:/R4/RIAA/PVA/north_sea_run.csv', quote=F, row.names=F)

## ##################################################################
# # read back in pva run data and extract required metrics and results
## ##################################################################
library(dplyr)

pva_res<-read.csv('C:/R4/RIAA/PVA/north_sea_run.csv')

pva_res%>%filter(YRS_From_First_Imp==60  & mort_type=='site')%>%select(1:7, 42:51)

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
