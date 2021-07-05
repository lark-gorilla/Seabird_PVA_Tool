## ##################################################################
## Created 16 Jan 2019, last modified 12 August 2019
## Version 1.9: various miscellaneous bug fixes, so code for all examples
##  now runs without error messages -- still needs checking for sanity
## Version 2.1 on 14 Feb 2019: various bugs fixed following feedback from Deena
## Version 2.2 on 15 Feb 2019: fixed inputs so there is compatability between
##   the structure of all high-level functions, so they all call the same 
##   internal function ("nepva.calcs"), and so all calculations are relegated
##   to that function (and child functions of it)
## ###########################################################
## Version 4.2 on 12 Aug 2019: added "nburn" to each example
## Version 4.11: added in Example 7, previously in a separate file
## ##################################################################

rm(list=objects())

setwd("C:/niras_rproj/Seabird_PVA_Tool/R/Rpackage/")

ff <- list.files(".", pattern="functions") ; for(k in 1:length(ff)){ source(ff[k])} ## Automated Version 4.8

modeoptions <- read.csv("ModeOptions.csv") ## Added Version 2.1

## ##################################################

library(popbio)

## ##################################################################
## Implement NIRAS/NE parameter choices
## lookup values
## ##################################################################


lookup.dir <- "C://Users//adam//Work//AlreadyBackedUp//CEH-NE-PVA-Tool//Current//Data//Lookup//Version3.4//Tables//"
## lookup.dir <- "C://Users//adam//Work//Active//NE-PVA//Data//Lookup//Version3.4//Tables//"

example7 <- nepva.calcdefaults(Species = "Herring Gull", 
                               poolregtype.BS = "Global", poolregion.BS = "Global", 
                               sourcepop.surv = "National", lookup.dir = 'C:/niras_rproj/Seabird_PVA_Tool/R/Rpackage/')

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
                                  demobase.prod = data.frame(Mean=0.5, SD=0.05), 
                                  demobase.survadult = data.frame(Mean =  0.88, SD = 0.03), # baseline survival mean(s) and SD(s)
                                  inipop.years = c(2012), # year(s) when inital count was made
                                  inipop.inputformat = "breeding.adults", # initial population size entered as breeding adults 
                                  inipop.vals = c(791), # initial population value(s)
                                  impacts.relative = TRUE, # specify relative impacts (% change in population from impact)
                                  impacts.splitimmat = FALSE, # different impacts specified for immatures? (no - NIRAS)
                                  impacts.year.start = 2013, # when to start impact (inipop.years+1 - NIRAS)
                                  impacts.year.end = 2020,# when does impact end? (impacts.year.start+60 - NIRAS) 
                                  impacts.splitpops = FALSE,  # different demographic rates for each subpopulation, ignored if npop=1
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

## ##################################################################
## Created Version 1.5 on 14 Jan 2018, last modified 16 Jan 2018
## Version 2.1: "demo.splitpops = TRUE" > 
##      "demobase.splitpops = TRUE" and "impacts.splitpops = FALSE"
## ##################################################################

example1 <- nepva.simplescenarios(model.envstoch = "betagamma",
                  model.demostoch = TRUE,
                  model.dd = "nodd", model.prodmax = TRUE,
                  mbs = 4, afb = 5, npop = 2, nscen = 3,
                  sim.n = 10, sim.seed = 43576, nburn = 3,
                  demobase.specify.as.params = FALSE,
                  demobase.splitpops = TRUE, demobase.splitimmat = FALSE, 
                  demobase.prod = data.frame(Mean = c(0.3, 0.4), SD = c(0.11, 0.13)),
                  demobase.survadult = data.frame(Mean = c(0.92, 0.88), SD = c(0.02, 0.03)),
                  inipop.years = c(2012, 2015), inipop.inputformat = "breeding.pairs", 
                  inipop.vals = c(791, 113),
                  impacts.splitimmat = FALSE,
                  impacts.provideses = FALSE,
                  impacts.year.start = 2035, impacts.year.end = 2055,
                  impacts.splitpops = FALSE,
                  impacts.scennames = c("ice", "fish", "bob"),
                  impacts.prod.mean = c(+0.03, 0, 0),
                  impacts.survadult.mean = c(+0.05, +0.12, +0.17),
                  output.agetype = "age.separated",
                  output.year.end = 2070, output.year.start = 2015,
                  output.popsize.target = 150, output.popsize.qe = 10, silent = FALSE, 
                  changetablenames = TRUE)

## ############################################

example2 <- nepva.simplescenarios(model.envstoch = "betagamma",
                   model.demostoch = TRUE,
                   model.dd = "dduloglin", model.prodmax = TRUE,
                   mbs = 4, afb = 5, npop = 2, nscen = 3,
                   sim.n = 10, sim.seed = 43576, nburn = 5,
                   demobase.specify.as.params = FALSE, 
                   demobase.splitpops = TRUE, demobase.splitimmat = FALSE,
                   demobase.prod = data.frame(Mean = c(0.4, 0.36), SD = c(0.11, 0.13), DD = c(-0.01, -0.02)),
                   demobase.survadult = data.frame(Mean = c(0.86, 0.89), SD = c(0.02, 0.03), DD = c(-0.03, -0.04)),
                   inipop.years = c(2012, 2015), inipop.vals = c(791, 113),
                   impacts.relative = TRUE, impacts.splitpops = FALSE,
                   impacts.splitimmat = FALSE,
                   impacts.provideses = TRUE,
                   impacts.year.start = 2020, impacts.year.end = 2070,
                   impacts.scennames = c("ice", "fish", "bob"),
                   impacts.prod.mean = c(-0.03, 0, 0),
                   impacts.prod.se = c(0.01, 0, 0),
                   impacts.survadult.mean = c(-0.05, -0.12, -0.17),
                   impacts.survadult.se = c(0.02, 0.02, 0.02),
                   output.agetype = "breeding.pairs",
                   output.year.end = 2070, output.year.start = 2019,
                   output.popsize.target = 150, output.popsize.qe = 10, silent = FALSE,
                   changetablenames = TRUE)

## ############################################

example3 <- nepva.validation(model.envstoch = "betagamma",
                  model.demostoch = TRUE,
                  model.dd = "nodd", model.prodmax = TRUE,
                  mbs = 4, afb = 5, nburn = 0,
                  sim.n = 10, sim.seed = 43576,
                  demobase.specify.as.params = FALSE,
                  demobase.prod = data.frame(Mean = 0.5, SD = 0.11),
                  demobase.survadult = data.frame(Mean = 0.82, SD = 0.02),
                  inipop.years = 2001, inipop.vals = 791,
                  output.agetype = "breeding.pairs",
                  output.year.end = 2018, 
                  output.popsize.target = 150, output.popsize.qe = 10,
                  output.validation.counts = c(501, 589, 612, 698),
                  output.validation.years = c(2004, 2008, 2009, 2016),
                  silent = FALSE, changetablenames = TRUE)
          
## ############################################

example4 <- nepva.sensitivity.local(model.envstoch = "betagamma",
                                 model.demostoch = TRUE,
                                 model.prodmax = TRUE,
                                 mbs = 4, afb = 5, nburn = 11,
                                 sim.n = 10, sim.seed = 43576,
                                 demobase.prod = data.frame(mean = 0.7, sd = 0.11),
                                 demobase.survadult = data.frame(mean = 0.92, sd = 0.02),
                                 inipop.years = 2012, inipop.vals = 791,
                                 impacts.year.start = 2030, impacts.year.end = 2050,
                                 impacts.relative = TRUE,
                                 impacts.prod.mean = -0.02, impacts.survadult.mean = -0.03,
                                 output.year.end = 2070,
                                 output.popsize.target = 250, output.popsize.qe = 10,
                                 sens.npvlocal = 3, sens.pcr = c(10, 5, 5, 50, 50), silent = FALSE,
                                 changetablenames = TRUE)

## ############################################

example5 <- nepva.sensitivity.global(model.envstoch = "betagamma",
                                    model.demostoch = TRUE,
                                    model.prodmax = TRUE,
                                    mbs = 4, afb = 5, nburn = 20,
                                    sim.n = 2, sim.seed = 43576,
                                    demobase.prod = data.frame(mean = 0.7, sd = 0.11),
                                    demobase.survadult = data.frame(mean = 0.92, sd = 0.02),
                                    inipop.years = 2012, inipop.vals = 791,
                                    impacts.year.start = 2070, impacts.year.end = 2070,
                                    impacts.relative = TRUE,
                                    impacts.prod.mean = -0.02, impacts.survadult.mean = -0.03,
                                    output.year.end = 2070,
                                    output.popsize.target = 250, output.popsize.qe = 10,
                                    sens.npvglobal = 2, sens.pcr = c(10, 5, 5, 50, 50), silent = TRUE,
                                    changetablenames = TRUE)

## ############################################

example6 <- nepva.simplescenarios(model.envstoch = "deterministic",
                                  model.demostoch = FALSE,
                                  model.dd = "nodd", model.prodmax = FALSE,
                                  mbs = NULL, afb = 5, npop = 1, nscen = 1,
                                  sim.n = 100, sim.seed = 43576, nburn = 6,
                                  demobase.specify.as.params = FALSE,
                                  demobase.splitpops = FALSE, demobase.splitimmat =
                                    FALSE,
                                  demobase.prod = data.frame(Mean = c(0.9), SD = c(0)),
                                  demobase.survadult = data.frame(Mean = c(0.82), SD = c(0)),
                                  inipop.years = c(2010), inipop.inputformat =
                                    "breeding.pairs",
                                  inipop.vals = c(10),
                                  impacts.splitimmat = FALSE,
                                  impacts.provideses = FALSE,
                                  impacts.year.start = 2035, impacts.year.end = 2055,
                                  impacts.splitpops = FALSE,
                                  impacts.scennames = c("bob"),
                                  impacts.prod.mean = c(0),
                                  impacts.survadult.mean = c(0),
                                  output.agetype = "breeding.adults",
                                  output.year.end = 2070, output.year.start = 2012,
                                  output.popsize.target = 150, output.popsize.qe = 10, silent = TRUE,
                                  changetablenames = TRUE)

## ############################################
## "Extra" example, previously in separate file
## "run-extra-example-for-defaults.R"

lookup.dir <- "C://Users//adam//Work//AlreadyBackedUp//CEH-NE-PVA-Tool//Current//Data//Lookup//Version3.4//Tables//"
## lookup.dir <- "C://Users//adam//Work//Active//NE-PVA//Data//Lookup//Version3.4//Tables//"

example7 <- nepva.calcdefaults(Species = "Herring Gull", 
                          poolregtype.BS = "Global", poolregion.BS = "Global", 
                          sourcepop.surv = "Skomer (1978-2010)", lookup.dir = lookup.dir)

## ############################################
## "Extra" examples, added Version 4.12

ina <- list(model.envstoch = "deterministic", model.demostoch = TRUE,
            model.dd = "nodd", model.prodmax = TRUE, mbs = 1, afb = 2,
            npop = 2, nscen = 1, sim.n = 2, sim.seed = 1, nburn = 0,
            demobase.specify.as.params = FALSE, demobase.splitpops = TRUE,
            demobase.splitimmat = TRUE, demobase.prod = array(dim=c(2,2), data = c(0.5,0.5,0,0)),
            demobase.survadult = array(dim=c(2,2), data = c(0.9,0.9,0,0)),
            demobase.survimmat = array(dim=c(2,2,2), data = c(rep(0.9,4),rep(0,4))),
            inipop.years = c(2010, 2010), inipop.inputformat = "breeding.adults",
            inipop.vals = c(50, 50), impacts.relative = TRUE, impacts.splitimmat = FALSE,
            impacts.provideses = FALSE, impacts.year.start = 2020, impacts.year.end = 2030,
            impacts.splitpops = FALSE, impacts.scennnames = "bob", impacts.prod.mean = 0.1,
            impacts.survadult.mean = 0.1, output.agetype = "breeding.pairs",
            output.year.start = 2010, output.year.end = 2040, output.popsize.target = NULL,
            output.popsize.qe = NULL, silent = TRUE, changetablenames = TRUE)

ota <- nepva.batchmode(ina, runtype = "simplescenarios") ; ota$tab[1:6,1:10]

inb <- ina ; inb$model.demostoch <- FALSE

otb <- nepva.batchmode(inb, runtype = "simplescenarios") ; otb$tab[1:6,1:10]

inc <- ina ; inc$impacts.splitimmat <- TRUE ; inc$impacts.survimmat.mean <- c(0.1, 0.1)

otc <- nepva.batchmode(inc, runtype = "simplescenarios") ; otc$tab[1:6,1:10]

## ###############################################################################################
