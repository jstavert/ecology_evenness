---
title: "Plant species dominance increases pollination complementarity and plant reproductive function"
author: "Jamie R. Stavert"
date: "4 April 2019"
---
  
#############################################################
#load required packages----
#############################################################

library(emmeans)
library(car)
library(lme4)
library(nlme)
library(plyr)
library(dplyr)
library(MuMIn)

#############################################################
#read in network metric data----
#############################################################

seed_production <- read.csv(file = "data/seed_production.csv", header=TRUE, row.names=NULL)

#############################################################
#read in seed set data----
#############################################################

network_metrics <- read.csv(file = "data/network_metrics.csv", header=TRUE, row.names=NULL)

#############################################################
#model seed production----
#############################################################

#run seed production model
seed.prod.mod <- lme(scaled_seed_production ~ treatment*plant_species,
                       random = ~1|cage_number,
                       data = seed_production,
                       contrasts=list(treatment=contr.sum, plant_species=contr.sum),
                       na.action = na.omit)

#run wald test
car::Anova(seed.prod.mod, type = "III")

#run emmeans pairwise comparisons
seed.prod.comp <- emmeans(seed.prod.mod, pairwise ~ plant_species|treatment)

#pairwise contrasts among treatments within plant species
CLD(seed.prod.comp, Letters = letters, level = .95, adjust = "fdr", by = "plant_species")

#pairwise contrasts amoung plant species within treatments
CLD(seed.prod.comp, Letters = letters, level = .95, adjust = "fdr", by = "treatment")

#############################################################
#model seed production as function of network metrics----
#############################################################

#run seed production - network metric model
network.seedprod.mod <- lme(scaled_seed_production ~ visitation_evenness+pollinator_sharing+no_links+selectivity+visitation_rate,
                           random = list(~1|cage_number,~1|plant_species),
                           data = seed_production)

#run dredge for model selection
network.seedprod.mod.dredge <- dredge(network.seedprod.mod)

#run model averaging for 95% model confidence set
confset.95p <- get.models(network.seedprod.mod.dredge, cumsum(weight) <= .95)
avg <- model.avg(confset.95p, beta)
summary(avg)

#################################################################
#run models for network-level interaction evenness----
#################################################################

#extract unique values
interaction.evenness <- distinct(network_metrics, cage_number, treatment, interaction_evenness)

#run interaction evenness model
interaction.evenness.mod <- lm(interaction_evenness ~ treatment, data = interaction.evenness)

#run lsmeans pairwise comparisons
interaction.evenness.comp <- emmeans(interaction.evenness.mod, pairwise ~ treatment)
CLD(interaction.evenness.comp, Letters = letters, level = .95, adjust = "fdr")

#################################################################
#run models for species-level metrics----
#################################################################
#note that these models generate the results presented in the main text only
#models can be altered to generate results presented in Appendix S1 by adding the 
#treatment*species interaction and removing the species random effect.

#run model for interaction evenness
selected <- c("visitation_evenness")
visitation.evenness <- network_metrics[network_metrics$metric %in% selected,]
visitation.evenness.model <- glmer(value ~ treatment + (1|cage_number) + (1|species),
                        data = visitation.evenness,
                        family = Gamma(link = "log"),
                        na.action = na.omit,
                        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
#run emmeans pairwise comparisons
visitation.evenness.comp <- emmeans(visitation.evenness.model, pairwise ~ treatment)
CLD(visitation.evenness.comp, Letters = letters, level = .95, type = "response", adjust = "fdr")

#run model for vistation rate
selected <- c("visitation_rate")
visitation.rate <- network_metrics[network_metrics$metric %in% selected,]
visitation.rate.model <- glmer(value ~ treatment + (1|cage_number) + (1|species),
                          data = visitation.rate,
                          family = Gamma(link = "log"),
                          na.action = na.omit,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
#run emmeans pairwise comparisons
visitation.rate.comp <- emmeans(visitation.rate.model, pairwise ~ treatment)
CLD(visitation.rate.comp, Letters = letters, level = .95, type = "response", adjust = "fdr")

#run model for pollinator selectivity
selectivity.pol <- network_metrics[(network_metrics$trophic_level == "pollinator" & network_metrics$metric=="selectivity"), ]
selectivity.pol.mod <- glmer(value ~ treatment + (1|cage_number) + (1|species),
                                  data = selectivity.pol,
                                  family = Gamma(link = "log"),
                                  na.action=na.exclude,
                                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
#run pairwise comparisons
selectivity.pol.comp <- emmeans(selectivity.pol.mod, pairwise ~ treatment)
CLD(selectivity.pol.comp, Letters = letters, level = .95, type = "response", adjust = "fdr")

#run model for plant selectivity
selectivity.pla <- network_metrics[(network_metrics$trophic_level == "plant" & network_metrics$metric=="selectivity"), ]
selectivity.pla.mod <- glmer(value ~ treatment + (1|cage_number) + (1|species),
                                  data = selectivity.pla,
                                  family = gaussian(link = "log"),
                                  na.action=na.exclude,
                                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
#run pairwise comparisons
selectivity.pla.comp <- emmeans(selectivity.pla.mod, pairwise ~ treatment)
CLD(selectivity.pla.comp, Letters = letters, level = .95, type="response", adjust = "fdr")

#run model for pollinator sharing
selected <- c("pollinator_sharing")
pol.sharing <- network_metrics[network_metrics$metric %in% selected,]
pol.sharing.mod <- lme(value ~ treatment,
                       random = list(~1|cage_number, ~1|species),
                       weights=varIdent(form=~1|species),
                       data = pol.sharing,
                       na.action=na.exclude)
#run lsmeans pairwise comparisons
pol.sharing.comp <- emmeans(pol.sharing.mod, pairwise ~ treatment)
CLD(pol.sharing.comp, Letters = letters, level = .95, adjust = "fdr")

#run model for pollinator niche overlap
selected <- c("niche_overlap")
pol.niche.overlap <- network_metrics[network_metrics$metric %in% selected,]
pol.niche.overlap.mod <- lme(value ~ treatment,
                             random = list(~1|cage_number, ~1|species),
                             data = pol.niche.overlap,
                             na.action=na.exclude)
#run lsmeans pairwise comparisons
pol.niche.overlap.comp <- emmeans(pol.niche.overlap.mod, pairwise ~ treatment)
CLD(pol.niche.overlap.comp, Letters = letters, level = .95, adjust = "fdr")

#run model for number of pollinator links
links.pol <- network_metrics[(network_metrics$trophic_level == "pollinator" & network_metrics$metric=="no_links"), ]
links.pol.mod <- glmer(value ~ treatment + (1|cage_number) + (1|species),
                         data = links.pol,
                         family = poisson(link = "log"),
                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
#run lsmeans pairwise comparisons
links.pol.comp <- emmeans(links.pol.mod, pairwise ~ treatment)
CLD(links.pol.comp, Letters = letters, level = .95, type="response", adjust = "fdr")

#run model for number of plant links
links.plant <- network_metrics[(network_metrics$trophic_level == "plant" & network_metrics$metric=="no_links"), ]
links.plant.mod <- glmer(value ~ treatment + (1|cage_number) + (1|species),
                         data = links.plant,
                         family = poisson(link = "log"),
                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000)))
#run lsmeans pairwise comparisons
links.plant.comp <- emmeans(links.plant.mod, pairwise ~ treatment)
CLD(links.plant.comp, Letters = letters, level = .95, type = "response", adjust = "fdr")

#############################################################
#END
#############################################################