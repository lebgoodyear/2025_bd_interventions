################################################################################
############################## Run full models #################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Feb 2025


# clear workspace
rm(list = ls())


################################################################################
################################### Set up #####################################


# set up path to required data
output_options <- "b1_0.25_b2_0.25_minw_0.6_sigmoid/"
path <- paste0("~/Documents/scripts/2025_bd_interventions/outputs/", output_options)

# load data
df <- as.data.frame(read.csv(paste0(path, "success_df.csv")))

# load packages
library("tidyverse")
library("brms")
library("tidybayes")

# reset factors so that reference levels are linked to highest success scores,
# based on bivairate models, plots and observation of means
# this provides more meaningful comparisons and makes the model more interpretable
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual, 
                                    levels = c("I", "H", "B"))
df$Intervention.category.itra.multi <- factor(df$Intervention.category.itra.multi, 
                                    levels = c(
                                      "Itraconazole", 
                                      "Bioaugmentation", 
                                      "Climate",
                                      "Multiple",
                                      "Other chemical", 
                                      "Population demographic"))        
df$Life.Stage <- factor(df$Life.Stage, 
                                    levels = c("Larvae", "Metamorph", "Juvenile", "Adult"))   

# zero inflated beta family cannot have values equal to 1 so set 1 as nearly 1
df$Success_brms <- df$Success
df$Success_brms[df$Success_brms == 1] <- 0.9999999     

# increase max size to enable reloo option for large sizes in loo()
options(future.globals.maxSize=8000* 1024^3)

# create output folder for brms outputs
path_out_brms <- paste0(path, "brms/full_models/")
ifelse(!dir.exists(file.path(path_out_brms)), 
        dir.create(file.path(path_out_brms), recursive=T), 
        FALSE)


################################################################################
################################ NULL MODEL ####################################


# run model with intercept only for comparison of fit of subsequent models
modb_null <- brm(Success_brms | weights(Uncertainty_weights) ~
                       1,
                data = df[!(is.na(df$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_null)
saveRDS(modb_null, paste0(path_out_brms, "modb_null.rds"))
modb_null <- readRDS(paste0(path_out_brms, "modb_null.rds"))
# diagnostics
#pp_check(modb_null)
#plot(modb_null)
loo_null <- loo(modb_null, moment_match=TRUE)
saveRDS(loo_null, paste0(path_out_brms, "loo_null.rds"))
loo_null <- readRDS(paste0(path_out_brms, "loo_null.rds"))


################################################################################
########################## FULL MODEL - NO CONTROLS ############################


# run full model with no controls
modb_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                Intervention.category.itra.multi+
                Habitat.or.Individual+
                Life.Stage+
                Therapeutic.or.Prophylactic+
                In.situ.or.Ex.situ+
                log1p(ClutchMn)+
                log1p(SVLMx)+
                Climate+
                Habitat,
                data = df[!(is.na(df$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_full)
saveRDS(modb_full, paste0(path_out_brms, "modb_full.rds"))
modb_full <- readRDS(paste0(path_out_brms, "modb_full.rds"))
# diagnostics
#pp_check(modb_full)
#plot(modb_full)
loo_full <- loo(modb_full, moment_match=TRUE)
saveRDS(loo_full, paste0(path_out_brms, "loo_full.rds"))
loo_full <- readRDS(paste0(path_out_brms, "loo_full.rds"))


################################################################################
############################### FULL MODEL #####################################


# run full model controlling for efficacy matrix
modb_full_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       Therapeutic.or.Prophylactic+
                       In.situ.or.Ex.situ+
                       log1p(ClutchMn)+
                       log1p(SVLMx)+
                       Climate+
                       Habitat+
                       (1 | Efficacy.Matrix),
                data = df[!(is.na(df$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_full_effmat)
saveRDS(modb_full_effmat, paste0(path_out_brms, "modb_full_effmat.rds"))
modb_full_effmat <- readRDS(paste0(path_out_brms, "modb_full_effmat.rds"))
# diagnostics
#pp_check(modb_full_effmat)
#plot(modb_full_effmat)
loo_full_effmat <- loo(modb_full_effmat, moment_match=TRUE)
loo_full_effmat
saveRDS(loo_full_effmat, paste0(path_out_brms, "loo_full_effmat.rds"))
loo_full_effmat <- readRDS(paste0(path_out_brms, "loo_full_effmat.rds"))


################################################################################
##################### MODEL WITH PERINTENT VAIRABLES ONLY ######################


# run full model controlling for efficacy matrix
modb_pert_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       #Therapeutic.or.Prophylactic+
                       #In.situ.or.Ex.situ+
                       #log1p(ClutchMn)+
                       #log1p(SVLMx)+
                       #Climate+
                       #Habitat+
                       (1 | Efficacy.Matrix),
                data = df[!(is.na(df$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_pert_effmat)
saveRDS(modb_pert_effmat, paste0(path_out_brms, "modb_pert_effmat.rds"))
modb_pert_effmat <- readRDS(paste0(path_out_brms, "modb_pert_effmat.rds"))
# diagnostics
#pp_check(modb_pert_effmat)
#plot(modb_pert_effmat)
loo_pert_effmat <- loo(modb_pert_effmat, moment_match=TRUE)
loo_pert_effmat
saveRDS(loo_pert_effmat, paste0(path_out_brms, "loo_pert_effmat.rds"))
loo_pert_effmat <- readRDS(paste0(path_out_brms, "loo_pert_effmat.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_pert_effmat, 
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha=0.05
)

hypothesis(
  modb_pert_effmat, 
  c("Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiOtherchemical = Intervention.category.itra.multiPopulationdemographic"),
  alpha=0.05/10
)

hypothesis(
  modb_pert_effmat, 
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha=0.05/3
)


# double check against model controlling for publication
modb_pert_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       #Therapeutic.or.Prophylactic+
                       #In.situ.or.Ex.situ+
                       #log1p(ClutchMn)+
                       #log1p(SVLMx)+
                       #Climate+
                       #Habitat+
                       (1 | Publication.Date.Authorship),
                data = df[!(is.na(df$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_pert_pub)
saveRDS(modb_pert_pub, paste0(path_out_brms, "modb_pert_pub.rds"))
modb_pert_pub <- readRDS(paste0(path_out_brms, "modb_pert_pub.rds"))
# diagnostics
#pp_check(modb_pert_pub)
#plot(modb_pert_pub)
loo_pert_pub <- loo(modb_pert_pub, moment_match=TRUE, reloo=TRUE)
loo_pert_pub
saveRDS(loo_pert_pub, paste0(path_out_brms, "loo_pert_pub.rds"))
loo_pert_pub <- readRDS(paste0(path_out_brms, "loo_pert_pub.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_pert_pub, 
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha=0.05
)

hypothesis(
  modb_pert_pub, 
  c("Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiOtherchemical = Intervention.category.itra.multiPopulationdemographic"),
  alpha=0.05/10
)

hypothesis(
  modb_pert_pub, 
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha=0.05/3
)


################################################################################
########################### MODEL COMPARISON ###################################


# load random effects model
loo_effmat <- readRDS(paste0(path_out_brms, "../grouping_variables/loo_effmat.rds"))

loo_compare(
  loo_effmat,
  loo_null,
  loo_full,
  loo_full_effmat,
  loo_pert_effmat,
  loo_pert_pub
)
# error is comparatively very large for all apart from null model
# model controlling for publication with pertinent variables is the best model within
# error that still contains variables of interest
# efficacy matrix was originally chosen as control over publication due to possible
# collinearity with publication and predictors. However, pertinent associations are
# the same in both models, but slightly stronger in the model controlling for publication


################################################################################
########################### RUN BEST MODEL AGAIN ###############################


# run best model again with full dataset (now including the datapoints with NA for 
# life history traits since no life history traits are in the best model)
modb_best <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       (1 |Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_best)
saveRDS(modb_best, paste0(path_out_brms, "../modb_best.rds"))
modb_best <- readRDS(paste0(path_out_brms, "../modb_best.rds"))
# diagnostics
#pp_check(modb_pert_effmat)
#plot(modb_pert_effmat)
loo_best <- loo(modb_best, moment_match=TRUE)
loo_best
saveRDS(loo_best, paste0(path_out_brms, "../loo_best.rds"))
loo_best <- readRDS(paste0(path_out_brms, "../loo_best.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_best, 
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha=0.05
)

hypothesis(
  modb_best, 
  c("Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiOtherchemical = Intervention.category.itra.multiPopulationdemographic"),
  alpha=0.05/10
)

hypothesis(
  modb_best, 
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha=0.05/3
)

####################### SAVE BEST MODEL OUTPUTS ################################


sink(file=paste0(path_out_brms, "../best_model_results.txt"))
print(summary(modb_best))
cat("\n\n")
hypothesis(
  modb_best, 
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha=0.05
)
cat("\n\n")
hypothesis(
  modb_best, 
  c("Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiOtherchemical = Intervention.category.itra.multiPopulationdemographic"),
  alpha=0.05/10
)
cat("\n\n")
hypothesis(
  modb_best, 
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha=0.05/3
)
sink()


################################################################################
################ MODEL CHECK USING DIFFERENT REFERENCE LEVELS ##################


############################# WORST COMBINATION ################################


# reset factors so that reference levels are different
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual, 
                                    levels = c("H", "I", "B"))
df$Intervention.category.itra.multi <- factor(df$Intervention.category.itra.multi, 
                                    levels = c(
                                      "Population demographic",                                      
                                      "Itraconazole", 
                                      "Bioaugmentation", 
                                      "Climate",
                                      "Multiple",
                                      "Other chemical"))        
df$Life.Stage <- factor(df$Life.Stage, 
                                    levels = c("Adult", "Larvae", "Metamorph", "Juvenile")) 

# run pertinent model controlling for publication
modb_wref <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       (1 |Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_wref)
saveRDS(modb_wref, paste0(path_out_brms, "check_modb_worst_reflevs.rds"))
modb_wref <- readRDS(paste0(path_out_brms, "check_modb_worst_reflevs.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_wref, 
  c("Habitat.or.IndividualI = Habitat.or.IndividualB"),
  alpha=0.05
)

hypothesis(
  modb_wref, 
  c("Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiBioaugmentation",
    "Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical"),
  alpha=0.05/10
)

hypothesis(
  modb_wref, 
  c("Life.StageLarvae = Life.StageMetamorph",
    "Life.StageLarvae = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageJuvenile"),
  alpha=0.05/3
)

# cross-validation
loo_wref <- loo(modb_wref, moment_match=TRUE)
saveRDS(loo_wref, paste0(path_out_brms, "check_loo_worst_reflevs.rds"))
loo_wref <- readRDS(paste0(path_out_brms, "check_loo_worst_reflevs.rds"))

# Itraconazole and individual perform better. Larvae is borderline better.
# In agreement with best model.


############################# RANDOM COMBINATION ###############################


# reset factors so that reference levels are different
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual, 
                                    levels = c("H", "I", "B"))
df$Intervention.category.itra.multi <- factor(df$Intervention.category.itra.multi, 
                                    levels = c(
                                      "Multiple",                                    
                                      "Itraconazole", 
                                      "Bioaugmentation", 
                                      "Climate",
                                      "Other chemical",
                                      "Population demographic"))        
df$Life.Stage <- factor(df$Life.Stage, 
                                    levels = c("Juvenile", "Adult", "Larvae", "Metamorph")) 

# run pertinent model controlling for publication
modb_rref <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       (1 |Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_rref)
saveRDS(modb_rref, paste0(path_out_brms, "check_modb_rand_reflevs.rds"))
modb_rref <- readRDS(paste0(path_out_brms, "check_modb_rand_reflevs.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_rref, 
  c("Habitat.or.IndividualI = Habitat.or.IndividualB"),
  alpha=0.05
)

hypothesis(
  modb_rref, 
  c("Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiBioaugmentation",
    "Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiItraconazole = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiPopulationdemographic",
    "Intervention.category.itra.multiOtherchemical = Intervention.category.itra.multiPopulationdemographic"),
  alpha=0.05/10
)

hypothesis(
  modb_rref, 
  c("Life.StageAdult = Life.StageLarvae",
    "Life.StageAdult = Life.StageMetamorph",
    "Life.StageLarvae = Life.StageMetamorph"),
  alpha=0.05/3
)

# cross-validation
loo_rref <- loo(modb_rref, moment_match=TRUE)
saveRDS(loo_rref, paste0(path_out_brms, "check_loo_rand_reflevs.rds"))
loo_rref <- readRDS(paste0(path_out_brms, "check_loo_rand_reflevs.rds"))

# All effects disappear apart from individual being better than habitat
# Itraconazole > population demographic reappears when alpha is not adjusted
# so this relationship is not as strong as the habitat < individual
# The larave > adult is borerline when alpha is not adjusted
# Note that these results both come from hypothesis testing, they are not
# model estimates


########################### COMPARE MODEL CHECKS ###############################


# compare models
loo_compare(
  loo_best,
  loo_wref,
  loo_rref
)

# best model is just outside of standard error of other two in being best model,
# however difference is very small.


## end of script
