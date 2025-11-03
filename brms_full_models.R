################################################################################
############################## RUN FULL MODELS #################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Oct 2025


# clear workspace
rm(list = ls())


################################################################################
################################### SET UP #####################################
################################################################################


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
df$Intervention.category.itra.multi <- factor(df$Intervention.category.itra.multi)
df$Intervention.category.itra.multi <- relevel(df$Intervention.category.itra.multi, ref = "Itraconazole")
df$TreatmentType <- as.factor(df$TreatmentType)
df$TreatmentType <- relevel(df$TreatmentType, ref = "Itraconazole")
df$Life.Stage <- factor(df$Life.Stage,
                        levels = c("Larvae", "Metamorph", "Juvenile", "Adult"))

# zero inflated beta family cannot have values equal to 1 so set 1 as nearly 1
df$Success_brms <- df$Success
df$Success_brms[df$Success_brms == 1] <- 0.9999999

# increase max size to enable reloo option for large sizes in loo()
options(future.globals.maxSize = 8000 * 1024^3)

# create output folder for brms outputs
path_out_brms <- paste0(path, "brms/full_models/")
ifelse(!dir.exists(file.path(path_out_brms)),
       dir.create(file.path(path_out_brms), recursive = TRUE),
       FALSE)

# load random effects models to compare
modb_effmat <- readRDS(paste0(path_out_brms, "../grouping_variables/modb_effmat.rds"))
modb_pub <- readRDS(paste0(path_out_brms, "../grouping_variables/modb_pub.rds"))
loo_effmat <- readRDS(paste0(path_out_brms, "../grouping_variables/loo_effmat.rds"))
loo_pub <- readRDS(paste0(path_out_brms, "../grouping_variables/loo_pub.rds"))


################################################################################
######################### INTERVENTION CATEGORY MODELS #########################
################################################################################


################################################################################
############################### FULL MODEL #####################################


# run full model controlling for efficacy matrix
modb_full_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~
                          Intervention.category.itra.multi +
                            Habitat.or.Individual +
                            Life.Stage +
                            Therapeutic.or.Prophylactic +
                            In.situ.or.Ex.situ +
                            log1p(ClutchMn) +
                            log1p(SVLMx) +
                            Climate +
                            Habitat +
                            (1 | Efficacy.Matrix),
                        data = df[!(is.na(df$SVLMx)), ],
                        family = zero_inflated_beta(),
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.99),
                        iter = 1e4,
                        cores = 4,
                        chains = 4)
summary(modb_full_effmat)
saveRDS(modb_full_effmat, paste0(path_out_brms, "intcat/modb_full_effmat.rds"))
modb_full_effmat <- readRDS(paste0(path_out_brms, "intcat/modb_full_effmat.rds"))
# leave one out cross-validation for model comparison
loo_full_effmat <- loo(modb_full_effmat, moment_match = TRUE)
saveRDS(loo_full_effmat, paste0(path_out_brms, "incat/loo_full_effmat.rds"))
loo_full_effmat <- readRDS(paste0(path_out_brms, "intcat/loo_full_effmat.rds"))


# run full model controlling for publication
modb_full_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi +
                         Habitat.or.Individual +
                         Life.Stage +
                         Therapeutic.or.Prophylactic +
                         In.situ.or.Ex.situ +
                         log1p(ClutchMn) +
                         log1p(SVLMx) +
                         Climate +
                         Habitat +
                         (1 | Publication.Date.Authorship),
                     data = df[!(is.na(df$SVLMx)), ],
                     family = zero_inflated_beta(),
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.99),
                     iter = 1e4,
                     cores = 4,
                     chains = 4)
summary(modb_full_pub)
saveRDS(modb_full_pub, paste0(path_out_brms, "intcat/modb_full_pub.rds"))
modb_full_pub <- readRDS(paste0(path_out_brms, "intcat/modb_full_pub.rds"))
# loo for comparison
loo_full_pub <- loo(modb_full_pub, moment_match = TRUE)
saveRDS(loo_full_pub, paste0(path_out_brms, "intcat/loo_full_pub.rds"))
loo_full_pub <- readRDS(paste0(path_out_brms, "intcat/loo_full_pub.rds"))


################################################################################
##################### MODEL WITH PERINTENT VAIRABLES ONLY ######################


# run full model controlling for efficacy matrix
modb_pert_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~
                          Intervention.category.itra.multi +
                            Habitat.or.Individual +
                            Life.Stage +
                            (1 | Efficacy.Matrix),
                        data = df[!(is.na(df$SVLMx)), ],
                        family = zero_inflated_beta(),
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.99),
                        iter = 1e4,
                        cores = 4,
                        chains = 4)
summary(modb_pert_effmat)
saveRDS(modb_pert_effmat, paste0(path_out_brms, "intcat/modb_pert_effmat.rds"))
modb_pert_effmat <- readRDS(paste0(path_out_brms, "intcat/modb_pert_effmat.rds"))
# loo for comparison
loo_pert_effmat <- loo(modb_pert_effmat, moment_match = TRUE)
saveRDS(loo_pert_effmat, paste0(path_out_brms, "intcat/loo_pert_effmat.rds"))
loo_pert_effmat <- readRDS(paste0(path_out_brms, "intcat/loo_pert_effmat.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_pert_effmat,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

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
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_pert_effmat,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis


# run full model controlling for publication
modb_pert_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi +
                         Habitat.or.Individual +
                         Life.Stage +
                         (1 | Publication.Date.Authorship),
                     data = df[!(is.na(df$SVLMx)), ],
                     family = zero_inflated_beta(),
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.99),
                     iter = 1e4,
                     cores = 4,
                     chains = 4)
summary(modb_pert_pub)
saveRDS(modb_pert_pub, paste0(path_out_brms, "intcat/modb_pert_pub.rds"))
modb_pert_pub <- readRDS(paste0(path_out_brms, "intcat/modb_pert_pub.rds"))
# loo for comparison
loo_pert_pub <- loo(modb_pert_pub, moment_match = TRUE, reloo = TRUE)
saveRDS(loo_pert_pub, paste0(path_out_brms, "intcat/loo_pert_pub.rds"))
loo_pert_pub <- readRDS(paste0(path_out_brms, "intcat/loo_pert_pub.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_pert_pub,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

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
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_pert_pub,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis


################################################################################
########################### MODEL COMPARISON ###################################


loo_compare(
  loo_effmat,
  loo_pub,
  loo_full_effmat,
  loo_full_pub,
  loo_pert_effmat,
  loo_pert_pub
)

# error is comparatively very large for all models
# model controlling for publication with pertinent variables is the best model
# within error that still contains variables of interest


################################################################################
############################ TREATMENT TYPE MODELS #############################
################################################################################


table(df$TreatmentType)
# remove low sample sized groups and groups that cause convergence issues
df_treatsub <- df[!(df$TreatmentType %in% c("BMP-NTf2", "General Tonic", "Medicinal antiparasitic")), ]
df_treatsub <- droplevels(df_treatsub)


################################################################################
############################### RANDOM MODELS ##################################


# run random models on treat type subset to enable loo comparison

# publication
modb_pubt <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                  (1 | Publication.Date.Authorship),
                data = df_treatsub,
                family = zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter = 1e4,
                cores = 4,
                chains = 4)
summary(modb_pubt)
saveRDS(modb_pubt, paste0(path_out_brms, "modb_pubt.rds"))
modb_pubt <- readRDS(paste0(path_out_brms, "modb_pubt.rds"))
# loo for comparison
loo_pubt <- loo(modb_pubt, moment_match = TRUE)
saveRDS(loo_pubt, paste0(path_out_brms, "loo_pubt.rds"))
loo_pubt <- readRDS(paste0(path_out_brms, "loo_pubt.rds"))


# efficacy matrix
modb_effmatt <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                  (1 | Efficacy.Matrix),
                data = df_treatsub,
                family = zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter = 1e4,
                cores = 4,
                chains = 4)
summary(modb_effmatt)
saveRDS(modb_effmatt, paste0(path_out_brms, "modb_effmatt.rds"))
modb_effmatt <- readRDS(paste0(path_out_brms, "modb_effmatt.rds"))
# loo for comparison
loo_effmatt <- loo(modb_effmatt, moment_match = TRUE)
saveRDS(loo_effmatt, paste0(path_out_brms, "loo_effmatt.rds"))
loo_effmatt <- readRDS(paste0(path_out_brms, "loo_effmatt.rds"))


################################################################################
################################### FULL MODEL #################################


# run full model controlling for efficacy matrix
modb_full_effmatt <- brm(Success_brms | weights(Uncertainty_weights) ~
                           TreatmentType +
                             Habitat.or.Individual +
                             Life.Stage +
                             Therapeutic.or.Prophylactic +
                             In.situ.or.Ex.situ +
                             (1 | Efficacy.Matrix),
                         data = df_treatsub,
                         family = zero_inflated_beta(),
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.99),
                         iter = 1e4,
                         cores = 4,
                         chains = 4)
summary(modb_full_effmatt)
saveRDS(modb_full_effmatt, paste0(path_out_brms, "modb_full_effmatt.rds"))
modb_full_effmatt <- readRDS(paste0(path_out_brms, "modb_full_effmatt.rds"))
# run loo for comparison
loo_full_effmatt <- loo(modb_full_effmatt, moment_match = TRUE)
saveRDS(loo_full_effmatt, paste0(path_out_brms, "loo_full_effmatt.rds"))
loo_full_effmatt <- readRDS(paste0(path_out_brms, "loo_full_effmatt.rds"))


# run full model controlling for publication
modb_full_pubt <- brm(Success_brms | weights(Uncertainty_weights) ~
                        TreatmentType +
                          Habitat.or.Individual +
                          Life.Stage +
                          Therapeutic.or.Prophylactic +
                          In.situ.or.Ex.situ +
                          (1 | Publication.Date.Authorship),
                      data = df_treatsub,
                      family = zero_inflated_beta(),
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.99),
                      iter = 1e4,
                      cores = 4,
                      chains = 4)
summary(modb_full_pubt)
saveRDS(modb_full_pubt, paste0(path_out_brms, "modb_full_pubt.rds"))
modb_full_pubt <- readRDS(paste0(path_out_brms, "modb_full_pubt.rds"))
# run loo for comparison
loo_full_pubt <- loo(modb_full_pubt, moment_match = TRUE)
saveRDS(loo_full_pubt, paste0(path_out_brms, "loo_full_pubt.rds"))
loo_full_pubt <- readRDS(paste0(path_out_brms, "loo_full_pubt.rds"))


################################################################################
##################### MODEL WITH PERINTENT VAIRABLES ONLY ######################


# run full model controlling for efficacy matrix
modb_pert_effmatt <- brm(Success_brms | weights(Uncertainty_weights) ~
                           TreatmentType +
                             Habitat.or.Individual +
                             Life.Stage +
                             (1 | Efficacy.Matrix),
                         data = df_treatsub,
                         family = zero_inflated_beta(),
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.99),
                         iter = 1e4,
                         cores = 4,
                         chains = 4)
summary(modb_pert_effmatt)
saveRDS(modb_pert_effmatt, paste0(path_out_brms, "modb_pert_effmatt.rds"))
modb_pert_effmatt <- readRDS(paste0(path_out_brms, "modb_pert_effmatt.rds"))
# loo for comparison
loo_pert_effmatt <- loo(modb_pert_effmatt, moment_match = TRUE)
saveRDS(loo_pert_effmatt, paste0(path_out_brms, "loo_pert_effmatt.rds"))
loo_pert_effmatt <- readRDS(paste0(path_out_brms, "loo_pert_effmatt.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_pert_effmatt,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_pert_effmatt,
  c("TreatmentTypeMultiple = TreatmentTypeClimate",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeMultiple = TreatmentTypePopulationdemographic",
    "TreatmentTypeMultiple = TreatmentTypeBioaugmentation",
    "TreatmentTypeMultiple = TreatmentTypeFungicide",
    "TreatmentTypeMultiple = TreatmentTypeDisinfectant",
    "TreatmentTypeMultiple = TreatmentTypeHerbicide",
    "TreatmentTypeMultiple = TreatmentTypeInsecticide",
    "TreatmentTypeMultiple = TreatmentTypeSalt",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeClimate = TreatmentTypePopulationdemographic",
    "TreatmentTypeClimate = TreatmentTypeBioaugmentation",
    "TreatmentTypeClimate = TreatmentTypeFungicide",
    "TreatmentTypeClimate = TreatmentTypeDisinfectant",
    "TreatmentTypeClimate = TreatmentTypeHerbicide",
    "TreatmentTypeClimate = TreatmentTypeInsecticide",
    "TreatmentTypeClimate = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeMedicinalantifungal = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantifungal= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantibiotic= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeSalt",
    "TreatmentTypePopulationdemographic= TreatmentTypeBioaugmentation",
    "TreatmentTypePopulationdemographic = TreatmentTypeFungicide",
    "TreatmentTypePopulationdemographic = TreatmentTypeDisinfectant",
    "TreatmentTypePopulationdemographic = TreatmentTypeHerbicide",
    "TreatmentTypePopulationdemographic = TreatmentTypeInsecticide",
    "TreatmentTypePopulationdemographic = TreatmentTypeSalt",
    "TreatmentTypeBioaugmentation = TreatmentTypeFungicide",
    "TreatmentTypeBioaugmentation = TreatmentTypeDisinfectant",
    "TreatmentTypeBioaugmentation = TreatmentTypeHerbicide",
    "TreatmentTypeBioaugmentation = TreatmentTypeInsecticide",
    "TreatmentTypeBioaugmentation = TreatmentTypeSalt",
    "TreatmentTypeFungicide = TreatmentTypeDisinfectant",
    "TreatmentTypeFungicide = TreatmentTypeHerbicide",
    "TreatmentTypeFungicide = TreatmentTypeInsecticide",
    "TreatmentTypeFungicide = TreatmentTypeSalt",
    "TreatmentTypeHerbicide = TreatmentTypeInsecticide",
    "TreatmentTypeHerbicide = TreatmentTypeSalt",
    "TreatmentTypeInsecticide = TreatmentTypeSalt"),
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_pert_effmatt,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis


# run full model controlling for publication
modb_pert_pubt <- brm(Success_brms | weights(Uncertainty_weights) ~
                        TreatmentType +
                          Habitat.or.Individual +
                          Life.Stage +
                          (1 | Publication.Date.Authorship),
                      data = df_treatsub,
                      family = zero_inflated_beta(),
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.99),
                      iter = 1e4,
                      cores = 4,
                      chains = 4)
summary(modb_pert_pubt)
saveRDS(modb_pert_pubt, paste0(path_out_brms, "modb_pert_pubt.rds"))
modb_pert_pubt <- readRDS(paste0(path_out_brms, "modb_pert_pubt.rds"))
# loo for comparison
loo_pert_pubt <- loo(modb_pert_pubt, moment_match = TRUE, reloo = TRUE)
saveRDS(loo_pert_pubt, paste0(path_out_brms, "loo_pert_pubt.rds"))
loo_pert_pubt <- readRDS(paste0(path_out_brms, "loo_pert_pubt.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_pert_pubt,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_pert_pubt,
  c("TreatmentTypeMultiple = TreatmentTypeClimate",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeMultiple = TreatmentTypePopulationdemographic",
    "TreatmentTypeMultiple = TreatmentTypeBioaugmentation",
    "TreatmentTypeMultiple = TreatmentTypeFungicide",
    "TreatmentTypeMultiple = TreatmentTypeDisinfectant",
    "TreatmentTypeMultiple = TreatmentTypeHerbicide",
    "TreatmentTypeMultiple = TreatmentTypeInsecticide",
    "TreatmentTypeMultiple = TreatmentTypeSalt",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeClimate = TreatmentTypePopulationdemographic",
    "TreatmentTypeClimate = TreatmentTypeBioaugmentation",
    "TreatmentTypeClimate = TreatmentTypeFungicide",
    "TreatmentTypeClimate = TreatmentTypeDisinfectant",
    "TreatmentTypeClimate = TreatmentTypeHerbicide",
    "TreatmentTypeClimate = TreatmentTypeInsecticide",
    "TreatmentTypeClimate = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeMedicinalantifungal = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantifungal= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantibiotic= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeSalt",
    "TreatmentTypePopulationdemographic= TreatmentTypeBioaugmentation",
    "TreatmentTypePopulationdemographic = TreatmentTypeFungicide",
    "TreatmentTypePopulationdemographic = TreatmentTypeDisinfectant",
    "TreatmentTypePopulationdemographic = TreatmentTypeHerbicide",
    "TreatmentTypePopulationdemographic = TreatmentTypeInsecticide",
    "TreatmentTypePopulationdemographic = TreatmentTypeSalt",
    "TreatmentTypeBioaugmentation = TreatmentTypeFungicide",
    "TreatmentTypeBioaugmentation = TreatmentTypeDisinfectant",
    "TreatmentTypeBioaugmentation = TreatmentTypeHerbicide",
    "TreatmentTypeBioaugmentation = TreatmentTypeInsecticide",
    "TreatmentTypeBioaugmentation = TreatmentTypeSalt",
    "TreatmentTypeFungicide = TreatmentTypeDisinfectant",
    "TreatmentTypeFungicide = TreatmentTypeHerbicide",
    "TreatmentTypeFungicide = TreatmentTypeInsecticide",
    "TreatmentTypeFungicide = TreatmentTypeSalt",
    "TreatmentTypeHerbicide = TreatmentTypeInsecticide",
    "TreatmentTypeHerbicide = TreatmentTypeSalt",
    "TreatmentTypeInsecticide = TreatmentTypeSalt"),
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_pert_pubt,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis


################################################################################
########################### MODEL COMPARISON ###################################


loo_compare(
  loo_effmatt,
  loo_pubt,
  loo_full_effmatt,
  loo_full_pubt,
  loo_pert_effmatt,
  loo_pert_pubt
)

# random effects models are better than models with predictors (even with the
# large standard errors). this is likely because some of the groupings involving
# treatment type are too small, leaving the data with high sparsity, along with
# more parameters to estimate than for the intervention category models

# the intervention category models are, then, more informative, depsite being
# more general and so will be used for the post hoc analysis

# note we can't compare across the two types of models beause treatment type
# had to have some observations removed due to low sample sizes and so the
# number of observations is not equal


################################################################################
########################### RUN BEST MODEL AGAIN ###############################
################################################################################


# treatment type groupings are really too sparse to provide useful information
# so use intervention category model intsead (even though it is more general)

# run best model again with full dataset (now including the datapoints with NA
# for life history traits since no life history traits are in the best model)

modb_best <- brm(Success_brms | weights(Uncertainty_weights) ~
                   Intervention.category.itra.multi +
                     Habitat.or.Individual +
                     Life.Stage +
                     (1 | Publication.Date.Authorship),
                 data = df,
                 family = zero_inflated_beta(),
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99),
                 iter = 1e4,
                 cores = 4,
                 chains = 4)
summary(modb_best)
saveRDS(modb_best, paste0(path_out_brms, "../modb_best.rds"))
modb_best <- readRDS(paste0(path_out_brms, "../modb_best.rds"))
# loo for comparison
loo_best <- loo(modb_best, moment_match = TRUE, reloo = TRUE)
loo_best
saveRDS(loo_best, paste0(path_out_brms, "../loo_best.rds"))
loo_best <- readRDS(paste0(path_out_brms, "../loo_best.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_best,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

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
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_best,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis


################################################################################
####################### SAVE BEST MODEL OUTPUTS ################################


sink(file = paste0(path_out_brms, "../best_model_results.txt"))
print(summary(modb_best))
cat("\n\n")
hypothesis(
  modb_best,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis
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
  alpha = 0.05
)$hypothesis
cat("\n\n")
hypothesis(
  modb_best,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis
sink()



################################################################################
###################### SAVE BEST TREATMENT MODEL OUTPUTS #######################


# also save outputs for best treatment model
saveRDS(modb_pert_pubt, paste0(path_out_brms, "../modb_bestt.rds"))
modb_pert_pubt <- readRDS(paste0(path_out_brms, "../modb_bestt.rds"))
saveRDS(loo_pert_pubt, paste0(path_out_brms, "../loo_bestt.rds"))
loo_pert_pubt <- readRDS(paste0(path_out_brms, "../loo_bestt.rds"))

# use publication rather than efficacy matrix for consistency
sink(file = paste0(path_out_brms, "../best_treatment_model_results.txt"))
print(summary(modb_pert_pubt))
cat("\n\n")
hypothesis(
  modb_pert_pubt,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis
cat("\n\n")
hypothesis(
  modb_pert_pubt,
  c("TreatmentTypeMultiple = TreatmentTypeClimate",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeMultiple = TreatmentTypePopulationdemographic",
    "TreatmentTypeMultiple = TreatmentTypeBioaugmentation",
    "TreatmentTypeMultiple = TreatmentTypeFungicide",
    "TreatmentTypeMultiple = TreatmentTypeDisinfectant",
    "TreatmentTypeMultiple = TreatmentTypeHerbicide",
    "TreatmentTypeMultiple = TreatmentTypeInsecticide",
    "TreatmentTypeMultiple = TreatmentTypeSalt",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeClimate = TreatmentTypePopulationdemographic",
    "TreatmentTypeClimate = TreatmentTypeBioaugmentation",
    "TreatmentTypeClimate = TreatmentTypeFungicide",
    "TreatmentTypeClimate = TreatmentTypeDisinfectant",
    "TreatmentTypeClimate = TreatmentTypeHerbicide",
    "TreatmentTypeClimate = TreatmentTypeInsecticide",
    "TreatmentTypeClimate = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeMedicinalantibiotic",
    "TreatmentTypeMedicinalantifungal = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantifungal= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantibiotic= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantibiotic = TreatmentTypeSalt",
    "TreatmentTypePopulationdemographic= TreatmentTypeBioaugmentation",
    "TreatmentTypePopulationdemographic = TreatmentTypeFungicide",
    "TreatmentTypePopulationdemographic = TreatmentTypeDisinfectant",
    "TreatmentTypePopulationdemographic = TreatmentTypeHerbicide",
    "TreatmentTypePopulationdemographic = TreatmentTypeInsecticide",
    "TreatmentTypePopulationdemographic = TreatmentTypeSalt",
    "TreatmentTypeBioaugmentation = TreatmentTypeFungicide",
    "TreatmentTypeBioaugmentation = TreatmentTypeDisinfectant",
    "TreatmentTypeBioaugmentation = TreatmentTypeHerbicide",
    "TreatmentTypeBioaugmentation = TreatmentTypeInsecticide",
    "TreatmentTypeBioaugmentation = TreatmentTypeSalt",
    "TreatmentTypeFungicide = TreatmentTypeDisinfectant",
    "TreatmentTypeFungicide = TreatmentTypeHerbicide",
    "TreatmentTypeFungicide = TreatmentTypeInsecticide",
    "TreatmentTypeFungicide = TreatmentTypeSalt",
    "TreatmentTypeHerbicide = TreatmentTypeInsecticide",
    "TreatmentTypeHerbicide = TreatmentTypeSalt",
    "TreatmentTypeInsecticide = TreatmentTypeSalt"),
  alpha = 0.05
)$hypothesis
cat("\n\n")
hypothesis(
  modb_pert_pubt,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)$hypothesis
sink()


################################################################################
################ MODEL CHECK USING DIFFERENT REFERENCE LEVELS ##################
################################################################################


################################################################################
############################# WORST COMBINATION ################################


# reset factors so that reference levels are different
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual,
                                   levels = c("H", "I", "B"))
df$Intervention.category.itra.multi <- factor(
  df$Intervention.category.itra.multi,
  levels = c(
    "Population demographic",
    "Itraconazole",
    "Bioaugmentation",
    "Climate",
    "Multiple",
    "Other chemical"
  )
)
df$Life.Stage <- factor(df$Life.Stage,
                        levels = c("Adult", "Larvae", "Metamorph", "Juvenile"))

# run pertinent model controlling for publication
modb_wref <- brm(Success_brms | weights(Uncertainty_weights) ~
                   Intervention.category.itra.multi +
                     Habitat.or.Individual +
                     Life.Stage +
                     (1 | Publication.Date.Authorship),
                 data = df,
                 family = zero_inflated_beta(),
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99),
                 iter = 1e4,
                 cores = 4,
                 chains = 4)
summary(modb_wref)
saveRDS(modb_wref, paste0(path_out_brms, "check_modb_worst_reflevs.rds"))
modb_wref <- readRDS(paste0(path_out_brms, "check_modb_worst_reflevs.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_wref,
  c("Habitat.or.IndividualI = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

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
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_wref,
  c("Life.StageLarvae = Life.StageMetamorph",
    "Life.StageLarvae = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageJuvenile"),
  alpha = 0.05
)$hypothesis

# cross-validation
loo_wref <- loo(modb_wref, moment_match = TRUE)
saveRDS(loo_wref, paste0(path_out_brms, "check_loo_worst_reflevs.rds"))
loo_wref <- readRDS(paste0(path_out_brms, "check_loo_worst_reflevs.rds"))

# Itraconazole and individual perform better. Larvae is borderline better.
# In agreement with best model.


################################################################################
############################# RANDOM COMBINATION ###############################


# reset factors so that reference levels are different
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual,
                                   levels = c("H", "I", "B"))
df$Intervention.category.itra.multi <- factor(
  df$Intervention.category.itra.multi,
  levels = c(
    "Multiple",
    "Itraconazole",
    "Bioaugmentation",
    "Climate",
    "Other chemical",
    "Population demographic"
  )
)
df$Life.Stage <- factor(df$Life.Stage,
                        levels = c("Juvenile", "Adult", "Larvae", "Metamorph"))

# run pertinent model controlling for publication
modb_rref <- brm(Success_brms | weights(Uncertainty_weights) ~
                   Intervention.category.itra.multi +
                     Habitat.or.Individual +
                     Life.Stage +
                     (1 | Publication.Date.Authorship),
                 data = df,
                 family = zero_inflated_beta(),
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99),
                 iter = 1e4,
                 cores = 4,
                 chains = 4)
summary(modb_rref)
saveRDS(modb_rref, paste0(path_out_brms, "check_modb_rand_reflevs.rds"))
modb_rref <- readRDS(paste0(path_out_brms, "check_modb_rand_reflevs.rds"))

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_rref,
  c("Habitat.or.IndividualI = Habitat.or.IndividualB"),
  alpha = 0.05
)$hypothesis

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
  alpha = 0.05
)$hypothesis

hypothesis(
  modb_rref,
  c("Life.StageAdult = Life.StageLarvae",
    "Life.StageAdult = Life.StageMetamorph",
    "Life.StageLarvae = Life.StageMetamorph"),
  alpha = 0.05
)$hypothesis

# cross-validation
loo_rref <- loo(modb_rref, moment_match = TRUE)
saveRDS(loo_rref, paste0(path_out_brms, "check_loo_rand_reflevs.rds"))
loo_rref <- readRDS(paste0(path_out_brms, "check_loo_rand_reflevs.rds"))

# Itraconazole and individual perform better. Larvae is borderline better.
# In agreement with best model.


################################################################################
########################### COMPARE MODEL CHECKS ###############################


# compare models
loo_compare(
  loo_best,
  loo_wref,
  loo_rref
)

# best model is just outside of standard error of random reference model
# however difference is very small, and best model uses model estimates rather
# than hypothesis testing


## end of script
