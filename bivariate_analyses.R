################################################################################
############################ Bivariate analyses ################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Sept 2025


# clear workspace
rm(list = ls())


################################################################################
################################### Set up #####################################


# set up path to required data
output_options <- "b1_0.25_b2_0.25_minw_0.6_sigmoid/"
path <- paste0(
  "~/Documents/scripts/2025_bd_interventions/outputs/",
  output_options
)

# load data
df <- as.data.frame(read.csv(paste0(path, "success_df.csv")))

# load packages
library("tidyverse")
library("FSA") # for dunnTest()
library("brms")
# load custom plotting functions
source(paste0(path, "../../functions_for_plotting.R"))


################################################################################
############################# Box plots ########################################


# create output folder for bxplots
path_out_box <- paste0(path, "box/")
ifelse(!dir.exists(file.path(path_out_box)),
       dir.create(file.path(path_out_box), recursive = TRUE),
       FALSE)

# set up dataset without 0 success
df_wozero <- df[df$Success != 0, ]

# plot and save box plots for intervention variables
# all success values
labs_life <- makeLabs(df, "Life.Stage")
lifebox <- plotBox(df, "Life.Stage", labs_life, zero = TRUE)
labs_situ <- makeLabs(df, "In.situ.or.Ex.situ")
situbox <- plotBox(df, "In.situ.or.Ex.situ", labs_situ, zero = TRUE)
labs_int <- makeLabs(df, "Intervention.category.itra.multi")
intbox <- plotBox(df, "Intervention.category.itra.multi", labs_int, zero = TRUE)
labs_treattype <- makeLabs(df, "TreatmentType")
treattypebox <- plotBox(df, "TreatmentType", labs_treattype, zero = TRUE)
labs_hab <- makeLabs(df, "Habitat.or.Individual")
habbox <- plotBox(df, "Habitat.or.Individual", labs_hab, zero = TRUE)
labs_thera <- makeLabs(df, "Therapeutic.or.Prophylactic")
therabox <- plotBox(df, "Therapeutic.or.Prophylactic", labs_thera, zero = TRUE)
# plot and save box plots for life history traits
labs_taxg <- makeLabs(df, "TaxaGroup")
taxagbox <- plotBox(df, "TaxaGroup", labs_taxg, zero = TRUE)
labs_mhab <- makeLabs(df, "Habitat")
mhabbox <- plotBox(df, "Habitat", labs_mhab, zero = TRUE)
labs_clim <- makeLabs(df, "Climate")
climbox <- plotBox(df, "Climate", labs_clim, zero = TRUE)
# rough scatterplot for the two continuous variables
ggplot(data = df, aes(x = log1p(SVLMx), y = Success)) +
  geom_jitter() +
  theme_bw()
ggsave(paste0(path_out_box, "scatter_bodysize.png"))
ggplot(data = df, aes(x = log1p(ClutchMn), y = Success)) +
  geom_jitter() +
  theme_bw()
ggsave(paste0(path_out_box, "scatter_clutchsize.png"))

# plot violins for intervention type variables
lifevio <- plotViolin(df, "Life.Stage", labs_life, path_out_box)
situvio <- plotViolin(df, "In.situ.or.Ex.situ", labs_situ, path_out_box)
intvio <- plotViolin(df, "Intervention.category.itra.multi", labs_int, path_out_box)
treattypevio <- plotViolin(df, "TreatmentType", labs_treattype, path_out_box)
habvio <- plotViolin(df, "Habitat.or.Individual", labs_hab, path_out_box)
theravio <- plotViolin(df, "Therapeutic.or.Prophylactic", labs_thera, path_out_box)

# without zero succes values
labs_life_wozero <- makeLabs(df_wozero, "Life.Stage")
lifebox_wozero <- plotBox(df_wozero, "Life.Stage", labs_life_wozero, zero = FALSE)
labs_situ_wozero <- makeLabs(df_wozero, "In.situ.or.Ex.situ")
situbox_wozero <- plotBox(df_wozero, "In.situ.or.Ex.situ", labs_situ_wozero, zero = FALSE)
labs_int_wozero <- makeLabs(df_wozero, "Intervention.category.itra.multi")
intbox_wozero <- plotBox(df_wozero, "Intervention.category.itra.multi", labs_int_wozero, zero = FALSE)
labs_treattype_wozero <- makeLabs(df_wozero, "TreatmentType")
treattypebox_wozero <- plotBox(df_wozero, "TreatmentType", labs_treattype_wozero, zero = FALSE)
labs_hab_wozero <- makeLabs(df_wozero, "Habitat.or.Individual")
habbox_wozero <- plotBox(df_wozero, "Habitat.or.Individual", labs_hab_wozero, zero = FALSE)
labs_multi_wozero <- makeLabs(df_wozero, "MultipleTreatments")
multibox_wozero <- plotBox(df_wozero, "MultipleTreatments", labs_multi_wozero, zero = FALSE)
labs_thera_wozero <- makeLabs(df_wozero, "Therapeutic.or.Prophylactic")
therabox_wozero <- plotBox(df_wozero, "Therapeutic.or.Prophylactic", labs_thera_wozero, zero = FALSE)
# plot and save box plots for life history traits
labs_taxg_wozero <- makeLabs(df_wozero, "TaxaGroup")
taxagbox_wozero <- plotBox(df_wozero, "TaxaGroup", labs_taxg_wozero, zero = FALSE)
labs_mhab_wozero <- makeLabs(df_wozero, "Habitat")
mhabbox_wozero <- plotBox(df_wozero, "Habitat", labs_mhab_wozero, zero = FALSE)
labs_clim_wozero <- makeLabs(df_wozero, "Climate")
climbox_wozero <- plotBox(df_wozero, "Climate", labs_clim_wozero, zero = FALSE)
# rough scatterplot for the two continuous variables
ggplot(data = df_wozero, aes(x = log1p(SVLMx), y = Success)) +
  geom_jitter() +
  theme_bw()
ggsave(paste0(path_out_box, "scatter_bodysize_wozero.png"))
ggplot(data = df_wozero, aes(x = log1p(ClutchMn), y = Success)) +
  geom_jitter() +
  theme_bw()
ggsave(paste0(path_out_box, "scatter_clutchsize_wozero.png"))


################################################################################
########################## Bayesian bivariate models ###########################


# reset factors so that reference levels are linked to highest success scores
# based on plots and observation of means as this provides more meaningful
# comparisons and makes the model more interpretable
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual,
                                   levels = c("I", "H", "B"))
df$Intervention.category.itra.multi <- factor(df$Intervention.category.itra.multi, 
                                              levels = c(
                                                "Itraconazole",
                                                "Bioaugmentation",
                                                "Climate",
                                                "Multiple",
                                                "Other chemical",
                                                "Population demographic"
                                              ))
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
path_out_brms <- paste0(path, "brms/bivariate/")
ifelse(!dir.exists(file.path(path_out_brms)),
       dir.create(file.path(path_out_brms), recursive = TRUE),
       FALSE)


####################### Amphibian-related predictors ###########################


# CLIMATE
modb_clim <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                   Climate,
                 data = df[!(is.na(df$SVLMx)), ],
                 family = zero_inflated_beta(),
                 save_pars = save_pars(all = TRUE),
                 iter = 1e4,
                 cores = 4,
                 chains = 4)
summary(modb_clim)
saveRDS(modb_clim, paste0(path_out_brms, "modb_clim.rds"))
modb_clim <- readRDS(paste0(path_out_brms, "modb_clim.rds"))
# Credible interval contains 0 and estimate is close to 0

modb_clim_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                       Climate +
                       (1 | Publication.Date.Authorship),
                     data = df[!(is.na(df$SVLMx)), ],
                     family = zero_inflated_beta(),
                     save_pars = save_pars(all = TRUE),
                     iter = 1e4,
                     cores = 4,
                     chains = 4)
summary(modb_clim_pub)
saveRDS(modb_clim_pub, paste0(path_out_brms, "modb_clim_pub.rds"))
modb_clim_pub <- readRDS(paste0(path_out_brms, "modb_clim_pub.rds"))
# Credible interval contains 0 and estimate is close to 0


# SVL
modb_svl <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                  log1p(SVLMx),
                data = df[!(is.na(df$SVLMx)), ],
                family = zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                iter = 1e4,
                cores = 4,
                chains = 4)
summary(modb_svl)
saveRDS(modb_svl, paste0(path_out_brms, "modb_svl.rds"))
modb_svl <- readRDS(paste0(path_out_brms, "modb_svl.rds"))
# Credible interval contains 0 and estimate is close to 0

modb_svl_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                      log1p(SVLMx) +
                      (1 | Publication.Date.Authorship),
                    data = df[!(is.na(df$SVLMx)), ],
                    family = zero_inflated_beta(),
                    save_pars = save_pars(all = TRUE),
                    iter = 1e4,
                    cores = 4,
                    chains = 4)
summary(modb_svl_pub)
saveRDS(modb_svl_pub, paste0(path_out_brms, "modb_svl_pub.rds"))
modb_svl_pub <- readRDS(paste0(path_out_brms, "modb_svl_pub.rds"))
# Credible interval contains 0 and estimate is close to 0


# CLUTCH SIZE
modb_clutch <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                     log1p(ClutchMn),
                   data = df[!(is.na(df$SVLMx)), ],
                   family = zero_inflated_beta(),
                   save_pars = save_pars(all = TRUE),
                   iter = 1e4,
                   cores = 4,
                   chains = 4)
summary(modb_clutch)
saveRDS(modb_clutch, paste0(path_out_brms, "modb_clutch.rds"))
modb_clutch <- readRDS(paste0(path_out_brms, "modb_clutch.rds"))
# Credible interval contains 0 and estimate is close to 0

modb_clutch_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                         log1p(ClutchMn) +
                         (1 | Publication.Date.Authorship),
                       data = df[!(is.na(df$SVLMx)), ],
                       family = zero_inflated_beta(),
                       save_pars = save_pars(all = TRUE),
                       iter = 1e4,
                       cores = 4,
                       chains = 4)
summary(modb_clutch_pub)
saveRDS(modb_clutch_pub, paste0(path_out_brms, "modb_clutch_pub.rds"))
modb_clutch_pub <- readRDS(paste0(path_out_brms, "modb_clutch_pub.rds"))
# Credible interval contains 0 and estimate 0


# HABITAT
modb_hab <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                  Habitat,
                data = df[!(is.na(df$SVLMx)), ],
                family = zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                iter = 1e4,
                cores = 4,
                chains = 4)
summary(modb_hab)
saveRDS(modb_hab, paste0(path_out_brms, "modb_hab.rds"))
modb_hab <- readRDS(paste0(path_out_brms, "modb_hab.rds"))
# Credible interval contains 0 and estimate is close to 0

modb_hab_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                      Habitat +
                      (1 | Publication.Date.Authorship),
                    data = df[!(is.na(df$SVLMx)), ],
                    family = zero_inflated_beta(),
                    save_pars = save_pars(all = TRUE),
                    iter = 1e4,
                    cores = 4,
                    chains = 4)
summary(modb_hab_pub)
saveRDS(modb_hab_pub, paste0(path_out_brms, "modb_hab_pub.rds"))
modb_hab_pub <- readRDS(paste0(path_out_brms, "modb_hab_pub.rds"))
# Credible interval contains 0 and estimate is close to 0


############################ Treatment predictors ##############################


# INTERVENTION CATEGORY
modb_intcat <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                     Intervention.category.itra.multi,
                   data = df,
                   family = zero_inflated_beta(),
                   save_pars = save_pars(all = TRUE),
                   iter = 1e4,
                   cores = 4,
                   chains = 4)
summary(modb_intcat)
saveRDS(modb_intcat, paste0(path_out_brms, "modb_intcat.rds"))
modb_intcat <- readRDS(paste0(path_out_brms, "modb_intcat.rds"))
# check all contrasts
hypothesis(
  modb_intcat,
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
)
# Credible interval does NOT contain 0 for multiple levels

modb_intcat_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                         Intervention.category.itra.multi +
                         (1 | Publication.Date.Authorship),
                       data = df,
                       family = zero_inflated_beta(),
                       save_pars = save_pars(all = TRUE),
                       iter = 1e4,
                       cores = 4,
                       chains = 4)
summary(modb_intcat_pub)
saveRDS(modb_intcat_pub, paste0(path_out_brms, "modb_intcat_pub.rds"))
modb_intcat_pub <- readRDS(paste0(path_out_brms, "modb_intcat_pub.rds"))
# check all contrasts
hypothesis(
  modb_intcat_pub,
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
)
# Credible interval does NOT contain 0 for multiple levels

modb_intcat_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                            Intervention.category.itra.multi +
                            (1 | Efficacy.Matrix),
                          data = df,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.95),
                          iter = 1e4,
                          cores = 4,
                          chains = 4)
summary(modb_intcat_effmat)
saveRDS(modb_intcat_effmat, paste0(path_out_brms, "modb_intcat_effmat.rds"))
modb_intcat_effmat <- readRDS(paste0(path_out_brms, "modb_intcat_effmat.rds"))
# check all contrasts
hypothesis(
  modb_intcat_effmat,
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
)
# Credible interval does NOT contain 0 for multiple levels


# TREATMENT TYPE
df_treatsub <- df[!(df$TreatmentType %in% c("BMP-NTf2", "General Tonic", "Medicinal antiparasitic")), ]
df_treatsub <- droplevels(df_treatsub)
modb_treat <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                    TreatmentType,
                  data = df_treatsub,
                  family = zero_inflated_beta(),
                  save_pars = save_pars(all = TRUE),
                  iter = 1e4,
                  cores = 4,
                  chains = 4)
summary(modb_treat)
saveRDS(modb_treat, paste0(path_out_brms, "modb_treat.rds"))
modb_treat <- readRDS(paste0(path_out_brms, "modb_treat.rds"))
# check all contrasts
hypothesis(
  modb_treat,
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
)
# Credible interval does NOT contain 0 for multiple levels

modb_treat_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                        TreatmentType +
                        (1 | Publication.Date.Authorship),
                      data = df_treatsub,
                      family = zero_inflated_beta(),
                      save_pars = save_pars(all = TRUE),
                      iter = 1e4,
                      cores = 4,
                      chains = 4)
summary(modb_treat_pub)
saveRDS(modb_treat_pub, paste0(path_out_brms, "modb_treat_pub.rds"))
modb_treat_pub <- readRDS(paste0(path_out_brms, "modb_treat_pub.rds"))
# check all contrasts
hypothesis(
  modb_treat_pub,
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
)
# Credible interval does NOT contain 0 for multiple levels

# remove medicinal antibiotics interventions (n=5) to aid with convergence
df_treatsubb <- df_treatsub[df_treatsub$TreatmentType != "Medicinal antibiotic", ]
modb_treat_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                           TreatmentType +
                           (1 | Efficacy.Matrix),
                         data = df_treatsubb,
                         family = zero_inflated_beta(),
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.99),
                         iter = 1e5,
                         cores = 4,
                         chains = 4)
summary(modb_treat_effmat)
saveRDS(modb_treat_effmat, paste0(path_out_brms, "modb_treat_effmat.rds"))
modb_treat_effmat <- readRDS(paste0(path_out_brms, "modb_treat_effmat.rds"))
# check all contrasts
hypothesis(
  modb_treat_effmat,
  c("TreatmentTypeMultiple = TreatmentTypeClimate",
    "TreatmentTypeMultiple = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeMultiple = TreatmentTypePopulationdemographic",
    "TreatmentTypeMultiple = TreatmentTypeBioaugmentation",
    "TreatmentTypeMultiple = TreatmentTypeFungicide",
    "TreatmentTypeMultiple = TreatmentTypeDisinfectant",
    "TreatmentTypeMultiple = TreatmentTypeHerbicide",
    "TreatmentTypeMultiple = TreatmentTypeInsecticide",
    "TreatmentTypeMultiple = TreatmentTypeSalt",
    "TreatmentTypeClimate = TreatmentTypeMedicinalantifungal",
    "TreatmentTypeClimate = TreatmentTypePopulationdemographic",
    "TreatmentTypeClimate = TreatmentTypeBioaugmentation",
    "TreatmentTypeClimate = TreatmentTypeFungicide",
    "TreatmentTypeClimate = TreatmentTypeDisinfectant",
    "TreatmentTypeClimate = TreatmentTypeHerbicide",
    "TreatmentTypeClimate = TreatmentTypeInsecticide",
    "TreatmentTypeClimate = TreatmentTypeSalt",
    "TreatmentTypeMedicinalantifungal = TreatmentTypePopulationdemographic",
    "TreatmentTypeMedicinalantifungal= TreatmentTypeBioaugmentation",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeFungicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeDisinfectant",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeHerbicide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeInsecticide",
    "TreatmentTypeMedicinalantifungal = TreatmentTypeSalt",
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
)
# Credible interval does NOT contain 0 for multiple levels


# HABITAT OR INDIVIDUAL
modb_habindv <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                      Habitat.or.Individual,
                    data = df,
                    family = zero_inflated_beta(),
                    save_pars = save_pars(all = TRUE),
                    iter = 1e4,
                    cores = 4,
                    chains = 4)
summary(modb_habindv)
saveRDS(modb_habindv, paste0(path_out_brms, "modb_habindv.rds"))
modb_habindv <- readRDS(paste0(path_out_brms, "modb_habindv.rds"))
# check all contrasts
hypothesis(
  modb_habindv,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)
# Credible interval contains 0 but some estimates are non-negligible
# and I-H upper CI is almost -ve

modb_habindv_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                          Habitat.or.Individual +
                          (1 | Publication.Date.Authorship),
                        data = df,
                        family = zero_inflated_beta(),
                        save_pars = save_pars(all = TRUE),
                        iter = 1e4,
                        cores = 4,
                        chains = 4)
summary(modb_habindv_pub)
saveRDS(modb_habindv_pub, paste0(path_out_brms, "modb_habindv_pub.rds"))
modb_habindv_pub <- readRDS(paste0(path_out_brms, "modb_habindv_pub.rds"))
# check all contrasts
hypothesis(
  modb_habindv_pub,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)
# I-H relationship is negative


# LIFE STAGE
modb_lifestage <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                        Life.Stage,
                      data = df,
                      family = zero_inflated_beta(),
                      save_pars = save_pars(all = TRUE),
                      iter = 1e4,
                      cores = 4,
                      chains = 4)
summary(modb_lifestage)
saveRDS(modb_lifestage, paste0(path_out_brms, "modb_lifestage.rds"))
modb_lifestage <- readRDS(paste0(path_out_brms, "modb_lifestage.rds"))
# check all contrasts
hypothesis(
  modb_lifestage,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)
# Credible intervals contain 0

modb_lifestage_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                            Life.Stage +
                            (1 | Publication.Date.Authorship),
                          data = df,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          iter = 1e4,
                          cores = 4,
                          chains = 4)
summary(modb_lifestage_pub)
saveRDS(modb_lifestage_pub, paste0(path_out_brms, "modb_lifestage_pub.rds"))
modb_lifestage_pub <- readRDS(paste0(path_out_brms, "modb_lifestage_pub.rds"))
# check all contrasts
hypothesis(
  modb_lifestage_pub,
  c("Life.StageMetamorph = Life.StageJuvenile",
    "Life.StageMetamorph = Life.StageAdult",
    "Life.StageJuvenile = Life.StageAdult"),
  alpha = 0.05
)
# Credible intervals contain 0


# THERAPEUTIC OR PROPHYLATIC
modb_thera <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                    Therapeutic.or.Prophylactic,
                  data = df,
                  family = zero_inflated_beta(),
                  save_pars = save_pars(all = TRUE),
                  iter = 1e4,
                  cores = 4,
                  chains = 4)
summary(modb_thera)
saveRDS(modb_thera, paste0(path_out_brms, "modb_thera.rds"))
modb_thera <- readRDS(paste0(path_out_brms, "modb_thera.rds"))
# Credible interval does NOT contain 0

modb_thera_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                        Therapeutic.or.Prophylactic +
                        (1 | Publication.Date.Authorship),
                      data = df,
                      family = zero_inflated_beta(),
                      save_pars = save_pars(all = TRUE),
                      iter = 1e4,
                      cores = 4,
                      chains = 4)
summary(modb_thera_pub)
saveRDS(modb_thera_pub, paste0(path_out_brms, "modb_thera_pub.rds"))
modb_thera_pub <- readRDS(paste0(path_out_brms, "modb_thera_pub.rds"))
# Credible interval contains 0 but is close to +ve

modb_thera_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                           Therapeutic.or.Prophylactic +
                           (1 | Efficacy.Matrix),
                         data = df,
                         family = zero_inflated_beta(),
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.95),
                         iter = 1e4,
                         cores = 4,
                         chains = 4)
summary(modb_thera_effmat)
saveRDS(modb_thera_effmat, paste0(path_out_brms, "modb_thera_effmat.rds"))
modb_thera_effmat <- readRDS(paste0(path_out_brms, "modb_thera_effmat.rds"))
# Credible interval contains 0
table(df$Therapeutic.or.Prophylactic, df$Efficacy.Matrix)
# Therapeutic appears to better than prophylatcic because prevelance matrix
# has a bias towards more success and there is a larger proportion of
# therapeutic methods using the prevelance matrix then prophylactic using
# the prevelance matrix (91.5%). So the variation in success due to
# therapeutic vs prophylactic is actually explained by the
# variation in success with efficacy matrix type.


# IN SITU OR EX SITU
modb_situ <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                   In.situ.or.Ex.situ,
                 data = df,
                 family = zero_inflated_beta(),
                 save_pars = save_pars(all = TRUE),
                 iter = 1e4,
                 cores = 4,
                 chains = 4)
summary(modb_situ)
saveRDS(modb_situ, paste0(path_out_brms, "modb_situ.rds"))
modb_situ <- readRDS(paste0(path_out_brms, "modb_situ.rds"))
# Credible interval contains 0 and estimate is close to 0

modb_situ_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                       In.situ.or.Ex.situ +
                       (1 | Publication.Date.Authorship),
                     data = df,
                     family = zero_inflated_beta(),
                     save_pars = save_pars(all = TRUE),
                     iter = 1e4,
                     cores = 4,
                     chains = 4)
summary(modb_situ_pub)
saveRDS(modb_situ_pub, paste0(path_out_brms, "modb_situ_pub.rds"))
modb_situ_pub <- readRDS(paste0(path_out_brms, "modb_situ_pub.rds"))
# Credible interval contains 0 and estimate is close to 0


####################### SAVE PERTINENT OUTPUTS #################################


sink(file = paste0(path_out_brms, "brms_results.txt"))

cat("No amphibian-related factors had estimates of interest.")

cat("\n\nTreatment predictors\n\n")

cat("INTERVENTION CATEGORY: Some credible intervals do NOT contain 0\n\n")
print(summary(modb_intcat))
cat("\n\nINTERVENTION CATEGORY (PUBLICATION): Some credible intervals do NOT 
contain 0\n\n")
print(summary(modb_intcat_pub))
cat("\n\nINTERVENTION CATEGORY (EFFICACY MATRIX): Some credible intervals do 
NOT contain 0\n\n")
print(summary(modb_intcat_effmat))
cat("\nNo hypothesis tests had estimates of interest")

cat("\n\nTREATMENT TYPE: Some credible intervals do NOT contain 0\n\n")
print(summary(modb_treat))
cat("\n\nTREATMENT TYPE (PUBLICATION): Some credible intervals do NOT 
contain 0\n\n")
print(summary(modb_treat_pub))
cat("\n\nTREATMENT TYPE (EFFICACY MATRIX): Some credible intervals do NOT 
contain 0\n\n")
print(summary(modb_treat_effmat))
cat("\nNo hypothesis tests had estimates of interest")

cat("\n\nHABITAT/INDIVIDUAL: Credible interval contains 0 but some estimates 
are non-negligible and I-H lower CI is almost -ve\n\n")
print(summary(modb_habindv))
cat("\n\nHABITAT/INDIVIDUAL (PUBLICATION): I-H lower CI is -ve\n\n")
print(summary(modb_habindv_pub))
cat("\nHypothesis test did not have estimate of interest")

cat("\n\nTHERAPEUTIC/PROPHYLATIC: Credible interval does NOT contain 0.\n\n")
cat("However, therapeutic only appears better than prophylactic because 
prevelance matrix has a bias towards more success and there is a larger 
proportion of therapeutic methods using the prevelance matrix compared to  
prophylactic methods using the prevelance matrix (91.5%). 
So the variation in success due to therapeutic vs prophylactic is actually 
mostly explained by the variation in success with efficacy matrix type.
We see this when the previously significant credible interval changes to 
inlcude zero when publication/efficacy matrix is included in the model 
(note, publication and efficacy matrix are highly correlated)\n")

cat("\n\nNeither LIFE STAGE nor IN SITU/EX SITU had estimates of interest.")
sink()


## end of script
