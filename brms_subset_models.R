################################################################################
############################ Run subset models #################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Oct 2025


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
path_out_brms <- paste0(path, "brms/full_models/subsets/")
ifelse(!dir.exists(file.path(path_out_brms)),
       dir.create(file.path(path_out_brms), recursive = T),
       FALSE)


################################################################################
############################### ADULT ##########################################


df_adult <- df[which(df$Life.Stage == "Adult"), ]
table(df_adult$Habitat.or.Individual)
table(df_adult$Intervention.category.itra.multi)
# remove single population demographic intervention
df_adult <- df_adult[-which(df_adult$Intervention.category.itra.multi == "Population demographic"), ]
df_adult <- droplevels(df_adult)
table(df_adult$Intervention.category.itra.multi)
table(df_adult$TreatmentType)
# group sizes are already small so don't split further into treatment type

modb_adult <- brm(Success_brms | weights(Uncertainty_weights) ~
                    Intervention.category.itra.multi +
                      Habitat.or.Individual +
                      Therapeutic.or.Prophylactic +
                      In.situ.or.Ex.situ +
                      (1 | Publication.Date.Authorship),
                  data = df_adult,
                  family = zero_inflated_beta(),
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.9999),
                  iter = 1e4,
                  cores = 4,
                  chains = 4)
summary(modb_adult)

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_adult,
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha = 0.05
)

hypothesis(
  modb_adult,
  c("Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical"),
  alpha = 0.05
)
# no predictors of interest in adult subset
# there is an inclination for itraconazole/individual to perform better
# (unequal credible intervals) but intervals are all quite large and contain 0
# (note that sample size very small)


################################################################################
############################### LARVAE #########################################


df_larva <- df[which(df$Life.Stage == "Larvae"), ]
table(df_larva$Habitat.or.Individual)
# remove 'both' interventions because it is effecting model building
df_larva <- df_larva[-which(df_larva$Habitat.or.Individual == "B"), ]
df_larva <- droplevels(df_larva)
table(df_larva$Habitat.or.Individual)
table(df_larva$Intervention.category.itra.multi)
table(df_adult$TreatmentType)
# group sizes are already small so don't split further into treatment type

# full model
modb_larva <- brm(Success_brms | weights(Uncertainty_weights) ~
                    Intervention.category.itra.multi +
                      Habitat.or.Individual +
                      Therapeutic.or.Prophylactic +
                      In.situ.or.Ex.situ +
                      (1 | Publication.Date.Authorship),
                  data = df_larva,
                  family = zero_inflated_beta(),
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99),
                  iter = 1e4,
                  cores = 4,
                  chains = 4)
summary(modb_larva)

# run with pertinent variables only to compare with best model in post hoc
modb_larva_pert <- brm(Success_brms | weights(Uncertainty_weights) ~
                    Intervention.category.itra.multi +
                      Habitat.or.Individual +
                      (1 | Publication.Date.Authorship),
                  data = df_larva,
                  family = zero_inflated_beta(),
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99),
                  iter = 1e4,
                  cores = 4,
                  chains = 4)
summary(modb_larva_pert)
saveRDS(modb_larva_pert, paste0(path_out_brms, "modb_larva_only.rds"))
modb_larva_pert <- readRDS(paste0(path_out_brms, "modb_larva_only.rds"))

hypothesis(
  modb_larva,
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
# no effect of habitat but all interventions are inclining towards worse than
# itraconazole with population demographic and other chemical not containing 0
# in credible interval


################################################################################
############################### IN SITU ########################################


# df_insitu <- df[which(df$In.situ.or.Ex.situ == "In situ"), ]
# table(df_insitu$Habitat.or.Individual)
# # remove 'both' interventions because it is effecting model building
# df_insitu <- df_insitu[-which(df_insitu$Habitat.or.Individual == "B"), ]
# df_insitu <- droplevels(df_insitu)
# table(df_insitu$Habitat.or.Individual)
# table(df_insitu$Intervention.category.itra.multi)
# table(df_insitu$Life.Stage)

# modb_insitu <- brm(Success_brms | weights(Uncertainty_weights) ~
#                      Intervention.category.itra.multi +
#                      Habitat.or.Individual +
#                      Life.Stage +
#                      (1 | Publication.Date.Authorship),
#                    data = df_insitu,
#                    family = zero_inflated_beta(),
#                    save_pars = save_pars(all = TRUE),
#                    control = list(adapt_delta = 0.999999),
#                    iter = 1e5,
#                    cores = 4,
#                    chains = 4)
# summary(modb_insitu)
# DOES NOT CONVERGE (likely due to small size of dataset)


################################################################################
############################### HABITAT ########################################


df_hab <- df[which(df$Habitat.or.Individual == "H"), ]
table(df_hab$Intervention.category.itra.multi)
df_hab <- df_hab[which(!(df_hab$Intervention.category.itra.multi %in% c("Itraconazole", "Population demographic"))), ]
table(df_hab$Intervention.category.itra.multi)
table(df_hab$Life.Stage)
table(df_hab$In.situ.or.Ex.situ)
# too few in situ so don't include in model
df_hab <- droplevels(df_hab)
modb_hab <- brm(Success_brms | weights(Uncertainty_weights) ~
                  Intervention.category.itra.multi +
                    Life.Stage +
                    (1 | Publication.Date.Authorship),
                data = df_hab,
                family = zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999),
                iter = 2e4,
                cores = 4,
                chains = 4)
summary(modb_hab)
# no evidence that intervention category or life stage impact success for
# habitat interventions


################################################################################
################################# CHEMICAL #####################################


# subset data by chemical only
chem <- df[which((df$Intervention.category.1 == "Chemical") | (df$Intervention.category.2 == "Chemical")), ]
chem$Specific.treatment.used.1 <- trimws(chem$Specific.treatment.used.1)
chem$Specific.treatment.used.2 <- trimws(chem$Specific.treatment.used.2)

table(chem$Habitat.or.Individual)
# remove both interventions to help reduce divergent transitions
chem <- chem[chem$Habitat.or.Individual != "B", ]
chem <- droplevels(chem)
table(chem$MultipleTreatments)
table(chem$Life.Stage)
table(chem$Therapeutic.or.Prophylactic)
table(chem$In.situ.or.Ex.situ)
table(chem$TreatmentType)

# remove single antiparasitic intervention
chemsub <- chem[chem$TreatmentType != "Medicinal antiparasitic", ]
# remove two BMP-NTf2 interventions
chemsub <- chemsub[chemsub$TreatmentType != "BMP-NTf2", ]
# remove general tonic interventions as they are causing issues in the model
chemsub <- chemsub[chemsub$TreatmentType != "General Tonic", ]
chemsub$TreatmentType <- as.factor(chemsub$TreatmentType)
chemsub <- droplevels(chemsub)
table(chemsub$TreatmentType)

# run full model
modb_chem_full_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                            TreatmentType +
                              Habitat.or.Individual +
                              Life.Stage +
                              Therapeutic.or.Prophylactic +
                              In.situ.or.Ex.situ +
                              (1 | Publication.Date.Authorship),
                          data = chemsub,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.99),
                          iter = 1e4,
                          cores = 4,
                          chains = 4)
summary(modb_chem_full_pub)
# no strong evidence of any predictors influencing success
# however, the estimates for habitat, metamorph, juvenile and adult
# are all noticeably negatively skewed

# run with pertinent variables only to compare with best model in post hoc
modb_chem_pert_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                            TreatmentType +
                              Habitat.or.Individual +
                              Life.Stage +
                              (1 | Publication.Date.Authorship),
                          data = chemsub,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.99),
                          iter = 1e4,
                          cores = 4,
                          chains = 4)
summary(modb_chem_pert_pub)
saveRDS(modb_chem_pert_pub, paste0(path_out_brms, "modb_chem_only.rds"))
modb_chem_pert_pub <- readRDS(paste0(path_out_brms, "modb_chem_only.rds"))
# no strong evidence of any predictors influencing success
# however, the estimates for habitat, metamorph, juvenile and adult
# are all noticeably negatively skewed


################################################################################
################################ ITRACONAZOLE ##################################


# subset data by chemical only
itra <- df[which((df$Specific.treatment.used.1 == "Itraconazole") | (df$Specific.treatment.used.2 == "Itraconazole")), ]

table(itra$Habitat.or.Individual)
# not enough habitat or both to use this predictor in model
table(itra$MultipleTreatments)
table(itra$Life.Stage)
# remove single metamorph intervention
itra <- itra[itra$Life.Stage != "Metamorph", ]
table(itra$Life.Stage)
table(itra$Therapeutic.or.Prophylactic)
table(itra$In.situ.or.Ex.situ)
itra <- droplevels(itra)


################################# Dosage #######################################


table(itra$Dosage)
# function to convert dosages to have the same units
convert_dosage <- function(itra) {
  itra$Dosage[itra$Dosage == "100 mg l−1 for the first 3 d, 5 mg l−1 for 6 d,  50 mg l−1 for the last day"] <- NA
  itra$Dosage <- gsub(" mg l -1", "", itra$Dosage)
  itra$Dosage <- gsub(" mg l–1", "", itra$Dosage)
  itra$Dosage <- gsub(" mg l−1", "", itra$Dosage)
  itra$Dosage <- gsub(" mg/L", "", itra$Dosage)

  for (i in seq_len(nrow(itra))) {
    if ((grepl("Itrafungol", itra$Dosage[i])) ||
          (grepl("Sporanox", itra$Dosage[i])) ||
          (grepl("vikron", itra$Dosage[i]))) {
      itra$Dosage[i] <- NA
    } else if (grepl("%", itra$Dosage[i])) {
      temp <- as.numeric(gsub("%", "", itra$Dosage[i]))
      itra$Dosage[i] <- temp * 10000
    } else if (grepl("μg/mL", itra$Dosage[i])) {
      temp <- as.numeric(gsub(" μg/mL", "", itra$Dosage[i]))
      itra$Dosage[i] <- temp / 1000
    }
  }
  return(itra)
}

itra <- convert_dosage(itra)
table(itra$Dosage)
itra$Dosage <- as.numeric(itra$Dosage)

# run model with just dosage and publication
modb_itra_dosage_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                              Dosage +
                                (1 | Publication.Date.Authorship),
                            data = itra,
                            family = zero_inflated_beta(),
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.99),
                            iter = 1e4,
                            cores = 4,
                            chains = 4)
summary(modb_itra_dosage_pub)
# dosage has no bivariate assoiciation with success
# estimate is 0 with very narrow credible intervals


############################### Full model #####################################


# run full model controlling for publication
modb_itra_full_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                            MultipleTreatments +
                              Life.Stage +
                              Therapeutic.or.Prophylactic +
                              In.situ.or.Ex.situ +
                              Dosage +
                              (1 | Publication.Date.Authorship),
                          data = itra,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.999999),
                          iter = 2e4,
                          cores = 4,
                          chains = 4)
summary(modb_itra_full_pub)
# no evidence of associations
# but note a few divergent transitions and large CIs


################################################################################
########################## WITHOUT ITRACONAZOLE ################################


# remove itraconazole
df_sub <- df[-which(df$Intervention.category.itra.multi == "Itraconazole"), ]
table(df_sub$Intervention.category.itra.multi)
df_sub <- droplevels(df_sub)
# reset factors now itraconazole has been removed
df_sub$Intervention.category.itra.multi <- factor(
  df_sub$Intervention.category.itra.multi,
  levels = c(
    "Population demographic",
    "Bioaugmentation",
    "Climate",
    "Multiple",
    "Other chemical"
  )
)
table(df_sub$Habitat.or.Individual)
table(df_sub$Life.Stage)
table(df_sub$In.situ.or.Ex.situ)
table(df_sub$Therapeutic.or.Prophylactic)

# run intervention category model controlling for publication
modb_subintcat_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                            Intervention.category.itra.multi +
                              Habitat.or.Individual +
                              Life.Stage +
                              Therapeutic.or.Prophylactic +
                              In.situ.or.Ex.situ +
                              (1 | Publication.Date.Authorship),
                          data = df_sub,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.99),
                          iter = 1e4,
                          cores = 4,
                          chains = 4)
summary(modb_subintcat_pub)
# no strong evidence but habitat and both are negatively skewed with
# respect to individual and intervention categories are slightly
# positively skewed with respect to population demographic (ref level)
# but all have 0 in the credible interval and some have large
# CIs. So trends match main model, which means the pattern that habitat may
# be worse than individual and population demographic performs worse is
# still (faintly) present

# run with pertinent variables only to compare with best model in post hoc
modb_subintcat_pert_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                            Intervention.category.itra.multi +
                              Habitat.or.Individual +
                              Life.Stage +
                              (1 | Publication.Date.Authorship),
                          data = df_sub,
                          family = zero_inflated_beta(),
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.99),
                          iter = 1e4,
                          cores = 4,
                          chains = 4)
summary(modb_subintcat_pert_pub)
saveRDS(modb_subintcat_pert_pub, paste0(path_out_brms, "modb_without_itra.rds"))
modb_subintcat_pert_pub <- readRDS(paste0(path_out_brms, "modb_without_itra.rds"))
# no strong evidence but habitat is are negatively skewed with
# respect to individual and intervention categories are slightly
# positively skewed with respect to population demographic (ref level)
# but all have 0 in the credible interval and some have large
# CIs. So trends match main model, which means the pattern that habitat may
# be worse than individual and population demographic performs worse is
# still (faintly) present

# model with treatment type
table(df_sub$TreatmentType)
# remove medicinal antiparasitic (only one observation)
# remove general tonic (causes divergent transitions)
df_sub <- df_sub[which(!(
  df_sub$TreatmentType %in% c(
    "Medicinal antiparasitic",
    "BMP-NTf2",
    "General Tonic"
  )
)), ]
df_sub <- droplevels(df_sub)
# run full model controlling for publication
modb_subt_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       TreatmentType +
                         Habitat.or.Individual +
                         Life.Stage +
                         Therapeutic.or.Prophylactic +
                         In.situ.or.Ex.situ +
                         (1 | Publication.Date.Authorship),
                     data = df_sub,
                     family = zero_inflated_beta(),
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.99),
                     iter = 1e4,
                     cores = 4,
                     chains = 4)
summary(modb_subt_pub)
# treatment type model adds no further insight to intervention model


################################################################################
#################### BINARY ITRACONAZOLE OR NOT CATEGORY #######################


# set new column with 0/1 depending on if treatment is itraconazole or not
df$Itra <- ifelse(df$Intervention.category.itra.multi == "Itraconazole", 1, 0)

modb_full_itra_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                            Itra +
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
summary(modb_full_itra_pub)
# habitat still worse than individual
# life stages are all negatively skewed compared to larvae
# binary itraconazole is positively skewed but 0 is still in credible interval


################################################################################
############## ITRACONAZOLE IN BOTH SINGLE AND MULTIPLE TREATMENTS #############


# make a new 1/0 column for if itraconazole has been used at all in
# interventions (including in interventions with multiple treatment types)
df$Itra_all <- ifelse((!is.na(df$Specific.treatment.used.1) & df$Specific.treatment.used.1 == "Itraconazole") | 
                      (!is.na(df$Specific.treatment.used.2) & df$Specific.treatment.used.2 == "Itraconazole"), 1, 0)

# run model with this as treatment type predictor
modb_itra_all_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                           Itra_all +
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
summary(modb_itra_all_pub)
# same as above
# habitat still worse than individual
# life stages are all negatively skewed compared to larvae
# binary itraconazole is positively skewed but 0 is still in credible interval


## end of script
