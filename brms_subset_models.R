################################################################################
############################ Run subset models #################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Mar 2025


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
df$Life.Stage <- factor(df$Life.Stage, 
                                    levels = c("Larvae", "Metamorph", "Juvenile", "Adult"))   

# zero inflated beta family cannot have values equal to 1 so set 1 as nearly 1
df$Success_brms <- df$Success
df$Success_brms[df$Success_brms == 1] <- 0.9999999     

# increase max size to enable reloo option for large sizes in loo()
options(future.globals.maxSize=8000* 1024^3)

# create output folder for brms outputs
path_out_brms <- paste0(path, "brms/full_models/subsets/")
ifelse(!dir.exists(file.path(path_out_brms)), 
        dir.create(file.path(path_out_brms), recursive=T), 
        FALSE)


################################################################################
####################### THERAPEUTIC VS PROPHYLACTIC ############################


# modb_thera <- brm(Success_brms | weights(Uncertainty_weights) ~
#                         Therapeutic.or.Prophylactic,
#                 data = df,
#                 family=zero_inflated_beta(),
#                 save_pars = save_pars(all = TRUE),
#                 control = list(adapt_delta = 0.99),
#                 iter=1e4,
#                 cores=4,
#                 chains=4)
# summary(modb_thera)

# modb_thera_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
#                         Therapeutic.or.Prophylactic+
#                        (1 |Publication.Date.Authorship),
#                 data = df,
#                 family=zero_inflated_beta(),
#                 save_pars = save_pars(all = TRUE),
#                 control = list(adapt_delta = 0.99),
#                 iter=1e4,
#                 cores=4,
#                 chains=4)
# summary(modb_thera_pub)

# modb_thera_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~
#                         Therapeutic.or.Prophylactic+
#                         (1 | Efficacy.Matrix),
#                 data = df,
#                 family=zero_inflated_beta(),
#                 save_pars = save_pars(all = TRUE),
#                 control = list(adapt_delta = 0.99),
#                 iter=1e4,
#                 cores=4,
#                 chains=4)
# summary(modb_thera_eff)

# strong effect of therapeutic/prophylactic disapears in mixed models when
# publication (or efficiacy matrix) is controlled for


################################################################################
############################### ADULT ##########################################


df_adult <- df[which(df$Life.Stage == "Adult"),]
table(df_adult$Habitat.or.Individual)
table(df_adult$Intervention.category.itra.multi)
# remove single population demographic intervention
df_adult <- df_adult[-which(df_adult$Intervention.category.itra.multi == "Population demographic"),]
df_adult <- droplevels(df_adult)
table(df_adult$Intervention.category.itra.multi)

modb_adult_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                Intervention.category.itra.multi+
                Habitat.or.Individual+
                Therapeutic.or.Prophylactic+
                In.situ.or.Ex.situ+
                log1p(ClutchMn)+
                log1p(SVLMx)+
                Climate+
                Habitat+
                (1 |Publication.Date.Authorship),
                data = df_adult,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_adult_full)

modb_adult_pert <- brm(Success_brms | weights(Uncertainty_weights) ~
                Intervention.category.itra.multi+
                Habitat.or.Individual+
                In.situ.or.Ex.situ+
                (1 |Publication.Date.Authorship),
                data = df_adult,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_adult_pert)

# hypothesis tests to check differences between all factor levels
hypothesis(
  modb_adult_pert, 
  c("Habitat.or.IndividualH = Habitat.or.IndividualB"),
  alpha=0.05
)

hypothesis(
  modb_adult_pert, 
  c("Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiClimate",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiBioaugmentation = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiMultiple",
    "Intervention.category.itra.multiClimate = Intervention.category.itra.multiOtherchemical",
    "Intervention.category.itra.multiMultiple = Intervention.category.itra.multiOtherchemical"),
  alpha=0.05/6
)

# no predictors of interest in adult subset
# but note that sample size very small


################################################################################
############################### LARVAE #########################################


df_larva <- df[which(df$Life.Stage == "Larvae"),]
table(df_larva$Habitat.or.Individual)
# remove 'both' interventions because it is effecting model building
df_larva <- df_larva[-which(df_larva$Habitat.or.Individual == "B"),]
df_larva <- droplevels(df_larva)
table(df_larva$Habitat.or.Individual)
table(df_larva$Intervention.category.itra.multi)

modb_larva_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                Intervention.category.itra.multi+
                Habitat.or.Individual+
                Therapeutic.or.Prophylactic+
                In.situ.or.Ex.situ+
                log1p(ClutchMn)+
                log1p(SVLMx)+
                Climate+
                Habitat+
                (1 |Publication.Date.Authorship),
                data = df_larva,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_larva_full)

modb_larva_pert <- brm(Success_brms | weights(Uncertainty_weights) ~
                Intervention.category.itra.multi+
                Habitat.or.Individual+
                In.situ.or.Ex.situ+
                (1 |Publication.Date.Authorship),
                data = df_larva,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_larva_pert)
saveRDS(modb_larva_pert, paste0(path_out_brms, "../modb_larva.rds"))
modb_larva_pert <- readRDS(paste0(path_out_brms, "../modb_larva.rds"))

hypothesis(
  modb_larva_pert, 
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

# no effect of habitat but all interventions are worse than itraconazole
# with other chemical and population demographic not containing 0 in credible intervals

sink(file=paste0(path_out_brms, "../larva_model_results.txt"))
print(summary(modb_larva_pert))
cat("\n\n")
hypothesis(
  modb_larva_pert, 
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
sink()


################################################################################
############################### IN SITU ########################################


# df_insitu <- df[which(df$In.situ.or.Ex.situ == "In situ"),]
# table(df_insitu$Habitat.or.Individual)
# # remove 'both' interventions because it is effecting model building
# #df_insitu <- df_insitu[-which(df_insitu$Habitat.or.Individual == "B"),]
# df_insitu <- droplevels(df_insitu)
# #table(df_insitu$Habitat.or.Individual)
# table(df_insitu$Intervention.category.itra.multi)
# table(df_insitu$Life.Stage)


# modb_insitu <- brm(Success_brms | weights(Uncertainty_weights) ~
#                 Intervention.category.itra.multi+
#                 Habitat.or.Individual+
#                 Therapeutic.or.Prophylactic+
#                 Life.Stage+
#                 (1 |Publication.Date.Authorship),
#                 data = df_insitu,
#                 family=zero_inflated_beta(),
#                 save_pars = save_pars(all = TRUE),
#                 control = list(adapt_delta = 0.99),
#                 iter=1e4,
#                 cores=4,
#                 chains=4)
# summary(modb_insitu)

# DOES NOT CONVERGE (likely due to small size of dataset)


################################################################################
############################### HABITAT ########################################


# df_hab <- df[which(df$Habitat.or.Individual == "H"),]
# modb_hab <- brm(Success_brms | weights(Uncertainty_weights) ~
#                 Intervention.category.itra.multi+
#                 Life.Stage+
#                 Therapeutic.or.Prophylactic+
#                 In.situ.or.Ex.situ+
#                 (1 |Publication.Date.Authorship),
#                 data = df_hab,
#                 family=zero_inflated_beta(),
#                 save_pars = save_pars(all = TRUE),
#                 control = list(adapt_delta = 0.99),
#                 iter=2e5,
#                 cores=4,
#                 chains=4)
# summary(modb_hab_full)

# DOES NOT CONVERGE


################################################################################
################################# CHEMICAL #####################################


# subset data by chemical only
chem <- df[which((df$Intervention.category.1 == "Chemical") | (df$Intervention.category.2 == "Chemical")),]

table(chem$Habitat.or.Individual)
# remove both interventions to help reduce divergent transitions
chem <- chem[chem$Habitat.or.Individual != "B",]
table(chem$MultipleTreatments)
table(chem$Life.Stage)
table(chem$Therapeutic.or.Prophylactic)
table(chem$In.situ.or.Ex.situ) 

chem$Treatment.used.multi <- chem$Treatment.used.1
for (row in 1:nrow(chem)) {
  if (chem$Intervention.category.1.itra[row] == "Itraconazole") {
      chem$Treatment.used.multi[row] <- "Itraconazole"
  }
  if (chem$MultipleTreatments[row] == "Yes") {
    chem$Treatment.used.multi[row] <- "Multiple"
  }
}

table(chem$Treatment.used.multi)
# remove single antiparasitic intervention
chem <- chem[chem$Treatment.used.multi != "Antiparisitic",]
table(chem$Treatment.used.multi)
chem <- droplevels(chem)
# reset factor levels to compare against itraconazole
chem$Treatment.used.multi <- factor(chem$Treatment.used.multi, 
                                    levels = c(
                                      "Itraconazole",
                                      "Antibiotic", 
                                      "Antifungal", 
                                      "Disinfectant", 
                                      "Herbicide", 
                                      "Insecticide",
                                      "Multiple",
                                      "Other",
                                      "Sodium"))     


########################## Full model - no controls ############################


# run full model with no controls
modb_chem_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                Treatment.used.multi+
                Habitat.or.Individual+
                Life.Stage+
                Therapeutic.or.Prophylactic+
                In.situ.or.Ex.situ,
                data = chem,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_chem_full)

# Adult is worse than larvae


############################### Full model #####################################


# run full model with no controls
modb_chem_full_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                Treatment.used.multi+
                Habitat.or.Individual+
                Life.Stage+
                Therapeutic.or.Prophylactic+
                In.situ.or.Ex.situ+
                (1 | Publication.Date.Authorship),
                data = chem,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_chem_full_pub)

# Note pertinent associations (note < 10 divergent transitions)
# No pertinent variables


################################################################################
################################ ITRACONAZOLE ##################################


# subset data by chemical only
itra <- df[which((df$Specific.treatment.used.1 == "Itraconazole") | (df$Specific.treatment.used.2 == "Itraconazole")),]

table(itra$Habitat.or.Individual)
# not enough habitat or both to use this preidtcor in model
table(itra$MultipleTreatments)
table(itra$Life.Stage)
# remove single metamorph intervention
itra <- itra[itra$Life.Stage != "Metamorph",]
table(itra$Life.Stage)
table(itra$Therapeutic.or.Prophylactic)
table(itra$In.situ.or.Ex.situ)
itra <- droplevels(itra)


################################# Dosage #######################################


table(itra$Dosage)

convert_dosage <- function(itra){
  itra$Dosage[itra$Dosage == "100 mg l−1 for the first 3 d, 5 mg l−1 for 6 d,  50 mg l−1 for the last day"] <- NA
  itra$Dosage <- gsub(" mg l -1", "", itra$Dosage)
  itra$Dosage <- gsub(" mg l–1", "", itra$Dosage)
  itra$Dosage <- gsub(" mg l−1", "", itra$Dosage)
  itra$Dosage <- gsub(" mg/L", "", itra$Dosage)

  for (i in 1:nrow(itra)) {
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

modb_itra_dosage <- brm(Success_brms | weights(Uncertainty_weights) ~
                Dosage,
                data = itra,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_dosage)

modb_itra_dosage_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                Dosage+
                (1 | Publication.Date.Authorship),
                data = itra,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_dosage_pub)

# Issues with starting values
# Dosage has no assoiciation with success


########################## Full model - no controls ############################


# run full model with no controls
modb_itra_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                MultipleTreatments+
                Life.Stage+
                Therapeutic.or.Prophylactic+
                In.situ.or.Ex.situ,
                data = itra,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_full)

# Adult is worse than larvae


############################### Full model #####################################


# run full model controlling for publication
modb_itra_full_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       MultipleTreatments+
                       Life.Stage+
                       Therapeutic.or.Prophylactic+
                       In.situ.or.Ex.situ+
                       (1 | Publication.Date.Authorship),
                data = itra,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_full_pub)

# No pertinent variables but adult is close to performing worse than larvae


################################################################################
########################## WITHOUT ITRACONAZOLE ################################


# remove itraconazole
df_sub <- df[-which(df$Intervention.category.itra.multi == "Itraconazole"),]
table(df_sub$Intervention.category.itra.multi)
df_sub <- droplevels(df_sub)

# reset factors now itraconazole has been removed
df_sub$Intervention.category.itra.multi <- factor(df_sub$Intervention.category.itra.multi, 
                                    levels = c(
                                      "Population demographic",                                      
                                      "Bioaugmentation", 
                                      "Climate",
                                      "Multiple",
                                      "Other chemical"))        

# run full model controlling for efficacy matrix
modb_sub_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Intervention.category.itra.multi+
                       Habitat.or.Individual+
                       Life.Stage+
                       (1 | Publication.Date.Authorship),
                data = df_sub[!(is.na(df_sub$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_sub_pub)

# No associations of interest (habitat/individual association disappears)


################################################################################
#################### BINARY ITRACONAZOLE OR NOT CATEGORY #######################


# set new column with 0/1 depending on if treatment is itraconazole or not
df$Itra <- ifelse(df$Intervention.category.itra.multi == "Itraconazole", 1, 0)
df$Intervention.category.itra.multi

modb_full_itra_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Itra+
                       Habitat.or.Individual+
                       Life.Stage+
                       (1 | Publication.Date.Authorship),
                data = df[!(is.na(df$SVLMx)),],
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_full_itra_pub)

# No associations of interest (habitat/individual and adult associations disappear
# but are borderline)


################################################################################
################ ITRACONAZOLE IN SINGLE AND MULTIPLE TREATMENTS ################


# make a new 1/0 column for if itraconazole has been used at all in intervention 
# (including in interventions with multiple treatment types)
df$Itra_all <- ifelse((!is.na(df$Specific.treatment.used.1) & df$Specific.treatment.used.1 == "Itraconazole") | 
                      (!is.na(df$Specific.treatment.used.2) & df$Specific.treatment.used.2 == "Itraconazole"), 1, 0)

# run model with this as treatment type predictor
modb_itra_all_pub <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Itra_all+
                       Habitat.or.Individual+
                       Life.Stage+
                       (1 | Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_all_pub)
loo_itra_full_pub <- loo(modb_itra_all_pub, moment_match=TRUE, reloo=TRUE)

# No associations of interest (habitat/individual and adult associations disappear
# but are borderline)


############### Compare itraconazole only with habitat only ####################


# try and tease apart effect of individual/habitat from that of itraconazole 
# (since most itraconazole treatments are individual so there is confoudning)

# run model on just itraconazole
modb_itra_all <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Itra_all+
                       (1 | Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_all)
loo_itra_all <- loo(modb_itra_all, moment_match=TRUE, reloo=TRUE)
# itraconazole is positively associated with success

# run model on itracoazole and other interventions predictors but 
# without habitat/individual
modb_itra_nohabindv_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Life.Stage+
                       Itra_all+
                       (1 | Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_itra_nohabindv_full)
loo_itra_nohabindv <- loo(modb_itra_nohabindv_full, moment_match=TRUE, reloo=TRUE)
# itraconazole is positively associated with success but no effect of life stage

# run model on just habitat/individual
modb_hab <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Habitat.or.Individual+
                       (1 | Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_hab)
loo_hab <- loo(modb_hab, moment_match=TRUE, reloo=TRUE)
# individual is positively associated with success

# run model on habitat/individual and other interventions predictors but 
# without itraconazole (or any treatment type predictor)
modb_noitra_habindv_full <- brm(Success_brms | weights(Uncertainty_weights) ~
                       Life.Stage+
                       Habitat.or.Individual+
                       (1 | Publication.Date.Authorship),
                data = df,
                family=zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99),
                iter=1e4,
                cores=4,
                chains=4)
summary(modb_noitra_habindv_full)
loo_habindv_noitra <- loo(modb_noitra_habindv_full, moment_match=TRUE, reloo=TRUE)
# itraconazole is positively associated with success but no effect of life stage


############################### Compare models #################################


loo_compare(
  loo_itra_all,
  loo_itra_nohabindv,
  loo_hab, 
  loo_habindv_noitra,
  loo_itra_full_pub
)
# worst model is modb_itra_nohabindv_full (itraconazole and lifestage 
# without habitat/individual). But all others have ELPD difference within std err.


## end of script
