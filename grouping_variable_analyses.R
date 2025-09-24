################################################################################
########################## Grouping variable analyses ##########################
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
path <- paste0("~/Documents/scripts/2025_bd_interventions/outputs/", output_options)

# load data
df <- as.data.frame(read.csv(paste0(path, "success_df.csv")))

# load packages
library("tidyverse")
library("brms")
# load custom plotting functions
source(paste0(path, "../../functions_for_plotting.R"))


# zero inflated beta family cannot have values equal to 1 so set 1 as nearly 1
df$Success_brms <- df$Success
df$Success_brms[df$Success_brms == 1] <- 0.9999999

# increase max size to enable reloo option for large sizes in loo()
options(future.globals.maxSize = 8000 * 1024^3)

# create output folder for brms outputs
path_out_brms_gv <- paste0(path, "brms/grouping_variables/")
ifelse(!dir.exists(file.path(path_out_brms_gv)),
       dir.create(file.path(path_out_brms_gv), recursive = TRUE),
       FALSE)


################################################################################
############################## Efficacy matrix #################################


# split dataset by intervention category
table(df$Efficacy.Matrix)

# infection matrix
infm <- df[df$Efficacy.Matrix == "Infection matrix", ]
mean(infm$Success)
table(infm$Intervention.category.1)
table(infm$Specific.treatment.used.1[infm$Intervention.category.1 == "Chemical"])

# prevelance matrix
prevm <- df[df$Efficacy.Matrix == "Prevelance matrix", ]
mean(prevm$Success)
table(prevm$Intervention.category.1)
table(prevm$Specific.treatment.used.1[prevm$Intervention.category.1 == "Chemical"])

# survivability matrix
survm <- df[df$Efficacy.Matrix == "Survivability matrix", ]
mean(survm$Success)
table(survm$Intervention.category.1)
table(survm$Specific.treatment.used.1[survm$Intervention.category.1 == "Chemical"])

# make vioin plot
eff_labs <- makeLabs(df, "Efficacy.Matrix")
eff_violin <- plotViolin(df, "Efficacy.Matrix", eff_labs, path_out_brms_gv)
eff_violin


################################################################################
############################### Bayesian models ################################


############################## PUBLICATION #####################################


modb_pub <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                  (1 | Publication.Date.Authorship),
                data = df[!(is.na(df$SVLMx)), ],
                family = zero_inflated_beta(),
                save_pars = save_pars(all = TRUE),
                iter = 1e4,
                cores = 4,
                chains = 4)
summary(modb_pub)
saveRDS(modb_pub, paste0(path_out_brms_gv, "modb_pub.rds"))
modb_pub <- readRDS(paste0(path_out_brms_gv, "modb_pub.rds"))
# diagnsotics
#pp_check(modb_pub)
#plot(modb_pub)
loo_pub <- loo(modb_pub, moment_match = TRUE, reloo = TRUE)
loo_pub
saveRDS(loo_pub, paste0(path_out_brms_gv, "loo_pub.rds"))
loo_pub <- readRDS(paste0(path_out_brms_gv, "loo_pub.rds"))


########################### EFFICACY MATRIX ####################################


modb_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                     (1 | Efficacy.Matrix),
                   data = df[!(is.na(df$SVLMx)), ],
                   family = zero_inflated_beta(),
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99),
                   iter = 1e5,
                   cores = 4,
                   chains = 4)
summary(modb_effmat)
saveRDS(modb_effmat, paste0(path_out_brms_gv, "modb_effmat.rds"))
modb_effmat <- readRDS(paste0(path_out_brms_gv, "modb_effmat.rds"))
# diagnostics
#pp_check(modb_effmat)
#plot(modb_effmat)
loo_effmat <- loo(modb_effmat, moment_match = TRUE, reloo = TRUE)
loo_effmat
saveRDS(loo_effmat, paste0(path_out_brms_gv, "loo_effmat.rds"))
loo_effmat <- readRDS(paste0(path_out_brms_gv, "loo_effmat.rds"))


################################# GENUS ########################################


modb_genus <- brm(Success_brms | weights(Uncertainty_weights) ~ 1 +
                    (1 | Genus),
                  data = df[!(is.na(df$SVLMx)), ],
                  family = zero_inflated_beta(),
                  save_pars = save_pars(all = TRUE),
                  iter = 1e4,
                  cores = 4,
                  chains = 4)
summary(modb_genus)
saveRDS(modb_genus, paste0(path_out_brms_gv, "modb_genus.rds"))
modb_genus <- readRDS(paste0(path_out_brms_gv, "modb_genus.rds"))
# diagnostics
#pp_check(modb_genus)
#plot(modb_genus)
loo_genus <- loo(modb_genus, moment_match = TRUE, reloo = TRUE)
loo_genus
saveRDS(loo_genus, paste0(path_out_brms_gv, "loo_genus.rds"))
loo_genus <- readRDS(paste0(path_out_brms_gv, "loo_genus.rds"))


################## PUBLICATION AND EFFICACY AND GENUS ##########################


modb_pub_effmat_genus <- brm(Success_brms | weights(Uncertainty_weights) ~
                               (1 | Publication.Date.Authorship) +
                                 (1 | Efficacy.Matrix) +
                                 (1 | Genus),
                             data = df[!(is.na(df$SVLMx)), ],
                             family = zero_inflated_beta(),
                             save_pars = save_pars(all = TRUE),
                             control = list(adapt_delta = 0.95),
                             iter = 1e4,
                             cores = 4,
                             chains = 4)
summary(modb_pub_effmat_genus)
saveRDS(modb_pub_effmat_genus, paste0(path_out_brms_gv, "modb_effmat_pub_genus.rds"))
modb_pub_effmat_genus <- readRDS(paste0(path_out_brms_gv, "modb_effmat_pub_genus.rds"))
# diagnostics
#pp_check(modb_pub_effmat_genus)
#plot(modb_pub_effmat_genus)
loo_pub_effmat_genus <- loo(modb_pub_effmat_genus, moment_match = TRUE, reloo = TRUE)
loo_pub_effmat_genus
saveRDS(loo_pub_effmat_genus, paste0(path_out_brms_gv, "loo_pub_effmat_genus.rds"))
loo_pub_effmat_genus <- readRDS(paste0(path_out_brms_gv, "loo_pub_effmat_genus.rds"))


####################### PUBLICATION AND EFFICACY ###############################


modb_pub_effmat <- brm(Success_brms | weights(Uncertainty_weights) ~
                         (1 | Publication.Date.Authorship) +
                           (1 | Efficacy.Matrix),
                       data = df[!(is.na(df$SVLMx)), ],
                       family = zero_inflated_beta(),
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.95),
                       iter = 1e4,
                       cores = 4,
                       chains = 4)
summary(modb_pub_effmat)
saveRDS(modb_pub_effmat, paste0(path_out_brms_gv, "modb_effmat_pub.rds"))
modb_pub_effmat <- readRDS(paste0(path_out_brms_gv, "modb_effmat_pub.rds"))
# diagnostics
#pp_check(modb_pub_effmat)
#plot(modb_pub_effmat)
loo_pub_effmat <- loo(modb_pub_effmat, moment_match = TRUE, reloo = TRUE)
loo_pub_effmat
saveRDS(loo_pub_effmat, paste0(path_out_brms_gv, "loo_pub_effmat.rds"))
loo_pub_effmat <- readRDS(paste0(path_out_brms_gv, "loo_pub_effmat.rds"))


################################################################################
############################## Comparison ######################################


# compare grouping models with loo
loo_compare(
  loo_pub,
  loo_effmat,
  loo_genus,
  loo_pub_effmat,
  loo_pub_effmat_genus
)

# compare grouping models with r2
models <- list(modb_pub_effmat, modb_pub_effmat_genus, modb_effmat, modb_pub, modb_genus)
names(models) <- c("pub_effmat", "pub_effmat_genus", "effmat", "pub", "genus")

# find r2 distributions
r2_comparison <- sapply(models, bayes_R2)

# print results in a friendly format
r2_results <- data.frame(
  R2_mean = r2_comparison[1, ],
  R2_sd = r2_comparison[2, ],
  R2_lower = r2_comparison[3, ],
  R2_upper = r2_comparison[4, ]
)
print(r2_results)


# loo and r2 show that genus doesn't affect success, and publication and
# efficacy matrix are probably capturing the same signal (correlated and have
# similar explanatory power). The combined model is likely suffering from
# overfitting because of this. Publication has more levels and is more
# likely to explain more of the variation so use only publication for
# random effect.


## end of script
