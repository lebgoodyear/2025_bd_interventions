################################################################################
######################## Calculate success and weights #########################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Sept 2025


# clear workspace
rm(list = ls())


################################################################################
################################# Prep #########################################


# load required packages
library("roam")
library("ggplot2")

# load data
path <- "~/Documents/scripts/2025_bd_interventions/"
df <- read.csv(paste0(path, "data/interventions_gabip_dataset.csv"))


############################## Data wrangling ##################################


# create binary longevity variable based on long-term efficacy
table(df$Long.term.monitoring)
table(df$Trend.in.efficacy.score)
# set all positive ratings to 1 and negatives to 0
df$Trend <- NA
df$Monitoring <- NA
df$Longevity <- NA
for (r in seq_len(nrow(df))) {
  if (df$Long.term.monitoring[r] == "Y") {
    df$Monitoring[r] <- 1
    if (df$Trend.in.efficacy.score[r] %in%
          c("Increase in score", "No change")) {
      df$Trend[r] <- 1
    } else {
      df$Trend[r] <- 0
    }
  } else {
    df$Monitoring[r] <- 0
    df$Trend[r] <- 0
  }
  df$Longevity[r] <- df$Monitoring[r] * df$Trend[r]
}
table(df$Longevity)

# make scalability numeric
table(df$Scalability)
df$Scalability[df$Scalability == "N"] <- 0
df$Scalability[df$Scalability == "Y"] <- 1
df$Scalability <- as.numeric(df$Scalability)
table(df$Scalability)

# remove any interventions with NA efficacy
df <- df[which(!is.na(df$Efficacy.of.intervention)), ]
# scale efficacy to between 0 and 1
table(df$Efficacy.of.intervention)
df$Efficacy <- df$Efficacy.of.intervention / max(df$Efficacy.of.intervention)
table(df$Efficacy)

# scale adverse effects to between 0 and 1, where 0 is worst
table(df$Adverse.Effects)
df$Adverse <- 1 - (df$Adverse.Effects / max(df$Adverse.Effects))
table(df$Adverse)

# check usability
table(df$Usability)
# remove any interventions with NA usability
df <- df[which(!is.na(df$Usability)), ]
# remove all interventions with a usability of 3
df <- df[-which(df$Usability == 3), ]

# remove leading and trailing whitespace
df$Specific.treatment.used.1 <- trimws(df$Specific.treatment.used.1)
df$Specific.treatment.used.2 <- trimws(df$Specific.treatment.used.2)
# group chemical treatments by use to reduce last group size of "Other Chemical"
df$TreatmentType <- as.character(df$Intervention.category.1.itra)

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Antibiotic mix",
                                                     "Florfenicol",
                                                     "Chloramphenicol",
                                                     "Trimethoprim-sulfadiazine")] <- "Medicinal antibiotic"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Itraconazole",
                                                     "Voriconazole",
                                                     "Terbinafine hydrochloride",
                                                     "Miconazole",
                                                     "Amphotericin B",
                                                     "Fluconazole")] <- "Medicinal antifungal"

df$TreatmentType[df$Specific.treatment.used.1 == "Itraconazole"] <- "Itraconazole"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Ivermectin",
                                                     "Malachite green")] <- "Medicinal antiparasitic"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Mandipropamid",
                                                     "Thiophanate-methyl",
                                                     "Chlorothalonil")] <- "Fungicide"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Atrazine",
                                                     "Herbicide mix",
                                                     "Glyphosate")] <- "Herbicide"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Insecticide mix",
                                                     "Malathion",
                                                     "Carbaryl")] <- "Insecticide"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Virkon Aquatic®",
                                                     "Virkon S",
                                                     "Steriplant N",
                                                     "F10SC",
                                                     "Bleach",
                                                     "Benzalkonium chloride",
                                                     "Hydrogen peroxide",
                                                     "Ethanol",
                                                     "Formalin")] <- "Disinfectant"

df$TreatmentType[df$Specific.treatment.used.1 %in% c("Sodium chloride",
                                                     "Natural sea salt",
                                                     "Increased salinity")] <- "Salt"

df$TreatmentType[df$Specific.treatment.used.1 == "Electrolyte"] <- "Electrolyte"
df$TreatmentType[df$Specific.treatment.used.1 == "BMP-NTf2"] <- "BMP-NTf2"
df$TreatmentType[df$Specific.treatment.used.1 == "General Tonic®"] <- "General Tonic"

table(df$TreatmentType)

# create new binary column stating whether more than one treatment was used
# and new intervention category column where multiple treatment is regarded
# as its own category
df$MultipleTreatments <- "No"
df$Intervention.category.itra.multi <- df$Intervention.category.1.itra
for (row in seq_len(nrow(df))) {
  if (!is.na(df$Intervention.category.2[row])) {
    df$MultipleTreatments[row] <- "Yes"
    df$TreatmentType[row] <- "Multiple"
    df$Intervention.category.itra.multi[row] <- "Multiple"
  }
}

table(df$TreatmentType)

# view distribution of sample sizes
table(df$Treatment.group.size)
# set maximum sample size for weight calculation
# i.e. what is the cut off sample size whereby any datapoint with a sample
# size larger receives full weight
max_sample_size <- 600


################################################################################
########################## Calculate utilities #################################


# weight for scalability
beta1 <- 0.25
# weight for longevity
beta2 <- 0.25
# minimum weight if scalability and longevity are 0
beta3 <- 0.5
# check beats sum to 1
beta1 + beta2 + beta3

# calculate success using utility function from dare package
df$Success <- sapply(
  seq_len(nrow(df)),
  function(i) {
    calc_metric(
      c(df$Efficacy[i], df$Adverse[i]),
      c(df$Scalability[i], df$Longevity[i]),
      c(beta1, beta2, beta3)
    )
  }
)

# view distribution of success scores
table(df$Success)
hist(df$Success, breaks = 10)


################################################################################
######################### Calculate uncertainties ##############################


############################# Usability weights ################################


# interventions with a usability of 0 are weighted at 100%, 1 at 80%, 2 at 60%
min_grade_weight <- 0.6
df$Usability_weight <- sapply(
  seq_len(nrow(df)),
  function(i) {
    calc_grade_weight(
      df$Usability[i],
      0,
      2,
      min_grade_weight
    )
  }
)


############################ Scaled sample size ################################


# initialise new sample size column to account for capped size and NAs
df$Sample_size <- df$Treatment.group.size
df$Sample_size[is.na(df$Sample_size)] <- 1

df$Sample_size_capped <- df$Sample_size

for (i in seq_along(df$Sample_size_capped)) {
  if (df$Sample_size_capped[i] > max_sample_size){
    df$Sample_size_capped[i] <- max_sample_size
  }
}

table(df$Sample_size_capped)
hist(df$Sample_size_capped, breaks = 20)

# scale samplesize to between 0 and 1
scale_method <- "sigmoid"
df$Scaled_sample_size <- scale_sample_size(df$Sample_size_capped, method = scale_method)

# view scaling graphically
ggplot(df, aes(x = log10(Sample_size_capped))) +
  geom_line(aes(y = Scaled_sample_size)) +
  labs(y = "scaled sample size") +
  scale_x_continuous("log10(sample size)", labels = c(0, 10, 100, 1000)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
#ggsave(paste0(path, "sigmoid_curve.png"))


############################ Uncertainty weights ###############################


# calculate weights to be used in brms meta regression
df$Uncertainty_weights <- calc_uncertainty_weight(df$Scaled_sample_size, df$Usability_weight)


######################### Standard error and CIs ###############################


# estimate standard error for Artificially Constructucted from Aggregate (ACA)
# beta distribution with success as mean and  using effective sample size
# calculated from sample size and usability weight
# use sample_size (where NA treatment size is set as 1) so that every
# datapoint has a standard error estimation
df$Std_error <- calc_std_error(df$Success, df$Sample_size, df$Usability_weight)

# use standard error to calculate confidence interval lower bound
df$Lower_ci <- sapply(seq_len(nrow(df)), function(i) {
  calc_ci(df$Success[i], df$Std_error[i])[1]
})

# use standard error to calculate confidence interval upper bound
df$Upper_ci <- sapply(seq_len(nrow(df)), function(i) {
  calc_ci(df$Success[i], df$Std_error[i])[2]
})

# view results
results <- data.frame(cbind(df$Success,
                            df$Treatment.group.size,
                            df$Lower_ci,
                            df$Upper_ci,
                            df$Uncertainty_weights))
names(results) <- c("Success", "Sample size", "LCI", "UCI", "Weights")
results


################################################################################
############################ Save result #######################################


# set date for output folder
set_up <- paste0("b1_", beta1, "_b2_", beta2, "_minw_", min_grade_weight, "_", scale_method)
path_out <- paste0(path, "outputs/", set_up, "/")

# create output folder
ifelse(!dir.exists(file.path(path_out)),
       dir.create(file.path(path_out), recursive = TRUE),
       FALSE)

# save dataset
write.csv(df, paste0(path_out, "success_df.csv"))


## end of script