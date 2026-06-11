################################################################################
############################ Chemical safety metric ############################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Jun 2026
# Last edited: Jun 2026


# clear workspace
rm(list = ls())


################################################################################
################################### Set up #####################################


# set up path to required data
path <- "~/Documents/scripts/2025_bd_interventions/"
pathout <- paste0(path, "outputs/b1_0.25_b2_0.25_minw_0.6_sigmoid/")

# load packages
library("tidyverse")
library("roam")
library("ggrepel")
library("ggthemes")

# load data
# intervention dataset
df <- read.csv(paste0(pathout, "success_df.csv"))
# chemical info dataset
chems <- read.csv(paste0(path, "data/chemical_info.csv"))


################################################################################
############################## Data wrangling ##################################


# only keep required columns
colnames(chems)
chems <- chems[,-c(21:25)]
colnames(chems)

# add column for binning number of interventions (for plotting)
table(chems$NoInterventions)
chems <- chems |>
  mutate(
    NoInterventionsBinned  = case_when(
        NoInterventions == 1 ~ "1",
        NoInterventions > 1 & NoInterventions < 5 ~ "2–4",
        NoInterventions >= 5 & NoInterventions <= 10 ~ "5–10",
        NoInterventions > 10 & NoInterventions < 20 ~ "11–20",
        NoInterventions > 20 ~ "21+",
    )
  )
chems$NoInterventionsBinned <- factor(chems$NoInterventionsBinned,
                                      levels = c("1", "2–4", "5–10", "11–20", "21+"))

# initiliase columns for success scores
chems$SuccessMn <- NA # mean
chems$SuccessMd <- NA # median
chems$SuccessMin <- NA # minimum
chems$SuccessMax <- NA # maximum

# extract treatments
df_treatments1 <- data.frame(cbind(df$Specific.treatment.used.1, df$Success))
colnames(df_treatments1) <- c("Treatment", "Success")
df_treatments2 <- data.frame(cbind(df$Specific.treatment.used.2, df$Success))
colnames(df_treatments2) <- c("Treatment", "Success")
df_treatments <- data.frame(rbind(df_treatments1, df_treatments2))
head(df_treatments)
nrow(df_treatments)
# remove NA treatments (resulting from second treatment column)
df_treatments <- df_treatments[!is.na(df_treatments$Treatment),]
# set success as numeric
df_treatments$Success <- as.numeric(df_treatments$Success)

# group by treatment and summarise over success
df_treatments_grouped <- df_treatments %>%
  group_by(Treatment) %>%
  summarise(SuccessMn = mean(Success, na.rm = TRUE),
            SuccessMd = median(Success, na.rm = TRUE),
            SuccessMin = min(Success, na.rm = TRUE),
            SuccessMax = max(Success, na.rm = TRUE))
data.frame(df_treatments_grouped)

# add success scores to chemical info dataset
for (t in 1:nrow(df_treatments_grouped)) {
    if (df_treatments_grouped$Treatment[t] %in% chems$Chemical) {
        chems$SuccessMn[chems$Chemical == df_treatments_grouped$Treatment[t]] <- df_treatments_grouped$SuccessMn[t]
        chems$SuccessMd[chems$Chemical == df_treatments_grouped$Treatment[t]] <- df_treatments_grouped$SuccessMd[t]
        chems$SuccessMin[chems$Chemical == df_treatments_grouped$Treatment[t]] <- df_treatments_grouped$SuccessMin[t]
        chems$SuccessMax[chems$Chemical == df_treatments_grouped$Treatment[t]] <- df_treatments_grouped$SuccessMax[t]        
    }
}

# remove any chemicals with NA success since they were not used in the analysis
nrow(chems)
chems <- chems[!is.na(chems$SuccessMn),]
nrow(chems)


################################################################################
############################# RoAM variable prep ###############################


# equation is Accute Aquatic Hazard x Chronic Aquatic Hazard x 
# (beta1 x Human Toxicity + beta2 x Human Damage + beta3 x Human Irritant + beta4 x Human Long-term Risks + beta5)

# reset Aquatic Hazards so worst score equates to 0, best score equates to 1
# initialise columns
chems$AcAqH <- NA
chems$ChAqH <- NA
chems$AcAqH[chems$AcuteAquaticHazard == 1] <- 0
chems$AcAqH[chems$AcuteAquaticHazard == 2] <- 0.25
chems$AcAqH[chems$AcuteAquaticHazard == 3] <- 0.5
chems$AcAqH[chems$AcuteAquaticHazard == "No"] <- 1
chems$ChAqH[chems$ChronicAquaticHazard == 1] <- 0
chems$ChAqH[chems$ChronicAquaticHazard == 2] <- 0.25
chems$ChAqH[chems$ChronicAquaticHazard == 3] <- 0.5
chems$ChAqH[chems$ChronicAquaticHazard == 4] <- 0.75
chems$ChAqH[chems$ChronicAquaticHazard == "No"] <- 1

# create Human Toxicity variable
chems <- chems |>
  mutate(
    WorstTox = pmin(AcuteToxicityDermal, AcuteToxicityInhalation, AcuteToxicityOral, na.rm = TRUE), # find worst category
    HumTox  = case_when(
      WorstTox == 1 ~ 0,
      WorstTox == 2 ~ 0.25,
      WorstTox == 3 ~ 0.5,
      WorstTox == 4 ~ 0.75,
      .default = 1
    )
  ) |>
  select(-WorstTox) # remove temporary column

# create Human Damage variable
chems$HumDam <- ifelse(chems$SkinIrritant == 1 | 
             chems$EyeIrritation == 1 |
             chems$OrganToxicitySingleExp == "Yes" |
             chems$OrganToxicityRepeatedEx == "Yes", 0, 1) # set to 1 if condition not met

# create Human Irritant variable
chems$HumIrr <- ifelse(chems$SkinIrritant %in% c(2, 3) | 
                       chems$EyeIrritation %in% c(2, 3) |
                       chems$SkinSensitisation == "Yes" |
                       chems$RespiratoryTractIrritation == "Yes", 0, 1)

# create Human Long-term Risks variable
chems$HumLTR <- ifelse(chems$Carcinogenicity == "Yes" | 
             chems$GermCellMutagenicity == "Yes" |
             chems$ReproductiveToxicity == "Yes", 0, 1)

# view data
chems


###############################################################################
############################ RoAM construction ################################


# equation is Accute Aquatic Hazard x Chronic Aquatic Hazard x 
# (beta1 x Human Toxicity + beta2 x Human Damage + beta3 x Human Irritant + beta4 x Human Long-term Risks + beta5)

# set beta
beta1 <- 0.1
beta2 <- 0.1
beta3 <- 0.05
beta4 <- 0.1
beta5 <- 0.65
# combine into vector
betas <- c(beta1, beta2, beta3, beta4, beta5)
# check betas sum to 1
sum(betas)

# calculate metric
chems$Metric <- sapply(
    1:nrow(chems), 
    function(i) {
        calc_metric(
            c(chems$AcAqH[i], chems$ChAqH[i]), 
            c(chems$HumTox[i], chems$HumDam[i], chems$HumIrr[i], chems$HumLTR[i]),
            betas
        )
    }
)

# view dataframe
table(chems$Metric[!is.na(chems$SuccessMn)])
chems$Metric[chems$Chemical == "Itraconazole"]
chems$SuccessMn[chems$Chemical == "Itraconazole"]
chems$Metric[chems$Chemical == "Amphotericin B"]
chems$SuccessMn[chems$Chemical == "Amphotericin B"]


################################################################################
########################### Visualise metric values ############################


# set colours for plotting
colstab <- as.vector(palette.colors(palette = "Tableau 10"))
cols <- colstab[c(7,1,4,2,3)]


ggplot(chems, aes(x = Metric, y = SuccessMn, colour = NoInterventionsBinned)) +
    geom_jitter() +
    geom_text_repel(
        data = filter(chems, Metric != 0),
        aes(label = Chemical),
        show.legend = FALSE,
        max.iter = 10000, # increase iterations for finding positions
        max.time = 1, # set max seconds ggrepel spends searching
        point.padding = 0.5, # space between label and point
        box.padding = 0.5 # space between labels
    ) +
    labs(x = "Chemical Safety Score", y = "Mean Success Score", colour = "Number of Interventions") +
    scale_colour_manual(values = cols) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")
    )
# save plot
ggsave(paste0(pathout, "chemical_safety_vs_success_", beta1, "_", beta2, "_", beta3, "_", beta4, "_", beta5, "_mean.png"))


## end of script