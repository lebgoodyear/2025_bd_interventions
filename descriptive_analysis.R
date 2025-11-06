################################################################################
############################ Descriptive analysis ##############################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Aug 2023
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


# set colours for plotting
cols <- as.vector(palette.colors(palette = "Classic Tableau"))
# assign colours for intervention categories for consistency
colours_intcat <- c(
  "Bioaugmentation" = cols[5],
  "Other chemical" = cols[2],
  "Climate" = cols[7],
  "Population demographic" = cols[1],
  "Multiple" = "#000000",
  "Itraconazole" = cols[3]
)
# assign colours for habitat/individual for consistency
cols_hab <- as.vector(palette.colors(palette = "Tableau 10"))
colours_hab <- c(
  "H" = cols_hab[4],
  "I" = cols_hab[6],
  "B" = cols_hab[3]
)

# load packages
library("tidyverse")
library("scales") # greyscale colour palette
source(paste0(path, "../../functions_for_plotting.R"))


################################################################################
############################### Pie charts #####################################


# pie chart for category breakdown
savePie(df, "Intervention.category.itra.multi", "Intervention category", colours_intcat)


# pie chart for in situ/ex situ breakdown
savePie(df, "In.situ.or.Ex.situ", "In situ or ex situ", c("#000000", "#ffffff"))


# pie chart for life stage breakdown
savePie(df, "Life.Stage", "Life stage", grey_pal()(length(unique(df$Life.Stage))))


# pie chart for habitat/individual breakdown
savePie(df, "Habitat.or.Individual", "Habitat or individual", colours_hab)#grey_pal()(length(unique(df$Habitat.or.Individual))))


# pie chart for multiple or single treatments breakdown
savePie(df, "MultipleTreatments", "Multiple treatments", c("#000000", "#ffffff"))


# pie chart for therapeutic/prophylactic eakdown
savePie(df, "Therapeutic.or.Prophylactic", "Therapeutic or Prophylactic", c("#000000", "#ffffff"))


# # pie chart for species breakdown set up
# subtax <- df
# # remove any species with less than 5 recrods
# subtax$Taxa[subtax$Taxa %in% tabtax$Taxa[tabtax$Freq < 5]] <- NA
# subtax <- droplevels(subtax)
# subtax <- subtax[!is.na(subtax$Taxa),]
# subtax$Taxa <- as.factor(subtax$Taxa)
# # pie chart for taxa breakdown
# savePie(subtax, "Taxa", "Species", grey_pal()(length(unique(subtax$Taxa))))


################################################################################
############################ Quantative breakdowns #############################


# store descriptive stats
sink(paste0(path, "descriptive_stats.txt"), append = TRUE)
cat(
  "\n\nTotal number of interventions:\n",
  nrow(df),
  "\nTotal number of papers:",
  length(unique(df$Publication.Date.Authorship)),
  "\nNumber of different species used in interventions:",
  length(unique(df$Taxa)),
  "\n\n\nEx situ/in situ breakdown"
)
print(table(df$In.situ.or.Ex.situ))
cat("\n\nLife stage breakdown")
print(table(df$Life.Stage))
cat("\n\nMultiple treatments breakdown")
table(df$MultipleTreatments)
cat("\n\nTherapeutic or prophylactic breakdown")
table(df$Therapeutic.or.Prophylactic)
cat("\n\nHabitat or individual breakdown")
table(df$Habitat.or.Individual)
cat("\n\nInterventions category breakdown")
table(df$Intervention.category.itra.multi)
cat("\n\nTreatment type breakdown")
table(df$TreatmentType)
sink()


################################################################################
##################### Publication/Treatment breakdown ##########################


# save publication vs treatment type breakdown
pub_treat_breakdown <- table(df$Publication.Date.Authorship, df$TreatmentType)
write.csv(pub_treat_breakdown, file = paste0(path, "pub_treat_breakdown.csv"))

# how many publications tried each treatment type
num_pubs <- data.frame(colSums(pub_treat_breakdown > 0))
# how many of each treatment type were trialled
num_treats <- data.frame(table(df$TreatmentType))
# store as one dataframe
treat_count_dat <- merge(num_treats, num_pubs, by.x = "Var1", by.y = "row.names")
colnames(treat_count_dat) <- c("Treatment", "NumTrials", "NumPublications")
# what is the average number of papers per treatment type
treat_count_dat$AvPubsPerTrial <- treat_count_dat$NumPublications / treat_count_dat$NumTrials
treat_count_dat$AvPubsPerTrial <- sprintf("%.2f", treat_count_dat$AvPubsPerTrial)
write.csv(treat_count_dat, file = paste0(path, "pub_treat_freq_breakdown.csv"))


## end of script
