################################################################################
############################ Descriptive analysis ##############################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Aug 2023
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


# set colours for plotting
cols <- palette.colors(palette = "Classic Tableau")
# assign colours for intervention categories for consistency
colours_intcat <- c(
  "Bioaugmentation"=cols[5], 
  "Other chemical"=cols[2], 
  "Climate"=cols[7], 
  "Population demographic"=cols[1],
  "Multiple"="#000000",
  "Itraconazole"=cols[3]
)
colours_intcat_heatmap <- c(
  "Bioaugmentation"=cols[5], 
  "Chemical"=cols[2], 
  "Climate"=cols[7], 
  "Population demographic"=cols[1],
  "Translocation"="#000000",
  "Habitat"="grey"
)
# assign colours for habitat/individual for consistency
cols_hab <- palette.colors(palette = "Tableau 10")
colours_hab <- c(
  "H"=cols_hab[4], 
  "I"=cols_hab[6], 
  "B"=cols_hab[3]
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
########################### Multiple treatments ################################


# view specific treatments used in multiple treatments
df_multi <- df[df$MultipleTreatments == "Yes",]

# plot these as a heat map
df_multi_plot_treat <- df_multi %>%
  group_by(Treatment.used.1, Treatment.used.2, Intervention.category.1, Intervention.category.2) %>%
  summarise(SuccessMn = mean(Success),
            count = n())
df_multi_plot_treat$Treatment.used.2[which(df_multi_plot_treat$Treatment.used.2 == "Other")][1] <- "Other.1"
df_multi_plot_treat$Treatment.used.2[which(df_multi_plot_treat$Treatment.used.2 == "Other")][2] <- "Other.2"
df_multi_plot_treat$Treatment.used.2[which(df_multi_plot_treat$Treatment.used.2 == "Other")][3] <- "Other.3"

df_multi_plot_treat$colour.1 <- NA
df_multi_plot_treat$colour.2 <- NA
for (i in 1:nrow(df_multi_plot_treat)) {
  for (j in 1:length(colours_intcat_heatmap)) {
    if (df_multi_plot_treat$Intervention.category.1[i] == names(colours_intcat_heatmap)[j]) {
      df_multi_plot_treat$colour.1[i] <- colours_intcat_heatmap[j]
    }
    if (df_multi_plot_treat$Intervention.category.2[i] == names(colours_intcat_heatmap)[j]) {
      df_multi_plot_treat$colour.2[i] <- colours_intcat_heatmap[j]
    }
  }
}
plot_cols_multi_treat.1 <- df_multi_plot_treat %>% 
  group_by(Treatment.used.1, colour.1) %>%
  summarise(count = n())
plot_cols_multi_treat.2 <- df_multi_plot_treat %>% 
  group_by(Treatment.used.2, colour.2) %>%
  summarise(count = n())
# create plot
ggplot(df_multi_plot_treat, aes(x = Treatment.used.1, y = Treatment.used.2, fill = SuccessMn)) +
  geom_tile() +
  geom_text(aes(label = ifelse(count == 0, "", count)), # Conditional labels
            vjust = 0.5, hjust = 0.5, size = 3) +
  scale_fill_gradient(low = "#F5F2D0", high = "darkorange", na.value="white") +
  labs(x = "Treatment 1", y = "Treatment 2", fill = "Success") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 320, hjust = 0, colour=plot_cols_multi_treat.1$colour.1),
        axis.text.y = element_text(colour=plot_cols_multi_treat.2$colour.2))
ggsave(paste0(path, "multiple_treatment_heatmap.png"))


# plot these as a heat map
df_multi_plot_spectreat <- df_multi %>%
  group_by(Specific.treatment.used.1, Specific.treatment.used.2, Intervention.category.1, Intervention.category.2) %>%
  summarise(SuccessMn = mean(Success),
            count = n())
df_multi_plot_spectreat$colour.1 <- NA
df_multi_plot_spectreat$colour.2 <- NA
for (i in 1:nrow(df_multi_plot_spectreat)) {
  for (j in 1:length(colours_intcat_heatmap)) {
    if (df_multi_plot_spectreat$Intervention.category.1[i] == names(colours_intcat_heatmap)[j]) {
      df_multi_plot_spectreat$colour.1[i] <- colours_intcat_heatmap[j]
    }
    if (df_multi_plot_spectreat$Intervention.category.2[i] == names(colours_intcat_heatmap)[j]) {
      df_multi_plot_spectreat$colour.2[i] <- colours_intcat_heatmap[j]
    }
  }
}
plot_cols_multi_spectreat.1 <- df_multi_plot_spectreat %>% 
  group_by(Specific.treatment.used.1, colour.1) %>%
  summarise(count = n())
plot_cols_multi_spectreat.2 <- df_multi_plot_spectreat %>% 
  group_by(Specific.treatment.used.2, colour.2) %>%
  summarise(count = n())
# create plot
ggplot(df_multi_plot_spectreat, aes(x = Specific.treatment.used.1, y = Specific.treatment.used.2, fill = SuccessMn)) +
  geom_tile() +
  geom_text(aes(label = ifelse(count == 0, "", count)), # Conditional labels
            vjust = 0.5, hjust = 0.5, size = 3) +
  scale_fill_gradient(low = "#F5F2D0", high = "darkorange", na.value="white") +
  labs(x = "Specific treatment 1", y = "Specific treatment 2", fill = "Success") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 320, hjust = 0, colour=plot_cols_multi_spectreat.1$colour.1),
        axis.text.y = element_text(colour=plot_cols_multi_spectreat.2$colour.2))
ggsave(paste0(path, "multiple_specific_treatment_heatmap.png"))

# quick way to get legend
savePie(df, "Intervention.category.2", "Intervention category 2", colours_intcat_heatmap)


################################################################################
############################ Quantative breakdowns #############################


# store descriptive stats
sink(paste0(path, "descriptive_stats.txt"), append=TRUE)
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
sink()


## end of script
