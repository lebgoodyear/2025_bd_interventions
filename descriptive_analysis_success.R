################################################################################
####################### Success descriptive analysis ###########################
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
library("ggpattern")
library("scales")
library("patchwork")
source(paste0(path, "../../functions_for_plotting.R"))

# set colours for plotting
#colours <- c("#009e73", "#000000", "#d55e00", "#0037a6", "#cc79a7", "#e69f00") ##56b4e9

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
colours_intcat_heatmap <- c(
  "Bioaugmentation" = cols[5],
  "Chemical" = cols[2],
  "Climate" = cols[7],
  "Population demographic" = cols[1],
  "Translocation" = "#000000",
  "Habitat" = "grey"
)
# assign colours for habitat/individual for consistency
cols_hab <- as.vector(palette.colors(palette = "Tableau 10"))
colours_hab <- c(
  "H" = cols_hab[4],
  "I" = cols_hab[6],
  "B" = cols_hab[3]
)


################################################################################
######################### DESCRIPTIVE STATS ####################################


sink(paste0(path, "descriptive_success_stats.txt"), append = TRUE)
cat(
  "\n\nMax success score:",
  max(df$Success),
  "\nMean success score:",
  mean(df$Success),
  "\nEarliest published intervention:",
  min(df$Year.of.publication),
  "\nLatest published interventions:",
  max(df$Year.of.publication),
  "\nNumber of experiments with reported study date:",
  nrow(df[which(!is.na(df$Year.of.treatment)), ]),
  "\nEariest study date:",
  min(df$Year.of.treatment[!is.na(df$Year.of.treatment)]),
  "\nLatest study date:",
  max(df$Year.of.treatment[!is.na(df$Year.of.treatment)])
)
sink()


################################################################################
########## SUCCESS BY TREATMENT TYPE / MULTIPLE TREATMENT BREAKDOWN  ###########


# break success down into bins
df$Success_bins <- cut(df$Success,
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       include.lowest = TRUE)

# set colours for consistency
cols_success <- brewer_pal()(length(unique(df$Success_bins)))
names(cols_success) <- levels(df$Success_bins)

create_treatment_plot <- function(df, cat) {
  cat_treat <- c(
    df$Specific.treatment.used.1[which(df$Intervention.category.1 == cat)],
    df$Specific.treatment.used.2[which(df$Intervention.category.2 == cat)]
  )
  cat_success <- c(
    as.character(df$Success_bins[which(df$Intervention.category.1 == cat)]),
    as.character(df$Success_bins[which(df$Intervention.category.2 == cat)])
  )
  cat_multiple <- c(
    as.character(df$MultipleTreatments[which(df$Intervention.category.1 == cat)]),
    as.character(df$MultipleTreatments[which(df$Intervention.category.2 == cat)])
  )
  cat_treat <- str_trim(cat_treat)
  cat_treat_levels <- levels(fct_rev(fct_infreq(as.factor(cat_treat))))

  cat_treat_success_tab <- as.data.frame(table(cat_treat, cat_multiple, cat_success))

  cat_treat_success_tab$cat_treat <- factor(cat_treat_success_tab$cat_treat,
                                            levels = cat_treat_levels)

  cat_treat_success_tab$cat_success <- factor(cat_treat_success_tab$cat_success, 
                                              levels = c(
                                                "[0,0.2]",
                                                "(0.2,0.4]",
                                                "(0.4,0.6]",
                                                "(0.6,0.8]",
                                                "(0.8,1]"
                                              ))

  return(list(cat_treat_success_tab, max(table(cat_treat))))
}

popdem_treat_success_tab_full <- create_treatment_plot(df, "Population demographic")
popdem_treat_success_tab <- popdem_treat_success_tab_full[[1]]
chem_treat_success_tab_full <- create_treatment_plot(df, "Chemical")
chem_treat_success_tab <- chem_treat_success_tab_full[[1]]
clim_treat_success_tab_full <- create_treatment_plot(df, "Climate")
clim_treat_success_tab <- clim_treat_success_tab_full[[1]]
bioaug_treat_success_tab_full <- create_treatment_plot(df, "Bioaugmentation")
bioaug_treat_success_tab <- bioaug_treat_success_tab_full[[1]]

# find maximum count for setting scale in plot
y_max <- max(
  chem_treat_success_tab_full[[2]],
  clim_treat_success_tab_full[[2]],
  popdem_treat_success_tab_full[[2]],
  bioaug_treat_success_tab_full[[2]]
)

# set size of text for plot
base_text_size <- 14

# function to generate stacked bar plot
plot_stacked_multiple <- function(df, axis_name, pattern_spacing, show_x) {
  ggplot(df[[1]]) +
    geom_bar_pattern(
      aes(x = cat_treat,
          y = Freq,
          fill = cat_success,
          pattern = cat_multiple,
          group = interaction(cat_multiple, cat_success)),
      position = position_stack(reverse = TRUE),
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_density = pattern_spacing,
      pattern_spacing = pattern_spacing,
      pattern_key_scale_factor = 0.6,
      stat = "identity",
      colour = "black",
      linewidth = 0.3,
      width = 0.8
    ) +
    scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
    scale_fill_manual(values = cols_success) +
    labs(
      x = axis_name,
      y = "Count",
      fill = "Success score",
      pattern = "Part of multiple treatment intervention"
    ) +
    theme_bw(base_family = "serif", base_size = base_text_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.line.x = if (!show_x) element_blank() else element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.text = element_text(colour = "black", size = base_text_size),
      axis.text.x = if (!show_x) element_blank() else element_text(size = base_text_size),
      axis.text.y = element_text(size = base_text_size),
      axis.ticks = element_line(colour = "black"),
      axis.ticks.x = if (!show_x) element_blank() else element_line(colour = "black"),
      axis.ticks.y = element_line(colour = "black"),
      axis.title.x = if (!show_x) element_blank() else element_text(size = base_text_size),
      axis.title.y = element_text(vjust = 1, angle = 360, size = base_text_size),
      legend.position = "none",
    ) +
    guides(
      fill = guide_legend(override.aes = list(pattern = "none")),
      pattern = guide_legend(
        override.aes = list(
          fill = "white",
          colour = "black",
          pattern_spacing = 0.01,
          pattern_density = 0.01,
          pattern_colour = "black"
        )
      )
    ) +
    coord_flip(clip = "off") +
    scale_x_discrete(expand = c(0.03, 0)) +
    scale_y_continuous(
      expand = c(0.01, 0),
      limits = c(0, y_max),
      breaks = seq(
        0,
        ceiling(y_max/10) * 10,
        by = ifelse((ceiling(y_max / 10) * 10) > 10, 10, 2)
      )
    )
}

# plot for chemical
p1 <- plot_stacked_multiple(chem_treat_success_tab_full, "(a)", pattern_spacing = 0.01, show_x=FALSE)
# add legend
p1 <- p1 +
  theme(
    legend.position = c(0.85, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.box = "vertical",
    legend.text = element_text(size = base_text_size),
    legend.title = element_text(size = base_text_size)
  )

# plot for climate
p2 <- plot_stacked_multiple(clim_treat_success_tab_full, "(b)", pattern_spacing = 0.05, show_x = FALSE)

# plot for population demographic
p3 <- plot_stacked_multiple(popdem_treat_success_tab_full, "(c)", pattern_spacing = 0.045,show_x = FALSE)

# plot for bioagumentation
p4 <- plot_stacked_multiple(bioaug_treat_success_tab_full, "(d)", pattern_spacing = 0.025, show_x = TRUE)

# calculate number of unique treatments for each panel to find relative heights
# combine data into one data frame
combined_df <- bind_rows(
  chem_treat_success_tab_full[[1]] %>% mutate(panel = "(a) Chemical"),
  clim_treat_success_tab_full[[1]] %>% mutate(panel = "(b) Climate"),
  popdem_treat_success_tab_full[[1]] %>% mutate(panel = "(c) Population/Demographic"),
  bioaug_treat_success_tab_full[[1]] %>% mutate(panel = "(d) Bioaugmentation")
)
# calculate number of unique categories for each panel
heights <- combined_df %>%
  group_by(panel) %>%
  summarise(n_categories = n_distinct(cat_treat)) %>%
  pull(n_categories)

# merge panels into one plot
final_plot <- (p1 / p2 / p3 / p4) +
  plot_layout(heights = heights) &
  theme(plot.margin = margin(20, 20, 20, 20))

# save plot
ggsave(paste0(path, "combined_treatment_success.png"),
       final_plot,
       width = 16, height = 16, unit = "in",
       dpi = 300)


################################################################################
################### FOREST PLOT FOR DIFFERENT TREATMENT TYPES ##################


# here, we are looking at success of specific treatment types so any interventions
# that use two treatments types are included twice, once for each treatment type
# (this means counts of each treatment type include all instances of treatments use)

# remove all preceding whitespace for treatment type
df$Specific.treatment.used.1 <- str_trim(df$Specific.treatment.used.1)
df$Specific.treatment.used.2 <- str_trim(df$Specific.treatment.used.2)

# group data by treatment type
grouped_data <- df %>%
  pivot_longer(
    cols = c(Specific.treatment.used.1, Specific.treatment.used.2),
    names_to = "treatment_type",
    values_to = "treatment"
  ) %>%
  filter(!is.na(treatment)) %>%
  group_by(treatment) %>%
  summarise(
    effect = mean(Success),
    lower = min(Success),
    upper = max(Success),
    .groups = "drop",
    count=n()
  ) %>%
  mutate(
    treatment = fct_reorder(treatment, effect),
    treatment_label = as.factor(paste0(treatment, " (", count, ")")),
    treatment_label = fct_reorder(treatment_label, effect)
  ) %>%
  arrange(desc(effect))

# create forest plot for treatment types
ggplot(grouped_data, aes(y = treatment_label, x = effect,
                         xmin = lower, xmax = upper)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(height = 0.2, color = "#D55E00", linewidth = 1) +
  theme_bw(base_family = "serif", base_size = 16) +
  # Define a consistent theme with larger fonts
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_blank(),
    strip.text = element_text(size = 18, hjust = 0),
    strip.background = element_blank()
  )
ggsave(paste0(path, "treatment_forest.png"), width = 11, height = 14)


################################################################################
########################### Multiple treatments ################################


# view specific treatments used in multiple treatments
df_multi <- df[df$MultipleTreatments == "Yes", ]

table(df$In.situ.or.Ex.situ[df$MultipleTreatments == "Yes"])
table(df$MultipleTreatments[df$Specific.treatment.used.1 == "Itraconazole"])
table(df$MultipleTreatments[df$Specific.treatment.used.2 == "Itraconazole"])
table(df$In.situ.or.Ex.situ[df$Intervention.category.itra.multi == "Multiple"])
table(df$In.situ.or.Ex.situ[df$Specific.treatment.used.2 == "Itraconazole"])
table(df$In.situ.or.Ex.situ[df$Specific.treatment.used.1 == "Itraconazole"])
table(df$Efficacy.Matrix[df$Specific.treatment.used.1 == "Itraconazole"])
table(df$Efficacy.Matrix[df$Specific.treatment.used.2 == "Itraconazole"])

# plot these as a heat map
df_multi_plot_spectreat <- df_multi %>%
  group_by(Specific.treatment.used.1, Specific.treatment.used.2, Intervention.category.1, Intervention.category.2) %>%
  summarise(SuccessMn = mean(Success),
            count = n())
df_multi_plot_spectreat$colour.1 <- NA
df_multi_plot_spectreat$colour.2 <- NA
for (i in seq_len(nrow(df_multi_plot_spectreat))) {
  for (j in seq_along(colours_intcat_heatmap)) {
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

# set font size
text_base_size_multi <- 16

# create plot
ggplot(df_multi_plot_spectreat, aes(x = Specific.treatment.used.1, y = Specific.treatment.used.2, fill = SuccessMn)) +
  geom_tile() +
  geom_text(aes(label = ifelse(count == 0, "", count)), # Conditional labels
            vjust = 0.5, hjust = 0.5, size = 5) +
  scale_fill_gradient(low = "#F5F2D0", high = "darkorange", na.value = "white") +
  labs(x = "Specific treatment 1", y = "Specific treatment 2", fill = "Success") +
  theme_bw(base_family = "serif", base_size = text_base_size_multi) +
  theme(
    axis.text.x = element_text(angle = 320, hjust = 0, colour=plot_cols_multi_spectreat.1$colour.1, size = text_base_size_multi),
    axis.text.y = element_text(colour = plot_cols_multi_spectreat.2$colour.2, size = text_base_size_multi),
    axis.title = element_text(colour = "black", size = text_base_size_multi + 1),
    legend.text = element_text(size = text_base_size_multi - 1),
    legend.title = element_text(size = text_base_size_multi)
  )
ggsave(paste0(path, "multiple_specific_treatment_heatmap.png"), height = 10, width = 17)

# quick way to get legend
savePie(df, "Intervention.category.1", "Intervention category 1", colours_intcat_heatmap)


################################################################################
############################### STACKED BAR ####################################


# create table for stacked bars
tab_stacked <- data.frame(table(df$Habitat.or.Individual, df$Intervention.category.itra.multi))
names(tab_stacked) <- c("Habitat.or.Individual", "Intervention.category.itra.multi", "Count")


# remove the two treatment types with less than 5 trials
df_treatsub <- df[!(df$TreatmentType %in% c("BMP-NTf2", "Medicinal antiparasitic")),]
tab_stacked_treattype <- data.frame(table(df_treatsub$Habitat.or.Individual, df_treatsub$TreatmentType))
names(tab_stacked_treattype) <- c("Habitat.or.Individual", "TreatmentType", "Count")

# make labels for stacked bar
habindv_labs <- makeLabs(df, "Habitat.or.Individual")
intcat_labs <- makeLabs(df, "Intervention.category.itra.multi")
treattype_labs <- makeLabs(df_treatsub, "TreatmentType")

# set font size
text_base_size_bar <- 14

# make plot for intervention category
ggplot(data=tab_stacked, aes(Intervention.category.itra.multi, y=Count, fill=Habitat.or.Individual, colour=Habitat.or.Individual)) +
  geom_bar(
    position = position_stack(reverse = TRUE),
    stat = "identity",
    linewidth = 0.3,
    width = 0.5,
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_hab, labels = habindv_labs) +
  scale_colour_manual(values = colours_hab, labels = habindv_labs) +
  scale_x_discrete(labels = intcat_labs) +
  labs(
    x = "Intervention category",
    y = "Count",
    fill = "Habitat or Indiviudual",
    colour = "Habitat or Indiviudual"
  ) +
  guides(fill = guide_legend(title = "Habitat or individual"),
         colour = guide_legend(title = "Habitat or individual")) +
  theme_bw(base_family = "serif", base_size = text_base_size_bar) +
  theme(
    panel.grid.major.y = element_blank(),  # Remove vertical gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor vertical gridlines
    panel.grid.minor.x = element_blank(),  # Remove minor horizontal gridlines
    panel.border = element_blank(),        # Remove panel border
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 70, hjust = 1),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )
ggsave(paste0(path, "stacked_bar_intcat_habindv.png"), width = 12, height = 6, dpi = 300)

# make plot for treatment type
ggplot(data=tab_stacked_treattype, aes(TreatmentType, y=Count, fill=Habitat.or.Individual, colour=Habitat.or.Individual)) +
  geom_bar(
    position = position_stack(reverse = TRUE),
    stat = "identity",
    linewidth = 0.3,
    width = 0.5,
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_hab, labels = habindv_labs) +
  scale_colour_manual(values = colours_hab, labels = habindv_labs) +
  scale_x_discrete(labels = treattype_labs) +
  labs(
    x = "Treatment Type",
    y = "Count",
    fill = "Habitat or Indiviudual",
    colour = "Habitat or Indiviudual"
  ) +
  guides(fill = guide_legend(title = "Habitat or individual"),
         colour = guide_legend(title = "Habitat or individual")) +
  theme_bw(base_family = "serif", base_size = text_base_size_bar) +
  theme(
    panel.grid.major.y = element_blank(),  # Remove vertical gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor vertical gridlines
    panel.grid.minor.x = element_blank(),  # Remove minor horizontal gridlines
    panel.border = element_blank(),        # Remove panel border
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 70, hjust = 1),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )
ggsave(paste0(path, "stacked_bar_treattype_habindv.png"), width = 10, height = 5)


################################################################################
################################ VIOLINS #######################################


# ggplot(df_plot, aes(x=Habitat.or.Individual, y=Success)) +
# geom_violin(alpha=0.5,colour="grey", width=0.8) +
# geom_violin(aes(fill=as.factor(Intervention.category.1.itra), colour=as.factor(Intervention.category.1.itra)), 
#   alpha = 0.5,
#   position = position_dodge(width=0.5), 
#   width = 0.7,
#   drop=TRUE,
#   scale="area") +
# scale_fill_manual(values=colours_intcat) +
# scale_colour_manual(values=colours_intcat) +
# labs(x="Habitat or individual") +
# guides(fill=guide_legend(title="Intervention category"),
#   colour=guide_legend(title="Intervention category")) +
# #geom_boxplot(aes(fill=as.factor(Intervention.category.1.itra)), alpha = 0.5, position = position_dodge(width=0.6)) +
# theme_minimal()

# intervention category
ggplot(df, aes(x = Intervention.category.itra.multi, y = Success)) +
  geom_violin(alpha = 0.5, colour = "grey", width = 0.8) +
  geom_violin(aes(fill = as.factor(Habitat.or.Individual), colour = as.factor(Habitat.or.Individual)), 
              linewidth = 0.3,
              position = position_dodge(width = 0.6),
              width = 0.7,
              alpha = 0.5,
              drop = TRUE,
              scale = "area") +
  scale_fill_manual(values = colours_hab, labels = habindv_labs) +
  scale_colour_manual(values = colours_hab, labels = habindv_labs) +
  scale_x_discrete(labels = intcat_labs) +
  labs(x = "Intervention category") +
  guides(fill = guide_legend(title = "Habitat or individual"),
         colour = guide_legend(title = "Habitat or individual")) +
  #geom_boxplot(aes(fill=as.factor(Intervention.category.1.itra)), alpha = 0.5, position = position_dodge(width=0.6)) +
  theme_bw(base_family = "serif", base_size = text_base_size_bar) +
  theme(
    panel.grid.major.y = element_blank(),  # Remove vertical gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor vertical gridlines
    panel.grid.minor.x = element_blank(),  # Remove minor horizontal gridlines
    panel.border = element_blank(),        # Remove panel border
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 70, hjust = 1),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )
ggsave(paste0(path, "violin_intcat_habindv.png"), width = 12, height = 8, dpi = 300)

# treatment type
ggplot(df_treatsub, aes(x = TreatmentType, y = Success)) +
  geom_violin(alpha = 0.5, colour = "grey", width = 0.8) +
  geom_violin(aes(fill = as.factor(Habitat.or.Individual), colour = as.factor(Habitat.or.Individual)), 
              linewidth = 0.3,
              position = position_dodge(width = 0.6),
              width = 0.7,
              alpha = 0.5,
              drop = TRUE,
              scale = "area") +
  scale_fill_manual(values = colours_hab, labels = habindv_labs) +
  scale_colour_manual(values = colours_hab, labels = habindv_labs) +
  scale_x_discrete(labels = treattype_labs) +
  labs(x = "Intervention category") +
  guides(fill = guide_legend(title = "Habitat or individual"),
         colour = guide_legend(title = "Habitat or individual")) +
  #geom_boxplot(aes(fill=as.factor(Intervention.category.1.itra)), alpha = 0.5, position = position_dodge(width=0.6)) +
  theme_bw(base_family = "serif", base_size = text_base_size_bar) +
  theme(
    panel.grid.major.y = element_blank(),  # Remove vertical gridlines
    panel.grid.minor.y = element_blank(),  # Remove minor vertical gridlines
    panel.grid.minor.x = element_blank(),  # Remove minor horizontal gridlines
    panel.border = element_blank(),        # Remove panel border
    axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 70, hjust = 1),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )
ggsave(paste0(path, "violin_treattype_habindv.png"), width = 15, height = 8)

# ggplot(df_plot, aes(x=Intervention.category.1.itra, y=Success)) +
# geom_violin() +
# theme_minimal()

# ggplot(df, aes(x=Habitat.or.Individual, y=Success)) +
# geom_violin() +
# theme_minimal()

#ggplot(df_plot, aes(x=Life.Stage, y=Success)) +
#geom_violin() +
#theme_minimal()

# ggplot(df_plot, aes(x=In.situ.or.Ex.situ, y=Success)) +
# geom_violin() +
# theme_minimal()

#ggplot(df_plot, aes(x=TaxaGroup, y=Success)) +
#geom_violin() +
#theme_minimal()

ggplot(df, aes(x = Efficacy.Matrix, y = Success)) +
  geom_violin() +
  theme_minimal()


################################################################################
############################ MEANS SCATTER PLOT ################################


dfs <- df
# get mean of success for different groups
int_means <- dfs %>%
  group_by(Intervention.category.itra.multi) %>%
  summarise(mean_int = mean(Success))

hori_means <- dfs %>%
  group_by(Habitat.or.Individual) %>%
  summarise(mean_hori = mean(Success))

tax_means <- dfs %>%
  group_by(TaxaGroup) %>%
  summarise(mean_tax = mean(Success))

ggplot(data  = dfs,
       aes(x = Intervention.category.itra.multi,
           y = Success,
           col = Intervention.category.itra.multi,
           shape = Habitat.or.Individual)) +
  #        size = Exsitu)) +
  geom_point(position = "jitter", size = 5) +
  geom_hline(data = int_means,
             aes(yintercept = mean_int,
                 col = Intervention.category.itra.multi)) +
  geom_hline(data = hori_means,
             aes(yintercept = mean_hori,
                 linetype = Habitat.or.Individual),
             col = "darkgrey") +
  labs(x = "Intervention category",
       colour = "Intervention category",
       shape = "Taxonomic group",
       linetype = "Habitat or individual") +
  scale_colour_manual(values = colours_intcat, name = "Intervention category") +
  #scale_linetype_discrete(labels=c("Both", "Habitat", "Individual")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank())
#ggsave(paste0(path_out, "plot_coloured_int.png"), width = 16, height = 8, dpi = 150, units = "in", device='png')


################################################################################
###################### SUCCESS AND WEIGHTS HISTOGRAMS ##########################


# success
histwidth <- 0.1
ggplot(data = df, aes(x = Success)) +
  geom_histogram(col = "black", fill = "white", binwidth = histwidth) +
  geom_density(aes(y = histwidth * after_stat(count))) +
  labs(y = "Count") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(paste0(path, "success_hist.png"),
       width = 16, height = 8, unit = "in")

# uncertainty weights
histwidth <- 0.1
ggplot(data = df, aes(x = Uncertainty_weights)) +
  geom_histogram(col = "black", fill = "white", binwidth = histwidth) +
  labs(y = "Count", x = "Uncertainty weights") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(paste0(path, "weights_hist.png"),
       width = 16, height = 8, unit = "in")


################################################################################
########################### TEMPORAL PLOTS #####################################


# cumulative interventions by year and intervention
df_cum_pubyear <- df %>%
  group_by(Year.of.publication, Intervention.category.itra.multi) %>%
  summarise(Count = n()) %>%
  group_by(Intervention.category.itra.multi) %>%
  mutate(CumCount = cumsum(Count))
# plot
ggplot(data = df_cum_pubyear,
       aes(x = Year.of.publication,
           y = CumCount,
           colour = Intervention.category.itra.multi)) +
  geom_jitter() +
  scale_colour_manual(values = colours_intcat) +
  theme_bw()
ggsave(paste0(path, "plot_cum_intervention_pubyear.png"))


# interventions count and success by treatment year
df_count_treatyear_mean_suc <- df %>%
  group_by(Year.of.treatment, Intervention.category.itra.multi) %>%
  mutate(SuccessMn = mean(Success)) %>%
  group_by(Year.of.treatment, Intervention.category.itra.multi, SuccessMn) %>%
  summarise(Count = n())
# number with year of treatment not NA
sum(df_count_treatyear_mean_suc$Count[!is.na(df_count_treatyear_mean_suc$Year.of.treatment)])
# scaling factor for axes by rounding up to nearest 10
df_count_treatyear <- df %>%
  group_by(Year.of.treatment) %>%
  summarise(Count = n())
top_limit_treatment <- ceiling((max(df_count_treatyear$Count)) / 10) * 10
# plot
ggplot() +
  geom_bar(
    data = df_count_treatyear_mean_suc,
    aes(x = Year.of.treatment,
        y = Count,
        fill = Intervention.category.itra.multi),
    stat = "identity",
    alpha = 0.3
  ) +
  geom_jitter(
    data = df,
    aes(x = Year.of.treatment,
        y = Success * top_limit_treatment,
        colour = Intervention.category.itra.multi),
    shape = 4,
    size = 2
  ) +
  geom_point(
    data = df_count_treatyear_mean_suc,
    aes(x = Year.of.treatment,
        y = SuccessMn * top_limit_treatment,
        colour = Intervention.category.itra.multi),
    shape = 19,
    size = 4
  ) +
  #ylim(0, top_limit_treatment) +
  #xlim(1999, 2023) +
  scale_y_continuous(
    # name of first axis
    name = "Number of interventions trialed",
    # add second axis with name and relation to first axis
    sec.axis = sec_axis(~ . / top_limit_treatment, name = "Success")
  ) +
  labs(x = "Year of treatment") +
  scale_fill_manual(values = colours_intcat) +
  scale_colour_manual(values = colours_intcat) +
  theme_bw() +
  theme(panel.border = element_blank())
# save plot
ggsave(paste0(path, "plot_intervention_treatyear.png"))


# interventions count and success by publication year
df_count_pubyear_mean_suc <- df %>%
  group_by(Year.of.publication, Intervention.category.itra.multi) %>%
  mutate(SuccessMn = mean(Success)) %>%
  group_by(Year.of.publication, Intervention.category.itra.multi, SuccessMn) %>%
  summarise(Count = n())
# number with year of treatment not NA
sum(df_count_pubyear_mean_suc$Count[!is.na(df_count_pubyear_mean_suc$Year.of.publication)])
# scaling factor for axes by rounding up to nearest 10
df_count_pubyear <- df %>%
  group_by(Year.of.publication) %>%
  summarise(Count = n())
top_limit_publication <- ceiling((max(df_count_pubyear$Count)) / 10) * 10
# plot
ggplot() +
  geom_bar(
    data = df_count_pubyear_mean_suc,
    aes(x = Year.of.publication,
        y = Count,
        fill = Intervention.category.itra.multi),
    stat = "identity",
    alpha = 0.3
  ) +
  geom_jitter(
    data = df,
    aes(x = Year.of.publication,
        y = Success * top_limit_publication,
        colour = Intervention.category.itra.multi),
    shape = 4,
    size = 2
  ) +
  geom_point(
    data = df_count_pubyear_mean_suc,
    aes(x = Year.of.publication,
        y = SuccessMn * top_limit_publication,
        colour = Intervention.category.itra.multi),
    shape = 19,
    size = 4
  ) +
  #ylim(0, top_limit_publication) +
  #xlim(1999, 2023) +
  scale_y_continuous(
    # name of first axis
    name = "Number of interventions trialed",
    # add second axis with name and relation to first axis
    sec.axis = sec_axis(~ . / top_limit_publication, name="Success")
  ) +
  labs(x = "Year of publication", colour = "Intervention category", fill = "Intervention category") +
  scale_fill_manual(values = colours_intcat) +
  scale_colour_manual(values = colours_intcat) +
  theme_bw(base_family = "serif")
  theme(panel.border = element_blank())
# save plot
ggsave(paste0(path, "plot_intervention_pubyear.png"), width = 10,height = 6)


################################################################################
############################# TOP SCORERS ######################################


# subset by top scorers
top <- df[df$Success > 0.7, ]
top$Success
table(top$Intervention.category.itra.multi)
table(top$TreatmentType)
table(top$Treatment.used.1)
table(top$Specific.treatment.used.1)
table(top$Habitat.or.Individual)
table(top$In.situ.or.Ex.situ)
table(top$Life.Stage)
table(top$Therapeutic.or.Prophylactic)
table(top$Scalability)
table(top$Longevity)
# number of different publications in top scorers
nrow(table(top$Publication.Date.Authorship))
# number of different species in top scorers
nrow(table(top$Taxa))
table(df$Specific.treatment.used.1)

# save csv of top scorers
write.csv(top, paste0(path, "top_scorers.csv"))

# histogram of uncertainty weights for top scorers
histwidth <- 0.05
ggplot(data = top, aes(x = Uncertainty_weights)) +
  geom_histogram(col = "black", fill = "white", binwidth = histwidth) +
  labs(y = "Count", x = "Uncertainty weights") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(paste0(path, "weights_hist_top.png"),
       width = 16, height = 8, unit = "in")

# treatment breakdown for top scorers
# remove all preceding whitespace for treatment type
top$Specific.treatment.used.1 <- str_trim(top$Specific.treatment.used.1)
top$Specific.treatment.used.2 <- str_trim(top$Specific.treatment.used.2)
# group data by treatment type
grouped_top <- top %>%
  pivot_longer(
    cols = c(Specific.treatment.used.1, Specific.treatment.used.2),
    names_to = "treatment_type",
    values_to = "treatment"
  ) %>%
  filter(!is.na(treatment)) %>%
  group_by(treatment) %>%
  summarise(
    effect = mean(Success),
    lower = min(Success),
    upper = max(Success),
    .groups = "drop",
    count = n()
  ) %>%
  mutate(
    treatment = fct_reorder(treatment, effect),
    treatment_label = as.factor(paste0(treatment, " (", count, ")")),
    treatment_label = fct_reorder(treatment_label, effect)
  ) %>%
  arrange(desc(effect))

# prep for plotting pie chart
tab_treat <- data.frame(grouped_top$treatment, grouped_top$count)
names(tab_treat) <- c("treatment", "Freq")
tab_treat <- tab_treat %>%
  mutate(
    treatment = fct_reorder(treatment, Freq, .desc = TRUE),
    labels = as.factor(paste0(treatment, " (", Freq, ")")),
    labels = fct_reorder(labels, Freq, .desc = TRUE)
  ) %>%
  arrange(desc(Freq))
labs_treat <- tab_treat$labels
# plot pie chart
pie_treat <- plotPie(tab_treat, "treatment", labs_treat)
pie_treat +
  scale_fill_manual(values = grey_pal()(nrow(tab_treat)), labels = labs_treat)


## end of script
