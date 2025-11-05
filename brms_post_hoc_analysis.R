################################################################################
########################## brms post hoc analysis ##############################
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
library("ggridges")
library("patchwork")
library("viridis")
source(paste0(path, "../../functions_for_plotting.R"))

# reset factors so that reference levels are linked to highest success scores,
# based on bivairate models, plots and observation of means
# this provides more meaningful comparisons and makes the model more interpretable
df$Habitat.or.Individual <- factor(df$Habitat.or.Individual,
                                   levels = c("I", "H", "B"))
df$Intervention.category.itra.multi <- factor(df$Intervention.category.itra.multi, 
                                              levels = c(
                                                "Itraconazole",
                                                "Bioaugmentation",
                                                "Climate",
                                                "Multiple",
                                                "Other chemical",
                                                "Population demographic"))
df$Life.Stage <- factor(df$Life.Stage,
                                    levels = c("Larvae", "Metamorph", "Juvenile", "Adult"))

# create output folder for brms outputs
path_out_brms <- paste0(path, "brms/")
ifelse(!dir.exists(file.path(path_out_brms)),
       dir.create(file.path(path_out_brms), recursive = TRUE),
       FALSE)

# read in best model for post hoc analysis
mod <- readRDS(paste0(path_out_brms, "modb_best.rds"))
# read in best treatment model
modt <- readRDS(paste0(path_out_brms, "modb_bestt.rds"))

# generate labels for plots
labs_life <- makeLabs(df, "Life.Stage")
labs_hab <- makeLabs(df, "Habitat.or.Individual")
labs_intcat <- makeLabs(df, "Intervention.category.itra.multi")


################################################################################
###################### PERTINENT VARIABLE ASSOCIATIONS #########################


df_pert_plot <- df %>%
  group_by(Intervention.category.itra.multi, Life.Stage, Habitat.or.Individual) %>%
  summarise(SuccessMn = mean(Success),
            count = n())

# view heatmap to see if some variable combinations are more represented than
# others or if there are any stand out interactions
ggplot(df_pert_plot, aes(
  x = Life.Stage,
  y = Intervention.category.itra.multi,
  fill = SuccessMn)) +
  geom_tile() +
  geom_text(aes(label = ifelse(count == 0, "", count)), # Conditional labels
            vjust = 0.5, hjust = 0.5, size = 3) +
  facet_wrap(~ Habitat.or.Individual) +
  scale_fill_gradient(low = "#F5F2D0", high = "darkorange", na.value = "white") +
  labs(x = "Life stage", y = "Intervention category", fill = "Success") +
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))
ggsave(paste0(path, "pert_heatmap.png"))


# view interaction plots
plotInteraction(
  df,
  "Habitat.or.Individual",
  "Life.Stage",
  labs_hab,
  labs_life
)
plotInteraction(
  df,
  "Intervention.category.itra.multi",
  "Habitat.or.Individual",
  labs_intcat,
  labs_hab
)
plotInteraction(
  df,
  "Intervention.category.itra.multi",
  "Life.Stage",
  labs_intcat,
  labs_life
)


################################################################################
################ PLOTS COMPARING ESTIMATE AND ERROR ACROSS LEVELS ##############


c_eff_habindv <- conditional_effects(mod, effects = "Habitat.or.Individual")
habindv_plot <- plot(c_eff_habindv, plot = FALSE)[[1]] +
  theme_minimal()
habindv_plot
c_eff_intcat <- conditional_effects(mod, effects = "Intervention.category.itra.multi")
intcat_plot <- plot(c_eff_intcat, plot = FALSE)[[1]] +
  theme_minimal()
intcat_plot
c_eff_life <- conditional_effects(mod, effects = "Life.Stage")
life_plot <- plot(c_eff_life, plot = FALSE)[[1]] +
  theme_minimal()
life_plot

# plot for treatment type from treatment model
c_eff_treat <- conditional_effects(modt, effects = "TreatmentType")
treat_plot <- plot(c_eff_treat, plot = FALSE)[[1]] +
  theme_minimal()
treat_plot
# as expected from model summary, treatment type only shows statistically
# variable success between itraconazole and population demographic so
# isn't any more informative than intervention category model


################################################################################
############# POSTERIOR DISTRIBUTIONS FOR COMPARISON ACROSS EFFECTS ############


# prepare data for density plots
posterior_draws <- tidy_draws(mod) %>%
  gather_variables() %>%
  filter(grepl("b_", .variable))
posterior_draws <- posterior_draws[-which(posterior_draws$.variable == "b_Intercept"), ]

# prepare data for point estimates and credible intervals
fixed_effects <- posterior_draws %>%
  group_by(.variable) %>%
  summarize(
    estimate = mean(.value),
    lower = quantile(.value, 0.025),
    upper = quantile(.value, 0.975)
  )

# combine data for faceting
combined_data <- bind_rows(
  posterior_draws %>% mutate(plot_type = "density"),
  fixed_effects %>% mutate(plot_type = "point_estimate")
)

# extract grouping variable
combined_data <- combined_data %>%
  mutate(group = str_remove(str_extract(.variable, "b_[^.]*"),"b_")) # Extract and clean base name

# set x-axis limits
x_limits <- c(-6, 3)

# legend labels
legend_labels <- c(
  "Intervention" = "Intervention category",
  "Habitat" = "Habitat/individual",
  "Life" = "Life stage"
)

# facet labels
facet_labels <- c(
  "b_Intervention.category.itra.multiBioaugmentation" = "Bioaugmentation",
  "b_Intervention.category.itra.multiClimate" = "Climate",
  "b_Intervention.category.itra.multiMultiple" = "Multiple",
  "b_Intervention.category.itra.multiOtherchemical" = "Other chemical",
  "b_Intervention.category.itra.multiPopulationdemographic" = "Population demographic",
  "b_Habitat.or.IndividualH" = "Habitat",
  "b_Habitat.or.IndividualB" = "Both individual and habitat",
  "b_Life.StageMetamorph" = "Metamorph",
  "b_Life.StageJuvenile" = "Juvenile",
  "b_Life.StageAdult" = "Adult"
)

# reorder factor levels to place intercept at the top
combined_data$.variable <- factor(combined_data$.variable, levels = names(facet_labels))

# create faceted plot
plot <- ggplot(combined_data, aes(x = .value, y = 1)) + # y=1 is a placeholder for faceting.
  geom_density_ridges(aes(fill = group),
                      data = combined_data %>% filter(plot_type == "density"),
                      scale = 1.5, rel_min_height = 0.01, alpha = 0.5) +
  geom_point(aes(x = estimate),
             data = combined_data %>% filter(plot_type == "point_estimate")) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 data = combined_data %>% filter(plot_type == "point_estimate"),
                 height = 0.2) +
  facet_wrap(~.variable, ncol = 1, labeller = labeller(.variable = facet_labels)) +
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  labs(x = "Estimate", fill = "Predictor type") +
  xlim(x_limits) +
  scale_fill_viridis_d(labels = legend_labels)
# add the continuous vline as a separate layer
plot + geom_vline(xintercept = 0, linetype = "dashed")
ggsave(paste0(path_out_brms, "fixed_eff_density_plots.png"), width = 7, height = 8)


################################################################################
########################## CALCULATE PROBABILITIES #############################


calc_probabilities <- function(df, var, reflev, mod) {
  # extract factor levels
  fctlevs <- levels(as.factor(df[[var]]))

  # initialise matrix to store pairwise comparisons
  n_fcts <- length(fctlevs)
  comparison_matrix <- matrix(NA, n_fcts, n_fcts)
  rownames(comparison_matrix) <- fctlevs
  colnames(comparison_matrix) <- fctlevs

  # extract posterior samples of fixed effects
  post_fixed <- as_draws_df(mod, pars = "^b_")

  # extract column names for all factor levels
  cols <- grep(var, colnames(post_fixed), value = TRUE)

  # loop through all pairs and store pairwise probabilities
  for (i in 1:n_fcts) {
    for (j in 1:n_fcts) {
      if (i != j) {
        # for reference level
        if (fctlevs[i] == reflev) {
          # Reference vs other: probability that reference > other
          other_col <- grep(gsub(" ", "", fctlevs[j]), cols, value = TRUE)[1]
          prob <- 1 - mean(post_fixed[[other_col]] > 0)
        } else if (j == which(fctlevs == reflev)) {
          # Other vs reference: probability that other > reference
          other_col <- grep(gsub(" ", "", fctlevs[i]), cols, value = TRUE)[1]
          prob <- mean(post_fixed[[other_col]] > 0)
        } else {
          # Neither is reference: direct comparison
          col_i <- grep(gsub(" ", "", fctlevs[i]), cols, value = TRUE)[1]
          col_j <- grep(gsub(" ", "", fctlevs[j]), cols, value = TRUE)[1]
          prob <- mean(post_fixed[[col_i]] > post_fixed[[col_j]])
        }
        comparison_matrix[i, j] <- prob
      }
    }
  }
  return(comparison_matrix)
}

# calculate comparison matrices for pertinent predictors
comp_mat_intcat <- calc_probabilities(df, "Intervention.category.itra.multi", "Itraconazole", mod)
write.csv(comp_mat_intcat, paste0(path_out_brms, "intcat_probability_comparison_matrix.csv"))
comp_mat_habindv <- calc_probabilities(df, "Habitat.or.Individual", "I", mod)
write.csv(comp_mat_habindv, paste0(path_out_brms, "hab_probability_comparison_matrix.csv"))
comp_mat_lifestage <- calc_probabilities(df, "Life.Stage", "Larvae", mod)
write.csv(comp_mat_lifestage, paste0(path_out_brms, "lifestage_probability_comparison_matrix.csv"))


# end of script
