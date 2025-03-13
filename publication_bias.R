################################################################################
############################# Publication bias #################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Feb 2025


# clear workspace
rm(list = ls())


################################################################################
################################# Prep #########################################


# set up path to required data
output_options <- "b1_0.25_b2_0.25_minw_0.6_sigmoid/"
path <- paste0("~/Documents/scripts/2025_bd_interventions/outputs/", output_options)

# load data
df <- as.data.frame(read.csv(paste0(path, "success_df.csv")))

# load required packages
library("tidyverse")
library("meta")
library("metafor")
library("patchwork") # for combining plots


################################################################################
############################## Functions #######################################


# function to create a form of funnel plot relevant for beta distributions
plot_beta_funnel <- function(successes, errors) {
  # create data frame of successes and errors
  df <- data.frame(
    success = successes,
    standard_error = errors
  )

  # calculate overall mean
  mean_effect <- mean(successes)
  
  # create theoretical boundary lines (for beta distributions)
  x_seq <- seq(0.001, 0.999, length.out = 1000)
  max_stderr <- sqrt(x_seq * (1 - x_seq))
  boundary_df <- data.frame(
    effect_size = x_seq,
    max_sd = max_stderr
  )
  
  # Create plot
  ggplot() +
    # Add theoretical boundary
    geom_line(data = boundary_df, aes(x = x_seq, y = max_stderr), 
              color = "red", size = 1, linetype = "dashed") +
    # Add observed points
    geom_point(data = df, aes(x = success, y = standard_error), 
               alpha = 0.5) +
    # Add mean line
    geom_vline(xintercept = mean_effect, linetype = "dotted") +
    # Customize appearance
    scale_y_reverse() +
    coord_cartesian(xlim = c(0, 1), ylim = c(max(errors), 0)) +
    labs(
      title = "Funnel Plot for Beta Distributed Variable",
      subtitle = "Red dashed line shows theoretical maximum standard deviation",
      x = "Mean Success",
      y = "Standard Error"
    ) +
    theme_minimal()
}


# view plots showing distribution of successes with respect to errors
# this gives further analysis to supplement funnel plot
perform_density_analysis <- function(successes, errors) {
  # create data frame of successes and errors
  df <- data.frame(
    success = successes,
    standard_error = errors
  )
  
  # plot for overall density of means
  p1 <- ggplot(df, aes(x = success)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    labs(title = "Density of Successes",
         x = "Success", 
         y = "Density") +
    theme_minimal()
  
  # 2D density plot
  p2 <- ggplot(df, aes(x = success, y = standard_error)) +
    geom_density_2d_filled(alpha = 0.8) +
    geom_point(alpha = 0.2) +
    scale_y_reverse() +
    labs(title = "2D Density of Successes and Standard Errors",
         x = "Success",
         y = "Standard Error") +
    theme_minimal()
  
  # plot stratified density
  # split standard errors into bins
  df <- df %>%
    mutate(precision_group = cut(standard_error,
                                breaks = quantile(standard_error, probs = c(0, 0.33, 0.66, 1)),
                                labels = c("High Precision", "Medium Precision", "Low Precision")))
  
  p3 <- ggplot(df, aes(x = success, fill = precision_group)) +
    geom_density(alpha = 0.5) +
    labs(title = "Success Density by Precision Level",
         x = "Success",
         y = "Density",
         fill = "Precision Group") +
    theme_minimal()
  
  # combine plots
  combined_plot <- (p1 / p2 / p3) +
    plot_layout(heights = c(1, 1.5, 1)) +
    plot_annotation(
      title = "Density Analysis for Publication Bias",
      subtitle = "Examining distribution patterns across precision levels"
    )
  
  return(combined_plot)
}


################################################################################
########################### Publication bias ###################################


# look at funnel plot for whole dataset
funnel_full <- plot_beta_funnel(successes=df$Success, errors=df$Std_error)
print(funnel_full)
density_analysis_full <- perform_density_analysis(df$Success, df$Std_error)
print(density_analysis_full)

# split dataset by habitat/individual
# habitat
hab <- df[df$Habitat.or.Individual == "H",]
funnel_hab <- plot_beta_funnel(hab$Success, hab$Std_error)
funnel_hab
density_analysis_hab <- perform_density_analysis(hab$Success, hab$Std_error)
print(density_analysis_hab)
# individual
indv <- df[df$Habitat.or.Individual == "I",]
funnel_indv <- plot_beta_funnel(indv$Success, indv$Std_error)
funnel_indv
density_analysis_indv <- perform_density_analysis(indv$Success, indv$Std_error)
print(density_analysis_indv)

# split dataset by intervention category
table(df$Intervention.category.1)
# chemical
chem <- df[df$Intervention.category.1 == "Chemical",]
funnel_chem <- plot_beta_funnel(chem$Success, chem$Std_error)
funnel_chem
density_analysis_chem <- perform_density_analysis(chem$Success, chem$Std_error)
print(density_analysis_chem)
# non-chemical
nonchem <- df[df$Intervention.category.1 != "Chemical",]
funnel_nonchem <- plot_beta_funnel(nonchem$Success, nonchem$Std_error)
funnel_nonchem
density_analysis_nonchem <- perform_density_analysis(nonchem$Success, nonchem$Std_error)
print(density_analysis_nonchem)


## end of script
