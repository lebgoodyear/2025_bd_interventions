# Bd Interventions

Date created: Feb 2025
Last updated : Feb 2025

Scripts pertaining to meta-analysis on Bd interventions

Workflow:  

1. calculate_success_and_weights.R
Uses utility calculation from the dare package to estimate success, and estimates
standard errors for each data point based on a beta distribution with a mean of success
and a standard deviation based on effective sample size (sample size weighted by 
usability) to give a measure of error. Also calculates weights for use in
meta-regressions.

2. explore_biases.R
Uses a form of funnel plot for beta distrbutions to check for publication bias.