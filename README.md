# Bd Interventions

Date created: Feb 2025
Last updated : Mar 2025

Scripts pertaining to meta-analysis on Bd interventions

Workflow:  

1. calculate_success_and_weights.R 
Uses utility calculation from the dare package to estimate success, and estimates
standard errors for each data point based on a beta distribution with a mean of success
and a standard deviation based on effective sample size (sample size weighted by 
usability) to give a measure of error. Also calculates weights for use in
meta-regressions.

2. publication_bias.R 
Uses a form of funnel plot for beta distrbutions to check for publication bias.

3. check_for_correlations.R 
Uses mixed method correlation calculations to evaluate correlations between all
pertinent variables to check for collinearity issues.

4. functions_for_plotting.R 

5. descriptive_analysis.R 
Plots of the different predictors.

6. mapping.R 

7. descriptive_analysis_success.R 
Plots of the different variables with respect to success.

8. grouping_variable_analyses.R 
Compare effects of efficacy matrix, publication, and amphibian taxa using random
effects models.

9. bivariate_analyses.R 
brms models with each predictor separately.

10. brms_full_models.R 
11. brms_subset_models.R 
12. brms_post_hoc_analysis.R 
Generation comparison of brms models with different variables. Produces post hoc 
probabilites and plots.