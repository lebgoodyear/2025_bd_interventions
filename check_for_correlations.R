################################################################################
########################### Check for correlations #############################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
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


# load packages
library("tidyverse")
library("rcompanion") # for CramerV()
library("polycor") # for polyserial()
library("corrplot") # for corrplot()


################################################################################
################################ Functions #####################################


# function to interpret correlation strength
interpret_correlations <- function(value, method, dof) {
  if (method == "cramers_v") {
    if (dof == 1) {
      if (value < 0.10) return("Weak")
      else if (value < 0.30) return("Moderate")
      else if (value < 0.50) return("Strong")
      else return("Very Strong")
    } else if (dof == 2) {
      if (value < 0.07) return("Weak")
      else if (value < 0.21) return("Moderate")
      else if (value < 0.35) return("Strong")
      else return("Very Strong")
    } else if (dof == 3) {
      if (value < 0.06) return("Weak")
      else if (value < 0.17) return("Moderate")
      else if (value < 0.29) return("Strong")
      else return("Very Strong")
    } else if (dof == 4) {
      if (value < 0.05) return("Weak")
      else if (value < 0.15) return("Moderate")
      else if (value < 0.25) return("Strong")
      else return("Very Strong")
    } else {
      if (value < 0.05) return("Weak")
      else if (value < 0.13) return("Moderate")
      else if (value < 0.22) return("Strong")
      else return("Very Strong")
    }
  } else {  # for Pearson's r and polyserial
    if (abs(value) < 0.30) return("Weak")
    else if (abs(value) < 0.50) return("Moderate")
    else if (abs(value) < 0.70) return("Strong")
    else return("Very Strong")
  }
}


# function to calculate correlation using either  pearson, cramers V  or polyserial methods
calculate_correlations <- function(df) {
  n <- ncol(df)
  cor_matrix <- matrix(1, n, n) # initialise correlation matrix
  method_matrix <- matrix("", n, n) # initialise method matrix
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (is.numeric(df[[i]]) && is.numeric(df[[j]])) { # pearson if both are numeric
        cor_matrix[i,j] <- cor_matrix[j,i] <- cor(df[[i]], df[[j]], use="pairwise.complete.obs")
        method_matrix[i,j] <- method_matrix[j,i] <- "pearson"
      } else if (is.factor(df[[i]]) && is.factor(df[[j]])) {  # use cramer's v if both are categorical
        cor_matrix[i,j] <- cor_matrix[j,i] <- cramerV(table(df[[i]], df[[j]]))
        method_matrix[i,j] <- method_matrix[j,i] <- "cramers_v"
      } else {
        # use polyserial is one is categorial and one numeric
        method_matrix[i,j] <- method_matrix[j,i] <- "polyserial"
        if (is.factor(df[[i]])) { # check they are the right way around for polyserial function
          cor_matrix[i,j] <- cor_matrix[j,i] <- polyserial(df[[j]], df[[i]])
        } else {
          cor_matrix[i,j] <- cor_matrix[j,i] <- polyserial(df[[i]], df[[j]])
        }
      }
    }
  }
  colnames(cor_matrix) <- rownames(cor_matrix) <- names(df)
  colnames(method_matrix) <- rownames(method_matrix) <- names(df)
  return(list(correlation = cor_matrix, method = method_matrix))
}


################################################################################
########################### Check for correlations #############################


# create a dataframe of predictors (as factors) only
dfggf <- as.data.frame(lapply(df[c("In.situ.or.Ex.situ", 
                                   "Life.Stage",
                                   "Intervention.category.itra.multi", 
                                   "Habitat.or.Individual",
                                   "Therapeutic.or.Prophylactic",
                                   "Climate",
                                   "TaxaGroup", 
                                   "Activity",
                                   "Habitat",
                                   "Efficacy.Matrix",
                                   "Publication.Date.Authorship")],
                              factor))
# combine predictors (factors and continuous) and response to test for correlations
dfgg <- as.data.frame(cbind(dfggf, log1p(df$SVLMx), log1p(df$ClutchMn), df$Success))
dfgg <- droplevels(dfgg) # drop any unused factor levels
str(dfgg)

# Calculate correlation matrix
cor_results <- calculate_correlations(dfgg)

# Create a custom color palette
colpal <- colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200)

# Prepare the correlation matrix for plotting
plot_cor <- cor_results$correlation

# create heatmap to show correlations
png(paste0(path, "correlations_heatmap.png"),
    width = 10, height = 6,
    units = "in",
    res = 300)
plot.new()
corrplot(plot_cor,
         method = "color",
         col = colpal,
         type = "upper",
         order = "original", # use original order to prevent clustering issue due to NaNs
         na.label = "NA",
         na.label.col = "grey50",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         mar = c(0,0,1,0),
         title = "Correlation Plot (Mixed Types)")
legend("bottom",
       legend = c("Pearson / Polyserial", "Cramer's V", "NA"),
       fill = c("#6D9EC1", "#E46726", "grey80"),
       border = "black",
       bty = "n",
       horiz = TRUE)
dev.off()

# set duplicates to NA
plot_cor[lower.tri(plot_cor)] <- NA
# save correlation matrix
write.csv(plot_cor, paste0(path, "correlations_matrix.csv"))

# make list of all correlations including method of calculation and degrees of freedom
weak_to_strong <- which(cor_results$correlation != 1, arr.ind = TRUE)
correlation_list_full <- data.frame(
  var1 = rownames(cor_results$correlation)[weak_to_strong[,1]],
  var2 = colnames(cor_results$correlation)[weak_to_strong[,2]],
  correlation = cor_results$correlation[weak_to_strong],
  method = cor_results$method[weak_to_strong],
  degrees_of_freedom = sapply(1:nrow(weak_to_strong), function(i) {
    if (cor_results$method[weak_to_strong][i] == "cramers_v") {
      num_dof_1 <- length(levels(dfgg[,rownames(cor_results$correlation)[weak_to_strong[i,1]]]))
      num_dof_2 <- length(levels(dfgg[,rownames(cor_results$correlation)[weak_to_strong[i,2]]]))
      dof <- min(num_dof_1, num_dof_2) - 1
    } else {
      dof <- nrow(dfgg) - 2
    }
    return(dof)
  })
)

# add strength to correlation list
correlation_list_full$strength <- sapply(1:nrow(weak_to_strong), function(i) 
                             interpret_correlations(
                               abs(cor_results$correlation[weak_to_strong[i,1], weak_to_strong[i,2]]), 
                               cor_results$method[weak_to_strong[i,1], weak_to_strong[i,2]],
                               correlation_list_full$degrees_of_freedom[i]
                             ))

# order list with strongest first
correlation_list_full_ord <- correlation_list_full[order(-abs(correlation_list_full$correlation)),]

# delete alternate entries since they are repeats (pairwise comparisons)
toDelete <- seq(1, nrow(correlation_list_full_ord), 2)
correlation_list <- correlation_list_full_ord[toDelete,]

# save correlation list with methods and degrees of freedom
write.csv(correlation_list, paste0(path, "correlations.csv"))

# only keep strong correlations
correlation_list_final <- correlation_list[!(correlation_list$strength %in% c("Weak", "Moderate")),]

# save strong correlations
write.csv(correlation_list_final, paste0(path, "sig_correlations.csv"))


## end of script