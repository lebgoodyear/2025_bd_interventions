################################################################################
########################## Functions for plotting ##############################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Feb 2025
# Last edited: Feb 2025


############################# MAKE TABLE #######################################


# function to generate table with specific column names
makeTable <- function(df, column) {
  tab <- as.data.frame(table(df[[column]]))
  names(tab) <- c(column, "Freq")
  return(tab)
}


############################# MAKE LABELS ######################################


# function to make labels for plots
makeLabs <- function(df, var) {
  df[[var]] <- as.factor(df[[var]]) # set categorical variable as factor
  labs <- c() # initialise empty vector for storing labels
  for (i in seq_len(length(unique(levels(df[[var]]))))) {
    # create labels with records of records per factor
    labs <- c(labs,
              paste0(levels(df[[var]])[i], " (n = ", table(df[[var]])[i], ")"))
  }
  return(labs)
}


############################### PLOT PIE #######################################


# plot pie chart
plotPie <- function(tabdf, column, labs) {
  pie <- ggplot(tabdf, aes(x = "", y = Freq, fill = .data[[column]])) +
    coord_polar("y", start = 0) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_discrete(labels = labs) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
  return(pie)
}


############## COMBINE MAKE TABLE, MAKE LABELS AND PLOT PIE ###################


# function that combines all of the above
savePie <- function(df, column, lab_title, cols) {
  tab <- makeTable(df, column)  
  labs <- makeLabs(df, column)
  plotPie(tab, column, labs) + 
    scale_fill_manual(values = cols, labels = labs) +
    labs(fill = lab_title)
  ggsave(paste0(path, "pie_", column, ".png")) 
}