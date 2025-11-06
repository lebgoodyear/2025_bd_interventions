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
    theme_bw(base_family = "serif") +
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


############################### BOX PLOT #######################################


# function to plot labeled boxplot
plotBox <- function(df, var, labs, zero) {

  ## INPUTS
  # df = dataframe containing at leasts columns var1 and var2
  # var1 = column name of categorical variable for boxplot (as string)
  # var2 = column name of continuous numeric variable for boxplot (as string)

  ## OUTPUTS
  # boxp = box plot ggplot2 object

  # create boxplot
  boxp <- ggplot(data = df,
                 aes(x = factor(.data[[var]]), y = Success)) + # nolint
    labs(y = "Success", x = var) +
    geom_boxplot() +
    scale_x_discrete(labels = labs) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270, hjust = 0))

  # save box plots
  if (zero == TRUE) {
    ggsave(paste0(path_out_box, "box_", var, ".png"))
  } else {
    ggsave(paste0(path_out_box, "box_", var, "wo_zero.png"))
  }

  # return boxplot to view in code editor
  return(boxp)
}


############################### PLOT VIOLIN ####################################


# function to plot labeled violin
plotViolin <- function(df, var, labs, pathout) {

  ## INPUTS
  # df = dataframe containing at leasts columns var1 and var2
  # var1 = column name of categorical variable for violin plot (as string)
  # var2 = column name of continuous numeric variable for violin plot(as string)
  # pathout = path location to save png

  ## OUTPUTS
  # violinp = violin plot ggplot2 object

  # create violin plot
  violinp <- ggplot(data = df,
                 aes(x = factor(.data[[var]]), y = Success)) + # nolint
    labs(y = "Success", x = var) +
    geom_violin() +
    scale_x_discrete(labels = labs) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270, hjust = 0))

  # save violin plot
  ggsave(paste0(pathout, "violin_", var, ".png"))

  # return vioin plot to view in code editor
  return(violinp)
}


############################ PLOT INTERACTION ##################################


plotInteraction <- function(df, x, trace, labs_x, labs_trace) {

    ## INPUTS
    # df = dataframe containing at leasts columns var1 and var2
    # x = column name of categorical variable to plot on x-axis
    # trace = column name of categorical variable to plot as trace
    # labs_x = labels for x-axis
    # labs_trace = labels for trace

    ## OUTPUTS
    # generates plot

    par(mar = c(5, 4, 4, 10),xpd = TRUE)
    plot.new()
    interaction.plot(
    x.factor = df[[x]],
    trace.factor = df[[trace]],
    response = df$Success,
    fun = median,
    ylab = "Success",
    xlab = "",
    trace.label = "",
    col = c("#0198f9", "#f95801", "#0198f9", "#f95801"),
    lwd = 4,
    lty = c(1, 1, 2, 2),
    xaxt = "n",
    bty = "l",
    legend = FALSE
    )
    axis(side=1, 1:length(labs_x), labels = labs_x)
    legend("topright", 
        legend = labs_trace, 
        col = c("#0198f9", "#f95801", "#0198f9", "#f95801"),
        lty = c(1, 1, 2, 2),
        bty = "n",
        inset = c(-0.25,0)
    )
}
