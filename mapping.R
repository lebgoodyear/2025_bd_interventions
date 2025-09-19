################################################################################
################################## Mapping #####################################
################################################################################


# Author: Luke Goodyear (lgoodyear01@qub.ac.uk)
# Date created: Oct 2024
# Last edited: Sept 2025


# clear workspace
rm(list = ls())


################################################################################
################################### Set up #####################################


# load packages
library("tidyverse")
library("sf")
library("rnaturalearth") # remember to also install r-rnaturalearthdata
#devtools::install_github("ropensci/rnaturalearthhires")
library("rnaturalearthhires")
theme_set(theme_bw())
theme_update(axis.line = element_blank(),
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

# set up path to required data
output_options <- "b1_0.25_b2_0.25_minw_0.6_sigmoid/"
path <- paste0("~/Documents/scripts/2025_bd_interventions/outputs/", output_options)

# load data
df <- as.data.frame(read.csv(paste0(path, "success_df.csv")))

# set colours for plotting points
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

# set lowest, highest and NA colours for plotting
cols_heatmap <- c("#873307", "#FCE7AA", "#ECECEC") #c("#0037a6", "#FF5608", "#A9A9A9")

# create output folder
path_out <- paste0(path, "map/")
ifelse(!dir.exists(file.path(path_out)),
       dir.create(file.path(path_out), recursive = TRUE),
       FALSE)

#set path to maps workflow
maps_path <- "~/Documents/scripts/BdSusceptibilityProject/mapping/maps_workflow/"
# set path
path_geo <- path_out

# create Mollweide projection coordinates
mollweide_crs <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"


################################################################################
################################ Create map ####################################


################################### Set up #####################################


print("Loading world map data...")
world <- ne_countries(scale = "large", returnclass = "sf")
#coast <- ne_coastline(scale = "large", returnclass = "sf")

area_sub <- world
# create Mollweide projection
#coast_projected <- st_transform(coast, crs = mollweide_crs)
# fix error by making sure object is valid
#coast_projected_fix <- st_make_valid(coast_projected)
# save coastline map to rds object for plotting world without country lines
#saveRDS(coast_projected_fix, paste0(path_geo, "moll_coast.rds"))

# create Mollweide projection
area_projected <- st_transform(area_sub, crs = mollweide_crs)
# fix error by making sure object is valid
area_projected_fix <- st_make_valid(area_projected)

# save map to rds object for plotting
saveRDS(area_projected_fix, paste0(path_geo, "moll.rds"))

# create grid over map
grid <- st_make_grid(area_sub, cellsize = 1, what = "polygons")

# transform grid to mollweide coordinates
grid_projected <- st_transform(grid, crs = mollweide_crs)

# combine grid and map
area_grid_projected <- st_intersection(grid_projected, area_projected_fix)

# create grid IDs for identifying each grid
area_grid_projected <- st_sf(gridId = seq_along(lengths(area_grid_projected)), area_grid_projected)
area_grid_projected <- st_cast(area_grid_projected)

# save combined grid and map to rds object
saveRDS(area_grid_projected, paste0(path_geo, "moll_grid.rds"))


################################# Prep shapefile ###############################


# path to shape file with location information (IUCN)
# column with taxa names must be called "BINOMIAL" (IUCN)
path_in_shape <- paste0(path, "../../data/shapefile/data_0.shp")
# import IUCN shapefile
shapefile <- read_sf(path_in_shape)

# subset shapefile by dataset
shapefile <- shapefile[which(shapefile$BINOMIAL %in% df$Taxa == TRUE), ]

# create mollweide projection of shapefile
shapefile_projected <- st_transform(shapefile, crs = mollweide_crs)


################################# Make grids ###################################


print("Joining data with map objects...")
# join with map
grid <- st_join(x = area_grid_projected,
                y = st_cast(shapefile_projected) %>% select(BINOMIAL),
                left = TRUE, join = st_intersects)

# save combined object
saveRDS(grid, paste0(path_geo, "moll_dataset_grid.rds"))


############################# Create map object ################################


print(paste0("Creating map oject for species richness..."))
# sum by number of species
grid_richness <- grid %>%
  group_by(gridId) %>%
  summarise(SpeciesRichness = n_distinct(BINOMIAL, na.rm = TRUE),
            spsList = paste(BINOMIAL, collapse = ";"),
            .groups = "drop")

# save object
saveRDS(grid_richness, paste0(path_geo, "grid_richness.rds"))
grid_richness <- readRDS(paste0(path_geo, "grid_richness.rds"))

#grid_richness <- readRDS(paste0(path_geo, "grid_richness.rds"))

# fix projection error due to different versions GDAL (if using both local and hpc)
st_crs(grid_richness) <- mollweide_crs


################################################################################
######################### Coordinates in situ interventions ####################


# function to convert dms to decimal
dms_to_dd <- function(dms) {
  # Extract degrees, minutes, seconds
  parts <- strsplit(dms, "[°'\"]")[[1]]
  degrees <- as.numeric(parts[1])
  minutes <- as.numeric(parts[2])
  seconds <- as.numeric(parts[3])

  # Convert to decimal degrees
  dd <- degrees + minutes / 60 + seconds / 3600

  # Adjust sign based on hemisphere
  if (grepl("[SWsw]", dms)) dd <- -dd

  return(dd)
}

# convert longs and lats to decimal
df$LatitudeSS_dd <- sapply(df$LatitudeSS, dms_to_dd)
df$LongitudeSS_dd <- sapply(df$LongitudeSS, dms_to_dd)


# subset by in situ only to mark sites on mark
df_insitu <- df[which(df$In.situ.or.Ex.situ == "In situ"), ]

# count number per intervention at each site to use count for size of point in plot
df_insitu_count <- df_insitu %>%
                    group_by(LatitudeSS_dd, LongitudeSS_dd, Intervention.category.itra.multi) %>%
                    summarise(count = n())

# convert to sf object for plotting
df_insitu_count_sf <- st_as_sf(df_insitu_count, coords = c("LongitudeSS_dd", "LatitudeSS_dd"), crs = 4326)

# create Mollweide projection
insitu_projected <- st_transform(df_insitu_count_sf, crs = mollweide_crs)


################################################################################
################################ Plot map ######################################


# plot richness
rich_plot <- ggplot() +
  geom_sf(data = grid_richness %>%
            mutate(SpeciesRichness = ifelse(SpeciesRichness == 0, NA, SpeciesRichness)),
          aes(fill = SpeciesRichness),
          colour = NA) +
  scale_fill_continuous(type = "gradient",
                        low=cols_heatmap[2], high = cols_heatmap[1], na.value = cols_heatmap[3],
                        limits = c(0, 14),
                        breaks = c(1, 4, 7, 10, 13)) +
  #geom_sf(data = coast_projected_fix, color='black', fill=NA, size=0.0001) +
  geom_sf(data = insitu_projected[which(insitu_projected$count != 1), ], shape = 19, aes(colour = factor(Intervention.category.itra.multi), size = count)) +
  geom_sf(data = insitu_projected[which(insitu_projected$count == 1),], shape = 4, size = 3, aes(colour = factor(Intervention.category.itra.multi))) +
  #geom_sf(data = insitu_projected, size = 1, shape = 19, aes(colour = factor(Intervention.category.1))) +
  scale_color_manual(values = colours_intcat) +
  scale_size(range = c(2, 4), breaks = c(2, 3, 4)) +
  labs(fill = "No. species", colour = "Intervention category", size = "No. interventions")

# save plot
ggsave(paste0(path_geo, "richness_plot_log.png"), rich_plot)


## end of script
