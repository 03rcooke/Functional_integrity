### Combine mammal data from each ecoregion

if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr)

# plyr: used to readadd filenames (ecoregion codes) to list of csvs # calls: ldply

csv_files <- dir(path = "~/ArcGIS/PhD/Ecoregion_CSVs", pattern = "csv$", full.names = TRUE)

empty_files <- csv_files[which(file.info(csv_files)$size<0.1)] # identify empty csvs

# function to add filename (ecoregion code) to new column named eco_code for each csv
read_csv_file <- function(filename){
  ret <- read.csv(filename, header = FALSE)
  ret$eco_code <- filename
  ret}

eco1 <- ldply(csv_files, read_csv_file) # runs function read_csv_file for each csv in list

eco1$eco_code <- basename(eco1$eco_code) # remove path name from ecoregion csv files
eco1$eco_code <- gsub("_M.csv", "", eco1$eco_code) # remove _M.csv from end of ecoregion csv files path

eco1 <- eco1[-c(1:2,31:32)] # drop geographic coordinates and empty columns

# Rename columns with shorter, lower case names - match column names in ArcGIS attribute tables
names(eco1) <- c("OBJECTid", "id_no", "binomial", "presence", "origin", "seasonal", "compiler", "year", "citation", "source", "dist_comm", "island", "subspecies", "subpop", "legend", "tax_comm", "kingdom_na", "phylum_nam", "class_name", "order_name", "family_nam", "genus_name", "species_na", "code", "shape_Leng", "Shape_Length", "Shape_Area", "eco_code")

