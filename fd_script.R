## --------------------------------------------------------------
## Name: fd_script.R
## Description: Code to run "fd_eco" function
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: March 2016 - 
## Inputs: trait data from function 'trait'
##         species data per ecoregion from eco_species.R script
## Outputs: 
## --------------------------------------------------------------

# check working directory
getwd()

# Set up required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, robustbase)

# dplyr: used to select columns, order data frames, create new columns # calls: select, arrange, mutate
# robustbase: used to calcultae row medians # calls: rowMedians

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

UK_trait <- read.csv("Trait_data_UK.csv", row.names = 2,
                     # add species names to rows
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
# assign simple names to variables

trait_P <- readRDS("trait_PanTHERIA.rds")
trait_A <- readRDS("trait_Amniote.rds")
trait_E <- readRDS("trait_Elton.rds")
trait_M <- readRDS("trait_MD.rds")

big <- Reduce(function(x, y) inner_join(x, y, by = c("id_no", "binomial", "presence", "origin", "shape_Area")), list(trait_P$trait_data, trait_A$trait_data, trait_E$trait_data, trait_M$trait_data)) # join trait data from all databases
big <- arrange(big, binomial) # order data by binomial A-Z
big <- big[,!duplicated(colnames(big))] # remove duplicated syn_name columns produced during inner_join
# nrow = 5235

big_IUCN_names <- Reduce(function(x, y) full_join(x, y, by = c("binomial_syn")), list(trait_P$IUCN_names, trait_A$IUCN_names, trait_E$IUCN_names, trait_M$IUCN_names)) # join unmatched names from all databases
# nrow = 673

### Taxonomic data frame ####

taxonomic <- select(big, id_no:origin, MSW05_Genus, MSW05_Species, genus, species, Scientific, Genus, Species) # select columns id_no to origin and then taxonomic related columns
# create data frame of taxonomic information across the databases to compare taxonomies used

# create new column with full scientific name
taxonomic <- mutate(taxonomic, binomial_PanTHERIA = paste(MSW05_Genus, MSW05_Species, sep = " "))
taxonomic <- mutate(taxonomic, binomial_Amniote = paste(genus, species, sep = " "))
taxonomic <- mutate(taxonomic, binomial_MD = paste(Genus, Species, sep = " "))
setnames(taxonomic, "Scientific", "binomial_Elton")

taxonomic <- select(taxonomic, id_no:origin, binomial_PanTHERIA, binomial_Amniote, binomial_MD, binomial_Elton)

### Trait data frames ####

body_mass <- select(big, id_no:origin, X5.1_AdultBodyMass_g, adult_body_mass_g, BodyMass.Value) # create data frame of body mass data to compare across databases
names(body_mass) <- c("id_no", "binomial", "presence", "origin", "body_mass_PanTHERIA", "body_mass_Amniote", "body_mass_Elton")
body_mass[body_mass==-999] <- NA # turn -999 to NAs
body_mass <- mutate(body_mass, body_mass_median = rowMedians(as.matrix(select(body_mass, starts_with("body_mass_"))), na.rm = TRUE))


litter_size <- select(big, id_no:origin, X15.1_LitterSize, litter_or_clutch_size_n) # create data frame of litter size data to compare across databases
names(litter_size) <- c("id_no", "binomial", "presence", "origin", "litter_size_PanTHERIA", "litter_size_Amniote")
litter_size[litter_size==-999] <- NA # turn -999 to NAs

## eco_trait

### FD function ####
out_fd <- fd_eco(UK_data, UK_trait, corr = "cailliez", spp_list = FALSE, tree = FALSE)

out_fd <- fd_eco(eco_trait, eco1, corr = "cailliez", spp_list = FALSE, tree = FALSE)

out_fd

#### Statistics
out_fd@stats

out_fd@stats$ecoregions # number of ecoregions assessed
out_fd@stats$spp_total # total number of species per ecoregion (including species with missing trait data)
out_fd@stats$spp_missing # number of species missing all trait data per ecoregion
out_fd@stats$spp_final # number of species per ecoregion after removing missing-data species
out_fd@stats$k # number of functional clusters as specified by the L method (Salvador & Chan, 2004)
out_fd@stats$CWM # community weighted mean
out_fd@stats$FRed # functional redundancy
out_fd@stats$FDis # functional dispersion

#### Plots
out_fd@plots

out_fd@plots$clus # Performance of clustering algorithms
out_fd@plots$eval # L method evaluation plot
out_fd@plots$dendro # Functional dendrogram: following best clustering algorithm and number of clusters as determined by L method

#### Session info
out_fd@session_info