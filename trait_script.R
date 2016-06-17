## --------------------------------------------------------------
## Name: trait_script.R
## Description: Code to run "trait" function
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: June 2016 - 
## Inputs: species data (list of species)
##         trait data (trait databases) set up with the species name column identifed as "binomial"
## Outputs: 
## --------------------------------------------------------------

# Set up required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr)

# dplyr: used to order data frames and duplicate binomial column # calls: arrange, mutate

## Read in species unique lists
IUCN_species <- read.csv("ALL_Mammals_1_2_3_5_unique.csv", stringsAsFactors = FALSE)
# All mammals for presence = extant; origin = native, reintroduced, introduced, origin uncertain
# nrow = 5235

### Perform code for each trait MAMMAL database: PanTHERIA, Amniote, EltonTraits, MammalDIET

######## PanTHERIA ##############

## read in PanTHERIA trait database
pan <- read.csv("PanTHERIA_1-0_WR05_Aug2008.csv", stringsAsFactors = FALSE)
# nrow = 5416

colnames(pan)[5] <- "binomial"
# set column MSW05_Binomial to binomial to match species data

pan <- arrange(pan, binomial)
# order data by binomial A-Z

# Have to name the necessary column binomial before using function

trait_PanTHERIA <- trait(trait = pan, species = IUCN_species)

IUCN_PanTHERIA <- trait_PanTHERIA$IUCN_names # IUCN species not matched by synonyms, i.e. species that need further matching efforts
trait_data_PanTHERIA <- trait_PanTHERIA$trait_data # Trait data for all IUCN species after trying synonyms

nrow(IUCN_PanTHERIA) # Number of species not matched to trait data
# 263
nrow(trait_data_PanTHERIA) - nrow(IUCN_PanTHERIA) # Number of species matched to trait data
# 4972
(nrow(trait_data_PanTHERIA) - nrow(IUCN_PanTHERIA))/nrow(trait_data_PanTHERIA)*100 # Percent of species macthed to trait data
# 95% matched

#saveRDS(trait_PanTHERIA, "trait_PanTHERIA.rds")
# t <- readRDS("trait_PanTHERIA.rds") # How to load the file under a different nameif needed

# export data as csvs
#write.csv(IUCN_PanTHERIA, "~/R/Functional_integrity/Missing_IUCN_PanTHERIA.csv")
#write.csv(trait_data_PanTHERIA, "~/R/Functional_integrity/Trait_PanTHERIA.csv")

#trait_data_PanTHERIA[trait_data_PanTHERIA==-999] <- NA # turn -999 to NAs
#apply(trait_data_PanTHERIA, 2, function(x) length(which(!is.na(x)))) # count per column number of values
# X5.1_AdultBodyMass_g = 3329 (64%) # % of IUCN species with data
# X15.1_LitterSize = 2366 (45%)

#subset data after function to remove species with just IUCN35

########### Amniote - mammals ####################

## read in Amniote trait database
amn <- read.csv("Amniote_Database_Aug_2015.csv", stringsAsFactors = FALSE)
# nrow = 21322

# create new column with full scientific name
amn <- mutate(amn, binomial = paste(genus, species, sep = " ")) # name as binomial to match species data
amn <- amn[c(1:7,37,8:36)] # reorder columns to move binomial from the end

# subset just for mammals
amn <- amn[amn$class == "Mammalia",]
# nrow = 4953

amn <- arrange(amn, binomial)
# order data by binomial A-Z

trait_Amniote <- trait(trait = amn, species = IUCN_species)

IUCN_Amniote <- trait_Amniote$IUCN_names # IUCN species not matched by synonyms, i.e. species that need further matching efforts
trait_data_Amniote <- trait_Amniote$trait_data # Trait data for all IUCN species after trying synonyms

nrow(IUCN_Amniote) # Number of IUCN species not matched to trait data
# 631
nrow(trait_data_Amniote) - nrow(IUCN_Amniote) # Number of species matched to trait data
# 4604
(nrow(trait_data_Amniote) - nrow(IUCN_Amniote))/nrow(trait_data_Amniote)*100 # Percent of species macthed to trait data
# 88% matched

#saveRDS(trait_Amniote, "trait_Amniote.rds")

# export data as csvs
#write.csv(IUCN_Amniote, "~/R/Functional_integrity/Missing_IUCN_Amniote.csv")
#write.csv(trait_data_Amniote, "~/R/Functional_integrity/Trait_Amniote.csv")

#trait_data_Amniote[trait_data_Amniote==-999] <- NA # turn -999 to NAs
#apply(trait_data_Amniote, 2, function(x) length(which(!is.na(x)))) # count per column number of values
# adult_body_mass_g = 4323 (83%)
# litter_or_clutch_size_n = 3262 (62%)

########### EltonTraits 1.0  - mammals ####################

## read in Elton Traits 1.0 database
et <- read.csv("MamFuncDat.csv", stringsAsFactors = FALSE)
# nrow = 5400

et <- mutate(et, binomial = Scientific)
# duplicate scientific column and name it binomial to match species data

et <- arrange(et, binomial)
# order data by binomial A-Z

trait_Elton <- trait(trait = et, species = IUCN_species)

IUCN_Elton <- trait_Elton$IUCN_names # IUCN species not matched by synonyms, i.e. species that need further matching efforts
trait_data_Elton <- trait_Elton$trait_data # Trait data for all IUCN species after trying synonyms

nrow(IUCN_Elton) # Number of species not matched to trait data
# 264
nrow(trait_data_Elton) - nrow(IUCN_Elton) # Number of species matched to trait data
# 4971
(nrow(trait_data_Elton) - nrow(IUCN_Elton))/nrow(trait_data_Elton)*100 # Percent of species macthed to trait data
# 95% matched

#saveRDS(trait_Elton, "trait_Elton.rds")

# export data as csvs
#write.csv(IUCN_Elton, "~/R/Functional_integrity/Missing_IUCN_Elton.csv")
#write.csv(trait_data_Elton, "~/R/Functional_integrity/Trait_Elton.csv")

#trait_data_Elton[trait_data_Elton==-999] <- NA # turn -999 to NAs
#apply(trait_data_Elton, 2, function(x) length(which(!is.na(x)))) # count per column number of values
# BodyMass.Value = 4971 (95%) # some interpolated
# Activity.Nocturnal = 4971 (95%)
# Diet = 4971 (95%)

########### MammalDiet ####################

## read in MammalDiet 1.0 trait database
md <- read.csv("MammalDIET_V1.0.csv")
# nrow = 5364

# create new column with full scientific name
md <- mutate(md, binomial = paste(Genus, Species, sep = " ")) # name as binomial to match species data
md <- md[c(1:5,31,6:30)] # reorder columns to move binomial from the end

md <- arrange(md, binomial)
# order data by binomial A-Z

trait_MD <- trait(trait = md, species = IUCN_species)

IUCN_MD <- trait_MD$IUCN_names # IUCN species not matched by synonyms, i.e. species that need further matching efforts
trait_data_MD <- trait_MD$trait_data # Trait data for all IUCN species after trying synonyms

nrow(IUCN_MD) # Number of species not matched to trait data
# 89
nrow(trait_data_MD) - nrow(IUCN_MD) # Number of species matched to trait data
# 5146
(nrow(trait_data_MD) - nrow(IUCN_MD))/nrow(trait_data_MD)*100 # Percent of species macthed to trait data
# 98% matched

#saveRDS(trait_MD, "trait_MD.rds")

# export data as csvs
#write.csv(IUCN_MD, "~/R/Functional_integrity/Missing_IUCN_MD.csv")
#write.csv(trait_data_MD, "~/R/Functional_integrity/Trait_MD.csv")

#trait_data_MD[trait_data_MD==-999] <- NA # turn -999 to NAs
#apply(trait_data_MD, 2, function(x) length(which(!is.na(x)))) # count per column number of values
# Diet = 5130-5146 (98%)

