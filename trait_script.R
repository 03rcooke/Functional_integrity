## Read in species unique lists
am12u <- read.csv("ALL_Mammals_1_2_unique.csv")
# All mammals for presence = extant; origin = native, reintroduced
# nrow = 5233

am1235u <- read.csv("ALL_Mammals_1_2_3_5_unique.csv")
# All mammals for presence = extant; origin = native, reintroduced, introduced, origin uncertain
# nrow = 5235

am12uT <- am1235u[am1235u$origin == 1 | am1235u$origin == 2,]

## SUBSET ALL MAMMALS BY ORIGIN FROM ORIGINAL EXCEL 
## OR JUST CARRY ORIGIN THROUGH TO THE END AND THEN SUBSET

### Perform code for each trait data set: PanTHERIA, Amniote, EltonTraits, MammalDIET

######## PanTHERIA ##############

## read in PanTHERIA trait database
pan <- read.csv("PanTHERIA_1-0_WR05_Aug2008.csv")
# nrow = 5416

colnames(pan)[5] <- "binomial"
# set column MSW05_Binomial to binomial to match species tables

pan <- arrange(pan, binomial)
# order data by binomial A-Z

trait_fun <- trait(trait = pan, species = am1235u)

trait_database <- trait_fun$trait_database
IUCN <- trait_fun$IUCN
trait_out <- trait_fun$trait_out
trait_out_final <- trait_fun$trait_out_final


#subset data after function to remove species with just IUCN35