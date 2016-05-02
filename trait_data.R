## Set up trait data for all mammal species

## read in species unique lists

am12u <- read.csv("ALL_Mammals_1_2_unique.csv")
# All mammals for presence = extant; origin = native, reintroduced

am1235u <- read.csv("ALL_Mammals_1_2_3_5_unique.csv")
# All mammals for presence = extant; origin = native, reintroduced, introduced, origin uncertain
# length = 5235

## read in PanTHERIA trait database
pan <- read.csv("PanTHERIA_1-0_WR05_Aug2008.csv")
# length = 5416

colnames(pan)[5] <- "binomial"
# set column MSW05_Binomial to binomial to match other tables

PanTHERIA <- anti_join(pan, am1235u, by = "binomial")
# Species listed by PanTHERIA but not listed by IUCN
# length = 612
IUCN <- anti_join(am1235u, pan, by = "binomial")
# IUCN species not listed in PanTHERIA
# length = 431

pan2 <- pan[!(pan$binomial %in% c(PanTHERIA)),]
# removed species listed by PanTHERIA but not listed by IUCN
# length = 4804

# create empty dataframe to match trait database columns
IUCN_empty <- data.frame(matrix(NA, nrow = length(IUCN$binomial), ncol = ncol(PanTHERIA)))
colnames(IUCN_empty) <- colnames(PanTHERIA)
IUCN_empty$binomial <- IUCN$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by PanTHERIA
pan2 <- rbind(pan2, IUCN_empty)

#write.table(PanTHERIA, "~/R/Functional_integrity/PanTHERIA.csv", sep=",")
#write.table(IUCN, "~/R/Functional_integrity/IUCN.csv", sep=",")
#write.table(pan2, "~/R/Functional_integrity/ALL_mammals_1235_traits.csv", sep=",")