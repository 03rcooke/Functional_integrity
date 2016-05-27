## Set up trait data for all mammal species

if (!require("pacman")) install.packages("pacman")
pacman::p_load(fossil, data.table, FD, clue, dplyr, qpcR, stats, cowplot, ggdendro)

# dplyr: used to compare two data frames, combine species names # calls: anti_join, mutate

## read in species unique lists

am12u <- read.csv("ALL_Mammals_1_2_unique.csv")
# All mammals for presence = extant; origin = native, reintroduced
# length = 5233

am1235u <- read.csv("ALL_Mammals_1_2_3_5_unique.csv")
# All mammals for presence = extant; origin = native, reintroduced, introduced, origin uncertain
# length = 5235

######## PanTHERIA ##############

## read in PanTHERIA trait database
pan <- read.csv("PanTHERIA_1-0_WR05_Aug2008.csv")
# nrow = 5416

colnames(pan)[5] <- "binomial"
# set column MSW05_Binomial to binomial to match species tables

## 12 All mammals for presence = extant; origin = native, reintroduced ##

PanTHERIA12 <- anti_join(pan, am12u, by = "binomial") # Fine to ignore Warning message: In anti_join_impl(x, y, by$x, by$y) : joining factors with different levels, coercing to character vector
# Species listed by PanTHERIA but not listed by IUCN
# nrow = 614
IUCN12P <- anti_join(am12u, pan, by = "binomial")
# IUCN species not listed in PanTHERIA
# nrow = 431

pan12 <- pan[!(pan$binomial %in% c(as.vector(PanTHERIA12$binomial))),]
# removed species listed by PanTHERIA but not listed by IUCN
# nrow = 4802

# create empty dataframe to match trait database columns
IUCN_empty12P <- data.frame(matrix(NA, nrow = length(IUCN12P$binomial), ncol = ncol(PanTHERIA12)))
colnames(IUCN_empty12P) <- colnames(PanTHERIA12)
IUCN_empty12P$binomial <- IUCN12P$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by PanTHERIA
pan12 <- rbind(pan12, IUCN_empty12P)

## 1235 All mammals for presence = extant; origin = native, reintroduced, introduced, origin ##

PanTHERIA1235 <- anti_join(pan, am1235u, by = "binomial")
# Species listed by PanTHERIA but not listed by IUCN
# nrow = 612
IUCN1235P <- anti_join(am1235u, pan, by = "binomial")
# IUCN species not listed in PanTHERIA
# nrow = 431

pan1235 <- pan[!(pan$binomial %in% c(as.vector(PanTHERIA1235$binomial))),]
# removed species listed by PanTHERIA but not listed by IUCN
# nrow = 4804

# create empty dataframe to match trait database columns
IUCN_empty1235P <- data.frame(matrix(NA, nrow = length(IUCN1235P$binomial), ncol = ncol(PanTHERIA1235)))
colnames(IUCN_empty1235P) <- colnames(PanTHERIA1235)
IUCN_empty1235P$binomial <- IUCN1235P$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by PanTHERIA
pan1235 <- rbind(pan1235, IUCN_empty1235P)

#write.table(PanTHERIA1235, "~/R/Functional_integrity/PanTHERIA1235.csv", sep=",")
#write.table(PanTHERIA12, "~/R/Functional_integrity/PanTHERIA12.csv", sep=",")
#write.table(IUCN1235P, "~/R/Functional_integrity/IUCN1235P.csv", sep=",")
#write.table(IUCN12P, "~/R/Functional_integrity/IUCN12P.csv", sep=",")
#write.table(pan1235, "~/R/Functional_integrity/ALL_Mammals_1235_traits_PanTHERIA.csv", sep=",")
#write.table(pan12, "~/R/Functional_integrity/ALL_Mammals_12_traits_PanTHERIA.csv", sep=",")

anti_join(PanTHERIA12, PanTHERIA1235, by = "binomial")
# species different between 12 and 1235

# Mus musculus - look at IUCN map

# Acomys nesiotes - Acomys nesiotes is endemic to Cyprus, where available evidence 
# suggests that it was possibly introduced by humans, and therefore may represent a 
# non-native population of Acomys cahirinus.








########### Amniote - mammals ####################

## read in Amniote trait database
amn <- read.csv("Amniote_Database_Aug_2015.csv")
# nrow = 21322

# create new column with full scientific name
amn <- mutate(amn, binomial = paste(genus, species, sep = " ")) # name as binomial to match species tables
amn <- amn2[c(1:7,37,8:36)] # reorder to move binomial from the end

# subset just for mammals
amn <- amn[amn$class == "Mammalia",]
# nrow = 4953

## 12 All mammals for presence = extant; origin = native, reintroduced ##

Amniote12 <- anti_join(amn, am12u, by = "binomial")
# Species listed by Amniote but not listed by IUCN
# nrow = 514
IUCN12A <- anti_join(am12u, amn, by = "binomial")
# IUCN species not listed in Amniote
# nrow = 794

amn12 <- amn[!(amn$binomial %in% c(as.vector(Amniote12$binomial))),]
# removed species listed by Amniote but not listed by IUCN
# nrow = 4439

# create empty dataframe to match trait database columns
IUCN_empty12A <- data.frame(matrix(NA, nrow = length(IUCN12A$binomial), ncol = ncol(Amniote12)))
colnames(IUCN_empty12A) <- colnames(Amniote12)
IUCN_empty12A$binomial <- IUCN12A$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by Amniote
amn12 <- rbind(amn12, IUCN_empty12A)

## 1235 All mammals for presence = extant; origin = native, reintroduced, introduced, origin ##

Amniote1235 <- anti_join(amn, am1235u, by = "binomial")
# Species listed by Amniote but not listed by IUCN
# length = 512
IUCN1235A <- anti_join(am1235u, amn, by = "binomial")
# IUCN species not listed in Amniote
# length = 794

amn1235 <- amn[!(amn$binomial %in% c(as.vector(Amniote1235$binomial))),]
# removed species listed by Amniote but not listed by IUCN
# length = 4441

# create empty dataframe to match trait database columns
IUCN_empty1235A <- data.frame(matrix(NA, nrow = length(IUCN1235A$binomial), ncol = ncol(Amniote1235)))
colnames(IUCN_empty1235A) <- colnames(Amniote1235)
IUCN_empty1235A$binomial <- IUCN1235A$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by Amniote
amn1235 <- rbind(amn1235, IUCN_empty1235A)

#write.table(Amniote1235, "~/R/Functional_integrity/Amniote1235_M.csv", sep=",")
#write.table(Amniote12, "~/R/Functional_integrity/Amniote12_M.csv", sep=",")
#write.table(IUCN1235A, "~/R/Functional_integrity/IUCN1235A_M.csv", sep=",")
#write.table(IUCN12A, "~/R/Functional_integrity/IUCN12A_M.csv", sep=",")
#write.table(amn1235, "~/R/Functional_integrity/ALL_Mammals_1235_traits_Amniote.csv", sep=",")
#write.table(amn12, "~/R/Functional_integrity/ALL_Mammals_12_traits_Amniote.csv", sep=",")





