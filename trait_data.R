## -------------------------------------------------------------------------------------
## Name: trait_data.R
## Description: Code to set up trait data for all mammal species, 
##              identifies taxonomic matches and mismatches for 
##              species data (IUCN) compared to the trait data 
##              (PanTHERIA, Amniote, EltonTraits, MammalDIET)
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: May 2016 - 
## Outputs: 
## -------------------------------------------------------------------------------------

# Set up required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, taxize, data.table, rlist, stats, tidyr)

# dplyr: used to compare two data frames, combine species names # calls: anti_join, mutate
# taxize: used to find taxonomic synonyms and subspecies for taxonomic mismatches # calls: synonyms
# data.table: used to concantenate a list of data frames # calls: rbindlist
# rlist: used to remove data frames from a list of data frames # calls: list.remove
# stats: used to aggregate synonyms # calls: aggregate
# tidyr: used to separate synonyms in to two columns # calls: separate

## Read in species unique lists
am12u <- read.csv("ALL_Mammals_1_2_unique.csv", stringsAsFactors = FALSE)
# All mammals for presence = extant; origin = native, reintroduced
# nrow = 5233

am1235u <- read.csv("ALL_Mammals_1_2_3_5_unique.csv", stringsAsFactors = FALSE)
# All mammals for presence = extant; origin = native, reintroduced, introduced, origin uncertain
# nrow = 5235

am12uT <- am1235u[am1235u$origin == 1 | am1235u$origin == 2,]

## SUBSET ALL MAMMALS BY ORIGIN FROM ORIGINAL EXCEL 
## OR JUST CARRY ORIGIN THROUGH TO THE END AND THEN SUBSET


### Perform code for each trait data set: PanTHERIA, Amniote, EltonTraits, MammalDIET

######## PanTHERIA ##############

## read in PanTHERIA trait database
pan <- read.csv("PanTHERIA_1-0_WR05_Aug2008.csv", stringsAsFactors = FALSE)
# nrow = 5416

colnames(pan)[5] <- "binomial"
# set column MSW05_Binomial to binomial to match species tables

pan <- arrange(pan, binomial)
# order data by binomial A-Z

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
IUCN_empty12P <- data.frame(matrix(NA, nrow = nrow(IUCN12P), ncol = ncol(PanTHERIA12)))
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
# Probably not needed - just inner join instead

hdg <- full_join(am1235u, pan, by = "binomial")

pan1235T <- pan[!(pan$binomial %in% c(as.vector(PanTHERIA1235$binomial))),]
# removed species listed by PanTHERIA but not listed by IUCN
# nrow = 4804

# create empty dataframe to match trait database columns
IUCN_empty1235P <- data.frame(matrix(NA, nrow = nrow(IUCN1235P), ncol = ncol(PanTHERIA1235)))
colnames(IUCN_empty1235P) <- colnames(PanTHERIA1235)
IUCN_empty1235P$binomial <- IUCN1235P$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by PanTHERIA
pan1235 <- rbind(pan1235T, IUCN_empty1235P)

#write.table(PanTHERIA1235, "~/R/Functional_integrity/PanTHERIA1235.csv", sep=",") # Species listed by PanTHERIA but not listed by IUCN 1235
#write.table(PanTHERIA12, "~/R/Functional_integrity/PanTHERIA12.csv", sep=",") # Species listed by PanTHERIA but not listed by IUCN 12
#write.table(IUCN1235P, "~/R/Functional_integrity/IUCN1235P.csv", sep=",") # Species listed by IUCN 1235 but not listed by PanTHERIA
#write.table(IUCN12P, "~/R/Functional_integrity/IUCN12P.csv", sep=",") # Species listed by IUCN 12 but not listed by PanTHERIA
#write.table(pan1235, "~/R/Functional_integrity/ALL_Mammals_1235_traits_PanTHERIA.csv", sep=",") # PanTHERIA trait data for all IUCN 1235 species
#write.table(pan12, "~/R/Functional_integrity/ALL_Mammals_12_traits_PanTHERIA.csv", sep=",") # PanTHERIA trait data for all IUCN 12 species

# Species different between 12 and 1235
anti_join(PanTHERIA12, PanTHERIA1235, by = "binomial")

# Mus musculus - look at IUCN map

# Acomys nesiotes - Acomys nesiotes is endemic to Cyprus, where available evidence 
# suggests that it was possibly introduced by humans, and therefore may represent a 
# non-native population of Acomys cahirinus.

## Find and try synonyms for trait data ##

spp1235P <- IUCN1235P$binomial[1:10] # species list to find synonyms for: species listed by IUCN1235 but not listed by PanTHERIA

syn <- synonyms(spp1235P, db = "itis") # find synonyms - uses taxize package
no_syn_na <- which(is.na(syn) == TRUE); no_syn_na <- names(no_syn_na) # list those species that have no synonym data
no_syn_col <- which(sapply(syn, NCOL) == 3); no_syn_col <- names(no_syn_col) # list those species that have data but no synonyms 
no_syn <- c(no_syn_na, no_syn_col)

spp1235P_syn <- setdiff(spp1235P, no_syn) # remove species from species list that have no synonyms

syndf <- list.remove(syn, no_syn) # remove species data frames from synonym list that are empty - uses rlist package

Syn1235P <- as.data.frame(data.table::rbindlist(syndf)) # collapse list of data frames for each species into new single data frame

seq_bi <- sapply(syndf, nrow); seq_bi <- unlist(seq_bi) # number of rows per species
seq_bi <- rep(spp1235P_syn, times = as.vector(seq_bi)) # repeat species names with synonyms by seq_bi
Syn1235P <- mutate(Syn1235P, binomial_IUCN = seq_bi)  # add sequence of species names to synonym data frame

acc_tsn_e <- ifelse(is.na(as.numeric(Syn1235P$acc_tsn)), Syn1235P$syn_name, Syn1235P$acc_tsn) # swap accepted tsn with syn names where needed
Syn1235P <- mutate(Syn1235P, acc_tsn_e = acc_tsn_e) # add edited acc_tsn column to data frame

syn_name_e <- ifelse(!is.na(as.numeric(Syn1235P$syn_name)), Syn1235P$acc_tsn, Syn1235P$syn_name) # swap syn names with accepted tsn where needed
Syn1235P <- mutate(Syn1235P, syn_name_e = syn_name_e) # add edited syn_name column to data frame

Syn1235P <- aggregate(Syn1235P, by = list(Syn1235P$binomial_IUCN), FUN = unique) # collapse species with more than one synonym into one row - uses stats package
Syn1235P <- separate(Syn1235P, syn_name_e, c("syn_name1", "syn_name2"), ", ", extra = "merge", fill = "right") # split multiple synonyms into separate columns - uses tidyr package

Syn1235P$syn_name1 <- gsub("c(", "", Syn1235P$syn_name1, fixed = TRUE) # remove leading "c(" introduced by separate
Syn1235P$syn_name1 <- gsub("[^[:alnum:][:space:]]", "", Syn1235P$syn_name1) # remove all punctuation introduced by separate
Syn1235P$syn_name2 <- gsub("[^[:alnum:][:space:]]", "", Syn1235P$syn_name2) # remove all punctuation introduced by separate

Syn1235P <- mutate(Syn1235P, syn_name3 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , Syn1235P$syn_name1))) # add column of trinomials reduced to binomials based on first synonyms - gsub removes last word of name and any spaces produced
Syn1235P <- mutate(Syn1235P, syn_name4 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , Syn1235P$syn_name2))) # add column of trinomial names reduced to binomial names based on second synonyms

Syn1235P <- Syn1235P[c(6,2,7:11)] # reorder columns to move binomial from the end and drop unneeded columns

## Synonyms

syn_a <- inner_join(Syn1235P, pan, by = c("syn_name1" = "binomial"))
# Species matched in trait database based on first synonyms

Syn1235P <- Syn1235P[!(Syn1235P$binomial_IUCN %in% c(as.vector(syn_a$binomial_IUCN))),]
# Remove species matched by first synonym - do not need to try second synonym

syn_b <- inner_join(Syn1235P, pan, by = c("syn_name2" = "binomial"))
# Species matched in trait database based on second synonyms

## Subspecies

Syn1235P <- Syn1235P[!(Syn1235P$binomial_IUCN %in% c(as.vector(syn_b$binomial_IUCN))),]
# Remove species matched by second synonym - do not need to try third synonym (subspecies)

syn_c <- inner_join(Syn1235P, pan, by = c("syn_name3" = "binomial"))
# Species matched in trait database based on third synonym (subspecies)

Syn1235P <- Syn1235P[!(Syn1235P$binomial_IUCN %in% c(as.vector(syn_c$binomial_IUCN))),]
# Remove species matched by third synonym (subspecies) - do not need to try fourth synonym (subspecies)

syn_d <- inner_join(Syn1235P, pan, by = c("syn_name4" = "binomial"))
# Species matched in trait database based on fourth synonym (subspecies)

syn_comb <- rbind(syn_a, syn_b, syn_c, syn_d); setnames(syn_comb, "binomial_IUCN", "binomial") # Combine trait data matched by synonyms and subspecies
syn_comb <- syn_comb[c(1,8:ncol(syn_comb),4:7)] # drop unneeded columns: sub_tsn and acc_tsn_e and reorder to put syn_names at the end
syn_comb <- inner_join(syn_comb, hdg, by = "binomial") # join trait data with species data for species matched by synonyms
syn_comb <- syn_comb[,!apply(syn_comb, 2, function(x) all(is.na(x)))] # remove columns of all NAs generated during inner_join
colnames(syn_comb) <- gsub(".x", "", colnames(syn_comb), fixed = TRUE) # remove .x from column names generated during inner_join

hdg[,c("syn_name1","syn_name2","syn_name3","syn_name4")] <- NA # add syn_names to hdg data frame

IUCN_syn <- anti_join(Syn1235P, pan, by = "binomial"); IUCN_syn <- arrange(IUCN_syn, binomial_IUCN) # order data by binomial A-Z
# IUCN species not listed in synonyms, i.e. species that need further matching efforts

hdg <- hdg[!(hdg$binomial %in% c(as.vector(syn_comb$binomial))),]
# Remove species with NA values that have now been matched by synonyms (subspecies)

top <- rbind(hdg, syn_comb)
# Combine matched species from synonyms with matched species from IUCN and trait database

########### Amniote - mammals ####################

## read in Amniote trait database
amn <- read.csv("Amniote_Database_Aug_2015.csv")
# nrow = 21322

# create new column with full scientific name
amn <- mutate(amn, binomial = paste(genus, species, sep = " ")) # name as binomial to match species tables
amn <- amn[c(1:7,37,8:36)] # reorder columns to move binomial from the end

# subset just for mammals
amn <- amn[amn$class == "Mammalia",]
# nrow = 4953

amn <- arrange(amn, binomial)
# order data by binomial A-Z

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
IUCN_empty12A <- data.frame(matrix(NA, nrow = nrow(IUCN12A), ncol = ncol(Amniote12)))
colnames(IUCN_empty12A) <- colnames(Amniote12)
IUCN_empty12A$binomial <- IUCN12A$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by Amniote
amn12 <- rbind(amn12, IUCN_empty12A)

## 1235 All mammals for presence = extant; origin = native, reintroduced, introduced, origin ##

Amniote1235 <- anti_join(amn, am1235u, by = "binomial")
# Species listed by Amniote but not listed by IUCN
# nrow = 512
IUCN1235A <- anti_join(am1235u, amn, by = "binomial")
# IUCN species not listed in Amniote
# nrow = 794

amn1235 <- amn[!(amn$binomial %in% c(as.vector(Amniote1235$binomial))),]
# removed species listed by Amniote but not listed by IUCN
# nrow = 4441

# create empty dataframe to match trait database columns
IUCN_empty1235A <- data.frame(matrix(NA, nrow = nrow(IUCN1235A), ncol = ncol(Amniote1235)))
colnames(IUCN_empty1235A) <- colnames(Amniote1235)
IUCN_empty1235A$binomial <- IUCN1235A$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by Amniote
amn1235 <- rbind(amn1235, IUCN_empty1235A)

#write.table(Amniote1235, "~/R/Functional_integrity/Amniote1235_M.csv", sep=",") # Species listed by Amniote[Mammals] but not listed by IUCN 1235
#write.table(Amniote12, "~/R/Functional_integrity/Amniote12_M.csv", sep=",") # Species listed by Amniote[Mammals] but not listed by IUCN 12
#write.table(IUCN1235A, "~/R/Functional_integrity/IUCN1235A_M.csv", sep=",") # Species listed by IUCN 1235 but not listed by Amniote[Mammals]
#write.table(IUCN12A, "~/R/Functional_integrity/IUCN12A_M.csv", sep=",") # Species listed by IUCN 12 but not listed by Amniote[Mammals]
#write.table(amn1235, "~/R/Functional_integrity/ALL_Mammals_1235_traits_Amniote.csv", sep=",") # Amniote[Mammals] trait data for all IUCN 1235 species
#write.table(amn12, "~/R/Functional_integrity/ALL_Mammals_12_traits_Amniote.csv", sep=",") # Amniote[Mammals] trait data for all IUCN 12 species







########### EltonTraits 1.0  - mammals ####################

## read in Elton Traits 1.0 database
et <- read.csv("MamFuncDat.csv")
# nrow = 5400

colnames(et)[3] <- "binomial"
# set column Scientific to binomial to match species tables

et <- arrange(et, binomial)
# order data by binomial A-Z

## 12 All mammals for presence = extant; origin = native, reintroduced ##

Elton12 <- anti_join(et, am12u, by = "binomial")
# Species listed by Elton but not listed by IUCN
# nrow = 599
IUCN12E <- anti_join(am12u, et, by = "binomial")
# IUCN species not listed in Elton
# nrow = 432

et12 <- et[!(et$binomial %in% c(as.vector(Elton12$binomial))),]
# removed species listed by Elton but not listed by IUCN
# nrow = 4801

# create empty dataframe to match trait database columns
IUCN_empty12E <- data.frame(matrix(NA, nrow = nrow(IUCN12E), ncol = ncol(Elton12)))
colnames(IUCN_empty12E) <- colnames(Elton12)
IUCN_empty12E$binomial <- IUCN12E$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by EltonTraits
et12 <- rbind(et12, IUCN_empty12E)

## 1235 All mammals for presence = extant; origin = native, reintroduced, introduced, origin ##

Elton1235 <- anti_join(et, am1235u, by = "binomial")
# Species listed by Elton but not listed by IUCN
# nrow = 597
IUCN1235E <- anti_join(am1235u, et, by = "binomial")
# IUCN species not listed in Elton
# nrow = 432

et1235 <- et[!(et$binomial %in% c(as.vector(Elton1235$binomial))),]
# removed species listed by Elton but not listed by IUCN
# nrow = 4803

# create empty dataframe to match trait database columns
IUCN_empty1235E <- data.frame(matrix(NA, nrow = nrow(IUCN1235E), ncol = ncol(Elton1235)))
colnames(IUCN_empty1235E) <- colnames(Elton1235)
IUCN_empty1235E$binomial <- IUCN1235E$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by EltonTraits
et1235 <- rbind(et1235, IUCN_empty1235E)

#write.table(Elton1235, "~/R/Functional_integrity/Elton1235_M.csv", sep=",") # Species listed by EltonTraits[Mammals] but not listed by IUCN 1235
#write.table(Elton12, "~/R/Functional_integrity/Elton12_M.csv", sep=",") # Species listed by EltonTraits[Mammals] but not listed by IUCN 12
#write.table(IUCN1235E, "~/R/Functional_integrity/IUCN1235E_M.csv", sep=",") # Species listed by IUCN 1235 but not listed by EltonTraits[Mammals]
#write.table(IUCN12E, "~/R/Functional_integrity/IUCN12E_M.csv", sep=",") # Species listed by IUCN 12 but not listed by EltonTraits[Mammals]
#write.table(et1235, "~/R/Functional_integrity/ALL_Mammals_1235_traits_Elton.csv", sep=",") # EltonTraits[Mammals] trait data for all IUCN 1235 species
#write.table(et12, "~/R/Functional_integrity/ALL_Mammals_12_traits_Elton.csv", sep=",") # EltonTraits[Mammals] trait data for all IUCN 12 species




########### MammalDiet ####################

## read in MammalDiet 1.0 trait database
md <- read.csv("MammalDIET_V1.0.csv")
# nrow = 5364

# create new column with full scientific name
md <- mutate(md, binomial = paste(Genus, Species, sep = " ")) # name as binomial to match species tables
md <- md[c(1:5,31,6:30)] # reorder columns to move binomial from the end

md <- arrange(md, binomial)
# order data by binomial A-Z

## 12 All mammals for presence = extant; origin = native, reintroduced ##

MDIET12 <- anti_join(md, am12u, by = "binomial")
# Species listed by MammalDIET but not listed by IUCN
# nrow = 246
IUCN12M <- anti_join(am12u, md, by = "binomial")
# IUCN species not listed in MammalDIET
# nrow = 115

md12 <- md[!(md$binomial %in% c(as.vector(MDIET12$binomial))),]
# removed species listed by MammalDIET but not listed by IUCN
# nrow = 5118

# create empty dataframe to match trait database columns
IUCN_empty12M <- data.frame(matrix(NA, nrow = nrow(IUCN12M), ncol = ncol(MDIET12)))
colnames(IUCN_empty12M) <- colnames(MDIET12)
IUCN_empty12M$binomial <- IUCN12M$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by MammalDIET
md12 <- rbind(md12, IUCN_empty12M)

## 1235 All mammals for presence = extant; origin = native, reintroduced, introduced, origin ##

MDIET1235 <- anti_join(md, am1235u, by = "binomial")
# Species listed by MammalDIET but not listed by IUCN
# nrow = 244
IUCN1235M <- anti_join(am1235u, md, by = "binomial")
# IUCN species not listed in MammalDIET
# nrow = 115

md1235 <- md[!(md$binomial %in% c(as.vector(MDIET1235$binomial))),]
# removed species listed by MammalDIET but not listed by IUCN
# nrow = 5120

# create empty dataframe to match trait database columns
IUCN_empty1235M <- data.frame(matrix(NA, nrow = nrow(IUCN1235M), ncol = ncol(MDIET1235)))
colnames(IUCN_empty1235M) <- colnames(MDIET1235)
IUCN_empty1235M$binomial <- IUCN1235M$binomial

# combine trait database with empty values for species listed by the IUCN but not listed by MammalDIET
md1235 <- rbind(md1235, IUCN_empty1235M)

#write.table(MDIET1235, "~/R/Functional_integrity/MDIET1235.csv", sep=",") # Species listed by MammalDIET but not listed by IUCN 1235
#write.table(MDIET12, "~/R/Functional_integrity/MDIET12.csv", sep=",") # Species listed by MammalDIET but not listed by IUCN 12
#write.table(IUCN1235M, "~/R/Functional_integrity/IUCN1235M.csv", sep=",") # Species listed by IUCN 1235 but not listed by MammalDIET
#write.table(IUCN12M, "~/R/Functional_integrity/IUCN12M.csv", sep=",") # Species listed by IUCN 12 but not listed by MammalDIET
#write.table(md1235, "~/R/Functional_integrity/ALL_Mammals_1235_traits_MammalDIET.csv", sep=",") # MammalDIET trait data for all IUCN 1235 species
#write.table(md12, "~/R/Functional_integrity/ALL_Mammals_12_traits_MammalDIET.csv", sep=",") # MammalDIET trait data for all IUCN 12 species


