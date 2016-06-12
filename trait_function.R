trait <- function(trait, species, binomial) 
{
  # Set up required packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(dplyr, taxize, data.table, rlist)
  
  # dplyr: used to compare two data frames, combine species names # calls: anti_join, mutate
  # taxize: used to find taxonomic synonyms and subspecies for taxonomic mismatches # calls: synonyms
  # data.table: used to concantenate a list of data frames # calls: rbindlist
  # rlist: used to remove data frames from a list of data frames # calls: list.remove

  ## 1235 All mammals for presence = extant; origin = native, reintroduced, introduced, origin ##
  
  trait_database <- anti_join(trait, species, by = "binomial")
  # Species listed by trait database but not listed by IUCN
  # nrow = 612

  trait_out <- full_join(am1235u, pan, by = "binomial")
  # nrow = 5847
  
  trait_out_final <- trait_out[!(trait_out$binomial %in% c(as.vector(trait_database$binomial))),]
  # removed species listed by trait database but not listed by IUCN
  # nrow = 5235
  
  
  ## Find and try synonyms for trait data ##
  
  spp <- IUCN$binomial[1:10] # species list to find synonyms for: species listed by IUCN1235 but not listed by trait database
  
  syno <- synonyms(spp, db = "itis") # find synonyms - uses taxize package
  no_syn <- which(is.na(syno) == TRUE); no_syn <- names(no_syn) # list those that have no synonyms
  
  spp_syn <- setdiff(spp, no_syn) # remove species from species list that have no synonyms
  syndf <- list.remove(syno, no_syn) # remove species data frames from synonym list that are empty - uses rlist package
  
  Syn <- as.data.frame(data.table::rbindlist(syndf)) # collapse list of data frames for each species into new single data frame
  
  seq_bi <- sapply(syndf, nrow) %>% unlist(seq_bi) # number of rows per species
  seq_bi <- rep(spp_syn, times = as.vector(seq_bi)) # repeat species names with synonyms by seq_bi
  Syn <- mutate(Syn, binomial = seq_bi)  # add sequence of species names to synonym data frame
  
  acc_tsn1 <- ifelse(is.na(as.numeric(Syn$acc_tsn)), Syn$syn_name, Syn$acc_tsn) # swap accepted tsn with syn names where needed
  Syn <- mutate(Syn, acc_tsn1 = acc_tsn1) # add edited acc_tsn column to data frame
  
  syn_name1 <- ifelse(!is.na(as.numeric(Syn$syn_name)), Syn$acc_tsn, Syn$syn_name) # swap syn names with accepted tsn where needed
  Syn <- mutate(Syn, syn_name1 = syn_name1) # add edited syn_name column to data frame
  
  Syn <- Syn[c(5,1,6:7,4)] # reorder columns to move binomial from the end
  
  # if else chained loop to split data in to first, second, third and fourth synonyms depending on the species with the most synonyms in data frame
  
  if(max(as.vector(table(Syn$binomial))) == 1) {
    Syn_1 <- Syn[!duplicated(Syn$binomial),] # create data frame of first synonyms
  } else {
    if(max(as.vector(table(Syn$binomial))) == 2) {
      Syn_1 <- Syn[!duplicated(Syn$binomial),] # create data frame of first synonyms
      Syn_2 <- Syn[duplicated(Syn$binomial),] # create data frame of second synonyms
    } else {
      if(max(as.vector(table(Syn$binomial))) == 3) {
        Syn_1 <- Syn[!duplicated(Syn$binomial),] # create data frame of first synonyms
        Syn_T <- Syn[duplicated(Syn$binomial),] # create data frame of second synonyms
        Syn_2 <- Syn_T[!duplicated(Syn_T$binomial),] # create temporary dataframe
        Syn_3 <- Syn_T[duplicated(Syn_T$binomial),] # create data frame of third synonyms
        rm(Syn_T)
      } else {
        if(max(as.vector(table(Syn$binomial))) == 4) {
          Syn_1 <- Syn[!duplicated(Syn$binomial),] # create data frame of first synonyms
          Syn_T <- Syn[duplicated(Syn$binomial),] # create temporary dataframe
          Syn_2 <- Syn_T[!duplicated(Syn_T$binomial),] # create data frame of second synonyms
          Syn_T <- Syn_T[duplicated(Syn_T$binomial),] # create temporary dataframe
          Syn_3 <- Syn_T[!duplicated(Syn_T$binomial),] # create data frame of third synonyms
          Syn_4 <- Syn_T[duplicated(Syn_T$binomial),] # create data frame of fourth synonyms
          rm(Syn_T)
        }}}}
  
  stats <- list(
    trait_database = trait_database,
    IUCN = IUCN,
    trait_out = trait_out,
    trait_out_final = trait_out_final)
  
  
  return(stats)
} # end of trait function
