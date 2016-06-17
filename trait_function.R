## --------------------------------------------------------------
## Name: trait_function.R
## Description: Function to match IUCN species data to trait databases
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: June 2016 - 
## Inputs: species data
##            species = the list of species to match the trait data against
##         trait database to be matched to the species data
##            trait = the trait data set up with the species name column identifed as "binomial"
## Outputs: Function named "trait"
## --------------------------------------------------------------

trait <- function(trait, species) 
{
  # Set up required packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(dplyr, taxize, data.table, rlist, stats, tidyr)
  
  # dplyr: used to compare two data frames, create new columns and order data # calls: anti_join, full_join, inner_join, mutate, arrange
  # taxize: used to find taxonomic synonyms and subspecies for taxonomic mismatches # calls: synonyms
  # data.table: used to concantenate a list of data frames # calls: rbindlist
  # rlist: used to remove data frames from a list of data frames # calls: list.remove
  # stats: used to aggregate synonyms # calls: aggregate
  # tidyr: used to separate synonyms in to two columns # calls: separate
  
  spp <- anti_join(species, trait, by = "binomial") # uses dplyr package
  # IUCN species not listed in trait database
  
  trait_database <- anti_join(trait, species, by = "binomial") # uses dplyr package
  # Species listed by trait database but not listed by IUCN

  trait_out <- full_join(species, trait, by = "binomial") # uses dplyr package
  
  
  ## Find and try synonyms for trait data ##
  
  spp <- spp$binomial # species list to find synonyms for: species listed by IUCN1235 but not listed by trait_outTHERIA
  
  syn <- synonyms(spp, db = "itis") # find synonyms - uses taxize package
  no_syn_na <- which(is.na(syn) == TRUE); no_syn_na <- names(no_syn_na) # list those species that have no synonym data
  no_syn_col <- which(sapply(syn, NCOL) == 3); no_syn_col <- names(no_syn_col) # list those species that have data but no synonyms 
  no_syn <- c(no_syn_na, no_syn_col)
  
  spp_syn <- setdiff(spp, no_syn) # remove species from species list that have no synonyms
  
  syndf <- list.remove(syn, no_syn) # remove species data frames from synonym list that are empty - uses rlist package
  
  SYN <- as.data.frame(data.table::rbindlist(syndf)) # collapse list of data frames for each species into new single data frame - uses data.table package
  
  seq_bi <- sapply(syndf, nrow); seq_bi <- unlist(seq_bi) # number of rows per species
  seq_bi <- rep(spp_syn, times = as.vector(seq_bi)) # repeat species names with synonyms by seq_bi
  SYN <- mutate(SYN, binomial_IUCN = seq_bi)  # add sequence of species names to synonym data frame - uses dplyr package
  
  acc_tsn_e <- suppressWarnings(ifelse(is.na(as.numeric(SYN$acc_tsn)), SYN$syn_name, SYN$acc_tsn)) # swap accepted tsn with syn names where needed - suppresses warning: "In ifelse(is.na(as.numeric(SYN$acc_tsn)), SYN$syn_name, SYN$acc_tsn) : NAs introduced by coercion"
  SYN <- mutate(SYN, acc_tsn_e = acc_tsn_e) # add edited acc_tsn column to data frame - uses dplyr package
  
  syn_name_e <- suppressWarnings(ifelse(!is.na(as.numeric(SYN$syn_name)), SYN$acc_tsn, SYN$syn_name)) # swap syn names with accepted tsn where needed - suppresses warning: "In ifelse(!is.na(as.numeric(SYN$syn_name)), SYN$acc_tsn, SYN$syn_name) : NAs introduced by coercion"
  SYN <- mutate(SYN, syn_name_e = syn_name_e) # add edited syn_name column to data frame - uses dplyr package
  
  SYN <- aggregate(SYN, by = list(SYN$binomial_IUCN), FUN = unique) # collapse species with more than one synonym into one row - uses stats package
  SYN <- separate(SYN, syn_name_e, c("syn_name1", "syn_name2", "syn_name3", "syn_name4", "syn_name5", "syn_name6"), ", ", extra = "merge", fill = "right") # split multiple synonyms into separate columns - uses tidyr package
  
  SYN$syn_name1 <- gsub("c(", "", SYN$syn_name1, fixed = TRUE) # remove leading "c(" introduced by separate
  SYN$syn_name1 <- gsub("[^[:alnum:][:space:]]", "", SYN$syn_name1) # remove all punctuation introduced by separate
  SYN$syn_name2 <- gsub("[^[:alnum:][:space:]]", "", SYN$syn_name2) # remove all punctuation introduced by separate
  SYN$syn_name3 <- gsub("[^[:alnum:][:space:]]", "", SYN$syn_name3) # remove all punctuation introduced by separate
  SYN$syn_name4 <- gsub("[^[:alnum:][:space:]]", "", SYN$syn_name4) # remove all punctuation introduced by separate
  SYN$syn_name5 <- gsub("[^[:alnum:][:space:]]", "", SYN$syn_name5) # remove all punctuation introduced by separate
  SYN$syn_name6 <- gsub("[^[:alnum:][:space:]]", "", SYN$syn_name6) # remove all punctuation introduced by separate
  
  SYN <- mutate(SYN, syn_name_sub1 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , SYN$syn_name1))) # add column of trinomials reduced to binomials based on first synonyms - gsub removes last word of name and any spaces produced - uses dplyr package
  SYN <- mutate(SYN, syn_name_sub2 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , SYN$syn_name2))) # add column of trinomials reduced to binomials based on second synonyms - gsub removes last word of name and any spaces produced - uses dplyr package
  SYN <- mutate(SYN, syn_name_sub3 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , SYN$syn_name3))) # add column of trinomials reduced to binomials based on third synonyms - gsub removes last word of name and any spaces produced - uses dplyr package
  SYN <- mutate(SYN, syn_name_sub4 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , SYN$syn_name4))) # add column of trinomials reduced to binomials based on fourth synonyms - gsub removes last word of name and any spaces produced - uses dplyr package
  SYN <- mutate(SYN, syn_name_sub5 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , SYN$syn_name5))) # add column of trinomials reduced to binomials based on fifth synonyms - gsub removes last word of name and any spaces produced - uses dplyr package
  SYN <- mutate(SYN, syn_name_sub6 = gsub("[[:space:]]+$", "", gsub("[[:alpha:]]+$", "" , SYN$syn_name6))) # add column of trinomials reduced to binomials based on sixth synonyms - gsub removes last word of name and any spaces produced - uses dplyr package
  
  SYN <- SYN[c(6,2,7:ncol(SYN))] # reorder columns to move binomial from the end and drop unneeded columns
  
  ## Synonyms
  
  syn_a <- inner_join(SYN, trait, by = c("syn_name1" = "binomial")) # uses dplyr package
  # Species matched in trait database based on first synonyms
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_a$binomial_IUCN))),]
  # Remove species matched by first synonym - do not need to try second synonym
  
  syn_b <- inner_join(SYN, trait, by = c("syn_name2" = "binomial")) # uses dplyr package
  # Species matched in trait database based on second synonyms
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_b$binomial_IUCN))),]
  # Remove species matched by second synonym - do not need to try third synonym
  
  syn_c <- inner_join(SYN, trait, by = c("syn_name3" = "binomial")) # uses dplyr package
  # Species matched in trait database based on third synonym
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_c$binomial_IUCN))),]
  # Remove species matched by third synonym - do not need to try fourth synonym
  
  syn_d <- inner_join(SYN, trait, by = c("syn_name4" = "binomial")) # uses dplyr package
  # Species matched in trait database based on fourth synonym
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_d$binomial_IUCN))),]
  # Remove species matched by fourth synonym - do not need to try fifth synonym
  
  syn_e <- inner_join(SYN, trait, by = c("syn_name5" = "binomial")) # uses dplyr package
  # Species matched in trait database based on fifth synonym
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_e$binomial_IUCN))),]
  # Remove species matched by fifth synonym - do not need to try sixth synonym
  
  syn_f <- inner_join(SYN, trait, by = c("syn_name6" = "binomial")) # uses dplyr package
  # Species matched in trait database based on sixth synonym
  
  ## Subspecies
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_f$binomial_IUCN))),]
  # Remove species matched by sixth synonym - do not need to try first collapsed (from trinomial into binomial) subspecies
  
  syn_g <- inner_join(SYN, trait, by = c("syn_name_sub1" = "binomial")) # uses dplyr package
  # Species matched in trait database based on first subspecies
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_g$binomial_IUCN))),]
  # Remove species matched by first subspecies - do not need to try second collapsed subspecies
  
  syn_h <- inner_join(SYN, trait, by = c("syn_name_sub2" = "binomial")) # uses dplyr package
  # Species matched in trait database based on second subspecies

  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_h$binomial_IUCN))),]
  # Remove species matched by second subspecies - do not need to try third collapsed subspecies  
  
  syn_i <- inner_join(SYN, trait, by = c("syn_name_sub3" = "binomial")) # uses dplyr package
  # Species matched in trait database based on third subspecies
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_i$binomial_IUCN))),]
  # Remove species matched by third subspecies - do not need to try fourth collapsed subspecies  

  syn_j <- inner_join(SYN, trait, by = c("syn_name_sub4" = "binomial")) # uses dplyr package
  # Species matched in trait database based on fourth subspecies
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_j$binomial_IUCN))),]
  # Remove species matched by fourth subspecies - do not need to try fifth collapsed subspecies  
  
  syn_k <- inner_join(SYN, trait, by = c("syn_name_sub5" = "binomial")) # uses dplyr package
  # Species matched in trait database based on fifth subspecies
  
  SYN <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_k$binomial_IUCN))),]
  # Remove species matched by fifth subspecies - do not need to try sixth collapsed subspecies 
  
  syn_l <- inner_join(SYN, trait, by = c("syn_name_sub6" = "binomial")) # uses dplyr package
  # Species matched in trait database based on sixth subspecies
  
  IUCN_more <- SYN[!(SYN$binomial_IUCN %in% c(as.vector(syn_l$binomial_IUCN))),]
  IUCN_more_names <- as.data.frame(c(no_syn, as.vector(IUCN_more$binomial_IUCN))); names(IUCN_more_names) <- c("binomial_syn")
  IUCN_more_names <- arrange(IUCN_more_names, binomial_syn) # order data A-Z - uses dplyr package
  # IUCN species not listed in synonyms, i.e. species that need further matching efforts
  
  syn_comb <- rbind(syn_a, syn_b, syn_c, syn_d, syn_e, syn_f, syn_g, syn_h, syn_i, syn_j, syn_k, syn_l); setnames(syn_comb, "binomial_IUCN", "binomial") # Combine trait data matched by synonyms and subspecies
  syn_comb <- syn_comb[c(1,4:ncol(syn_comb))] # drop unneeded columns: sub_tsn and acc_tsn_e
  syn_comb <- inner_join(syn_comb, trait_out, by = "binomial") # join trait data with species data for species matched by synonyms - uses dplyr package
  syn_comb <- syn_comb[,!apply(syn_comb, 2, function(x) all(is.na(x)))] # remove columns of all NAs generated during inner_join
  colnames(syn_comb) <- gsub(".x", "", colnames(syn_comb), fixed = TRUE) # remove .x from column names generated during inner_join
  
  syn_names <- c("syn_name1", "syn_name_sub1", "syn_name2", "syn_name_sub2", "syn_name3", "syn_name_sub3", "syn_name4", "syn_name_sub4", "syn_name5", "syn_name_sub5", "syn_name6", "syn_name_sub6")
  times <- max(match(colnames(syn_comb), syn_names), na.rm = TRUE) # how many synonym columns have been produced
  trait_out[,syn_names[0:times]] <- NA # add syn_names to trait_out data frame
  
  trait_out_temp <- trait_out[!(trait_out$binomial %in% c(as.vector(syn_comb$binomial))),]
  # Remove species with NA values that have now been matched by synonyms (subspecies)
  
  trait_out2 <- rbind(trait_out_temp, syn_comb)
  # Combine matched species from synonyms with matched species from IUCN and trait database
  # (nrow(trait_out) == nrow(trait_out2))
  
  trait_out_final <- trait_out2[!(trait_out2$binomial %in% c(as.vector(trait_database$binomial))),]
  # removed species listed by trait database but not listed by IUCN
  # nrow = 5235
  
  #### Outputs ####
  
  stats <- list(
    IUCN_names = IUCN_more_names, # $IUCN_names
    trait_data = trait_out_final) # $trait_data
  
  return(stats)
} # end of trait function
