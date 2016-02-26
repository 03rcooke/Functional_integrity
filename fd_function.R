FUN <- function(site, trait, plot_dendro, spp_list) 
{
  # load necessary packages: fossil, reshape2, data.table, FD, plyr
  
  ##################### fossil ##################
  if(require("fossil")){
    print("fossil is loaded correctly")
  } else {
    print("trying to install fossil")
    install.packages("fossil")
    if(require("fossil")){
      print("fossil installed and loaded")
    } else {
      stop("could not install fossil")
    }
  }
  
  ##################### reshape2 ######
  if(require("reshape2")){
    print("reshape2 is loaded correctly")
  } else {
    print("trying to install reshape2")
    install.packages("reshape2")
    if(require("reshape2")){
      print("reshape2 installed and loaded")
    } else {
      stop("could not install reshape2")
    }
  }
  
  ##################### data.table ######
  if(require("data.table")){
    print("data.table is loaded correctly")
  } else {
    print("trying to install data.table")
    install.packages("data.table")
    if(require("data.table")){
      print("data.table installed and loaded")
    } else {
      stop("could not install data.table")
    }
  }
  
  ##################### FD ######
  if(require("FD")){
    print("FD is loaded correctly")
  } else {
    print("trying to install FD")
    install.packages("FD")
    if(require("FD")){
      print("FD installed and loaded")
    } else {
      stop("could not install FD")
    }
  }
  
  ##################### plyr ######
  if(require("plyr")){
    print("plyr is loaded correctly")
  } else {
    print("trying to install plyr")
    install.packages("plyr")
    if(require("plyr")){
      print("plyr installed and loaded")
    } else {
      stop("could not install plyr")
    }
  }
  
  #### set up SITE data ###
  
  # sum the areas of each species per ecoregion
  UK_data_sum <- as.data.table(site[c(-3)])[, lapply(.SD, sum), by = list(id_no, binomial, eco_code)] # sum shape areas by species name and ecoregion
  UK_data_sum <- as.data.frame(UK_data_sum) # convert data table to data frame
  UK_data_sum <- with(UK_data_sum, UK_data_sum[order(binomial),]) # reorder by name of species
  
  # list of unique species
  Species <- unique(site[c("id_no", "binomial")])
  Species <- with(Species, Species[order(binomial),]) # reorder by name of species
  
  # create speciesxsite matrix
  site_m <- create.matrix(site, tax.name="binomial", locality="eco_code") # uses fossil package
  site_m <- t(site_m) # transpose
  
  ### set up TRAITS data ##
  
  trait$id_no <- NULL 
  # Turn off IUCN species id numbers (keep in data in case they become useful)
  
  trait[trait==-999] <- NA
  # convert -999's to NA
  
  trait$activity <- factor(trait$activity) ; trait$terrestriality <- factor(trait$terrestriality) ; trait$trophic <- factor(trait$trophic)
  # set activity, terrestriality and trophic level as factors
  
  ########### data cleaning for NAs
  species_remove <- trait[!rowSums(is.na(trait))<length(trait),] # returns the rows that have all NAs for traits
  spp_col <- unique(grep(paste(as.character(rownames(species_remove)),collapse="|"), colnames(site_m))) # finds column numbers in site data for species with all NAs for traits
  
  # total number of species per ecoregion (including species with missing trait data)
  spp_total <- apply(site_m, 1, sum)
  
  # missing data species
  id_missing <- site_m[,c(spp_col)] # matrix of missing data species per ecoregion
  spp_missing <- rowSums(id_missing) # number of missing data species per ecoregion
  
  # edit site data to match species in trait data
  site_m <- site_m[,-c(spp_col)] # removes NA species from site data
  trait <- trait[rowSums(is.na(trait))<length(trait),] # returns the rows that have at least one non-NA value for trait data
  
  isTRUE(nrow(trait) == ncol(site_m)) 
  # do the site and trait data sets contain the same number of species? - should be TRUE
  
  # number of species per ecoregion
  spp_final <- apply(site_m, 1, sum)
  
  # number of ecoregions assessed
  e <- dim(site_m)[1]
  
  ############### Functional indices ##############
  
  # combine trait and site data
  UK <- list(trait, site_m) ; names(UK) <- c("trait","site")
  
  # calculate species x species distance matrix based on effect traits
  gd <- gowdis(UK$trait)
  
  # Functional dispersion
  UK_dis <- fdisp(gd, UK$site) 
  
  # Community-weighted means
  UK_CWM <- functcomp(trait, UK$site) # CWM.type = "all" if I want frequencies of each ordinal class
  
  # plot dengrogram of species based on effect traits
  dendro <- hclust(gd, method = "average")

  
  
  # find number of groups and return species assignation to groups
  egroup <- cutree(dendro, k = 8)
  
  # gr = number of effect groups
  gr <- length(unique(egroup))
  
  if (spp_list == TRUE)
    cat(write.csv(Species, file = "Species_list.csv", row.names=FALSE)) # export species list as .csv
  
  if (plot_dendro == TRUE) 
    plot(dendro, main = "Cluster dengrogram based on effect traits", cex = 0.8)
  
  result <- list(
    eco = e[1],
    spp_total = spp_total,
    spp_missing = spp_missing,
    spp_final = spp_final,
    FunctDisp = UK_dis$FDis, 
    CWM = UK_CWM)
  return(result)
} # end of FUN