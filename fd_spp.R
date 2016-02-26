# check working directory
getwd()

# remove all current objects in environment
rm(list=ls())

# load necessary packages: fossil, reshape2, data.table, FD, plyr

##################### fossil ######
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

##### setting up data sheets ####

#### set up SITE data ###

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

# check the data loaded correctly:
head(UK_data)
str(UK_data)

# sum the areas of each species per ecoregion
UK_data_sum <- as.data.table(UK_data[c(-3)])[, lapply(.SD, sum), by = list(id_no, binomial, eco_code)] # sum shape areas by species name and ecoregion
UK_data_sum <- as.data.frame(UK_data_sum) # convert data table to data frame
UK_data_sum <- with(UK_data_sum, UK_data_sum[order(binomial),]) # reorder by name of species

# list of unique species
Species <- unique(UK_data[c("id_no", "binomial")])
Species <- with(Species, Species[order(binomial),]) # reorder by name of species
# write.csv(Species, file = "Species_list.csv", row.names=FALSE) # export species list as .csv

# create speciesxsite matrix
UK_site <- create.matrix(UK_data, tax.name="binomial", locality="eco_code") # uses fossil package
UK_site <- t(UK_site) # transpose

### set up TRAITS data ##

# read data into R:
UK_trait <- read.csv("Trait_data_UK.csv", row.names = 2,
                     # add species names to rows
                     # colClasses means -999 values stay as a factor level
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
                    # assign simple names to variables
UK_trait$id_no <- NULL 
  # Turn off IUCN species id numbers (keep in data in case they become useful)

UK_trait[UK_trait==-999] <- NA
  # convert -999's to NA

UK_trait$activity <- factor(UK_trait$activity) ; UK_trait$terrestriality <- factor(UK_trait$terrestriality) ; UK_trait$trophic <- factor(UK_trait$trophic)
  # set activity, terrestriality and trophic level as factors

# check the data has been set up correctly:
head(UK_trait)
str(UK_trait)

########### data cleaning for NAs
species_remove <- UK_trait[!rowSums(is.na(UK_trait))<length(UK_trait),] # returns the rows that have all NAs for traits
spp_col <- unique(grep(paste(as.character(rownames(species_remove)),collapse="|"), colnames(UK_site))) # finds column numbers in site data for species with all NAs for traits

# total number of species per ecoregion (including species with missing trait data)
spp_total <- apply(UK_site, 1, sum)

# missing data species
spp_missing <- UK_site[,c(spp_col)] # matrix of missing data species per ecoregion
rowSums(missing_spp) # number of missing data species per ecoregion

# edit site data to match species in trait data
UK_site <- UK_site[,-c(spp_col)] # removes NA species from site data
UK_trait <- UK_trait[rowSums(is.na(UK_trait))<length(UK_trait),] # returns the rows that have at least one non-NA value for trait data

isTRUE(nrow(UK_trait) == ncol(UK_site)) 
  # do the site and trait data sets contain the same number of species? - should be TRUE

# number of species per ecoregion
spp_final <- apply(UK_site, 1, sum)

# number of ecoregions assessed
e <- dim(UK_site)[1]
if (e > 1) print(paste("number of ecoregions assessed =", e[1]))
  # print number of ecoregions assessed

#UK_trait <- as.matrix(UK_trait)

############### Functional indices ##############

# combine trait and site data
UK <- list(UK_trait, UK_site) ; names(UK) <- c("trait","site")

# calculate species x species distance matrix based on effect traits
gd <- gowdis(UK$trait)


dbUK <- dbFD(UK$trait, UK$site, corr = "cailliez", m = "min") # need a m argument to get it to run
dbUK <- dbFD(UK$trait, UK$site, corr = "cailliez", calc.FRic = FALSE) # need FRic to be false to get it to run

#dbUK <- dbFD(UK$trait, UK$site, corr = "cailliez") # try running this on remote desktop

# Functional dispersion

UK_dis <- fdisp(gd, UK$site)
UK_dis$FDis

# Community-weighted means

UK_CWM <- functcomp(UK$trait, UK$site) # CWM.type = "all" if I want frequencies of each ordinal class
UK_CWM

# plot dengrogram of species based on effect traits
dendro <- hclust(gd, method = "average")
plot(dendro, main = "Cluster dengrogram based on effect traits", cex = 0.8)

# find number of groups and return species assignation to groups
egroup <- cutree(dendro, k = 8)

# p = number of plots
p <- nrow(UK$site)

# gr = number of effect groups
gr <- length(unique(egroup))




