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
                     colClasses = c("character","factor","factor","numeric","factor","factor","numeric","numeric","factor","factor"),
                     # assign data types to variables
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
                    # assign simple names to variables
UK_trait$id_no <- NULL 
  # Turn off IUCN species id numbers (keep in data in case they become useful)

for (i in (1:length(UK_trait))) {UK_trait[,i]=ifelse(UK_trait[,i]==-999, NA, UK_trait[,i])} 
# convert -999's to NA

UK_trait$diet <- ordered(UK_trait$diet) ; UK_trait$habitat <- ordered(UK_trait$habitat)
# set diet and habitat breadth as ordered factors
UK_trait$activity <- factor(UK_trait$activity) ; UK_trait$terrestriality <- factor(UK_trait$terrestriality) ; UK_trait$trophic <- factor(UK_trait$trophic)
# set activity, terrestriality and trophic level as factors

# check the data has been set up correctly:
head(UK_trait)
str(UK_trait)

# data cleaning for NAs

species_remove <- UK_trait[!rowSums(is.na(UK_trait))<length(UK_trait),]
  # returns the rows that have all NAs for traits
spp_col <- unique(grep(paste(as.character(rownames(species_remove)),collapse="|"), 
                        colnames(UK_site)))
  # finds column numbers in site data for species with all NAs for traits
UK_site <- UK_site[,-c(spp_col)]
  # removes NA species from site data

UK_trait <- UK_trait[rowSums(is.na(UK_trait))<length(UK_trait),] 
  # returns the rows that have at least one non-NA value

# combine trait and site data
UK <- list(UK_trait, UK_site) ; names(UK) <- c("trait","site")
dbUK <- dbFD(UK$trait, UK$site, corr = "cailliez")

gd <- gowdis(UK$trait)
UK_dis <- fdisp(gd, UK$site)
UK_dis$FDis



dendro = hclust(gd, method = "average")
plot(dendro)

cut_dendro <- cutree(dendro, k = 10)


