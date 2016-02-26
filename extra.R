# EXTRA

#### set up speciesxsite data frame ###
# UK_site <- as.data.frame(UK_ss) # convert matrix to data frame
# UK_site <- cbind(eco_code = rownames(UK_site), UK_site) # add species name column to data frame
# rownames(UK_site) <- NULL # turn off rownames for data frame

# birdTraitsr <- subset(birdTraits, select = c("Common", "logLen", "abun"))
# select traits required

# birdsr=merge(birdTraitsr, birds, by="Common")
# rownames(birdsr) <- birdsr[,"Common"]
# merge trait and site data

# my.dist.mat.2 = dist(as.matrix(traits[, 2]), method = "euclidean")
# distance matrix for trait 2 only
# dist(traits, method = "euclidean")
# distance matrix using all traits

UK_trait <- subset(UK_trait, select = -c(id_no, binomial))
# Turn off IUCN species id numbers and scientific names

names(UK_trait) <- c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic")
# rename columns with simple names

for (i in (1:length(UK_trait))) {UK_trait[,i]=ifelse(UK_trait[,i]==-999, NA, UK_trait[,i])} 
# convert -999's to NA

#UK_trait[UK_trait == -999] <- NA
# convert -999's to NA

final[rowSums(is.na(final))<(length(final)-1),] # return the rows that have at least TWO non-NA values

UK_trait <- read.csv("Trait_data_UK.csv", 
                     colClasses = c("character","factor","factor","numeric","factor","factor","numeric","numeric","factor","factor"),
                     # assign data types to variables
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
# assign simple names to variables

rownames(UK_trait) <- UK_trait$binomial
  # add species names to rows

UK_trait$activity <- factor(UK_trait$activity) ; UK_trait$terrestriality <- factor(UK_trait$terrestriality) ; UK_trait$trophic <- factor(UK_trait$trophic)
# set activity, terrestriality and trophic level as factors


################################## Attempts to deal with NAs ###########################

species_remove <- rownames(species_remove)

species_remove <- c(as.character("Myotis mystacinus"))

`%ni%` <- Negate(`%in%`)
subset(UK_site, select = names(UK_site) %ni% species_remove)

UK_site[,-which(names(UK_site) %in% species_remove)]

to.remove <- c("hp","drat","wt","qsec")
mtcars[,-which(names(mtcars) %in% to.remove)]

UK_sf <- as.data.frame(UK_site)
subset(UK_sf, select = -c(as.character("Myotis mystacinus")))

UK_site[ , as.character(c("Myotis mystacinus","Pipistrellus pygmaeus"))]

UK_site[ , !as.character(c(rownames(species_remove)))]

grep(as.character("Myotis mystacinus"), colnames(UK_site))

UK_site[ , -which(names(as.character(UK_site)) %in% as.character(c(rownames(species_remove))))]



# UK_trait[which(UK_trait==-999)] <- NA

for (i in (1:length(UK_trait))) {UK_trait[,i]=ifelse(UK_trait[,i]==-999, NA, UK_trait[,i])} 

UK_trait$diet <- ordered(UK_trait$diet) ; UK_trait$habitat <- ordered(UK_trait$habitat)
# set diet and habitat breadth as ordered factors

if (e >= 0) print(paste("number of ecoregions assessed =", e[1]))
# print number of ecoregions assessed
