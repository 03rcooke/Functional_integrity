# check working directory
getwd()

# load necessary packages: fossil, reshape2, data.table, plyr

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

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

# check the data loaded correctly:
head(UK_data)
str(UK_data)

# sum the areas of each species per ecoregion
UK_data_sum <- as.data.table(UK_data[c(-3)])[, lapply(.SD, sum), by = list(id_no, binomial, eco_code)] # sum shape areas by species name and ecoregion
UK_data_sum <- as.data.frame(UK_data_sum) # convert data table to data frame
with(UK_data_sum, UK_data_sum[order(binomial),]) # reorder by name of species

# create speciesxsite matrix
UK_ss <- create.matrix(UK_data, tax.name="binomial", locality="eco_code") # uses fossil package

UK_ss_t <- t(UK_ss) # transpose?

#### set up speciesxsite data frame ###
UK <- as.data.frame(UK_ss) # convert matrix to data frame
UK <- cbind(binomial = rownames(UK), UK) # add species name column to data frame
rownames(UK) <- NULL # turn off rownames for data frame

### set up traits data ##




# EXTRA

# birdTraitsr <- subset(birdTraits, select = c("Common", "logLen", "abun"))
  # select traits required

# for (i in (2:length(birds))) {birds[,i]=ifelse(birds[,i]==0, NA, birds[,i])} 
  # if I need to convert 0's to NA

# birdsr=merge(birdTraitsr, birds, by="Common")
# rownames(birdsr) <- birdsr[,"Common"]
  # merge trait and site data

# my.dist.mat.2 = dist(as.matrix(traits[, 2]), method = "euclidean")
  # distance matrix for trait 2 only
# dist(traits, method = "euclidean")
  # distance matrix using all traits


# library(FD)
# gowdis(traits)
  # Gower distance matrix
