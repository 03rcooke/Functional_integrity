# check working directory
getwd()

# load necessary packages: fossil, reshape2

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

##### setting up sitexspecies matrix ####

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

# check the data loaded correctly:
head(UK_data)
str(UK_data)

UK_ss <- create.matrix(UK_data,tax.name="binomial",locality="eco_code")

#### data frame

UK <- as.data.frame(UK_ss) # convert matrix to data frame
UK <- cbind(binomial = rownames(UK), UK) # add species name column to data frame
rownames(UK) <- NULL # turn off rownames for data frame

