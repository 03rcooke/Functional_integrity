# check working directory
getwd

# load necessary packages

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

##### setting up sitexspecies matrix ####

# read data into R:
UK <- read.csv("ALL_Species_Ecoregions.csv")

# check the data:
head(UK)
str(UK)

create.matrix(UK,tax.name="binomial",locality="eco_code")


