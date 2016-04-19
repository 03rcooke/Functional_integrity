###GBIF

library(dismo)  # check also the nice 'rgbif' package! 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dismo, rworldmap, foreign)


otter <- gbif("Lutra", "lutra")

# get data frame with spatial coordinates (points)
locs <- subset(otter, select = c("country", "lat", "lon"))
head(locs)  # a simple data frame with coordinates

# Discard data with errors in coordinates:
locs <- subset(locs, locs$lat < 90)

coordinates(locs) <- c("lon", "lat")  # set spatial coordinates
plot(locs)

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
proj4string(locs) <- crs.geo  # define projection system of our data
summary(locs)

plot(locs, pch = 20, col = "steelblue")

# library rworldmap provides different types of global maps, e.g:
data(coastsCoarse)
data(countriesLow)
plot(coastsCoarse, add = T)

locs_gb <- subset(locs, locs$country == "United Kingdom")  # select only locs in UK
plot(locs_gb, pch = 20, cex = 2, col = "steelblue")
title("Otter occurrences in the UK")
plot(countriesLow, add = T)
summary(locs_gb)




library(foreign)
data <- read.dbf("<Path to your file>")