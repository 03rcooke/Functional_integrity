# check working directory
getwd()

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

UK_trait <- read.csv("Trait_data_UK.csv", row.names = 2,
                     # add species names to rows
                     # colClasses means -999 values stay as a factor level
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
# assign simple names to variables

out <- FUN(UK_data, UK_trait, corr = "cailliez", spp_list = FALSE)

out

out$ecoregions # number of ecoregions assessed
out$spp_total # total number of species per ecoregion (including species with missing trait data)
out$spp_missing # number of missing data species per ecoregion
out$spp_final # number of species per ecoregion after removing missing-data species
out$CWM # community weighted mean
out$k # number of functional clusters as specified by the L method (Salvador & Chan, 2004)
out$FunctDisp # functional dispersion

print(out$evaluation) # L method evaluation plot

# plot_dendro is the functional dendrogram
# spp_list is a .csv of all species in the dataset