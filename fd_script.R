# check working directory
getwd()

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

UK_trait <- read.csv("Trait_data_UK.csv", row.names = 2,
                     # add species names to rows
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
# assign simple names to variables

All_trait <- read.csv("ALL_mammals_1235_traits.csv", row.names = 2)

out <- FUN(UK_data, UK_trait, corr = "cailliez", spp_list = FALSE, tree = FALSE)

out

#### Statistics
out@stats

out@stats$ecoregions # number of ecoregions assessed
out@stats$spp_total # total number of species per ecoregion (including species with missing trait data)
out@stats$spp_missing # number of species missing all trait data per ecoregion
out@stats$spp_final # number of species per ecoregion after removing missing-data species
out@stats$k # number of functional clusters as specified by the L method (Salvador & Chan, 2004)
out@stats$CWM # community weighted mean
out@stats$FR # functional redundancy
out@stats$FD # functional dispersion

#### Plots
out@plots

out@plots$clus # Performance of clustering algorithms
out@plots$eval # L method evaluation plot
out@plots$dendro # Functional dendrogram: following best clustering algorithm and number of clusters as determined by L method

# spp_list is a .csv of all species in the dataset
# tree is a Newick format functional dendrogram