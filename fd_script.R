## --------------------------------------------------------------
## Name: fd_script.R
## Description: Code to run "fd_eco" function
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: March 2016 - 
## Outputs: 
## --------------------------------------------------------------

# check working directory
getwd()

# read data into R:
UK_data <- read.csv("ALL_Species_Ecoregions.csv")

UK_trait <- read.csv("Trait_data_UK.csv", row.names = 2,
                     # add species names to rows
                     col.names = c("id_no","binomial","activity","mass","diet","habitat","litter","longevity","terrestriality","trophic"))
# assign simple names to variables

trait_P <- readRDS("trait_PanTHERIA.rds")
trait_A <- readRDS("trait_Amniote.rds")
trait_E <- readRDS("trait_Elton.rds")
trait_M <- readRDS("trait_MD.rds")

########################## Combine into one data frame with the columns I need




out_fd <- fd_eco(UK_data, UK_trait, corr = "cailliez", spp_list = FALSE, tree = FALSE)

out_fd

#### Statistics
out_fd@stats

out_fd@stats$ecoregions # number of ecoregions assessed
out_fd@stats$spp_total # total number of species per ecoregion (including species with missing trait data)
out_fd@stats$spp_missing # number of species missing all trait data per ecoregion
out_fd@stats$spp_final # number of species per ecoregion after removing missing-data species
out_fd@stats$k # number of functional clusters as specified by the L method (Salvador & Chan, 2004)
out_fd@stats$CWM # community weighted mean
out_fd@stats$FRed # functional redundancy
out_fd@stats$FDis # functional dispersion

#### Plots
out_fd@plots

out_fd@plots$clus # Performance of clustering algorithms
out_fd@plots$eval # L method evaluation plot
out_fd@plots$dendro # Functional dendrogram: following best clustering algorithm and number of clusters as determined by L method

#### Session info
out_fd@session_info