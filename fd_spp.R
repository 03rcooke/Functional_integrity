# check working directory
#getwd()

# remove all current objects in environment
#rm(list=ls())

# load necessary packages: fossil, data.table, FD, clue, dplyr, qpcR, stats, cowplot
# fossil makes matrix
# data.table used to make data.tables for range areas
# FD calculates fucntional indices and gower dissimilarity
# clue used to process dendrograms and clusters
# dplyr used ...
# qpcR calculates RMSE
# stats used for cophenetic distances
# cowplot for plots

if (!require("pacman")) install.packages("pacman")
pacman::p_load(fossil, data.table, FD, clue, dplyr, qpcR, stats, cowplot, ggdendro)

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
rowSums(spp_missing) # number of missing data species per ecoregion

# edit site data to match species in trait data
UK_site <- UK_site[,-c(spp_col)] # removes NA species from site data
UK_trait <- UK_trait[rowSums(is.na(UK_trait))<length(UK_trait),] # returns the rows that have at least one non-NA value for trait data

isTRUE(nrow(UK_trait) == ncol(UK_site)) 
  # do the site and trait data sets contain the same number of species? - should be TRUE

# number of species per ecoregion
spp_final <- apply(UK_site, 1, sum)

# number of ecoregions assessed
e <- dim(UK_site)[1]
if (e >= 0) print(paste("number of ecoregions assessed =", e[1]))
  # print number of ecoregions assessed

#UK_trait <- as.matrix(UK_trait)

############### Functional indices ##############

# combine trait and site data
UK <- list(UK_trait, UK_site) ; names(UK) <- c("trait","site")

# calculate species x species distance matrix based on effect traits
gd <- gowdis(UK$trait, ord = c("podani"))
# applied the Podani 1999 correction to account for ordered traits (Lefcheck et al., 2014)

######### Plot multiple dengrograms of species based on effect traits #########
hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty")
# average = UPGMA, mcquitty = WPGMA
# omitted UPGMC and WPGMC methods because they are not appropriate 
# for non-metric distances (Lefcheck et al., 2014)
hclust_results <- lapply(hclust_methods, function(m) hclust(gd, m))
names(hclust_results) <- hclust_methods

dendro_con <- cl_consensus(hclust_results)

######## Testing the performance of different clustering algorithms

######### co-phenetic correlations: similarities (Kreft & Jetz, 2010) ###############

co_cor <- lapply(hclust_results, function(m) cophenetic(m))
# cophenetic function - stats package
cpc <- lapply(co_cor, function(m) cor(gd, m))
# add consensus method
con_cophe <- cophenetic(dendro_con)
cor_con <- cor(gd, con_cophe)
all_cpc <- c(cpc, consensus = cor_con)
all_cpc2 <- all_cpc[order(unlist(all_cpc), decreasing = TRUE)]

barplot(unlist(all_cpc2), ylim = c(0,(unlist(all_cpc2[[1]])+0.2*unlist(all_cpc2[[1]]))), xlab = "Linkage function", ylab = "Co-phenetic correlation")
# just adjust 0.2 to change scaling of y acis

# dendrogram with highest co-phenetic correlation
t <- all_cpc2[1] # already ordered by correlation
x <- names(t)
if (is.null(hclust_results[[x]])){
  plot(dendro_con, hang = -1, main = "Functional dengrogram (based on effect traits) \n with the highest co-phenetic correlation", xlab = "method = consensus", cex = 0.8)
} else {
  plot(hclust_results[[x]], hang = -1, main = "Functional dengrogram (based on effect traits) \n with the highest co-phenetic correlation", xlab = "method = ", cex = 0.8)
}
# hang = -1 puts all species lables on same level

########## matrix 2-norm: dissimilarities (Lefcheck et al., 2014) ##########

ud <- lapply(hclust_results, function(m) cl_ultrametric(m)) 
# ultrametric function - clue package
ultra <- lapply(ud, function(x) cl_dissimilarity(x, gd, method = "spectral"))
# add consensus method
con_ultra <- cl_dissimilarity(dendro_con, gd, method = "spectral")
all_ultra <- c(ultra, consensus = con_ultra)
all_ultra2 <- all_ultra[order(unlist(all_ultra), decreasing = FALSE)]

barplot(1/unlist(all_ultra2), ylim = c(0,(1/unlist(all_ultra2[[1]])+0.1*1/unlist(all_ultra2[[1]]))), xlab = "Linkage function", ylab = "1/spectral norm (2-norm) dissimilarities")
# just adjust 0.1 to change scaling of y axis

# dendrogram with lowest 2-norm value
u <- all_ultra2[1] # already ordered
v <- names(u)
if (is.null(hclust_results[[v]])){
  plot(dendro_con, hang = -1, main = "Functional dengrogram (based on effect traits) \n with the lowest (2-norm) dissimilarity", xlab = "method = consensus", cex = 0.8)
} else {
  plot(hclust_results[[v]], hang = -1, main = "Functional dengrogram (based on effect traits) \n with the lowest (2-norm) dissimilarity", xlab = "method = ", cex = 0.8)
}
# hang = -1 puts all species lables on same level

####### Quantifying the number of clusters (k) ############

###### L method (Salvador & Chan, 2004)
mer <- sort(hclust_results[[1]]$height, decreasing = TRUE)
# merge heights of clusters
mer2 <- data.frame(c(NA, mer), 1:(length(hclust_results[[1]]$height)+1)); names(mer2) <- c("h","k")
# combine merge heights and number of clusters

b <- max(mer2$k) # n = b-1
c <- seq(from = 3, to = b-2, by = 1)

xl <- lapply(c, function(c) seq(from = min(c), to = max(c), by = 1))
xr <- lapply(c, function(c) seq(from = (min(c+1)), to = max((c+1)), by = 1))

lc <- lapply(xl, function(xl) do.call("lm",list(h ~ k, data = quote(mer2), subset = 2:xl)))
rc <- lapply(xr, function(xr) do.call("lm",list(h ~ k, data = quote(mer2), subset = xr:b)))

RMSEl <- lapply(lc, function(lc) RMSE(lc))
RMSEr <- lapply(rc, function(rc) RMSE(rc))

rmsed <- data.frame(unlist(RMSEl), unlist(RMSEr), c, rep(b, length(RMSEl))); names(rmsed) <- c("RMSEl", "RMSEr", "c", "b")
rmsed <- mutate(rmsed, RMSEc = c-1/b-1*RMSEl + b-c/b-1*RMSEr) ## equation Salvador & Chan: RMSEc = c-1/b-1*RMSEl + b-c/b-1*RMSEr
min_rmsed <- filter(rmsed, RMSEc == min(RMSEc))
c_ <- min_rmsed$c
row_c <- which(rmsed$c == c_)

plot(mer2$k, mer2$h, xlab = "Number of clusters", ylab = "Merging height", main = "Evaluation plot")
abline(lc[[row_c]])
abline(rc[[row_c]])

ggplot(data = mer2, aes(k, h)) + # sets up plot
  geom_point(size = 2.5, shape = 16) + # adds points
  labs(title = "Evaluation plot", x = "Number of clusters", y = "Merging height") + # adds axis titles and main title
  geom_vline(aes(xintercept = c_)) + # adds vertical line at c_
  annotate("text", x = 0.5*b, y = 0.7*max(mer2$h, na.rm = TRUE), label = paste("RMSE(min): k =", c_), size = 6) + # adds text to centre of plot
  theme(axis.text.y = element_text(size=15), # changes size of axis labels
        axis.text.x = element_text(size=15), # changes size of axis labels
        axis.title.x = element_text(size=17), # changes size of axis title
        axis.title.y = element_text(size=17), # changes size of axis title
        plot.title = element_text(size = 20)) # changes size of main title

################ Functional redundancy & Functional dispersion ##################

# find number of groups and return species assignation to groups
e_gr <- cutree(hclust_results[[v]], k = c_)

# p = number of plots
p <- nrow(UK$site)

# c_ = number of effect groups
gr = c_

gr2 <- unique(e_gr)
gr3 <- gr2[order(gr2)]
t2 <- rep(p,gr)
e_group <- rep(gr3, t2)
site <- rep(row.names(UK$site), gr)

e_gr1 <- rep(e_gr, p)
e_gr_m <- matrix(e_gr1, p, length(e_gr), byrow = T, dimnames = list(rownames(UK$site), colnames(UK$site)))
mats <- list()
FRed1 <- data.frame()
for (i in 1:gr){
  t <- ifelse(e_gr_m == i, 1, 0)
  mats[[i]] <- t * UK$site
  FRed1 <- rbind(FRed1, mats[[i]])
}

FRed <- FRed1 %>% transmute(nbsp_gr = rowSums(FRed1)) # sum number of species per ecoregion per group

results <- dbFD(UK$trait, UK$site, corr = "cailliez", calc.FRic = FALSE, calc.FDiv = FALSE)
results2 <- data.frame(site, e_group, FRed$nbsp, results$FDis); names(results2) <- c("site", "group", "FRed", "FDis")
res <- results2[order(results2$site),]; row.names(res) <- NULL

# Community-weighted means
UK_CWM <- functcomp(UK$trait, UK$site) # CWM.type = "all" if I want frequencies of each ordinal class

