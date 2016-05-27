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




dendro_a <- hclust(gd, method = "single")
dendro_b <- hclust(gd, method = "complete")
dendro_c <- hclust(gd, method = "ward.D")
dendro_d <- hclust(gd, method = "ward.D2")
dendro_e <- hclust(gd, method = "average") # UPGMA
dendro_f <- hclust(gd, method = "mcquitty") # WPGMA


#FRed3 <- cbind(site, rep(gr2, p), FRed2); names(FRed3) <- c("site","group","nbsp")
#FRed <- FRed3[order(FRed3$site),]

#e_gr_fac <- factor(e_gr)
#nbsp_e_gr <- tapply(e_gr_fac, e_gr_fac, length)

dbUK <- dbFD(UK$trait, UK$site, corr = , m = "min") # need a m argument to get it to run
dbUK <- dbFD(UK$trait, UK$site, corr = "cailliez", calc.FRic = FALSE) # need FRic to be false to get it to run

#dbUK <- dbFD(UK$trait, UK$site, corr = "cailliez") # try running this on remote desktop

# Functional dispersion

UK_dis <- fdisp(gd, UK$site)
UK_dis$FDis

# Community-weighted means

UK_CWM <- functcomp(UK$trait, UK$site) # CWM.type = "all" if I want frequencies of each ordinal class
UK_CWM



dendro_plot <- ggplot(dend$segments) + # sets up plot
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + # adds heights of branchs
  ylab("") + # remove label for y axis
  coord_flip() + # make dendrogram horizontal
  scale_y_reverse(expand = c(0, 0), breaks = seq(0, 1, by = 0.2)) + # reverse scale so that 0 is on the right
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), axis.line.y=element_blank(),
        plot.margin = unit(c(1,7,1,1), "lines"), axis.ticks.length = unit(0.5, "lines")) # remove y axis information

dp <- dendro_plot + geom_text(data = dend$labels, aes(x, y, label = label),
                              hjust = 0, size = 4)

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(dp))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)

########## matrix 2-norm: dissimilarities (Lefcheck et al., 2014) ############
ud <- lapply(hclust_results, function(m) cl_ultrametric(as.hclust(m)))
# convert dendrograms to ultrametric

dendro_en <- cl_ensemble(list = ud)
class(dendro_en)
dendro_con <- cl_consensus(dendro_en)

rownames(uld)[which(rownames(uld)=="average")]= "UPGMA"; rownames(uld)[which(rownames(uld)=="mcquitty")]= "WPGMA"; rownames(uld)[which(rownames(uld)=="single")]= "Single"; rownames(uld)[which(rownames(uld)=="complete")]= "Complete"; rownames(uld)[which(rownames(uld)=="consensus")]= "Consensus"; rownames(uld)[which(rownames(ul)=="ward.D2")]= "Ward D2"; rownames(ul)[which(rownames(ul)=="ward.D")]= "Ward D"


TM <- c(as.vector(PanTHERIA$binomial), as.vector(IUCN$binomial))
sort(TM)
# Vector of taxonomic mismatches between the PanTHERIA and IUCN mammal species lists
# length = 1043

# R.utils: used to identify empty csvs (no species within ecoregion) # calls: countLines
empty_files <- lapply(Filter(function(x) countLines(x)<=1, csv_files), unlink)
