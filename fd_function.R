## --------------------------------------------------------------
## Name: fd_function.R
## Description: Function to calculate functional diversity indices
## Author: R.S.C. Cooke, R.S.Cooke@soton.ac.uk
## Date: March 2016 - 
## Outputs: Function named 'FUN'
## --------------------------------------------------------------

FUN <- function(site, trait, corr = "cailliez", spp_list = FALSE, tree = FALSE) 
{
  # Set up required packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(fossil, data.table, FD, clue, dplyr, qpcR, stats, cowplot, ggdendro, ape, dendextend, methods)
  
  # fossil: used to create matrices # calls: create.matrix
  # data.table: used to make data.tables for range areas # calls: as.data.table
  # FD: calculates functional indices and gower dissimilarity # calls: gowdis, fdisp, functcomp
  # clue: used to process dendrograms and clusters # calls: cl_ensemble, cl_consensus, cl_ultrametric, cl_dissimilarity
  # dplyr: used to calculate across rows # calls: mutate, filter, transmute
  # qpcR: used to calculate RMSE # calls: RMSE
  # stats: used for clustering and data handling # calls: hclust, as.dendrogram, aggregate
  # cowplot: simplify ggplots # calls: ggplot, plot_grid
  # ggdendro: for functional dendrograms # calls: dendro_data
  # ape: used to write trees in Newick format # calls: write.tree
  # dendextend: used to cut a non-hclust dendrogram # calls: cutree
  # methods: used to create slots in output # calls: setClass
  
  #### Set up data ###
  
  # sum the areas of each species per ecoregion
  UK_data_sum <- as.data.table(site[c(-3)])[, lapply(.SD, sum), by = list(id_no, binomial, eco_code)] # sum shape areas by species name and ecoregion
  UK_data_sum <- as.data.frame(UK_data_sum) # convert data table to data frame
  UK_data_sum <- with(UK_data_sum, UK_data_sum[order(binomial),]) # reorder by name of species
  
  # All mammals for presence = extant; origin = native, reintroduced
  am12 <- read.csv("ALL_Mammals_1_2.csv")
  am12_sum <- aggregate(shape_Area~binomial, am12, FUN=sum)
  am12_merge <- merge(UK_data_sum, am12_sum, by.x = "binomial", by.y="binomial")
  names(am12_merge) <- c("binomial", "id_no", "eco_code", "eco_spp_area", "total_spp_area")
  am12 <- mutate(am12_merge, prop_spp_area = eco_spp_area/total_spp_area*100) # sum number of species per ecoregion per group
  
  # list of unique species
  Species <- unique(site[c("id_no", "binomial")])
  Species <- with(Species, Species[order(binomial),]) # reorder by name of species
  
  # create speciesxsite matrix
  site_m <- create.matrix(site, tax.name="binomial", locality="eco_code")
  site_m <- t(site_m) # transpose
  
  ### set up TRAITS data ##
  
  trait$id_no <- NULL 
  # Turn off IUCN species id numbers (keep in data in case they become useful)
  
  trait[trait==-999] <- NA
  # convert -999's to NA
  
  trait$activity <- factor(trait$activity) ; trait$terrestriality <- factor(trait$terrestriality) ; trait$trophic <- factor(trait$trophic)
  # set activity, terrestriality and trophic level as factors
  
  ########### data cleaning for NAs
  species_remove <- trait[!rowSums(is.na(trait))<length(trait),] # returns the rows that have all NAs for traits
  spp_col <- unique(grep(paste(as.character(rownames(species_remove)),collapse="|"), colnames(site_m))) # finds column numbers in site data for species with all NAs for traits
  
  # total number of species per ecoregion (including species with missing trait data)
  spp_total <- apply(site_m, 1, sum)
  
  # missing data species
  id_missing <- site_m[,c(spp_col)] # matrix of missing data species per ecoregion
  spp_missing <- rowSums(id_missing) # number of missing data species per ecoregion
  
  # edit site data to match species in trait data
  site_m <- site_m[,-c(spp_col)] # removes NA species from site data
  trait <- trait[rowSums(is.na(trait))<length(trait),] # returns the rows that have at least one non-NA value for trait data
  
  isTRUE(nrow(trait) == ncol(site_m)) 
  # do the site and trait data sets contain the same number of species? - should be TRUE
  
  # number of species per ecoregion
  spp_final <- apply(site_m, 1, sum)
  
  # number of ecoregions assessed
  e <- dim(site_m)[1]
  
  ############### Set up functional data ##############
  
  # combine trait and site data
  UK <- list(trait, site_m) ; names(UK) <- c("trait","site")
  
  # calculate species x species distance matrix based on effect traits
  gd <- gowdis(UK$trait, ord = c("podani"))
  # applied the Podani 1999 correction to account for ordered traits (Lefcheck et al., 2014)
  
  ###############################################################################
  ######### Plot multiple dengrograms of species based on effect traits #########
  ###############################################################################
  
  hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty")
  # average = UPGMA, mcquitty = WPGMA
  # omitted UPGMC and WPGMC methods because they are not appropriate for non-metric distances (Lefcheck et al., 2014)
  hclust_results <- lapply(hclust_methods, function(m) hclust(gd, method = m))
  names(hclust_results) <- hclust_methods
  
  ######## Testing the performance of different clustering algorithms ##########
  
  ########## matrix 2-norm: dissimilarities (Lefcheck et al., 2014) ############
  ud <- lapply(hclust_results, function(m) cl_ultrametric(m))
  # convert dendrograms to ultrametric
  
  dendro_en <- cl_ensemble(list = hclust_results)
  class(dendro_en)
  dendro_con <- cl_consensus(dendro_en)
  # build consensus dendrogram
  
  all_ultra <- c(ud, dendro_con[1]); names(all_ultra) <- c(hclust_methods, "consensus")
  (ul <- lapply(all_ultra, function(x) cl_dissimilarity(x, gd, method = "spectral")))
  # calculate dissimilarity values
  
  ul2 <- do.call(rbind, ul)
  min_dendro <- which.min(ul2); names(all_ultra)[min_dendro]
  # dendrogram with lowest 2-norm value
  v <- names(all_ultra)[min_dendro]
  # name of dendrogram with lowest 2-norm value
  min_dist <- all_ultra[names(all_ultra) == v][[1]]
  # distance matrix for dendrogram with lowest 2-norm value
  
  min_dist <- min_dist/max(min_dist)
  # scale between 0-1
  
  dend <- dendro_data(as.dendrogram(min_dist), type = "rectangle")
  # extract cluster data from min_dist dendrogram
  
  ####### Quantifying the number of clusters (k) ############
  
  ###### L method (Salvador & Chan, 2004)
  mer <- sort(unique(dend$segment[,"y"]), decreasing = TRUE)
  # merge heights of clusters
  mer2 <- data.frame(c(NA, mer), 1:(length(mer)+1)); names(mer2) <- c("h","k")
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
  
  ################ Functional redundancy ##################
  
  # find number of groups and return species assignation to groups
  e_gr <- cutree(as.dendrogram(min_dist), k = c_)
  
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
  results <- data.frame(site, e_group, FRed$nbsp); names(results) <- c("site", "group", "FRed")
  res <- results[order(results$site),]; row.names(res) <- NULL
  
########### Functional dispersion #########
  FD <- fdisp(gd, UK$site)
  
########### Community-weighted means ############
  CWM <- functcomp(UK$trait, UK$site) # CWM.type = "all" if frequencies of each ordinal class needed
  
########### Outputs and Plots ################
  if (spp_list == TRUE)
    cat(write.csv(Species, file = "Species_list.csv", row.names=FALSE)) # export species list as .csv
  
  ###### write newick tree
  if (tree == TRUE)
  write.tree(as.phylo(as.dendrogram(min_dist)),"Functional_dendrogram.new")
  
  ###### Plot: performance of clustering algorithms
  uld <- data.frame(unlist(ul)); names(uld) <- c("ultra"); rownames(uld)[which(rownames(uld)=="average")]= "UPGMA"; rownames(uld)[which(rownames(uld)=="mcquitty")]= "WPGMA"; rownames(uld)[which(rownames(uld)=="single")]= "Single"; rownames(uld)[which(rownames(uld)=="complete")]= "Complete"; rownames(uld)[which(rownames(uld)=="consensus")]= "Consensus"; rownames(uld)[which(rownames(uld)=="ward.D2")]= "Ward D2"; rownames(uld)[which(rownames(uld)=="ward.D")]= "Ward D"
  uld <- uld[order(uld$ultra), , drop = FALSE] # reorder data frame by increasing dissimilarity
  clus <- ggplot(data = uld, aes(x = row.names(uld), y = 1/ultra)) +
    geom_bar(stat = "identity", fill = "grey") + # adds bars
    scale_x_discrete(limits = row.names(uld)) + # maintains order of row names in ul data frame
    labs(title = "Performance of clustering algorithms", x = "Linkage function", y = "Spectral norm (2-norm) similarities") + # adds axis titles and main title
    theme(axis.text.y = element_text(size=15), # changes size of axis labels
          axis.text.x = element_text(size=15), # changes size of axis labels
          axis.title.x = element_text(size=17), # changes size of axis title
          axis.title.y = element_text(size=17), # changes size of axis title
          plot.title = element_text(size = 20)) # changes size of main title
  
  ###### Plot: Evaluation plot
  eval <- ggplot(data = mer2, aes(k, h)) + # sets up plot
    geom_point(size = 2.5, shape = 16) + # adds points
    labs(title = "Evaluation plot", x = "Number of clusters", y = "Merging height") + # adds axis titles and main title
    geom_vline(aes(xintercept = c_)) + # adds vertical line at c_
    annotate("text", x = 0.5*b, y = 0.7*max(mer2$h, na.rm = TRUE), label = paste("RMSE(min): k =", c_), size = 6) + # adds text to centre of plot
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme(axis.text.y = element_text(size=15), # changes size of axis labels
          axis.text.x = element_text(size=15), # changes size of axis labels
          axis.title.x = element_text(size=17), # changes size of axis title
          axis.title.y = element_text(size=17), # changes size of axis title
          plot.title = element_text(size = 20)) # changes size of main title
  
  ###### Plot: Functional dendrogram
  # Set up dendrogram
  dendro_plot <- ggplot(dend$segments) + # sets up plot
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + # adds heights of branchs
    ylab("") + # remove label for y axis
    coord_flip() + # make dendrogram horizontal
    scale_y_reverse(expand = c(0,0), breaks = seq(0, 1, by = 0.2)) + # reverse scale so that 0 is on the right
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), axis.line.y=element_blank(),
          axis.ticks.length = unit(0.5, "lines"),
          plot.margin = unit(c(1,0,1,1), "lines")) # remove y axis information
  
  # Add species labels separately to ensure axis is clipped to limits
  dendro_text <- ggplot(dend$labels) +
    scale_y_continuous(limits = c(0,1)) +
    geom_text(data = dend$labels, aes(x = x, y = y, label = label),
              hjust = 0, size = 4.5) + # adds species names
    coord_flip() + # make dendrogram horizontal
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), axis.line.y=element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), axis.line.x=element_blank(),
          plot.margin = unit(c(1,1,1,0), "lines")) # remove y axis information
  
  # Combine dendrogram and text plots
  dendro <- plot_grid(dendro_plot, dendro_text, align="h", rel_widths = c(3,1))
  
  stats <- list(
    ecoregions = e[1],
    spp_total = spp_total,
    spp_missing = spp_missing,
    spp_final = spp_final,
    k = c_,
    CWM = CWM,
    FRed = res,
    FDis = FD$FDis)
  
  plots <- list(
    clus = clus,
    eval = eval,
    dendro = dendro)
  
  result <- setClass("result", slots = c(stats = "list", plots = "list"))
  result <- result(stats = stats, plots = plots)
  
  return(result)
} # end of FUN function