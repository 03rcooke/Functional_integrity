FUN <- function(site, trait, corr = "cailliez", spp_list = FALSE) 
{
  # load necessary packages: fossil, data.table, FD, clue, dplyr, qpcR, stats, cowplot
  # fossil makes matrix
  # data.table used to make data.tables for range areas
  # FD calculates fucntional indices and gower dissimilarity
  # clue used to process dendrograms and clusters
  # dplyr used ...
  # qpcR calculates RMSE
  # stats used for cophenetic distances
  # cowplot used to simplify ggplots
  # ggdendro for functional dendrograms
  
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(fossil, data.table, FD, clue, dplyr, qpcR, stats, cowplot, ggdendro)
  
  #### set up data ###
  
  # sum the areas of each species per ecoregion
  UK_data_sum <- as.data.table(site[c(-3)])[, lapply(.SD, sum), by = list(id_no, binomial, eco_code)] # sum shape areas by species name and ecoregion
  UK_data_sum <- as.data.frame(UK_data_sum) # convert data table to data frame
  UK_data_sum <- with(UK_data_sum, UK_data_sum[order(binomial),]) # reorder by name of species
  
  # list of unique species
  Species <- unique(site[c("id_no", "binomial")])
  Species <- with(Species, Species[order(binomial),]) # reorder by name of species
  
  # create speciesxsite matrix
  site_m <- create.matrix(site, tax.name="binomial", locality="eco_code") # uses fossil package
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
  
  ######### Plot multiple dengrograms of species based on effect traits #########
  hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty")
  # average = UPGMA, mcquitty = WPGMA
  # omitted UPGMC and WPGMC methods because they are not appropriate 
  # for non-metric distances (Lefcheck et al., 2014)
  hclust_results <- lapply(hclust_methods, function(m) hclust(gd, m))
  names(hclust_results) <- hclust_methods
  
  dendro_con <- cl_consensus(hclust_results)
  
  ######## Testing the performance of different clustering algorithms
  
  ########## matrix 2-norm: dissimilarities (Lefcheck et al., 2014) ##########
  
  ud <- lapply(hclust_results, function(m) cl_ultrametric(m)) 
  # ultrametric function - clue package
  ultra <- lapply(ud, function(x) cl_dissimilarity(x, gd, method = "spectral"))
  # add consensus method
  con_ultra <- cl_dissimilarity(dendro_con, gd, method = "spectral")
  all_ultra <- c(ultra, consensus = con_ultra)
  all_ultra2 <- all_ultra[order(unlist(all_ultra), decreasing = FALSE)]
  
  # dendrogram with lowest 2-norm value
  u <- all_ultra2[1] # already ordered
  v <- names(u)
  
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
  row_c <- which(rmsed$c == c_) # what's this used for?
  
  ################ Functional redundancy ##################
  
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
  results <- data.frame(site, e_group, FRed$nbsp); names(results) <- c("site", "group", "FRed")
  res <- results[order(results$site),]; row.names(res) <- NULL
  
########### Functional dispersion #########
  FD <- fdisp(gd, UK$site)
  
########### Community-weighted means ############
  CWM <- functcomp(UK$trait, UK$site) # CWM.type = "all" if I want frequencies of each ordinal class
  
########### Outputs and Plots ################
  if (spp_list == TRUE)
    cat(write.csv(Species, file = "Species_list.csv", row.names=FALSE)) # export species list as .csv
  
    if (is.null(hclust_results[[v]])){
          
      
      plot(dendro_con, hang = -1, main = "Functional dengrogram (based on effect traits) \n with the lowest (2-norm) dissimilarity", xlab = "method = consensus", cex = 0.8)
      
      dendro_con <- as.dendrogram(dendro_con$.Data)
      ggdendrogram(as.dendrogram(dendro_con), leaf_labels = TRUE, size = 2) +
        theme_dendro()
      
      } else {
          plot(hclust_results[[v]], hang = -1, main = "Functional dengrogram (based on effect traits) \n with the lowest (2-norm) dissimilarity", xlab = "method = ", cex = 0.8)
    
        dend <- dendro_data(hclust_results[[v]], type = "rectangle")
        ggplot(dend$segments) + 
          geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
          geom_text(data = dend$labels, aes(x, y, label = label),
                    hjust = 1, angle = 90, size = 3) +
          theme(axis.line.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.x=element_blank(),
                axis.title.x=element_blank(),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank())
        }
  
    eval <- ggplot(data = mer2, aes(k, h)) + # sets up plot
              geom_point(size = 2.5, shape = 16) + # adds points
              labs(title = "Evaluation plot", x = "Number of clusters", y = "Merging height") + # adds axis titles and main title
              geom_vline(aes(xintercept = c_)) + # adds vertical line at c_
              annotate("text", x = 0.5*b, y = 0.7*max(mer2$h, na.rm = TRUE), label = paste("RMSE(min): k =", c_), size = 6) + # adds text to centre of plot
              theme(axis.text.y = element_text(size=15), # changes size of axis labels
                    axis.text.x = element_text(size=15), # changes size of axis labels
                    axis.title.x = element_text(size=17), # changes size of axis title
                    axis.title.y = element_text(size=17), # changes size of axis title
                    plot.title = element_text(size = 20)) # changes size of main title
  
    ul <- data.frame(unlist(all_ultra2)); names(ul) <- c("ultra"); rownames(ul)[which(rownames(ul)=="average")]= "UPGMA"; rownames(ul)[which(rownames(ul)=="mcquitty")]= "WPGMA"; rownames(ul)[which(rownames(ul)=="single")]= "Single"; rownames(ul)[which(rownames(ul)=="complete")]= "Complete"; rownames(ul)[which(rownames(ul)=="consensus")]= "Consensus"; rownames(ul)[which(rownames(ul)=="ward.D2")]= "Ward D2"; rownames(ul)[which(rownames(ul)=="ward.D")]= "Ward D"
    clus <- ggplot(data = ul, aes(x = row.names(ul), y = 1/ultra)) +
        geom_bar(stat = "identity", fill = "darkblue") + # adds bars
        scale_x_discrete(limits = row.names(ul)) + # maintains order of row names in ul data frame
        labs(title = "Performance of clustering algorithms", x = "Linkage function", y = "Spectral norm (2-norm) similarities") + # adds axis titles and main title
        theme(axis.text.y = element_text(size=15), # changes size of axis labels
              axis.text.x = element_text(size=15), # changes size of axis labels
              axis.title.x = element_text(size=17), # changes size of axis title
              axis.title.y = element_text(size=17), # changes size of axis title
              plot.title = element_text(size = 20)) # changes size of main title
  
  plots <- c("dendro", "eval", "clus")
  
  result <- list(
    ecoregions = e[1],
    spp_total = spp_total,
    spp_missing = spp_missing,
    spp_final = spp_final,
    plots = plots,
    CWM = CWM,
    k = c_,
    FR = res,
    FD = FD$FDis)
  return(result)
} # end of FUN