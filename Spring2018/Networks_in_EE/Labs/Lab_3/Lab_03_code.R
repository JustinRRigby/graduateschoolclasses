## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, out.width = '100%', out.height='40%')

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------
library(igraph)
library(bipartite)
library(RColorBrewer)

## ------------------------------------------------------------------------
chesapeake_nodes <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Chesapeake_bay_nodes.csv?token=ADd8OB99ouaZ72g32T037z_fAk33r_R_ks5a1IuYwA%3D%3D', header=F)
names(chesapeake_nodes) <- c('nodeId','species_name')
chesapeake_links <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Chesapeake_bay_links.csv?token=ADd8OMagGAfcXWC6rTyJV7gCv4OfclgXks5a1IurwA%3D%3D', header=F)
names(chesapeake_links) <- c('from','to','weight')
ches_web <- graph.data.frame(chesapeake_links, vertices = chesapeake_nodes, directed = T)

## ------------------------------------------------------------------------
ches_web_unweighted <- ches_web
ches_web_unweighted <- as.undirected(ches_web_unweighted)
E(ches_web_unweighted)$weight <- 1
cl <- cluster_louvain(ches_web_unweighted) # Can also use the weights = NULL argument
class(cl) # the result is of class communities
module_membership <- membership(cl)
cols <- data.frame(mem=unique(module_membership), col= brewer.pal(length(unique(module_membership)), 'Set2'))
V(ches_web_unweighted)$module_membership <- module_membership
V(ches_web_unweighted)$color <- cols$col[match(V(ches_web_unweighted)$module_membership, cols$mem)]
plot(ches_web_unweighted, vertex.color=V(ches_web_unweighted)$color, vertex.size=5, vertex.label=NA, edge.arrow.width=0.3, edge.arrow.curve=0.5)

## ------------------------------------------------------------------------
cl_wt <- cluster_louvain(as.undirected(ches_web), weights = E(ches_web)$weight) #Notice the `weights` argument.
class(cl_wt) # the result is of class communities
module_membership_wt <- membership(cl_wt)
cols <- data.frame(mem=unique(module_membership_wt), col= brewer.pal(length(unique(module_membership_wt)), 'Set1'))
V(ches_web)$module_membership_wt <- module_membership_wt
V(ches_web)$color_wt <- cols$col[match(V(ches_web)$module_membership_wt, cols$mem)]
plot(as.undirected(ches_web), vertex.color=V(ches_web)$color_wt, vertex.size=5, vertex.label=NA, edge.arrow.width=0.3, edge.arrow.curve=0.5)

## ------------------------------------------------------------------------
mod <- computeModules(memmott1999)
slotNames(mod) # see ?moduleWeb for details
mod@likelihood # This is the value of the modularity function Q. NOTICE THE @ SIGN (instead of $).
module_list <- listModuleInformation(mod) # The output is rather cumbersome...
plotModuleWeb(mod)
module_list <- module_list[[2]] # let's look at the modules. The first element in the list is the whole network, so start with 2
for (i in 1:2){ # Show the two first modules.
  message(paste('Module:',i))
  print(module_list[[i]])
}

## ------------------------------------------------------------------------
# Transform the lists to data frames
m <- length(module_list) # Number of modules
mod_plants <- unlist(module_list, recursive = F)[seq(1,2*m,2)] # Assignments of plants
names(mod_plants) <- 1:m
mod_pollinators <- unlist(module_list, recursive = F)[seq(2,2*m,2)] # Assignments of pollinators
names(mod_pollinators) <- 1:m
tmp_plants <- data.frame(module = rep(names(mod_plants), sapply(mod_plants, length)), species = unlist(mod_plants), type='Plants')
tmp_poillnators <- data.frame(module = rep(names(mod_pollinators), sapply(mod_pollinators, length)), species = unlist(mod_pollinators), type='Pollinators')
# Make one data frame
module_assignments <- rbind(tmp_plants,tmp_poillnators)

## ------------------------------------------------------------------------
library(vegan)
?commsim

## ------------------------------------------------------------------------
num_iterations <- 10
memmott1999_binary <- memmott1999
memmott1999_binary[memmott1999_binary>0] <- 1
null_model <- vegan::nullmodel(memmott1999_binary, method = 'r00')
shuffled_r00 <- simulate(null_model, nsim = num_iterations)
apply(shuffled_r00, 3, dim) # This produces an array with 10 matrices, each of them is a shuffled matrix.
# Calculate the density of each of these matrices:
connectance <- function(m){
  d <- sum(m>0)/(nrow(m)*ncol(m))
  return(d)
}
apply(shuffled_r00, MARGIN = 3, connectance)

## ------------------------------------------------------------------------
num_iterations <- 10
null_model_curveball <- vegan::nullmodel(memmott1999_binary, method = 'curveball') # Only preserve matrix fill (density)
shuffled_curveball <- simulate(null_model_curveball, nsim = num_iterations, burnin = 1000)
apply(shuffled_curveball, MARGIN = 3, connectance) # Calculate the density of each of these matrices
# Did the algorithm match the row  and column sums?
all(apply(shuffled_curveball, MARGIN = 3, rowSums)==rowSums(memmott1999_binary))
all(apply(shuffled_curveball, MARGIN = 3, colSums)==colSums(memmott1999_binary))
# Did it shuffle any interactions?
for (i in 1:num_iterations){
  if(!identical(shuffled_curveball[,,i], memmott1999_binary)){
    print(paste('Shuffled matrix',i,'is different than the original'))
  }
}

## ---- eval=FALSE---------------------------------------------------------
## web <- memmott1999_binary
## web[web>0] <- 1
## mod_observed <- computeModules(web)
## Q_obs <- mod_observed@likelihood # This is the value of the modularity function Q
## 
## # Now shuffle
## num_iterations <- 100
## null_model <- vegan::nullmodel(web, method = 'r00') # Only preserve matrix fill (density)
## shuffled_r00 <- simulate(null_model, nsim = num_iterations)
## # Calculate modularity of the shuffled networks
## Q_shuff <- NULL
## for (i in 1:100){
##   mod <- computeModules(shuffled_r00[,,i])
##   Q_shuff <- c(Q_shuff, mod@likelihood)
## }
## #Which proportion of the simulations have a Q value larger than the observed?
## P_value <- sum(Q_obs < Q_shuff)/length(Q_shuff)

