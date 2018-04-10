## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, out.width = '100%', out.height='40%')

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------
library(igraph)
library(bipartite)
library(ggplot2)

## ------------------------------------------------------------------------
data(olesen2002flores) # Check out the help for information on this data set!
olesen2002flores_binary <- 1 * (olesen2002flores > 0) # Make the data binary (unweighted)
I <- nrow(olesen2002flores_binary) # Number of lower level species (e.g., hosts, plants)
J <- ncol(olesen2002flores_binary) # Number of higher level species (e.g., parasites, pollinators)
S <- I + J # Total number of species, aka: Network size
L <- sum(olesen2002flores_binary > 0) # Number of edges in the network
A_i <- rowSums(olesen2002flores_binary) # The degree of hosts
A_j <- colSums(olesen2002flores_binary) # The degree of parasites
C <- L / (I * J) # Connectance
cc_high <- colSums(olesen2002flores_binary) / nrow(olesen2002flores_binary) # Clustering coefficient higher level (the number of realized links divided by the number of possible links for each species)
cc_low <- rowSums(olesen2002flores_binary) / ncol(olesen2002flores_binary) # Clustering coefficient lower level (the number of realized links divided by the number of possible links for each species)

## ------------------------------------------------------------------------
S_i <- rowSums(olesen2002flores) # Node strength of hosts
S_j <- colSums(olesen2002flores) # Node strength of parasites

## ------------------------------------------------------------------------
# Species level:
sl_metrics <- specieslevel(olesen2002flores)
names(sl_metrics$`higher level`)
sl_metrics$`higher level`$degree # The degree of flower visitors

# Group level:
gl_metrics <- grouplevel(olesen2002flores)
gl_metrics

# Network level
nl_metrics <- networklevel(olesen2002flores)
nl_metrics

## ------------------------------------------------------------------------
degreedistr(memmott1999)

## ------------------------------------------------------------------------
# Original
networklevel(memmott1999, index='connectance')

# Project plants
plants_projected <- tcrossprod(memmott1999>0) # Number of shared pollinators
diag(plants_projected) <- 0
g <- graph.adjacency(plants_projected, mode = 'undirected', weighted = T)
par(mar=c(0,0,0,0))
plot(g, vertex.size=6, vertex.label=NA, edge.color='black', edge.width=log(E(g)$weight), layout=layout.circle)
qplot(E(g)$weight)

## ------------------------------------------------------------------------
plants_projected <- as.one.mode(memmott1999, project = 'lower')
g <- graph.adjacency(plants_projected, mode = 'undirected', weighted = T)
qplot(E(g)$weight)

## ------------------------------------------------------------------------
# First, load the data 
otago_nodes <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Nodes.csv?token=ADd8ONa_lfqD2gKIE5h7vig2-1_xfezXks5azMwWwA%3D%3D')
otago_links <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Links.csv?token=ADd8ONRlDulEIlfx0L0qUfayAEPtQlfIks5azMwlwA%3D%3D')
otago_web <- graph.data.frame(otago_links, vertices = otago_nodes, directed = T)

# Also load a weighted food web
chesapeake_nodes <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Chesapeake_bay_nodes.csv?token=ADd8OP102bE0HzQNE_1hondsNVhvFiRjks5azMxFwA%3D%3D', header=F)
names(chesapeake_nodes) <- c('nodeId','species_name')
chesapeake_links <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Chesapeake_bay_links.csv?token=ADd8OFq6aj3oEZ84lFQ05iJ0Yu3cQ8qZks5azMxWwA%3D%3D', header=F)
names(chesapeake_links) <- c('from','to','weight')
ches_web <- graph.data.frame(chesapeake_links, vertices = chesapeake_nodes, directed = T)
plot(ches_web, edge.width=log(E(ches_web)$weight)/2, layout=layout.circle)

## ------------------------------------------------------------------------
deg_dist_out <- igraph::degree(otago_web, mode = 'out')
deg_dist_in <- igraph::degree(otago_web, mode = 'in')
df <- data.frame(deg=c(deg_dist_out,deg_dist_in), direction=c(rep('out',length(deg_dist_in)),rep('in',length(deg_dist_in))))
ggplot(df, aes(deg, fill=direction))+geom_histogram(alpha=0.3)

## ------------------------------------------------------------------------
deg_dist_out <- igraph::degree(ches_web, mode = 'out')
deg_dist_in <- igraph::degree(ches_web, mode = 'in')
df <- data.frame(deg=c(deg_dist_out,deg_dist_in), direction=c(rep('out',length(deg_dist_in)),rep('in',length(deg_dist_in))))
ggplot(df, aes(deg, fill=direction))+geom_histogram(alpha=0.3)

s_dist_out <- igraph::strength(ches_web, mode = 'out')
s_dist_in <- igraph::strength(ches_web, mode = 'in')
df <- data.frame(s=c(s_dist_out,s_dist_in),
                 direction=c(rep('out',length(s_dist_in)),rep('in',length(s_dist_in))))
ggplot(df, aes(s, fill=direction))+geom_histogram(alpha=0.3)

## ------------------------------------------------------------------------
transitivity(otago_web, type = 'local') # The ratio of the triangles connected to the vertex and the triples centered on the vertex.
transitivity(otago_web, type = 'global') # The ratio of the triangles and the connected triples in the graph

## ------------------------------------------------------------------------
geodesics <- distances(otago_web, mode = 'all') # Assume an undirected graph
geodesics[1:5,1:5]

## ------------------------------------------------------------------------
# Assume undirected networks
CC <- igraph::closeness(otago_web, mode = 'all', normalized = T) 
BC <- igraph::betweenness(otago_web, directed = F, normalized = T)
EC <- igraph::eigen_centrality(otago_web, directed = F, scale = T)
EC <- EC$vector

# Plot the web, re-sizing nodes by centrality.
V(otago_web)$BC <- BC
par(mar=c(0,0,0,0))
plot(otago_web, vertex.size=BC*100, vertex.label=NA, edge.arrow.width=0.5, edge.arrow.size=0.5, edge.curved=0.5, layout=layout.circle)

## ------------------------------------------------------------------------
par(mfrow=c(4,4), mar=c(.75,.75,.75,.75))
for (i in 0:15){ # note that counting starts at 0
  plot(graph_from_isomorphism_class(size = 3, number = i),
       edge.arrow.size = 0.4,
       edge.color='black',
       main = i + 1)
  box(col='red')
}

## ------------------------------------------------------------------------
motifs(otago_web) #absolute numbers
motifs(otago_web)/count_motifs(otago_web) # Proportion

