## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, out.width = '100%', out.height='40%')

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------
library(igraph)
library(Matrix)
library(rgl)
library(reshape2)
library(dplyr)

## ------------------------------------------------------------------------
df2matrix <- function(df,binary=F){
  rownames(df) <- df[,1]
  df <- df[,-1]
  df <- data.matrix(df)
  if (binary){df[df>0] <- 1}
  return(df)
}

# Network parameters
H <- 22 # Number of host species
P <- 56 # number of parasite species
S <- H+P # total number of species
L <- 6 # layers

# Get data
dat <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Pilosof_2017NEE_temporal_data.csv?token=ADd8OByYcPDFBzT2aHWq3wawxh8S9Xdqks5a7mKRwA%3D%3D')
attach(dat)

# Step 1: Build the six layers
hostAbundYear <- as.matrix(table(Host, YearCollected)) # abundance of hosts in different years
parasiteAbundanceYear <- aggregate(.~dat$YearCollected,data = dat[,3:ncol(dat)], FUN = sum) # abundance of parasites in each year (across hosts)
names(parasiteAbundanceYear)[1] <- 'YearCollected'
parasiteAbundanceYear <- df2matrix(parasiteAbundanceYear)

## Create network layers
data_matrices <- list()
years <- 1982:1987
for (y in years){
  idx <- which(years==y)
  d <- dat[dat$YearCollected==y,]
  d <- aggregate(.~d$Host, data=d[,2:ncol(d)], sum) # The total number of parasites found on a given host
  d <- d[-2]
  d <- df2matrix(d)
  d <- sweep(d, 1, hostAbundYear[rownames(hostAbundYear)%in%rownames(d), idx], '/') # Average parasite abundance per host
  missingHosts <- setdiff(rownames(hostAbundYear),rownames(d)) # All hosts have to appear in all matrices even if they were not present
  d <- rbind(d,matrix(0,length(missingHosts),ncol(d),dimnames =  list(missingHosts,colnames(d)))) # Add missing host species so all matrices will have the same size
  d <- d[sort(rownames(d)),] # sort by host so all matrices will have the same order
  data_matrices[[idx]] <- d
  # write.table(d, paste('../Network_data/host_parasite_abundance_weighted_layer_',idx,'.csv',sep=''), row.names = F, col.names = F,sep=',')
}
names(data_matrices) <- years
sapply(data_matrices,dim)

## Data frames with node data
physical_nodes <- data.frame(id=as.numeric(1:S),
                        type=c(rep('host',H),rep('paras',P)),
                        species=c(rownames(data_matrices[[1]]),colnames(data_matrices[[1]]))
                        )
state_nodes <- do.call("rbind", replicate(L, physical_nodes, simplify = FALSE))
state_nodes$state_id <- as.numeric(rownames(state_nodes))
state_nodes$layer <- rep(1:L, each=S)
nrow(state_nodes)

## ------------------------------------------------------------------------
## Transform the matrices to an extended edge list
EEL_intra <- c() # EEL stands for extended edge list
for (l in 1:L){
  x <- data_matrices[[l]]
  rownames(x) <- NULL
  colnames(x) <- NULL
  g <- graph.incidence(x, directed = F, weighted = T) # Note the graph.incidence instead of graph.adjacency!!
  x <- igraph::as_data_frame(g, what='edges')
  edge_list_layer <- cbind(x[,1], l, x[,2], l, x[,3]) # This builds the format: node_from layer_from node_to layer_to weight
  EEL_intra <- rbind(EEL_intra, edge_list_layer) # Add the layer to the bottom of the extended edge list
}
colnames(EEL_intra) <- c('node_from', 'layer_from', 'node_to', 'layer_to', 'weight')
EEL_intra <- as.data.frame(EEL_intra)
head(EEL_intra)

sort(unique(EEL_intra[,1])) # The hosts should be numbered 1-22
sort(unique(EEL_intra[,3])) # The parasites should be numbered 23-78


## ------------------------------------------------------------------------
# This is for the 'from' nodes
x <- EEL_intra[,c(1,2)]
names(x)[1:2] <- c('node','layer')
y <- state_nodes[,c(1,5,4)]
names(y)[1:2] <- c('node','layer')
z <- left_join(x, y)
# This is for the 'to' nodes
x <- EEL_intra[,c(3,4)]
names(x)[1:2] <- c('node','layer')
y <- state_nodes[,c(1,5,4)]
names(y)[1:2] <- c('node','layer')
w <- left_join(x, y)
EEL_intra_state_nodes <- data.frame(from=as.numeric(z[,3]), to=as.numeric(w[,3]), weight=EEL_intra$weight)
SAM <- Matrix::sparseMatrix(i=EEL_intra_state_nodes$from, j=EEL_intra_state_nodes$to, x = EEL_intra_state_nodes$weight, dims=c(S*L,S*L))
class(SAM)
SAM[80:85,134:139]

## ------------------------------------------------------------------------
## Hosts
### Use relative changes in species abundance as interlayer edge weights
mat <- hostAbundYear
tot <- ncol(mat)
interlayerEdgesHost <- matrix(0,nrow(mat),tot-1,dimnames=list(rownames(mat),colnames(mat)[-1]))
for (x in 1:(tot-1)){
  interlayerEdgesHost[,x] <- mat[,x+1]/mat[,x]
}
# By convention NA and Inf are 0 (no interlayer edges)
interlayerEdgesHost[is.nan(interlayerEdgesHost)] <- 0
interlayerEdgesHost[is.infinite(interlayerEdgesHost)] <- 0

## Transform to an extended edge list
x <- melt(interlayerEdgesHost)
x <- with(x, x[order(Var1,Var2),])
interlayerEdgesHost_EEL <- data.frame(node_from=physical_nodes$id[match(x$Var1, physical_nodes$species)],
                                      layer_from=rep(1:(L-1),H),
                                      node_to=physical_nodes$id[match(x$Var1, physical_nodes$species)],
                                      layer_to=rep(2:L,H),
                                      weight=x$value)
## Parasites
mat <- t(parasiteAbundanceYear)
tot <- ncol(mat)
interlayerEdgesParas <- matrix(0,nrow(mat),tot-1,dimnames=list(rownames(mat),colnames(mat)[-1]))
for (x in 1:(tot-1)){
  interlayerEdgesParas[,x] <- mat[,x+1]/mat[,x]
}
# By convention NA and Inf are 0 (no interlayer edges)
interlayerEdgesParas[is.nan(interlayerEdgesParas)] <- 0
interlayerEdgesParas[is.infinite(interlayerEdgesParas)] <- 0
## Transform to an extended edge list
x <- melt(interlayerEdgesParas)
x <- with(x, x[order(Var1,Var2),])
interlayerEdgesParas_EEL <- data.frame(node_from=physical_nodes$id[match(x$Var1, physical_nodes$species)],
                                      layer_from=rep(1:(L-1),P),
                                      node_to=physical_nodes$id[match(x$Var1, physical_nodes$species)],
                                      layer_to=rep(2:L,P),
                                      weight=x$value)

EEL_inter <- rbind(interlayerEdgesHost_EEL,interlayerEdgesParas_EEL)

## ------------------------------------------------------------------------
# This is for the 'from' nodes
x <- EEL_inter[,c(1,2)]
names(x)[1:2] <- c('node','layer')
y <- state_nodes[,c(1,5,4)]
names(y)[1:2] <- c('node','layer')
z <- left_join(x, y)
# This is for the 'to' nodes
x <- EEL_inter[,c(3,4)]
names(x)[1:2] <- c('node','layer')
y <- state_nodes[,c(1,5,4)]
names(y)[1:2] <- c('node','layer')
w <- left_join(x, y)
EEL_inter_state_nodes <- data.frame(from=as.numeric(z[,3]), to=as.numeric(w[,3]), weight=EEL_inter$weight)
SAM <- Matrix::sparseMatrix(i=EEL_inter_state_nodes$from, j=EEL_inter_state_nodes$to, x = EEL_inter_state_nodes$weight, dims=c(S*L,S*L))
class(SAM)
isSymmetric(SAM)
sum(SAM[lower.tri(SAM)])
SAM[1:10,S:85]

## ------------------------------------------------------------------------
#I adapted some of the code from http://www.blopig.com/blog/2016/10/plotting-and-storing-a-3d-network-in-r/).
circpos=function(n,r=1){#Coordinates on a circle
  rad=seq(0,2*pi,length.out=n+1)[-1];x=cos(rad)*r;y=sin(rad)*r
  return(cbind(x,y))
}
EEL_intra_state_nodes$col <- 'black'
EEL_intra_state_nodes$edge_type <- 'intra'
EEL_inter_state_nodes$edge_type <- 'inter'
EEL_inter_state_nodes$col <- 'gray'
edges <- rbind(EEL_intra_state_nodes, EEL_inter_state_nodes)

edges <- subset(edges, weight>0)

state_nodes$nodecol <- rep(c(rep('blue',H),rep('red',P)),L)
state_nodes <- with(state_nodes, state_nodes[order(layer,state_id),]) # This is important so the order of plotting will be by layer instead of by node!

# Make the network in igraph. For that the state_id needs to be the first column.
state_nodes <- state_nodes[,c('state_id','id','layer','species','type','nodecol')]
g <- graph.data.frame(edges, directed = T, vertices = state_nodes)

lay=rbind(cbind(circpos(table(V(g)$layer)[1],r=5), runif(n = table(V(g)$layer)[1],1,1.5)),
          cbind(circpos(table(V(g)$layer)[2],r=5), runif(n = table(V(g)$layer)[2],3,3.5)),
          cbind(circpos(table(V(g)$layer)[3],r=5), runif(n = table(V(g)$layer)[3],5,5.5)),
          cbind(circpos(table(V(g)$layer)[4],r=5), runif(n = table(V(g)$layer)[4],7,7.5)),
          cbind(circpos(table(V(g)$layer)[5],r=5), runif(n = table(V(g)$layer)[5],9,9.5)),
          cbind(circpos(table(V(g)$layer)[6],r=5), runif(n = table(V(g)$layer)[6],11,11.5))
)
# rgl.open()
# rglplot(g, layout=lay,
#         vertex.size=5,
#         # vertex.label.color='white',
#         vertex.label=NA,
#         vertex.color=V(g)$nodecol,
#         edge.color=E(g)$col,
#         edge.weight=E(g)$weight*50)


