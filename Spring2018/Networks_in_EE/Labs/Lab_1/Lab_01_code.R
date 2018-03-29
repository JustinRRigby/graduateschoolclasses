## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, out.width = '100%', out.height='40%')

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------
# install.packages('igraph')
library('igraph')
# install.packages('bipartite')
library('bipartite')

## ----echo=TRUE-----------------------------------------------------------
?igraph::betweenness()

## ------------------------------------------------------------------------
##Undirected matrix with symmteric form (both triangles are equal to one another)
A_u <- matrix(c(0,1,1,0,0, # An example input matrix
              1,0,0,1,1,
              1,0,0,0,0,
              0,1,0,0,0,
              0,1,0,0,0),5,5, byrow=F)
isSymmetric(A_u)
g <- igraph::graph.adjacency(A_u, mode = 'undirected')
par(mar=c(0,0,0,0))
plot(g)

## ------------------------------------------------------------------------
A_d <- matrix(c(0,1,1,0,1, # An example input matrix
              1,0,0,1,1,
              1,0,0,0,0,
              0,1,0,0,0,
              0,1,1,0,0),5,5, byrow=F)
isSymmetric(A_d)
g <- igraph::graph.adjacency(A_d, mode = 'directed')
par(mar=c(0,0,0,0)) #Sets margins of the igraph plot
plot(g)
##Question 1, 2, and 4
A_d <- matrix(c(0,0,1,0,1, # An example input matrix
                1,0,0,0,0,
                1,0,0,0,0,
                0,0,0,1,0,
                0,1,1,0,0),5,5, byrow=F)
isSymmetric(A_d)
g <- igraph::graph.adjacency(A_d, mode = 'directed')
par(mar=c(0,0,0,0)) #Sets margins of the igraph plot
plot(g)



## ------------------------------------------------------------------------
##Amount of biomass obtained. Pollinator networks (Amount of biomass lost by vistors or the probability fo visitation.)
A_w <- matrix(c(0,1,1,0,0, # An example input matrix
              1,0,0,1,1,
              1,0,0,0,0,
              0,1,0,0,0,
              0,1,0,0,0),5,5, byrow=F)
random_weights <- round(rnorm(10, 10, 4),2) # take weights from a normal distribution.
A_w[lower.tri(A_w)] <- A_w[lower.tri(A_w)]*random_weights # Fill the lower traiangle
A_w <- A_w+t(A_w) # This makes the matrix symmetric
isSymmetric(A_w)
g <- igraph::graph.adjacency(A_w, mode = 'undirected', weighted = T)
E(g)$weight ##E() is to measure the attribute of the matrix -- i.e. the weights (edge.attributes() gives all attributes.)
par(mar=c(0,0,0,0))
plot(g, edge.width=E(g)$weight)

##Question 4

A_w <- matrix(c(0,1,1,0,0, # An example input matrix
                1,0,0,1,1,
                1,0,0,0,0,
                0,1,0,0,0,
                0,1,0,0,0),5,5, byrow=F)
random_weights <- round(rnorm(10, 10, 4),2) # take weights from a normal distribution.
A_w[lower.tri(A_w)] <- A_w[lower.tri(A_w)]*random_weights # Fill the lower traiangle
#A_w <- A_w+t(A_w) # This makes the matrix symmetric
isSymmetric(A_w)
g <- igraph::graph.adjacency(A_w, mode = 'directed', weighted = T)
E(g)$weight ##E() is to measure the attribute of the matrix -- i.e. the weights (edge.attributes() gives all attributes.)
par(mar=c(0,0,0,0))
plot(g, edge.width=E(g)$weight, edge.arrow.width = 1.9)

## ------------------------------------------------------------------------
L_u <- data.frame(i=c(1,1,2,2),
                j=c(2,3,4,5))
L_u
g <- igraph::graph.data.frame(L_u, directed = F) ## instead of 'adjancency' we use a data frame
par(mar=c(0,0,0,0))
plot(g)

## ------------------------------------------------------------------------
L_u <- data.frame(i=c(1, 1, 2, 2, 2, 3, 3, 4, 5, 5), # i = from
                j=c(2, 3, 1, 4, 5, 1, 5, 2, 1, 2)) # j = to
g <- igraph::graph.data.frame(L_u, directed = T)
par(mar=c(0,0,0,0))
plot(g)

## ------------------------------------------------------------------------
L_w <- data.frame(i=c(1,1,2,2),
                j=c(2,3,4,5),
                weight=round(rnorm(4, 10, 4),2) # take weights from a normal distribution.
                )
L_w
g <- igraph::graph.data.frame(L_w, directed = F)
E(g)$weight
par(mar=c(0,0,0,0))
plot(g, edge.width=E(g)$weight)

## ------------------------------------------------------------------------
L_wd <- data.frame(from=c(1, 1, 2, 2, 2, 3, 3, 4, 5, 5),
                to=c(2, 3, 1, 4, 5, 1, 5, 2, 1, 2),
                weight=round(rnorm(10, 1, 0.2),2))
g <- igraph::graph.data.frame(L_wd, directed = T)
g
E(g)$weight
par(mar=c(0,0,0,0))
plot(g, edge.width=log(E(g)$weight)*10, # possible to rescale edge weights when plotting 
     vertex.label = NA,  
     edge.arrow.size=1.2,
       edge.curved=0.5,
     edge.color='black')

## ------------------------------------------------------------------------
A_w
g <- igraph::graph.adjacency(A_w, mode = 'directed', weighted = T)
L <- igraph::as_data_frame(g, what = 'edges')
L

##Question 5: Create a function to change a matrix into a list
A_w
A_w <- as.data.frame(A_w)
A_w



## ------------------------------------------------------------------------
L_wd
g <- igraph::graph.data.frame(L_wd, directed = T)
g
A <- igraph::as_adjacency_matrix(g, attr = 'weight', sparse=F)
A

## ------------------------------------------------------------------------
data("memmott1999") # load
class(memmott1999)
memmott1999[1:4,1:4] # view first lines

## ------------------------------------------------------------------------
visweb(memmott1999)

## ------------------------------------------------------------------------
memmott1999_binary <- 1*(memmott1999>0) ##Creating a binary matrix in order to see where there are edges.
visweb(memmott1999_binary)

## ------------------------------------------------------------------------
visweb(memmott1999,prednames = F, prey.lablength = 10)

## ------------------------------------------------------------------------
plotweb(memmott1999)

## ------------------------------------------------------------------------
ural_data <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Ural_valley_A_HP_048.csv?token=ADd8OH8_dgTW2e0V71XhqCX5EXsT0i2mks5awrkRwA%3D%3D')

## ------------------------------------------------------------------------
ural_data[1:4,1:4]

## ------------------------------------------------------------------------
rownames(ural_data) <- ural_data[,1] # Set row names
num_hosts_sampled <- ural_data[,2] # save in a variable
ural_data <- ural_data[,-2] # remove column
ural_data <- ural_data[,-1] # remove column
class(ural_data) # This is a data frame!
ural_data <- data.matrix(ural_data) # Transform to a matrix format
ural_data[1:4,1:4]

## ------------------------------------------------------------------------
plant_species <- rownames(memmott1999)
flower_visitor_species <- colnames(memmott1999)
head(plant_species, 3)
head(flower_visitor_species, 3)

## ------------------------------------------------------------------------
otago_nodes <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Nodes.csv?token=ADd8OGnaV0qMolv_YslyRlkYegPJVjLxks5awrkjwA%3D%3D')
otago_links <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Links.csv?token=ADd8OH_GXVK3HoD_9YVCFGkAXp0PrUAEks5awrk0wA%3D%3D')
otago_nodes[1:4,1:6]
otago_links[1:4,1:8]

## ------------------------------------------------------------------------
otago_web <- graph.data.frame(otago_links, vertices = otago_nodes, directed = T)

## ------------------------------------------------------------------------
names(edge.attributes(otago_web))
unique(E(otago_web)$LinkType)

## ------------------------------------------------------------------------
names(vertex.attributes(otago_web))
head(unique(V(otago_web)$name)) #Command 'V' gives all attributes. 

## ------------------------------------------------------------------------
par(mar=c(0,0,0,0)) #Reduce margin size
plot(otago_web)

## ------------------------------------------------------------------------
par(mar=c(0,0,0,0))
plot(otago_web, vertex.size=3, edge.arrow.size=0.4, vertex.label=NA, layout=layout.circle)

## ------------------------------------------------------------------------
E(otago_web)$color <- "grey" # First, we set a default color
E(otago_web)[otago_links$LinkType == 'Predation']$color <- "black"
E(otago_web)[otago_links$LinkType == 'Macroparasitism']$color <- "blue"
E(otago_web)[otago_links$LinkType == 'Trophic Transmission']$color <- "red"

# Now plot
par(mar=c(0,0,0,0))

plot(otago_web, vertex.size=2, edge.arrow.size=0.2, vertex.label=NA, layout=layout.circle)

## ------------------------------------------------------------------------
# Basal species (those that do not consume) -- do not have incoming links
basal <- which(igraph::degree(otago_web, mode = 'in') == 0)
# Top species do not have outgoing links
top <- which(igraph::degree(otago_web, mode = 'out') == 0)
# Intermediate are all the rest
interm <- V(otago_web)[which(!V(otago_web) %in% c(basal,top))]
# Are all the nodes included?
all(c(basal,top,interm) %in% V(otago_web))
all(V(otago_web) %in% c(basal,top,interm))

## ------------------------------------------------------------------------
V(otago_web)$troph_pos <- rep(0,length(V(otago_web)))
V(otago_web)$troph_pos[which(V(otago_web)$name %in% basal)] <- 1
V(otago_web)$troph_pos[which(V(otago_web)$name %in% top)] <- 3
V(otago_web)$troph_pos[which(V(otago_web)$name %in% interm)] <- 2
# create a matrix forthe layout coordinates.
coords <- matrix(nrow=length(V(otago_web)), ncol=2) #
# The x positions are randomly selected
coords[,1] <- runif(length(V(otago_web)))
# The y positions are the trophoc positions
coords[,2] <- V(otago_web)$troph_pos
par(mar=c(0,0,0,0))
plot(otago_web,layout=coords,
            vertex.color=V(otago_web)$troph_pos,
            vertex.label=NA,
            vertex.size=8,
            edge.color='black',
            edge.arrow.size=.3,
            edge.width=.5)

