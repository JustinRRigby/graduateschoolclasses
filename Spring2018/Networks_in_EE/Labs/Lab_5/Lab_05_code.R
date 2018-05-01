## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, message=FALSE, warning=FALSE)

## ---- eval=FALSE, echo=TRUE, message=FALSE, warning=FALSE----------------
## install.packages('devtools')
## library(devtools)
## install_github('quicklizard99/Basingstoke')

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------
library(Basingstoke)
library(igraph)
library(bipartite)
library(RColorBrewer)
library(vegan)
# library(cheddar)
library(gtools)
library(tidyverse)
library(reshape2)

# We will need this function later. It is based on the code from lab 1.
plot_web_trophic <- function(g){
  basal <- which(igraph::degree(g, mode = 'in') == 0) # Basal species do not have ingoing links
  top <- which(igraph::degree(g, mode = 'out') == 0) # Top species do not have outgoing links
  interm <- V(g)[which(!V(g) %in% c(basal,top))] # Intermediate are all the rest
  
  V(g)$troph_pos <- rep(0,length(V(g)))
  V(g)$troph_pos[which(V(g)$name %in% basal)] <- 1
  V(g)$troph_pos[which(V(g)$name %in% top)] <- 3
  V(g)$troph_pos[which(V(g)$name %in% interm)] <- 2
  # create a matrix forthe layout coordinates.
  coords <- matrix(nrow=length(V(g)), ncol=2) #
  # The x positions are randomly selected
  coords[,1] <- runif(length(V(g)))
  # The y positions are the trophoc positions
  coords[,2] <- V(g)$troph_pos
  par(mar=c(0,0,0,0))
  plot(g,layout=coords,
              vertex.color=V(g)$troph_pos,
              vertex.label=NA,
              vertex.size=8,
              edge.color='black',
              edge.arrow.size=.3,
              edge.width=.5)
}

## ------------------------------------------------------------------------
# Create 1 web with S species and a density of C
S <- 20 # Number of species
C <- 0.1 # Connectance (density)
web <- RandomLinks(1,S,C)
# Results is a list of "edge list" objects.
web <- web[[1]]
# Plot the web
g <- igraph::graph.data.frame(web)
graph.density(g)
plot(g)
# Now transform the edge-list to an adjacency matrix.
A <- igraph::get.adjacency(g, sparse = F)
A[mixedsort(rownames(A)),mixedsort(colnames(A))][1:5,1:5] # mixedsort comes from library gtools

## ------------------------------------------------------------------------
web <- matrix(data = rbinom(S*S,1,C),nrow = S, ncol = S) # rbinom here generates either 0 or 1 with a probability of success for 1 equals C.
g <- igraph::graph.adjacency(web)
graph.density(g)

## ------------------------------------------------------------------------
webs <- CommunityFactory(S=S, n=1, generator=RandomLinks, C=C) # use this function to generate an object of class "community"
class(webs) # This is a list
# The community class has these features
names(webs$`Artificial community 1`)
# Now plot the matrix
PlotPredationMatrix(webs[[1]])

## ------------------------------------------------------------------------
web <- CascadeModelLinks(1,S,C)
# Results is a list of "edge list" objects.
web <- web[[1]]
# Plot the web
g <- igraph::graph.data.frame(web)
graph.density(g)
plot(g)
# Now transform the edge-list to an adjacency matrix.
A <- igraph::get.adjacency(g, sparse = F)
#A[mixedsort(rownames(A)),mixedsort(colnames(A))] # mixedsort comes from library gtools

## ------------------------------------------------------------------------
web <- matrix(0, ncol=S, nrow=S)
num_links_allowed <- S*(S-1)/2
length(web[upper.tri(web)])==num_links_allowed 
web[upper.tri(web)] <- 2*C*S/(S - 1) > runif(num_links_allowed)
heatmap(web, Rowv = NA, Colv = NA, symm = F, scale = 'none', revC = T, col=brewer.pal(2,"RdBu"))

## ------------------------------------------------------------------------
webs <- CommunityFactory(S=S, n=1, generator=CascadeModelLinks, C=C) # use this function to generate an object of class "community"
class(webs) # This is a list
# The community class has these features
names(webs$`Artificial community 1`)
# Now plot the matrix
PlotPredationMatrix(webs[[1]])

## ------------------------------------------------------------------------
# Generate 4 communities with 3 basal species
webs <- CommunityFactory(S=S, n=4, C=C, generator=CascadeModelLinks,
                         accept=function(z) sum(IsBasalNode(z))==3)
# All have 3 basal species
sapply(webs, function(z) length(BasalNodes(z)))
par(mfrow=c(2,2))
for (i in 1:4){PlotPredationMatrix(webs[[i]])}

## ------------------------------------------------------------------------
g <- graph.data.frame(webs[[1]]$trophic.links, vertices = webs[[1]]$nodes)
plot(g)

## ------------------------------------------------------------------------
web <- NicheModelLinks(1,S,C)[[1]]
PlotPredationMatrix(webs[[1]])
# Can also plot with igraph
web <- as.data.frame(web)
web <- with(web, web[order(resource,consumer),])
nodes <- data.frame(name=1:S)
g <- graph.data.frame(web, vertices = nodes)
plot_web_trophic(g)

## ------------------------------------------------------------------------
S=100
C=0.15
n_i <- sort(runif(S, 0, 1)) # Positions on a niche axis between 0 and 1.
mean(n_i) # The expected mean of n_i should be 0.5
beta <- (S-1)/(2*S*C)-1
r_i <- rbeta(S, 1, beta) # Draw r_i values from a beta distribution with alpha=1. 
mean(r_i) # The expected mean of r_i should be 2*C
r_i <- r_i*n_i # The final size of r_i is obtained by multiplying by n_i to obtain C
mean(r_i)

r_i[1] <- 0 # The species with the smallest n_i (which is the first in the vector) has r_i = 0 to have at least one basal species

# uniformly draw the centre of the range (c_i) from [r_i/2, n_i]
# c_i <- t(runif(S, min=r_i/2, max=n_i)) # This would not work because the basal species would be always disconnected
c_i <- matrix(runif(S, min = r_i/2, max = pmin(n_i, 1 - r_i/2)), ncol = S)

# A species i consumes whatever is in the range [c_i - r_i/2, c_i + r_i/2] 
consume <- outer(n_i, c_i - r_i/2, ">=") & outer(n_i, c_i + r_i/2, "<=")
sum(consume)/S^2 # Do we recover C?
# Make the matrix
web <- matrix(0,S,S)
web[which(consume==T)] <- 1
# Plot
heatmap(web, Rowv = NA, Colv = NA, symm = F, scale = 'none', revC = T, col=brewer.pal(2,"RdBu"))

# Can also plot with igraph
colnames(web) <- rownames(web) <- 1:S
g <- graph.adjacency(web)
plot_web_trophic(g)

## ------------------------------------------------------------------------
generate_group_model <- function(groups, P, S){
    ## Pij is the probability of connecting a member of group i with a member of group j
    ## S is the total number of species
    ## groups is a vector of group memberships  
    M <- matrix(rbinom(S^2, size = 1, prob = P[groups, groups]), S, S)
    ## sort elements in M by groups
    M <- M[order(groups), order(groups)]
    return(M)
}

S <- 10
groups <- rbinom(S, size = 2, 0.5)
m <- max(groups)
P <- matrix(runif((m + 1)^2), m + 1, m + 1)

M <- generate_group_model(groups, P, S)
groups <- groups[order(groups)] + 1
## Make a plot by groups
Mplot <- tbl_df(melt(t(M))) %>% mutate(Var2 = S - Var2 + 1, group1 = groups[Var1], group2 = groups[Var2], tag = (group1-1) * (m + 1) + group2)

Mplot %>% ggplot(aes(x = Var1, y = Var2, fill = as.character(tag))) + geom_tile(alpha = 0.6) + theme_bw() + geom_point(size = 5, alpha = 0.5, aes(fill = as.factor(value), colour=  as.factor(value))) + scale_colour_manual(values = c("white", "black")) + theme(legend.position = "none")
    

## ------------------------------------------------------------------------
S=20
L=50
C=L/S^2
web <- matrix(0, S, S)
web[upper.tri(web)] <- c(rep(1, L), rep(0, (S^2-S)/2-L))[order(runif((S^2-S)/2))]
dimnames(web) <- list(1:length(web[,1]), 1:length(web[,1]))
heatmap(web, Rowv = NA, Colv = NA, symm = F, scale = 'none', revC = T, col=brewer.pal(2,"RdBu"))

## ------------------------------------------------------------------------
calculateProperties <- function(g){
  A <- as_adjacency_matrix(g, type = 'both', sparse = F)
  L <- sum(A!=0)
  S <- nrow(A)
  C <- L/S^2
  #C <- graph.density(g) # Also possible
  GenSD <- sd(1/(L/S)*rowSums(A)) # Generality of species i
  VulSD <- sd(1/(L/S)*colSums(A)) # Vulnerability of species i
  
  Top <- sum(rowSums(A)==0)/vcount(g) # Top species do not have outgoing links
  Bas <- sum(colSums(A)==0)/vcount(g) # Basal species do not have ingoing links
  Int <- 1-(Top+Bas) # Intermediate are all the rest
  
  basal_species <- which(colSums(A) == 0)
  chain_length_mean <- mean(distances(g, v=basal_species))
  chain_length_sd <- sd(distances(g, v=basal_species))
  Cannib <- sum(diag(A)>0)/S
  
  return(data.frame(L=L,C=C,GenSD=GenSD,VulSD=VulSD,Top=Top,Bas=Bas,Int=Int,ChnLg=chain_length_mean,ChnSd=chain_length_sd,Cannib=Cannib))
}

## ------------------------------------------------------------------------
otago_nodes <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Nodes.csv?token=ADd8OEIltgIeGLGKTneIrBNkopbsObJXks5a4hKrwA%3D%3D')
otago_links <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Links.csv?token=ADd8OAomGiQRFtpaN-AFt2ePFIEa8ZPQks5a4hPRwA%3D%3D')
otago_web <- graph.data.frame(otago_links, vertices = otago_nodes, directed = T)

plot_web_trophic(otago_web)
empirical <- calculateProperties(otago_web)

## ----message=FALSE, warning=FALSE----------------------------------------
S=vcount(otago_web)
C=graph.density(otago_web)
N <- 100

## Random model
webs <- NULL
for (i in 1:N){
  webs[[i]] <- matrix(data = rbinom(S*S,1,C),nrow = S, ncol = S)
}
graphs_list <- NULL
for (i in 1:N){
  graphs_list[[i]] <- graph.adjacency(webs[[i]], mode = 'directed')
}
generated_properties_random <- lapply(graphs_list, calculateProperties)
generated_properties_random <- do.call(rbind, generated_properties_random)
generated_properties_random <- reshape2::melt(generated_properties_random)


## Cascade model
webs <- CommunityFactory(S=S, n=N, generator=CascadeModelLinks, C=C) # use this function to generate an object of class "community"
graphs_list <- NULL
for (i in 1:N){
  graphs_list[[i]] <- graph.data.frame(webs[[i]]$trophic.links, directed = T, vertices = webs[[i]]$nodes)
}
generated_properties_cascade <- lapply(graphs_list, calculateProperties)
generated_properties_cascade <- do.call(rbind, generated_properties_cascade)
generated_properties_cascade <- reshape2::melt(generated_properties_cascade)

## Niche model
webs <- CommunityFactory(S=S, n=N, generator=NicheModelLinks, C=C) # use this function to generate an object of class "community"
graphs_list <- NULL
for (i in 1:N){
  graphs_list[[i]] <- graph.data.frame(webs[[i]]$trophic.links, directed = T, vertices = webs[[i]]$nodes)
}
generated_properties_niche <- lapply(graphs_list, calculateProperties)
generated_properties_niche <- do.call(rbind, generated_properties_niche)
generated_properties_niche <- reshape2::melt(generated_properties_niche)


## ------------------------------------------------------------------------
# Data frame from the empirical value, to plot vertical lines on the histograms
intercept_data <- as.data.frame(t(empirical)[,1])
names(intercept_data)[1] <- 'intercept'
intercept_data$variable <- rownames(intercept_data)

# Plot
generated_properties_random$model <- 'Random'
generated_properties_cascade$model <- 'Cascade'
generated_properties_niche$model <- 'Niche'
d <- rbind(generated_properties_random,generated_properties_cascade,generated_properties_niche)
ggplot(d, aes(value, fill=model))+geom_histogram(alpha=0.3)+facet_wrap(~variable, scales='free')+geom_vline(data = intercept_data, aes(xintercept = intercept), color='red')+theme_classic()

## ------------------------------------------------------------------------
S=100
C=0.15
web <- matrix(data = rbinom(S*S,1,C),nrow = S, ncol = S)
L=sum(web>0)
p_max=L/S^2
LL_vector=NULL # Store results of log-likelihood
p_range <- seq(p_max-0.1,p_max+0.1,0.001)
for (p in p_range){
  LL <- L*log(p)+(S^2-L)*log(1-p)
  LL_vector <- c(LL_vector, LL)
}
qplot(x=p_range, y=LL_vector)+geom_vline(xintercept = p_max, color='red')
p_range[which.max(LL_vector)] # Connectance corresponding to the maximum log-likelihood

## ------------------------------------------------------------------------
# Generate a random network with the cascade model
S=100
C=0.15
web <- matrix(0, ncol=S, nrow=S)
num_links_allowed <- S*(S-1)/2
web[upper.tri(web)] <- 2*C*S/(S - 1) > runif(num_links_allowed)
L=sum(web>0) # Number of links in the generated network
C=L/S^2 # Connectance of the generatd network
p_max=2*C*S/(S - 1) # The value which maximizes the likelihood
LL_vector=NULL # Store results of log-likelihood
p_range <- seq(p_max-0.1,p_max+0.1,0.001) # A range to calculate the likelihoods
for (p in p_range){
  LL <- L*log(p)+(S*(S-1)/2-L)*log(1-p)
  LL_vector <- c(LL_vector, LL)
}
# d <- data.frame(p_range=p_range, LL_vector=LL_vector)
# d <- d %>% arrange(LL_vector, p_range)
qplot(x=p_range, y=LL_vector)+geom_vline(xintercept = p_max, color='red')
p_range[which.max(LL_vector)]-p_max # Connectance corresponding to the maximum log-likelihood

## ------------------------------------------------------------------------
S <- vcount(otago_web)
L <- ecount(otago_web)

# random
p <- L/S^2
LL_random <- L*log(p)+(S^2-L)*log(1-p)

# Modified cascade
A <- as_adjacency_matrix(otago_web, sparse=F)
L1 <- sum(A[upper.tri(A)]>0)
L2 <- sum(A[lower.tri(A)]>0)+sum(diag(A))
q1 <- 2*L1/(S*(S-1))
q2 <- 2*L2/(S*(S+1))
LL_modified1 <- L1*log(q1)+(S*(S-1)/2-L1)*log(1-q1)
LL_modified2 <- L2*log(q2)+(S*(S+1)/2-L2)*log(1-q2)
LL_modified <- LL_modified1+LL_modified2
data.frame(model=c('R','M'), loglik=c(LL_random,LL_modified)) %>% arrange(desc(loglik))

