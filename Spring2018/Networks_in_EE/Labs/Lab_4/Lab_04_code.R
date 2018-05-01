## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, out.width = '100%', out.height='40%')

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------
library(bipartite) # will also load vegan
library(ggplot2)

## ------------------------------------------------------------------------
data("memmott1999")
memmott1999_binary <- 1*(memmott1999>0)
nodf <- nestednodf(memmott1999_binary)
nodf
plot(nodf)

## ------------------------------------------------------------------------
nodf_eval_1 <- oecosimu(memmott1999_binary, nestednodf, "r00", nsimul = 100) #ecological simulation
nodf_eval_1$statistic # This is the NODF of the observed network
names(nodf_eval_1$oecosimu) # look at the ?oecosimu help for details on what are these values.
nodf_eval_1 # The networks is significantly nested.

## ------------------------------------------------------------------------
data("memmott1999")
visweb(memmott1999)
# Follow secondary extinction of pollinators, while removing plants.
ex_random <- second.extinct(memmott1999, method = 'random', participant = 'lower', nrep=100, details=F) # Replicate simulations because it is random.
ex_random[1:3,]
ex_least <- second.extinct(memmott1999, method = 'abundance', participant = 'lower', details=F)
ex_least[1:3,]
ex_most <- second.extinct(memmott1999, method = 'degree', participant = 'lower', details=F)
ex_most[1:3,]









## ------------------------------------------------------------------------
par(mfrow=c(2,2))
slope.bipartite(ex_random, cex.axis=0.1) #fit a curve to the extinction data
slope.bipartite(ex_least) #fit a curve to the extinction data
slope.bipartite(ex_most) #fit a curve to the extinction data

## ------------------------------------------------------------------------
#The area curve is the measurement of the robustness of the network.
# Robustness
robustness(ex_random)
robustness(ex_least)
robustness(ex_most)



