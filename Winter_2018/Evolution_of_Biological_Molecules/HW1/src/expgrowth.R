file.out = F
output.type = 'pdf'

# Plot trajectories for exponential growth, given Poisson reproduction.

n0 <- 100 # initial number of individuals
max.n <- 1e6 # total (maximum) number of individuals before we quit (to guard against horrible slowdowns)
max.gen <- 100 #n0 # maximum number of generations we'll track
n.reps <- 100 # the number of replicate populations to run and accumulate into a single set of statistics
n.super.reps <- 1 # if you'd like to run multiple sets of n.reps replicates, and compute statistics on each one

# Growth rates
log.fitness <- 0 #1/n0
average.offspring <- exp(log.fitness) # average number of offspring per generation.

# Plot all the individual trajectories?
plot.raw.runs <- TRUE
show.statistics <- T

grow <- function(init.n, avg.offspring, max.gen, max.n=1e6) {
	n.array <- rep(NA,max.gen+1)
	n.array[1] <- init.n
    # Number of organisms
	n <- init.n
    # Generation number -- array index will be gen+1
	gen <- 0
    # While there are >0 organisms, fewer than max.gen generations, and fewer than max.n organisms...go.
	while (n > 0 & n < max.n & gen < max.gen) {
		gen <- gen+1
		# Sample a poisson for each individual, and sum to get the total population
		#n.next <- sum(rpois(n, avg.offspring))
		# This is mathematically equivalent, because the sum of poisson random variables is itself a poisson random variable, to:
		n.next <- rpois(1, n*avg.offspring)
		# ...and one should verify that these are precisely the same.
        # Add to gen+1 because gen is 0-based, but R is 1-based.
		n.array[gen+1] <- n.next
		n <- n.next
	}
	# If we bailed out early because we hit zero, fill in all the rest of the zeros.
	if (n==0) {
		n.array[(gen+1):(max.gen+1)] <- 0
	}
	n.array
}

# Run the simulation, creating all the replicates.
all.res <- lapply(1:n.super.reps, function(m) {replicate(n.reps, grow(n0, average.offspring, max.gen, max.n))})

# Merge the individual simulations into one giant list, which is useful for some purposes
res <- NULL
for (i in 1:n.super.reps) {
  res <- cbind(res, all.res[[i]])
}

# Compute statistics
stats <- data.frame(means=NULL, sds=NULL)
stats <- lapply(1:n.super.reps, function(rid) {
  r <- all.res[[rid]]
  means <- rowMeans(r, na.rm=F)
  sds <- apply(r, 1, sd, na.rm=F)
  list(means=means, sds=sds)
})

# Plot
# Want to show zeros on a log scale
bump <- 1

if (file.out) dev.out("exp.growth", output.type=output.type)
split.screen(c(1,2))
screen(1)
par(mar=c(4,4,2,1)) # Set plot margins
ylim=c(1, max(sapply(res,max)))
if (!plot.raw.runs) {
  max.fluct <- max(sapply(stats, function(m){max(m$means+m$sds, na.rm=TRUE)}), na.rm=TRUE)
  ylim = c(1, max.fluct)
}

# We want to distinguish replicates, so add a little noise to the data so that we can see them.
jitter <- replicate(ncol(res), rnorm(nrow(res), sd=ylim[1]/100))

# Time course over which we'll plot
t <- 0:max.gen

plot(1:10, 1:10, xlim=c(0,max.gen), ylim=ylim, type='n', log='y', las=1, 
	 xlab='Generations', ylab='Population size', yaxt='n', 
	 main='Population trajectories')
abline(v=n0)
abline(h=n0+bump)
axis(2, at=c(0,1,10,100,1000,10000)+bump, labels=c(0,1,10,100,1000,10000), las=1)
if (plot.raw.runs) {
  matplot(t, res+bump+jitter, type='l', add=TRUE, col='gray70', lty='solid', lwd=0.5)
}
# Statistics
if (show.statistics) {
	for (nrep in 1:n.super.reps) {
	  means <- stats[[nrep]]$means
	  sds <- stats[[nrep]]$sds
	  lines(t, means+bump, col='blue', lwd=2)
	  lines(t, means+sds+bump, col='forestgreen', lwd=2)
	  lines(t, means-sds+bump, col='forestgreen', lwd=2)
	}
}
## An exact expression
predicted.means <- n0*exp(log.fitness*t)
lines(t, predicted.means+bump, type='l', col='red', lwd=2)

screen(2)
par(mar=c(4,4,2,1)) # Set plot margins
## Fraction of trajectories that have died out after max.gen generations
end.res <- res[max.gen+1,]
barplot(c(length(end.res[end.res==0]), length(end.res[end.res>0])), names.arg=c('0', '>0'), main='Individuals remaining in population', las=1, ylab='Count')
cat("Fraction of dead populations after", max.gen, "generations =", length(end.res[end.res==0])/n.reps, '\n')
close.screen(all=TRUE)
if (file.out) dev.off()
