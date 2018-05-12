## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = TRUE, fig.height=6, fig.width=6)

## ----load libraries, echo=TRUE, message=FALSE, warning=FALSE-------------




library(tidyverse)
library(MASS)
library(RColorBrewer)

## theme for ggplot

mytheme <- theme_bw() + theme(
  legend.title  = element_text( size=17),
  #  legend.position = "bottom",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text( size=17),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=19),
  strip.text.x = element_text(family = "Helvetica", size = 10),
  strip.text.y = element_text(family = "Helvetica", size = 10),
  panel.border = element_rect( colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)


## ----Mays results--------------------------------------------------------
# Here sigma=1
S <- 500
M <- matrix(runif(S^2, min = -sqrt(3), max = sqrt(3)), S, S) / sqrt(S) # Note that the division by sqrt(S) normalizes the eigenvalues to fall on the range [-2,0]
diag(M) <- -1
e <- eigen(M, only.values = TRUE)$values
data_eigenvalues <- data_frame(x = Re(e), y = Im(e), S = S)

data_eigenvalues %>% ggplot(aes(x = x, y = y)) + geom_point() + mytheme

## ----put circle----------------------------------------------------------
normalized_expectation_girko <- data_frame(x = rep(seq(-2,0,0.001), 2) , type = rep(c(1, -1), each = length(x)/2), y = sqrt(1 - (x+1)^2) * type)
data_eigenvalues %>% ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_path(data = normalized_expectation_girko, colour = "red") +
  mytheme +
  xlab(expression(Re(lambda))) +
  ylab(expression(Im(lambda)))

## ----Role of sigma and d-------------------------------------------------
S <- 500
sigma_seq <- c(0.5, 1, 1.5)
diag_seq <- c(-3, 0, 3)

# This function creates a random matrix with a given sigma and d, and returns the ev.
get_eigenvalues_sigma_d <- function(sigma, d){
  M <- matrix(runif(S^2, min = -sqrt(3)*sigma, max = sqrt(3)*sigma), S, S)/sqrt(S)
  M <- M + diag(S) * d
  e <- eigen(M, only.values = TRUE)$values
  return(data_frame(x = Re(e), y = Im(e), S = S, d = d, sigma = sigma))
}

data_ev_sigma_d <- NULL
for (sigma in sigma_seq){
  for (d in diag_seq){
    data_ev_sigma_d <- rbind(data_ev_sigma_d, get_eigenvalues_sigma_d(sigma, d))
  }
}
data_ev_sigma_d %>% arrange(desc(sigma)) %>% ggplot(aes(x = x, y = y, colour = factor(sigma))) +
  geom_point(size = 2.5, alpha = 0.5) +
  mytheme +
  xlab(expression(Re(lambda))) +
  ylab(expression(Im(lambda))) +
  scale_colour_discrete(name = expression(sigma)) +
  coord_fixed() +
  facet_grid(d~.) +
  scale_x_continuous(breaks = c(-3, 0, 3)) +
  geom_vline(xintercept = c(-1.5,1.5), linetype='dashed')


## ----Probability of stability--------------------------------------------
S_seq <- c(25, 50, 100) # A curve per network size
nsim <- 200 # Generate multiple random matrices for each S
sigma_seq <- seq(0.8, 1.2, 0.01) # vary sigma
C <- 0.15

# Function to create a matrix and test if its righmost ev is negative.
stable_matrix <- function(i){
  M <- matrix(runif(S^2, min = -sqrt(3)*sigma, max = sqrt(3)*sigma) * rbinom(S^2, size = 1, prob = C), S, S)
  diag(M) <- 0
  M <- M + diag(S) * d
  e <- eigen(M, only.values = TRUE)$values
  return (max(Re(e)) < 0) # The rightmost eigenvalue should be negative for the matrix to be locally-stable
}

prob_stability_df <- NULL
for (S in S_seq){
  print(S)
  d <- -sqrt(S * C) # Set the diagonal to sqrt(SC) to make sigma=1
  for (sigma in sigma_seq){
    p <- 0
    for (i in 1:nsim){
      p <- p+stable_matrix()
    }
    p_stability <- p/nsim
    prob_stabili %>% ty_df <- rbind(prob_stability_df, data.frame(S=S, sigma=sigma, p_stability=p_stability))
  }
}

prob_stability_df %>% ggplot(aes(x = sigma, y = p_stability, colour = factor(S), shape = factor(S))) +
  geom_point(size = 3) +
  geom_line() +
  mytheme +
  xlab(expression(sigma)) +
  scale_colour_discrete(name = "S") +
  scale_shape_discrete(name = "S") +
  geom_vline(xintercept = 1, linetype='dashed')

## ----The role of mean----------------------------------------------------
S <- 500
mu_seq <- c(-0.1, 0.01, 0.1)

outliers <- data_frame(x = (S - 1) * mu_seq, y = 0, mu = mu_seq) # Estimation of where the ev will fall

data_eigenvalues_mean <- map(mu_seq,
                             function(mu){
                               M <- matrix(runif(S^2, min = -sqrt(3), max = sqrt(3)), S, S)
                               M <- M + mu - diag(S) * mu
                               e <- eigen(M, only.values = TRUE)$values
                               data_frame(x = Re(e), y = Im(e), S = S, mu = mu)
                             }) %>% bind_rows()


data_eigenvalues_mean %>% ggplot(aes(x = x, y = y)) + geom_point(size = 1.3, alpha = 0.5) + mytheme + xlab(expression(Re(lambda))) + ylab(expression(Im(lambda))) + facet_grid(mu ~.) + coord_fixed() + geom_point(data = outliers, size = 2, colour = "red")

## ----Role of correlation-------------------------------------------------
S <- 300
rho_seq <- c(-0.75, 0, 0.75)

generate_eigenvalues_cor <- function(rho){
  M <- matrix(0, S, S)
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  M <- diag(rnorm(S))
  for(i in seq(1, S - 1)) {
    for (j in seq(i + 1, S)){
      Z <- mvrnorm(1, c(0,0), Sigma) # Now sample from a bivariate normal distribution.
      M[i,j] <- Z[1]
      M[j,i] <- Z[2]
    }
  }
  M <- M / sqrt(S) # Rescale EV
  e <- eigen(M, only.values = TRUE)$values
  data_frame(x = Re(e), y = Im(e), S = S, rho = rho)
}

data_eigenvalues_cor <- map(rho_seq, generate_eigenvalues_cor) %>%bind_rows()

data_eigenvalues_cor %>% ggplot(aes(x = x, y = y, colour = factor(rho))) +
  geom_point(size = 2.5, alpha = 0.5) +
  mytheme +
  xlab(expression(Re(lambda))) +
  ylab(expression(Im(lambda))) +
  scale_colour_discrete(name = expression(rho)) +
  coord_fixed() +
  scale_x_continuous(breaks = c(1 + rho_seq, -1 - rho_seq)) +
  scale_y_continuous(breaks = c(1 - rho_seq, -1 + rho_seq))

## ----Role of structure---------------------------------------------------
S <- 200
C <- 0.4
alpha <- 1/4
m <- S * alpha  ## size of the top module
sigma <- 1
mu_seq <- c(-1/4, 0, 1/4)
#n <- length(mu_seq)
rho_seq <- c(-1/4, -3/4, -0.5)
Q_seq <- c(-0.5, 0, 0.35) # This is the modularity.

# K <- matrix(0, S, S) ## adjacency matrix
# M <- matrix(0, S, S) ## Community matrix


generate_K <- function(Q){
  K <- matrix(0, S, S) ## adjacency matrix
  Cw <- C * (1 + Q/(alpha^2 + (1-alpha)^2)) # within block connectance
  Cb <- C * (1 - Q/(2 * alpha * (1-alpha))) # between block conenctnve
  # Between-block interactions
  AB <- matrix(runif(m * (S-m)) < Cb, m , S - m) 
  # Within-block interactions
  K[1:m, 1:m] <- matrix(runif(m^2) < Cw, m, m)
  K[(m + 1):S, (m+1):S] <- matrix(runif((S - m)^2) < Cw, S - m, S - m)
  # Add the between-block interactions to K
  K[1:m, (m+1):S] <- AB
  K[(m+1):S, 1:m] <- t(AB)
  return(K)
}

generate_M <- function(mu, rho, K, Q){
  M <- matrix(0, S, S) ## Community matrix
  Sigma <- matrix(c(sigma^2, rho * sigma^2, rho * sigma^2, sigma^2), nrow = 2) # Covariance matrix for bivariate normal distribution
  # Fill in M with interaction values
  for(i in seq(1, S - 1)) {
    for (j in seq(i + 1, S)){
      if (K[i,j] > 0){
        Z <- mvrnorm(1, c(mu, mu), Sigma) #
        M[i,j] <- Z[1] 
        M[j,i] <- Z[2]
      }   
    }
  }
  return(M)
}

## ----plot heatmaps-------------------------------------------------------
# Plot the matrices
K <- generate_K(0.35)
heatmap(K, Rowv = NA, Colv = NA, symm = F, scale = 'none', revC = T, col=brewer.pal(2,"RdBu"))
M <- generate_M(0, 0.5, K)
heatmap(M, Rowv = NA, Colv = NA, symm = F, scale = 'none', revC = T)

## ----find stability, fig.height=4----------------------------------------
data_modules <- NULL
for (Q in Q_seq){
  K <- generate_K(Q)
  for (i in 1:length(mu_seq)){
    mu <- mu_seq[i]
    rho <- rho_seq[i]
    print (paste(Q, mu, rho))
    M <- generate_M(mu, rho, K, Q)
    # Get the eigenvalues of M
    e <- eigen(M, only.values = TRUE)$values
    data_modules <- rbind(data_modules, data_frame(x = Re(e), y = Im(e), S = S, mu = mu, rho = rho, Q = Q))
  }
}


data_modules %>% ggplot(aes(x = x, y = y, colour = factor(mu))) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  mytheme + xlab(expression(Re(lambda))) + 
  ylab(expression(Im(lambda))) + 
  scale_colour_discrete(name = expression(mu))+ 
  coord_fixed() +
  facet_grid(Q ~ mu, 
             labeller = label_bquote(cols = paste(mu," = ",.(mu), collapse = ""), 
                                     rows = paste(Q," = ",.(Q), collapse = "")))

## ----Empirical networks--------------------------------------------------
K <- as.matrix(read.table("../Network_data/littlerock_181.txt"))
K <- K - (t(K) + K == 2) ## remove loops of length 1 and 2
size <- nrow(K)
realizations <- 250
theta_seq <- runif(realizations, max = 2) ## get mean efficiency values

# This function comes from Stefano Allesina
RandomizeEmpirical <- function(M){
    S <- dim(M)[1]
    ## Shuffle rows and columns
    my_order <- sample(1:S)
    M <- M[my_order, my_order]
    ## Shuffle the pairs
    Mij <- M[upper.tri(M)]
    Mji <- t(M)[upper.tri(M)]
    my_order <- sample(1:(S * (S - 1) / 2))
    Mij <- Mij[my_order]
    Mji <- Mji[my_order]
    Mprime <- M
    Mprime[upper.tri(Mprime)] <- Mij
    Mprime <- t(Mprime)
    Mprime[upper.tri(Mprime)] <- Mji
    Mprime <- t(Mprime)
    ## Get eigenvalues
    ReL1 <- max(Re(eigen(M, only.values = TRUE)$values))
    ReL1tilda <- max(Re(eigen(Mprime, only.values = TRUE)$values))
    return(c(ReL1, ReL1tilda, ReL1 / ReL1tilda))
}

stability_ratio <- map_dbl (theta_seq, function(theta){
    M <- K * matrix(-abs(rnorm(size^2, -1, 0.2)), size, size) # This is negative becuae it is the effect of predator on prey.
    ## set positive part because it is the effect of prey on predator. The positive sign comes from the minus before the transpose.
    M <- M + -t(M) * matrix(runif(size^2, min=0, max=2*theta), size, size) # Sampling uniformly from 0 to 2\theta will have an expected value of theta (which is what we want for the mean efficeincy)
    L <- RandomizeEmpirical(M)
    return (L[3])
})

Data_empirical  <- data_frame(theta = theta_seq, stability_ratio = stability_ratio)
Data_empirical %>% ggplot(aes(x = theta, y = stability_ratio)) + geom_point(size = 3, alpha = 0.5) + mytheme + geom_hline(yintercept = 1, linetype = "dashed") + xlab(expression(Mean~efficiency~theta)) + ylab(expression(Re(lambda[M]^1)/Re(lambda[tilde(M)]^{1})))

## ------------------------------------------------------------------------
library(igraph)
otago_nodes <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Nodes.csv?token=ADd8OKe0yZq3ndjNT_M-TNSCgvpOpS9Yks5a9cNywA%3D%3D')
otago_links <- read.csv('https://raw.githubusercontent.com/pascualgroup/Networks_course/master/Network_data/Otago_Data_Links.csv?token=ADd8OHcXUTZl3fc_vxAI8au-GNen6wlwks5a9cNiwA%3D%3D')

