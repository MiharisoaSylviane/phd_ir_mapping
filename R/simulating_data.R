
## WHO bioassay data simulation
# ---- Beta-Binomial helpers ----

# Map (Q_star, rho) -> (alpha, beta), rho should be between 0 and 1, if rho = 0,
# it means overdispersion, if rho =1 it means extreme overdispersion
# same for Q* as it is a probability it should be between 0 and 1
# m <- alpha + beta and could be got also from m = 1/rho - 1,
# alpha (a) = Q_star * m, beta (b)= (1- Q*) * m

rho_to_ab <- function(Q_star, rho) {
  stopifnot(all(rho > 0 & rho < 1), all(Q_star > 0 & Q_star < 1))
  m <- 1 / rho - 1           # alpha + beta
  a <- Q_star * m
  b <- (1 - Q_star) * m
  list(a = a, b = b)
}


# Single draw from Beta-Binomial via Beta then Binomial
rbetabinom1 <- function(N, Q_star, rho) {
  ab <- rho_to_ab(Q_star, rho)
  theta <- stats::rbeta(1, ab$a, ab$b)
  stats::rbinom(1, size = N, prob = theta)
}

# Vectorised over Q_star (and optional vector N, rho)
rbetabinom <- function(N, Q_star, rho) {
  # recycle N and rho if scalar
  N    <- rep(N,    length.out = length(Q_star))
  rho  <- rep(rho,  length.out = length(Q_star))
  out <- numeric(length(Q_star))
  for (i in seq_along(Q_star)) out[i] <- rbetabinom1(N[i], Q_star[i], rho[i])
  out
}
# after your loop, you had (example variable names):
# new_allele_freq_rr: T x B matrix = posterior freq of R allele at each (t, locus)
Q_star_mat <- new_allele_freq_ss  # dimension T x B, values in (0,1)


# You built: RR_per_locus <- t(Z_post) %*% RR_locus  # T x B
Q_star_mat <- RR_per_locus          # dimension T x B




set.seed(123)

T <- nrow(Q_star_mat)
B <- ncol(Q_star_mat)

# Choose sample sizes per (t, locus) — can be scalar, vector per-locus, or full matrix
N_mat   <- matrix(100, nrow = T, ncol = B)   # e.g., 100 mosquitoes per time & locus
rho_val <- 0.05                              # modest overdispersion (ICC)

# Simulate
D_mat <- matrix(NA_integer_, nrow = T, ncol = B)
for (t in 1:T) {
  D_mat[t, ] <- rbetabinom(N = N_mat[t, ], Q_star = Q_star_mat[t, ], rho = rho_val)
}

# D_mat[t, i] is your simulated count of "resistant" (or whatever success means)
# N_mat[t, i] is the corresponding total
# Quick sanity check: empirical mean vs Q*
colMeans(D_mat / N_mat)  # should roughly track colMeans(Q_star_mat)

