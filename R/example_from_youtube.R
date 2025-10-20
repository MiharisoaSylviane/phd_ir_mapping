# https://www.sumsar.net/files/posts/2017-bayesian-tutorial-exercises/modeling_exercise1.html
# Number of random draws from the prior
n_draws <- 10000

prior <- runif(n_draws, min = 0, max = 1) # Here you sample n_draws draws from the prior  
hist(prior) # It's always good to eyeball the prior to make sure it looks ok.

#### 
# Here you define the ggenratinenerative model
generative_model <- function(rate) {
  subscribers <- rbinom(1, size = 16, prob = rate)
  subscribers
}

# Simulating the data
subscribers <- rep(NA, n_draws)
for(i in 1:n_draws) {
  subscribers[i] <- generative_model(prior[i])
}

# Filtering out those parameter values that didn't result in the
# data that we actually observed
post_rate <- prior[subscribers == 6]

# Checking that there enough samples left
length(post_rate)

hist(post_rate, xlim = c(0, 1))

mean(post_rate)


quantile(post_rate, c(0.025, 0.975))
## question 2 
sum(post_rate > 0.2) / length(post_rate)

# question 3
# This can be done with a for loop
singnups <- rep(NA, length(post_rate))
for(i in 1:length(post_rate)) {
  singnups[i] <- rbinom(n = 1, size = 100, prob = post_rate[i])
}

# But since rbinom is vectorized we can simply write it like this:
signups <- rbinom(n = length(post_rate), size = 100, prob = post_rate)

hist(signups, xlim = c(0, 100))


quantile(signups, c(0.025, 0.975))
