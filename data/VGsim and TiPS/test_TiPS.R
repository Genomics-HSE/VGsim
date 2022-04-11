library(TiPS)
library(ape)

reactions <- c("S [beta*S*I] -> I",
               "I [gamma*I] -> R")

sir_simu <- build_simulator(reactions)

N <- 10000000

initialStates <- c(I = 1, S = N - 1, R = 0)

time_start <- 0
time_finish <- 20

time <- c(time_start, time_finish)

theta <- list(gamma = 1, beta = 2.5 / N)

safe_run <- function(f, ...) {
  out <- list()
  while(! length(out)) {out <- f(...)}
  out
}

safe_sir_simu <- function(...) safe_run(sir_simu, ...)

time1 <- Sys.time()


traj_dm <- safe_sir_simu(
    paramValues = theta,
    initialStates = initialStates,
    method = "exact",
    times = time)

time2 <- Sys.time()
print("Forward time exact")
print(time2-time1)



number_of_samples <- 10000

time1 <- Sys.time()

sir_tree <- simulate_tree(
  simuResults = traj_dm,
  #dates = dates,
  dates = seq(from = 12, to = 12, length.out = number_of_samples),
  deme = c("I"), # the type of individuals that contribute to the phylogeny
  sampled = c(I = 1), # the type of individuals that are sampled and their proportion of sampling
  root = "I", # type of individual at the root of the tree
  isFullTrajectory = FALSE, # deads do not generate leaves
  nTrials = 1,
  addInfos = FALSE) # additional info for each node

time2 <- Sys.time()
print("Backward time")
print(time2-time1)

#ape::plot.phylo(sir_tree, cex = .5)
