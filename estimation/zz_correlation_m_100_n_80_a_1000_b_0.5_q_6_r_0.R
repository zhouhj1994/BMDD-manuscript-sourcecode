nsim <- 10

m <- 100
n <- 80
a <- 1000
b <- 0.5
q <- 6
r <- 0

###########################
type <- 'Correlation'

model <- 'gamma'
setup <- 'S6'
source('simu_fun.R')

model <- 'log-normal'
setup <- 'S7'
source('simu_fun.R')

model <- 'Poisson'
setup <- 'S8'
source('simu_fun.R')

model <- 'neg-binomial'
setup <- 'S9'
source('simu_fun.R')


