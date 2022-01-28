############### instance 2 paper 2016
# - blocking factor
# - same n° of levels
# - etas scalar
# - I, D, A
# - quadratic

library(multiDoE)

facts <- list(c(), 1:3)
units <- list(7, 4)
levels <- 3
etas <- list(1)
criteria <- c("I", "D")
model <- "quadratic"

msopt <- MSOpt(facts, units, levels, etas, criteria, model)









############### instance 4 paper 2016 ####
# - 3 strata
# - same n° of levels
# - etas vector
# - Id, Ds, As
# - interaction

library(multiDoE)

facts <- list(1, 2, 3)
units <- list(3, 3, 3)
levels <- 3
etas <- list(1, 1)
model2 <- "interaction"

############### main + livelli diversi ####

