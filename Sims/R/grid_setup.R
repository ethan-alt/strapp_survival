## Set up simulation scenarios
remove(list = ls())

library(dplyr)
library(survival)
library(posterior)

## Directory to save simulation scenarios 
save.dir <- '/proj/ibrahimlab/strapp_survival/Sims/R'

n1        <- seq(100, 325, by = 25)

grid <- data.frame(
  id = seq_len(length(n1)),
  n1 = n1,
  stringsAsFactors = FALSE
)

ndatasets <- 5000  ## total number of data sets
each.cl   <- 20 ## how many data sets to run on a single node
ncl       <- ceiling(ndatasets / each.cl) ## how many nodes per data set

## Repeat each row of the grid ncl times (want nrow(grid) <= 800 ideally)
grid <- grid[rep(1:nrow(grid), each = ncl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1)
grid = grid[,-1]

set.seed(1)
grid$seed <- sample(seq_len(1000 * nrow(grid)), nrow(grid), replace = FALSE)
saveRDS(grid, file = file.path(save.dir, "grid.rds"))
