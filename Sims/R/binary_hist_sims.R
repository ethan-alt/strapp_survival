## Simulation for (E1684, E1690) - binary historical and mixture cure rate for current
remove(list = ls())

library(dplyr)
library(survival)
library(posterior)
library(bayesplot)

# ## Set cmdstanr directory (necessary for sending code to cluster)
# library(cmdstanr)
# cmdstanr::set_cmdstan_path("/nas/longleaf/home/xinxinc/.cmdstan/cmdstan-2.30.1")

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_survival/R'
source(file.path(wrapper.dir, 'wrappers_logistic.R'))

## Sampling parameters
nburnin  = 2000
nsamples = 25000
nchains  = 1

## Obtain scenarios for computation: taking bootstrap samples from E1690
## with sample size (n1) varying from 100 to 325 
grid          <- readRDS(file = "/proj/ibrahimlab/strapp_survival/Sims/R/grid.rds")
each.cl       <- 20 ## how many data sets to run on a single node
priors        <- c('refprior', 'pp', 'strapp', 'genstrapp')
a0_vals       <- c(0, 1, 0.25, 0.25)
model.cur     <- 'curepwe'
histdata.type <- 'logistic'

## Get onyen of who is running job
grid$onyen <- ''
grid$onyen[1:625]     <- 'ethanalt'
grid$onyen[626:1250]  <- 'xinxinc'
grid$onyen[1251:1875] <- 'cococo'
grid$onyen[1875:2500] <- 'aniland'

## Get simulation situation based on cluster ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id) )
  id <- 1
grid.id <- grid[id, ]
n1.id   <- grid.id$n1

## Get save folder based on onyen
onyen    <- grid.id$onyen
o1       <- substr(onyen, 1, 1)
o2       <- substr(onyen, 2, 2)
save.dir <- file.path('/pine/scr', o1, o2, onyen, 'strapp_survival', 'Sims')

## Obtain file name based on id
filename <- file.path(save.dir, paste0('id_', id, '_', 'n1_', n1.id, '_rep_', grid.id$start, '-', grid.id$end, '.rds'))

## load data
data.dir   <- '/proj/ibrahimlab/strapp_survival/Data'
histdata  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
curdata   <- readRDS(file.path(data.dir, 'data_cur_1690.rds'))
## Replace 0 times with 0.50 (half day)
curdata <- curdata %>% mutate(failtime = if_else(failtime == 0, 0.50, failtime))

## Center and scale age
curdata$age_c <- as.numeric( scale(curdata$age) )
histdata$age_c <- as.numeric( scale(histdata$age) )

## formula for current and historical data (including intercept)
fmla.hist <- fail_2yr                ~ trt + age_c + sex + perform
fmla.cur  <- Surv(failtime, rfscens) ~ trt + age_c + sex + perform

## function to obtain post_mean, post_sd, CI for each prior
getStats <- function(fit, prior) {
  var_names <- names(fit$draws(format = 'draws_df'))
  var_names <- var_names[! var_names %in% c('lp__', '.chain', '.iteration', '.draw')]
  res = lapply(var_names, function(x){
    fit$summary(x, c(post_mean = mean, post_sd = sd, q=~quantile(.x, probs = c(0.025, 0.975))))
  })
  res       <- do.call(rbind, res)
  res$prior <- prior
  return(res)
}

## Simulation parameters
## By real data analysis results, 
## a0 = 0.25 leads to min DIC for straPP
## a0 = 0.25 leads to min DIC for Gen-straPP
## a0 = 1 leads to min DIC for power prior
 
## Set seed based on task ID
set.seed(grid.id$seed)
for(i in seq_len(each.cl)){
  idx = sample(1:nrow(curdata), n1.id, replace = TRUE)
  
  ## Take bootstrap samples from current data
  curdata.bootstrap = curdata[idx, ]
  nbreaks = 5
  probs   = 1:nbreaks / nbreaks
  breaks  = curdata.bootstrap %>%
    filter(rfscens == 1) %>%
    summarize(quant = quantile(failtime, probs = probs)) %>%
    unlist
  breaks = c(0, breaks)
  breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
  
  fit.list <- lapply(1:4, function(k){
    surv.sample(
      fmla.cur, fmla.hist, curdata.bootstrap, histdata
      , model.cur = model.cur
      , histdata.type = histdata.type, prior = priors[k]
      , breaks = breaks
      , iter_warmup = nburnin, iter_sampling = nsamples
      , chains = nchains
      , parallel_chains = 1
      , a0 = a0_vals[k]
      , f.tol = 1e-3, rel.tol = 1e-3, max.steps = 20
    )
  })
  
  simres.i <- lapply(1:4, function(k){
    getStats(fit.list[[k]], priors[k])
  })
  simres.i <- do.call(rbind, simres.i)
  
  if(i == 1) {
    simres <- simres.i
  } else{
    simres <- rbind(simres, simres.i)
  }
  print(paste0("#################### Completed iteration ", i, " ######################"))
}

colnames(simres)[c(4, 5)] <- c("lower", "upper")

res <- list(
  'id'       = id
  , 'n1'     = n1.id
  , 'simres' = simres
)

saveRDS(res, filename)

