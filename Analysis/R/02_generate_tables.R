remove(list = ls())

ndigits <- 2

library(dplyr)
library(posterior)
library(tables)

proj.dir <- '/proj/ibrahimlab/strapp_survival'
samples.dir <- file.path(proj.dir, 'Analysis', 'Results')
save.dir <- file.path(proj.dir, 'Analysis')
file.list   <- list.files(samples.dir, pattern = '.rds')


varnames <- c('trt', 'age_c', 'sex', 'perform')
for ( i in 1:length(file.list) ) {
  file_i <- readRDS( file.path(samples.dir, paste0('analysis_', i, '.rds' ) ) )
  samples <- file_i[[1]]
  samples <- samples %>% 
    subset_draws(variable = varnames) %>%
    summarize_draws(
      function(x) c(mean = mean(x), sd = sd(x), quantile(x, probs = c(0.025, 0.975)))
    )
  names(samples)[4:5] <- c('lower95', 'upper95')
  scen_i  <- data.frame(c(id = i, dic = file_i$dic, file_i$scen))
  scen_i  <- data.frame(scen_i)
  scen_i  <- scen_i[rep(1, times = nrow(samples)), ]
  samples <- cbind(scen_i, samples)
  if ( i == 1 ) {
    res <- samples
  } else {
    res <- rbind(res, samples)
  }
}
rownames(res) <- NULL
## Construct credible interval as string
textci <- function(lower, upper, digits = 3, mathmode = TRUE) {
  lower.ch <- formatC(lower, digits = digits, format = 'f')
  lower.ch[lower > 0] <- paste0('~', lower.ch[lower > 0])
  upper.ch <- formatC(upper, digits = digits, format = 'f')
  upper.ch[upper > 0] <- paste0('~', upper.ch[upper > 0])
  ci <- paste0('(', lower.ch, ', ', upper.ch, ')')
  if (mathmode)
    ci <- paste0('$', ci, '$')
  ci
}
res <- res %>%
  mutate(ci = textci(lower95, upper95, digits = ndigits))



## Recode values for table
res$variable <- recode(res$variable, trt = 'Treatment', age_c = 'Age', sex = 'Sex', perform = 'Performance')
res$priors <- recode(
  res$priors, refprior = 'Reference', pp = 'PP', strapp = 'straPP', genstrapp = 'Gen-straPP'
)
res$models <- recode(
  res$models, pwe = 'PH', 'curepwe' = "Mixture cure", 'promotiontime' = "Promotion time"
)


## Change a0 to factor variable with equal number of decimals
res$a0 <- formatC(res$a0, digits = 2, format = 'f')
# res$a0 <- factor(res$a0, levels = c('0.00', '0.25', '0.75', '1.00'), ordered = TRUE)

## Change prior to factor so it appears in consistent order
res$priors <- factor(
  res$priors, levels = c('Reference', 'straPP', 'Gen-straPP', 'PP'), ordered = TRUE
)

## Change variable to ordered factor
res$variable <- factor(
  res$variable, levels = c('Treatment', 'Age', 'Sex', 'Performance'), ordered = TRUE
)

## Create character version of DIC
res$dic_ch <- formatC(res$dic, digits = ndigits, format = 'f')


## Select best model for each prior based on DIC
best.models <- res %>%
  group_by(priors, models, hist.type, variable) %>%
  mutate(min_dic = min(dic)) %>%
  filter(dic == min_dic) %>%
  select(id, dic, dic_ch, a0, priors, models, hist.type, variable, mean, ci)



## Use booktabs
booktabs()

## Create table for normal historical data
tbl.normal <- tabular(
  RowFactor(models, 'Model') * Factor(priors, 'Prior') ~
    Heading() * identity * (Heading('DIC') * dic_ch + Heading("$a_0$") * a0) + 
      Heading() * identity * Factor(variable, 'Variable') * (
        Heading("Mean") * Format(digits = ndigits) * round(mean, ndigits) + Heading("$95\\%$ CI") * ci
      )
  , data = best.models %>% filter(hist.type == 'normal')
)
filename = file.path(save.dir, 'normal_results.tex')
fileConn <- file(filename)
writeLines(as.character(toLatex(tbl.normal)), fileConn)
close(fileConn)

## Create table for logistic historical data
tbl.logistic <- tabular(
  Factor(priors, 'Priors') ~
    Heading() * identity * (Heading('DIC') * dic_ch + Heading("$a_0$") * a0) + 
    Heading() * identity * Factor(variable, 'Variable') * (
      Heading("Mean") * Format(digits = ndigits) * round(mean, ndigits) + Heading("$95\\%$ CI") * ci
    )
  , data = best.models %>% 
    filter(hist.type == 'logistic', models == 'Mixture cure')
)
filename = file.path(save.dir, 'logistic_results.tex')
fileConn <- file(filename)
writeLines(as.character(toLatex(tbl.logistic)), fileConn)
close(fileConn)







## Create table for normal historical data with all a0
res$a0_f <- formatC(res$a0, digits = 2, format = 'f')
tbl.normal.full <- tabular(
  RowFactor(models, 'Model') * Factor(priors, 'Prior') * Factor(a0_f, '$a_0$') ~
    Heading() * identity * (Heading('DIC') * dic_ch) + 
    Heading() * identity * Factor(variable, 'Variable') * (
      Heading("Mean") * Format(digits = ndigits) * round(mean, ndigits) + Heading("$95\\%$ CI") * ci
    )
  , data = res %>% filter(hist.type == 'normal')
)
## Drop missing a0s
tbl.normal.full <- tbl.normal.full[!is.na(tbl.normal.full[, 1]), ]
filename = file.path(save.dir, 'normal_results_full.tex')
fileConn <- file(filename)
writeLines(as.character(toLatex(tbl.normal.full)), fileConn)
close(fileConn)

## Create table for logistic historical data with all a0
tbl.logistic.full <- tabular(
  Factor(priors, 'Priors') * Factor(a0_f, '$a_0$') ~
    Heading() * identity * (Heading('DIC') * dic_ch) + 
    Heading() * identity * Factor(variable, 'Variable') * (
      Heading("Mean") * Format(digits = ndigits) * round(mean, ndigits) + Heading("$95\\%$ CI") * ci
    )
  , data = res %>% 
    filter(hist.type == 'logistic', models == 'Mixture cure')
)
tbl.logistic.full <- tbl.logistic.full[!is.na(tbl.logistic.full[, 1]), ]
filename = file.path(save.dir, 'logistic_results_full.tex')
fileConn <- file(filename)
writeLines(as.character(toLatex(tbl.logistic.full)), fileConn)
close(fileConn)







