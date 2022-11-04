library(vtable)
library(gtsummary)
library(gt)
library(tidyverse)

## load data
data.dir   <- '/proj/ibrahimlab/strapp_survival/Data'
curdata  <- readRDS(file.path(data.dir, 'data_cur_1690.rds'))
e2696 <- readRDS( file.path(data.dir, 'data_hist_2696.rds') )
e1694 <- readRDS( file.path(data.dir, 'data_cur_1694.rds') )
e1684 <- readRDS( file.path(data.dir, 'data_hist_1684.rds') )
e1690 <- readRDS( file.path(data.dir, 'data_cur_1690.rds') )
e2696$trt <- e2696$ifn_gmk
e1694$trt <- e1694$gmk


e2696$Study <- 'E2696'
e1694$Study <- 'E1694'
e1684$Study <- 'E1684'
e1690$Study <- 'E1690'
pooled <- plyr::rbind.fill(e2696, e1694, e1684, e1690)
pooled$Study <- factor(pooled$Study)

xvars   <- c('perform', 'sex', 'age')
trtvar  <- 'trt'
pooled  <- pooled %>% select(c(xvars, trtvar, 'Study'))
pooled <- pooled %>%
  mutate(
    Performance = factor( ifelse(perform == 1, 'Ambulatory', 'Fully active'), levels = c('Fully active', 'Ambulatory') )
    , Sex = factor( ifelse(sex == 1, 'Female', 'Male'), levels = c('Male', 'Female') )
    , Treatment = factor( ifelse(trt == 0, 'Control', 'Treatment'), levels = c('Control', 'Treatment') )
    , Age = age
  )

# pooled %>%
#   sumtable(
#     vars = c('Treatment', 'Sex', 'Performance', 'Age')
#     , group = 'Study'
#     , digits = 1, fixed.digits = TRUE
#   )

gt <- pooled %>%
  select(c('Treatment', 'Sex', 'Performance', 'Age', 'Study')) %>%
  tbl_summary(
    by = 'Study'
    , digits = all_continuous() ~ 1
    ,  statistic = all_continuous() ~ "{mean} ({sd})"
  ) %>%
  modify_footnote(
    all_stat_cols() ~ "Mean (SD) or Frequency (%)"
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2", 'stat_3', 'stat_4') ~ "**Study**") %>%
  bold_labels()
gt <- as_gt(gt)
gt


gtsave(gt, '/proj/ibrahimlab/strapp_survival/Analysis/summarystats.tex')
