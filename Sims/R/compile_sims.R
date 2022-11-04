remove(list = ls())
library(dplyr)

sim.dir <- '/proj/ibrahimlab/strapp_survival/Sims/Results'
files   <- list.files(path = sim.dir, pattern = '.rds')
n1_vals <- sapply(files, function(x){as.numeric(strsplit(x, '_')[[1]][4])})
var_names = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]') # covariate names in simres

res_true  <- readRDS('/proj/ibrahimlab/strapp_survival/Analysis/Results/analysis_53.rds')
beta_true <- colMeans(res_true[[1]])[8:11]
names(beta_true) <- c("Treatment", "Age", "Gender", "Performance Status")
rm(res_true)

save.dir <- file.path('/proj/ibrahimlab/strapp_survival/Sims/R')
filename <- file.path(save.dir, 'combined_sims.rds')

combine_sims = function(n1){
  temp <- files[n1_vals == n1]
  for (i in seq_len(length(temp))) {
    simres.i <- readRDS(file.path(sim.dir, temp[i]))$simres
    simres.i <- simres.i %>% filter(variable %in% var_names)
    
    simres.i <- simres.i %>%
      rowwise() %>% 
      mutate(ci.ind  = ifelse(beta_true[which(var_names == variable)] >= lower &
                               beta_true[which(var_names == variable)] <= upper, 1, 0),
             log_var = ifelse(post_sd > 0, log(post_sd^2), NA), # replace post_sd = 0 by NA
             diff    = post_mean - beta_true[which(var_names == variable)],
             diff_sq = diff^2)
    
    if(i == 1) {
      simres <- simres.i
    } else{
      simres <- rbind(simres, simres.i)
    }
  }
  
  res <- simres %>% group_by(prior, variable) %>% summarize(
    avg_log_var = mean(log_var, na.rm = T)
    ,      bias = mean(diff)
    ,      mse = mean(diff_sq)
    ,      ci.cov = mean(ci.ind)
  )
  
  # change variable names
  res <- res %>%
    rowwise() %>% 
    mutate(variable = (names(beta_true))[which(var_names == variable)])
  res$n1 <- n1
  return(res)
}

res_df <- lapply(unique(n1_vals), combine_sims)
res_df <- do.call(rbind, res_df)
res_df <- as_tibble(res_df)
saveRDS(res_df, file = filename)



################################### Generate Plots ###############################
library(ggplot2)
library(ggpubr)
library(tidyr)
library(gridExtra)
library(grid)
library(dplyr)

res_df <- readRDS("/proj/ibrahimlab/strapp_survival/Sims/R/combined_sims.rds")
res_df$prior[res_df$prior == "genstrapp"] <- "Gen-straPP"
res_df$prior[res_df$prior == "strapp"]    <- "straPP"
res_df$prior[res_df$prior == "pp"]        <- "PP"
res_df$prior[res_df$prior == "refprior"]  <- "Reference"

#p3 = res_df %>%
#  ggplot(aes(x = n1, y = log(mse), linetype = prior, shape = prior, color = prior)) +
#  labs(y = "Log MSE", x = "(c)") +
#  geom_line() + 
#  geom_point() + 
#  facet_wrap(~variable)

p11 = res_df %>%
  filter(variable == "Age") %>%
  ggplot(aes(x = n1, y = avg_log_var, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(a)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p12 = res_df %>%
  filter(variable == "Gender") %>%
  ggplot(aes(x = n1, y = avg_log_var, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(b)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p13 = res_df %>%
  filter(variable == "Performance Status") %>%
  ggplot(aes(x = n1, y = avg_log_var, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(c)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p14 = res_df %>%
  filter(variable == "Treatment") %>%
  ggplot(aes(x = n1, y = avg_log_var, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(d)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p21 = res_df %>%
  filter(variable == "Age") %>%
  ggplot(aes(x = n1, y = bias, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(e)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p22 = res_df %>%
  filter(variable == "Gender") %>%
  ggplot(aes(x = n1, y = bias, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(f)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p23 = res_df %>%
  filter(variable == "Performance Status") %>%
  ggplot(aes(x = n1, y = bias, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(g)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p24 = res_df %>%
  filter(variable == "Treatment") %>%
  ggplot(aes(x = n1, y = bias, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(h)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p31 = res_df %>%
  filter(variable == "Age") %>%
  ggplot(aes(x = n1, y = log(mse), linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(i)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p32 = res_df %>%
  filter(variable == "Gender") %>%
  ggplot(aes(x = n1, y = log(mse), linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(j)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p33 = res_df %>%
  filter(variable == "Performance Status") %>%
  ggplot(aes(x = n1, y = log(mse), linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(k)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p34 = res_df %>%
  filter(variable == "Treatment") %>%
  ggplot(aes(x = n1, y = log(mse), linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(l)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p41 = res_df %>%
  filter(variable == "Age") %>%
  ggplot(aes(x = n1, y = ci.cov, linetype = prior, shape = prior, color = prior)) +
  ylim(0.8, 1) +
  labs(y = "", x = "(m)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p42 = res_df %>%
  filter(variable == "Gender") %>%
  ggplot(aes(x = n1, y = ci.cov, linetype = prior, shape = prior, color = prior)) +
  ylim(0.8, 1) +
  labs(y = "", x = "(n)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p43 = res_df %>%
  filter(variable == "Performance Status") %>%
  ggplot(aes(x = n1, y = ci.cov, linetype = prior, shape = prior, color = prior)) +
  ylim(0.8, 1) +
  labs(y = "", x = "(o)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

p44 = res_df %>%
  filter(variable == "Treatment") %>%
  ggplot(aes(x = n1, y = ci.cov, linetype = prior, shape = prior, color = prior)) +
  ylim(0.8, 1) +
  labs(y = "", x = "(p)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "none")

pl = list(p11, p21, p31, p41, p12, p22, p32, p42, p13, p23, p33, p43, p14, p24, p34, p44)

# Create row and column titles
col.titles = c("Age", "Gender", "Performance Status", "Treatment")
row.titles = c("Avg Log Var", "Bias", "Log MSE", "Coverage Prob")

# Add row titles
pl[1:4] = lapply(1:4, function(i) arrangeGrob(pl[[i]], 
                                              left=row.titles[i]
                                              #left = text_grob(row.titles[i], size = 10, face = "bold")
                                              )
                 )

# function to extract legend from plot
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# p11 with legend
p11_legend <- res_df %>%
  filter(variable == "Age") %>%
  ggplot(aes(x = n1, y = avg_log_var, linetype = prior, shape = prior, color = prior)) +
  labs(y = "", x = "(a)") +
  geom_line() + 
  geom_point() +
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.key.width= unit(2, 'cm'))
  
# extract legend from plot1 using above function
legend <- get_legend(p11_legend) 

grid.arrange(grobs = lapply(c(1, 5, 9, 13), function(i) {
  arrangeGrob(grobs=pl[i:(i+3)], 
              top=col.titles[i/4 + 1],
              #top=text_grob(col.titles[i/4 + 1], size = 10, face = "bold"),
              ncol=1)
}), bottom = legend, ncol=4)
