

## Input Onyen, e.g., ethanalt
onyen <- ''

## Create folders to save log, error, and rout files
o1         <- substr(onyen, 1, 1)
o2         <- substr(onyen, 2, 2)
scrdir     <- file.path('/pine/scr', o1, o2, onyen)
foldername <- 'strapp_survival/Sims'
dirs       <- c('Log', 'Error', 'Rout')

for ( x in dirs )
  dir.create(file.path(scrdir, foldername, x), recursive = TRUE)