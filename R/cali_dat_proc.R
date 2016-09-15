library(readxl)
library(dplyr)
library(tidyr)
library(missMDA)
library(FactoMineR)

######
# taxonomic abundance data
# output is list for each taxa level, in long format with abund, site name, wwtp, and treatment (cont)

# municipalities that are names of each spread sheet
munis <- c('ORco', 'SDci', 'LAco', 'LAci')

# import all taxonomic data for each site from each sheet
sht_ls <- list()
for(muni in munis){
  fl <- list.files('ignore', muni, full.names = T)
  shts <- excel_sheets(fl) %>% 
    grep('^CA', ., ignore.case = F, value = T)
  for(sht in shts){
    tmp <- read_excel(fl, sheet = sht)
    sht_ls[[sht]] <- tmp
  }
}

# add site meta (from ignore folder)
site_meta <- data.frame(
  site = paste0('CA', c(1, 2, 3, 4, 7, 8, 9, 11, 5, 10, 12, 6, 13, 14)),
  wwtp = c(rep('LAco', 4), rep('SDci', 4), rep('ORco', 3), rep('LAci', 3)), 
  cont = c('int', 'clo', 'int', 'far', 'clo', 'far', 'int', 'int', 'far', 'int', 'clo', 'clo', 'int', 'far'),  # int intermediate, clo closest to pip,  far farthest from pipe
  stringsAsFactors = F
)

# empty llist to fill
phylls <- c('DOMAIN', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES')
outs <- vector('list', length = length(phylls))
names(outs) <- phylls
  
# create separate list element for each taxanomic level, sites combined at each level
for(phyll in seq_along(phylls)){
  
  out_ls <- list()

  for(sht in names(sht_ls)){

    toget <- phylls[phyll]
    
    # find start row for phyll
    x <- sht_ls[[sht]]
    x <- as.data.frame(x, stringsAsFactors = F)
    strt <- grep(paste0('^', toget), x[, 1])
    
    # get NA values, min diff from strt is end for phyll
    if(phyll == length(phylls)){
      stp <- nrow(x)
    } else {
      stp <- grep(paste0('^', phylls[phyll + 1]), x[, 1]) - 1
    }

    # getcol
    getcol <- grep('^DOMAIN', x[, 1])
    getcol <- grep('abundance|count', x[getcol, ], ignore.case = T)[1]

    # get phyll data, format
    dat <- x[seq(strt, stp), ]
    dat <- dat[, c(1, getcol)]
    nms <- dat[1, ]
    dat <- na.omit(dat[-1,])
    names(dat) <- nms
    row.names(dat) <- 1:nrow(dat)
    names(dat)[1] <- tolower(names(dat)[1])
    names(dat)[2] <- 'abund'
    dat[, 1] <- gsub('[[:space:]]*$', '', dat[, 1])
    
    out_ls[[sht]] <- dat
    
  } 

  # combine sites for the phyll, add metadata
  comb <- reshape2::melt(out_ls) %>% 
    rename(site = L1) %>% 
    mutate(abund = as.numeric(abund)) %>% 
    left_join(., site_meta, by = 'site')
  
  # sometimes there are multiple abundance measures for the same taxa level
  # these need to be combined
  names(comb)[1] <- 'tmp'
  comb <- group_by(comb, tmp, site, wwtp, cont) %>% 
    summarize(abund = sum(abund)) %>% 
    ungroup %>% 
    data.frame
  names(comb)[1] <- tolower(phylls[phyll])

  # append to output
  outs[[phyll]] <- comb
  
}

abudat <- outs
save(abudat, file = 'data/abudat.RData')

##
# environmental data
# parameters with more than half of the values NA were removed, sites with more than half of the value NA were also removed

envdat <- read_excel('ignore/master environmental data.xlsx', sheet = 'master environmental') %>% 
  data.frame %>% 
  mutate(
    parameters = gsub('mg|kg|ug|/', '', parameters),
    parameters = gsub('[[:space:]]+$|[[:space:]]%$', '', parameters)
  )

rowkp <- apply(envdat[, -1], 1, function(x) sum(is.na(x)) < length(x)/2)
envdat <- envdat[rowkp, ]
colkp <- apply(envdat[, -1], 2, function(x) sum(is.na(x)) < length(x)/2)
envdat <- envdat[, c(T, colkp)] %>% 
  gather('site', 'val', -parameters) %>% 
  spread(parameters, val) 

# impute missing values
noaxes <- estim_ncpPCA(envdat[, -1])
tmp <- imputePCA(envdat[, -1], ncp = noaxes$ncp, scale = TRUE)
tmp <- as.data.frame(tmp$completeObs)
envdat <- cbind(envdat$site, tmp)
names(envdat)[1] <- 'site'

# combine with cats, split
envcat <- read.csv('ignore/envcat.csv', header = T, stringsAsFactors = F) 
envdat <- gather(envdat, 'parameters', 'val', -site) %>% 
  mutate(parameters = as.character(parameters)) %>% 
  left_join(., envcat, by = 'parameters')

# add site meta (from ignore folder)
site_meta <- data.frame(
  site = paste0('CA', c(1, 2, 3, 4, 7, 8, 9, 11, 5, 10, 12, 6, 13, 14)),
  wwtp = c(rep('LAco', 4), rep('SDci', 4), rep('ORco', 3), rep('LAci', 3)), 
  cont = c('int', 'clo', 'int', 'far', 'clo', 'far', 'int', 'int', 'far', 'int', 'clo', 'clo', 'int', 'far'),  # int intermediate, clo closest to pip,  far farthest from pipe
  stringsAsFactors = F
)
envdat <- right_join(site_meta, envdat, by = 'site')

envdat <- split(envdat, envdat$cat) %>% 
  lapply(., function(x) {
    out <- spread(x, parameters, val) %>% 
      mutate(site = as.character(site)) %>% 
      arrange(site) %>% 
      data.frame %>% 
      select(-cat)
    out
  })

save(envdat, file = 'data/envdat.RData')

##
# principal components analysis of envdat
# separate pca for envdat categories - contaminants, environmental, metals
# axes that explain at least 99% of the variance in each category are retained

data(envdat)

pcx <- lapply(envdat, function(x){

  tomod <- x[, -c(1:3)]
  row.names(tomod) <- x$site
  PCA(tomod, graph = FALSE)
  
})

# get scores for each pc analysis
scrs <- lapply(pcx, function(x){
  
  # axes to get
  toget <- which(x$eig[, 3] < 99)
  if(length(toget) == 0) toget <- 1
  out <- x$ind$coord[, toget, drop = F]
  out <- data.frame(site = envdat[[1]]$site, wwtp = envdat[[1]]$wwtp, cont = envdat[[1]]$cont, out)
  return(out) 

})
scrs <- reshape2::melt(scrs) %>% 
  mutate(
    variable = gsub('Dim\\.', '', variable),
    axis = paste(L1, variable, sep = '')
  ) %>% 
  select(axis, value, site, wwtp, cont) %>% 
  spread(axis, value)

envdatpca <- scrs
save(envdatpca, file = 'data/envdatpca.RData')

##
# envdatpca by municipality
rm(list = ls())

data(envdat)
muni <- unique(envdat[[1]]$wwtp)
categ <- names(envdat)
grd <- expand.grid(muni, categ)
names(grd) <- c('muni', 'categ')

txt <- 2
arrow <- 0.1

pcx <- vector('list', length = nrow(grd))
names(pcx) <- paste(grd[, 1], grd[, 2], sep = '_')

for(i in 1:nrow(grd)){
    
  muni <- grd[i, 'muni']
  categ <- grd[i, 'categ']
  
  tomod <- envdat[[categ]]
  tomod <- tomod[tomod$wwtp %in% muni, ]
  row.names(tomod) <- tomod$site
  mod <- PCA(tomod[, -c(1:3)], graph = FALSE)
  nms <- paste(muni, categ, sep = '_')
  
  pcx[[nms]] <- mod
  
}

# get scores for each pc analysis
scrs <- lapply(pcx, function(x){

  # axes to get
  toget <- which(x$eig[, 3] < 99)
  if(length(toget) == 0) toget <- 1
  out <- x$ind$coord[, toget, drop = F]
  out <- data.frame(out)
  out$site <- row.names(out)
  row.names(out) <- 1:nrow(out)
  return(out) 

})
scrs <- reshape2::melt(scrs) %>% 
  separate(L1, c('wwtp', 'cat'), sep = '_') %>% 
  unite(var, cat, variable, sep = '') %>% 
  mutate(var = gsub('Dim\\.', '', var)) %>% 
  rename(variable = var) %>% 
  spread(variable, value)

scrs <- split(scrs, scrs$wwtp) %>% 
  lapply(., function(x){
    
    torm <- apply(x, 2, function(v) sum(is.na(v)) == length(v))
    x[, !torm]
    
  })

envdatpca_muni <- scrs
save(envdatpca_muni, file = 'data/envdatpca_muni.RData')

##
# abudat with full taxonomic hierarchy

# # these links are broken 
# taxurl <- list(
#   CA2 = 'http://nbc.ece.drexel.edu/results/f22c0edd3f7545d426c17fb098ab3ad6',
#   CA3 = 'http://nbc.ece.drexel.edu/results/1f19ff149d601ac45d2ad391426b7eb1',
#   CA4 = 'http://nbc.ece.drexel.edu/results/1b7248a9c4a8f825ba71c34ce2794f53',
#   CA1 = 'http://nbc.ece.drexel.edu/results/39b5a9b78aefc7117e5f56ee6215a503',
#   CA7 = 'http://nbc.ece.drexel.edu/results/58f7dc39f300992900f5ac5d617921ec',
#   CA6 = 'http://nbc.ece.drexel.edu/results/328417aa0d49ab40a3f74956f336cb15',
#   CA5 = 'http://nbc.ece.drexel.edu/results/b008a5918798fcad483c7c21a30d67b8',
#   CA8 = 'http://nbc.ece.drexel.edu/results/06d109352045840929b6d45c4b516f28',
#   CA9 = 'http://nbc.ece.drexel.edu/results/7b6abf151b4d0357be5c3c6525711166',
#   CA10 = 'http://nbc.ece.drexel.edu/results/d54ea6dc293aa71259e1082a2f5f2d20',
#   CA11 = 'http://nbc.ece.drexel.edu/results/da303a4c95664d0d088f27998ac82184',
#   CA12 = 'http://nbc.ece.drexel.edu/results/95e6cb9e641f32a1a0618f1533ec805f',
#   CA14 = 'http://nbc.ece.drexel.edu/results/d02817c08ea4acad2d57d372fd2ee1b1'
# )
# 
# # have to read lines first, save as csv
# # direct read with read.table was not working
# for(site in names(taxurl)){
#   url <- paste0(taxurl[[site]], '/summarized_results.txt')
#   dat <- readLines(url)
#   dat <- gsub('\t', ',', dat)
#   flnm <- paste0('ignore/', site, '_tax.csv') 
#   write.csv(dat, flnm, quote = F, row.names = F)
# }

# reimport csv files
fls <- list.files('ignore', '^CA.*tax\\.csv', full.names = T)
tax <- vector('list', length = length(fls))
sites <- gsub('ignore/|_tax\\.csv', '', fls)
names(tax) <- fls
for(fl in fls) tax[[fl]] <- read.csv(fl, header = T, skip = 1, stringsAsFactors = FALSE, na.strings = '')

# combine by site, get uniques
tax <- do.call('rbind', tax) %>% 
  select(-QUERY_ID, -BEST_MATCH, -SCORE) %>% 
  unique
row.names(tax) <- 1:nrow(tax)
names(tax) <- tolower(names(tax))

# go through abudat, left merge with taxa given taxonomic level
data(abudat)
abudat <- lapply(abudat, function(x){
  
  # subset tax by level in abudat (unless domain)
  mrglev <- names(x)[1]
  tomrg <- which(names(tax) == mrglev)
  if(tomrg == ncol(tax)) return(x)
  tomrg <- tax[, tomrg:ncol(tax)] %>% 
    unique
  
  # join
  out <- left_join(x, tomrg, by = mrglev)
  
  return(out)
  
})

save(abudat, file = 'data/abudat.RData')

######
# get sites lat/long and euclid dist

# add site meta (from ignore folder)
site_meta <- data.frame(
  site = paste0('CA', c(1, 2, 3, 4, 7, 8, 9, 11, 5, 10, 12, 6, 13, 14)),
  wwtp = c(rep('LAco', 4), rep('SDci', 4), rep('ORco', 3), rep('LAci', 3)), 
  cont = c('int', 'clo', 'int', 'far', 'clo', 'far', 'int', 'int', 'far', 'int', 'clo', 'clo', 'int', 'far'),  # int intermediate, clo closest to pip,  far farthest from pipe
  stringsAsFactors = F
)

# locs
dists <- read_excel('ignore/Latitude longitude data 2Feb2016.xlsx') %>% 
  select(-wwtp) %>% 
  left_join(site_meta, ., by = 'site') %>% 
  group_by(wwtp) %>% 
  arrange(cont) %>% 
  nest() %>% 
  mutate(dists = 
      purrr::map(data, function(x){
        dists <- dist(x[, c('latitude', 'longitude')])
        c(0, dists[1:(nrow(x) - 1)])
        })
  ) %>% 
  unnest %>% 
  ungroup %>% 
  arrange(site)

save(dists, file = 'data/dists.RData')
  

  

