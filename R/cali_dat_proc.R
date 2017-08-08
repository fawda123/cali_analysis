library(readxl)
library(dplyr)
library(tidyr)
library(missMDA)
library(FactoMineR)

######
# taxonomic abundance data
# output is list for each taxa level, in long format with abund, site name, wwtp, and treatment (cont)

# add site meta (from ignore folder)
site_meta <- data.frame(
  site = paste0('CA', c(1, 2, 3, 4, 7, 8, 9, 11, 5, 10, 12, 6, 13, 14)),
  wwtp = c(rep('LAco', 4), rep('SDci', 4), rep('ORco', 3), rep('LAci', 3)), 
  cont = c('int', 'clo', 'int', 'far', 'clo', 'far', 'int', 'int', 'far', 'int', 'clo', 'clo', 'int', 'far'),  # int intermediate, clo closest to pip,  far farthest from pipe
  stringsAsFactors = F
)

# taxonomy data for the OTUs
taxdat <- read_xlsx('ignore/Run4_16S_R1.taxonomy.xlsx', sheet = 'PRIMER') %>% 
  rename(Domain = Kingdom)
names(taxdat) <- toupper(names(taxdat))

# OTU abundances by site
otuabu <- read_xlsx('ignore/Run4_16S_R1.27499.subsample.shared.xlsx', sheet = 'Run4_16S_R1.27499.subsample') %>% 
  select(-label, -numOtus) %>% 
  gather('otu', 'abund', -Group) %>% 
  rename(
    site = Group, 
    SPECIES = otu
    ) %>% 
  left_join(site_meta, by = 'site') %>% 
  left_join(taxdat, by = 'SPECIES')
  
# empty llist to fill
phylls <- c('DOMAIN', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES')
outs <- vector('list', length = length(phylls))
names(outs) <- phylls
  
# create separate list element for each taxanomic level, sites combined at each level
for(i in seq_along(phylls)){
  
  # taxa level to grab
  togrb <- phylls[i]
  
  # get relevant columns from otuabu, rename tax level for goruping
  comb <- otuabu[, c('site', 'wwtp', 'cont', 'abund', togrb)]
  names(comb)[names(comb) %in% togrb] <- 'grpvar'
  
  # group by tax level, sum abundances
  comb <- group_by(comb, site, wwtp, cont, grpvar) %>% 
    summarise(abund = sum(abund, na.rm = TRUE)) %>% 
    select(grpvar, everything()) %>% 
    arrange(grpvar)
  
  # get upper level taxa info for merging
  hitax <- taxdat[, phylls[1:i], drop = FALSE] %>% 
    unique
  names(hitax)[names(hitax) %in% togrb] <- 'grpvar'

  # merge upper level taxa info, format names
  comb <- left_join(comb, hitax, by = 'grpvar')
  names(comb)[names(comb) %in% 'grpvar'] <- togrb
  names(comb) <- tolower(names(comb))
  
  # append to output
  outs[[i]] <- comb
  
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
  

  

