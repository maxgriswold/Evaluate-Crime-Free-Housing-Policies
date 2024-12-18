# Calculate SFD estimators
# Max Griswold
# 11/14/21

library(tidygeocoder)
library(data.table)
library(plyr)
library(dplyr)
library(sf)
library(conleyreg)

sf::sf_use_s2(FALSE)

library(tidycensus)
options(tigris_use_cache = TRUE)

census_key <- ""

analysis_type <- "main_analysis"

main_dir  <- "./data"
model_dir <- paste0("./cfho_analysis/sfd/", analysis_type)

prep_cfho_locs <- F
save_tables    <- T
  
setwd(main_dir)

# Grab location hierarchy, tiger shapefiles, and census variables by block group
loc_hierarchy <- fread("helpful/census_cal_loc_hierarchy.csv")

census_vars <- fread('./census/acs5/california_block_group_2019_processed.csv')
setnames(census_vars, 'geoid', 'block_geoid')
census_vars[, block_geoid := as.numeric(block_geoid)]

tiger_block <- get_acs(geography = "block group", year = 2019, variables = "B01003_001", 
                       state = "CA", geometry = T, key = census_key, survey = "acs5") %>%
                        setDT %>%
                        .[, .(GEOID, NAME, geometry)] %>%
                        setnames(., c("NAME", "GEOID"), c("block_name", "block_geoid")) %>%
                        .[, block_geoid := as.numeric(block_geoid)] %>%
                        st_as_sf()

tiger_place <- get_acs(geography = "place", year = 2019, variables = "B01003_001", 
                       state = "CA", geometry = T, key = census_key, survey = "acs5") %>%
                        setDT %>%
                        .[, .(GEOID, NAME, geometry)] %>%
                        setnames(., c("NAME", "GEOID"), c("place_name", "place_geoid")) %>%
                        .[, place_geoid := as.numeric(place_geoid)] %>%
                        st_as_sf()

evict <- fread("./data/evict_block.csv")
evict <- evict[year == 2019,]
evict[, fips_code := as.numeric(fips_code)]
setnames(evict, "fips_code", "block_geoid")

evict <- evict[, .(location, evict_count, evict_rate)]

#Somehow, block group geoids in the eviction dataset got scrambled for 
# ~1k blocks. So, I'll merge back on geoid using the tiger files:
tiger_helper <- as.data.table(tiger_block)
tiger_helper <- tiger_helper[, .(block_geoid, block_name)]
setnames(tiger_helper, "block_name", "location")

evict <- join(evict, tiger_helper, type = 'left', by = "location")
evict <- join(evict, census_vars, type = 'left', by = "block_geoid")

#Sets up rolling batches to be processed by the geocoder, given input limits to
#the census API.

batch_geocodes <- function(d, batches = 2e3, add = ""){
  
  #Create splitting column
  d  <- d[, batch := ceiling(.I/batches)]
  ds <- split(d, by = "batch", keep.by = F)
  
  processed <- llply(ds, code_coordinates, added_text = add)
  processed <- rbindlist(processed)
  
  return(processed)
  
}

# Iterative function which applies geocoders to addresses. For each address,
# try a geocoder. If it doesn't work, try the next one.

code_coordinates <- function(df, added_text = ""){
  
  #Geocoders are super temperamental. Given this, we'll need to process geocodes
  #in batches. We'll split datasets into sets of 2k. For each batch, we'll try
  #out a cascade of open-source geocoders, starting w/ Census, then arc-gis, 
  #then finally open-street-maps. If all else fails, we can try Google at a later point (I have
  #an API key but it's set up to bill a different PTN)
  
  #Lots of unit tests in this function; I was running into  issues
  #from certain geocoders not producing new results. There's likely a better
  #way to do the below; too much rewritten code!
  
  dd <- df %>%
    geocode(address, method = "census") %>%
    setDT() %>%
    .[, method := "census"]
  
  #Add on additional text to address to create additional redundancies, which
  #tends to improve reliability of geocodes
  dd[, address := paste(address, added_text)]
  
  res <- dd[!is.na(lat)]
  dd  <- dd[is.na(lat),]
  
  #Remove lat/long column for uncoded addresses since keeping this column
  #in the dataset is creating issues for arcgis method:
  dd <- dd[, `:=`(lat = NULL, long = NULL)]
  
  if (dim(dd)[1] >= 1){
    dd <- dd %>%
      geocode(address, method = "arcgis") %>%
      setDT() %>%
      .[, method := "arcgis"]
    
    if (dim(dd[!is.na(lat)])[1] >= 1){
      res <-setDT(rbind(res, dd[!is.na(lat)]))
    }
    dd  <- dd[is.na(lat), ]
    dd <- dd[, `:=`(lat = NULL, long = NULL)]
  }
  
  if (dim(dd)[1] >= 1){
    dd <- dd %>%
      geocode(address, method = "osm") %>%
      setDT() %>%
      .[, method := "osm"]
    
    if (dim(dd[!is.na(lat)])[1] >= 1){
      res <- setDT(rbind(res, dd[!is.na(lat)]))
    }
  }
  
  print(paste0("Could not code ", dim(dd)[1], " addresses."))
  print(res)
  return(res)
  
}


# For a given spatial point, merge with census block shapefile

geo_in_block <- function(d, geo){
  
  d <- st_as_sf(d, coords = c("long", "lat"))
  d <- st_set_crs(d, st_crs(tiger_block))
  
  d <- st_join(d, geo, join = st_intersects, largest = T)
  setDT(d)
  
  d[, block_geoid := as.numeric(block_geoid)]
  d[, .(address, block_geoid, geometry)]
  
  return(d)
  
}

# Using a place shapefile, determine all blocks that touch that polygon
# Then, for each block, calculate the percentage of area within a census 
# place

block_in_place <- function(place, blocks){
  
  overlap               <- st_intersection(place, blocks)
  overlap$area_overlap  <- as.numeric(st_area(overlap))
  
  # Intersections can occur at the boundary between two the city & neighboring
  # city blocks. Drop these
  overlap <- overlap[overlap$area_overlap > 0,]
  
  setDT(overlap)
  overlap <- overlap[, .(block_geoid, area_overlap)]
  
  # Get original area of blocks and merge with intersected blocks
  old_blocks <- tiger_block[tiger_block$block_geoid %in% overlap$block_geoid,]
  old_blocks$area_original  <- as.numeric(st_area(old_blocks))
  
  setDT(old_blocks)
  old_blocks <- old_blocks[, .(block_geoid, area_original, geometry)]
  
  overlap <- join(overlap, old_blocks, by = "block_geoid", type = "left")
  overlap[, overlap_percent := area_overlap/area_original, ]
  overlap[, overlap_80 := ifelse(overlap_percent >= 0.8, 1, 0)]
  overlap[, overlap_50 := ifelse(overlap_percent >= 0.5, 1, 0)]
  
  return(overlap)
  
}

construct_analysis_df <- function(city_name){
  
  # Load up CFHO locations, geocode addresses, then merge with census blocks
  city_cfho <- fread(sprintf("cfho_locs/%s_2019.csv", city_name))
  
  if (city_name == "san_diego_unincorporated"){
    city_cfho[, address := paste0(address, sprintf(", %s, CA", city))]
  }else{
    city_cfho[, address := paste0(address, sprintf(", %s, CA", city_name))]
  }
  
  city_cfho <- batch_geocodes(city_cfho)
  city_cfho <- geo_in_block(city_cfho, tiger_block)
  
  city_cfho[, cfho_num := .N, by = "block_geoid"]
  city_cfho <- unique(city_cfho[, .(block_geoid, cfho_num)])
  
  if (city_name == "san_diego_unincorporated"){
    
    city_tiger <- get_acs(geography = "block group", year = 2019, variables = "B01003_001", 
                          state = "CA", geometry = T, county = "san diego", key = census_key, survey = "acs5") %>%
      setDT %>%
      .[, .(GEOID, NAME, geometry)] %>%
      setnames(., c("NAME", "GEOID"), c("block_name", "block_geoid")) %>%
      .[, block_geoid := as.numeric(block_geoid)]
    
  }else{
    
    tiger_place$place_name <- tolower(tiger_place$place_name)
    city <- tiger_place[tiger_place$place_name %like% paste(city_name, "city"),]
    city_tiger <- block_in_place(city, tiger_block)
    
  }
  
  # Calculate proportion of census block within city, then merge CFHO locations,
  # eviction data, and census data on canonical map
  df_city <- join(city_tiger, city_cfho, by = c("block_geoid"), type = "left")
  df_city <- join(df_city, evict, by = "block_geoid", type = "left")
  
  df_city[is.na(cfho_num), cfho_num := 0]
  df_city[, cfho_any := ifelse(cfho_num > 0, 1, 0)]
  

  df_city[is.na(evict_count), evict_count := 0]
  df_city[is.na(evict_rate)|evict_rate == Inf, evict_rate := 0]
  
  if (analysis_type == "reviewer_analysis"){
    df_city[, evict_count := asinh(evict_count)]
  }
  
  df_city <- df_city[, .(block_geoid, evict_count, cfho_num, cfho_any, number_rental_units,
                         per_capita_income, pop_white, pop_black, pop_asian, 
                         pop_native_american, pop_latin_hispanic,
                         pop_below_2_fpl, pop_male, number_rental_units, geometry)]
  
  # Change denominator for per-capita income and number of renting households so
  # that parameter values are more interpretable in table:
  
  df_city[, number_rental_units_100 := number_rental_units/100]
  df_city[, per_capita_income_10k := per_capita_income/10000]
  
  # Do not incorporate blocks without rental units as a comparators (since evictions
  # are 0 except commercial evictions)
  df_city <- df_city[df_city$number_rental_units_100 > 0,]
  
  # Hold onto minimum rows and columns for SFD models. SFD should only include
  # complete cases:
  df_city <- na.omit(df_city, cols = names(df_city))

  df_city <- st_as_sf(df_city) %>%
              st_transform(4326)
  
  sfd_save <- paste0("./cfho_analysis/sfd/", city_name, "_prepped_sfd.geojson")
  st_write(df_city, sfd_save, append = F, delete_dsn = T, delete_layer = T)
  
  return(df_city)
  
}

# Original algorithm: Hannah Druckenmiller, March 11 2018
# Some updates by Vincent Tanutama, April 20 2018
# Additional updates and adapted by Max Griswold Oct 1 2022

construct_sfd <- function (spatial_df, tol = 0, dependent_var, independent_vars, 
                           rotation = 0, needs_balanced = FALSE, balanced_df = NA, 
                           plot = TRUE, pngfilename) {
  
  library(magrittr, quietly = "true")
  library(dplyr, quietly = "true")
  library(sp, quietly = "true")
  library(sf, quietly = "true")
  library(geosphere, quietly = "true")
  
  # Rotation function
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  
  # Results holder
  robSDF <- data.frame()
  SP <- spatial_df # SP is now the original spatial_df
  SP$rowno <- 1:nrow(SP)
  
  # Here I defined rotation to take a numeric vector 
  # in degrees, e.g. seq(15, 180, 15)
  for (dg in rotation) {
    
    # Convert sp to sf object for more stablity
    spatial_df <- st_as_sf(SP) 
    st_crs(spatial_df) <- 4326
    # Identify row index of polygon i's neighbors
    neigh <- st_touches(spatial_df, spatial_df, sparse = T)
    # Produce centroid of polygon i for bearing calculation
    cent  <- lapply(st_centroid(st_geometry(spatial_df)), as.vector)
    
    # Rotate shapefile
    if (dg!=0) {
      spdf.geom <- st_geometry(spatial_df) 
      # need to rotate object that contains no crs information because
      # rotation can produce points outside the defined extent of the globe, 
      # i.e. bigger than [-180,180]x[-90,90]
      spdf2 <- spdf.geom * rot(dg*pi/180)
      spatial_df <- st_sf(cbind(as.data.frame(spatial_df),spdf2))
    }     
    
    # Compare centroid-based relative position of each polygon and its neighbors
    pone <- Map(rep, as.list(1:nrow(spatial_df)), lapply(neigh,length))
    p1 <- lapply(pone, as.list)
    p2 <- lapply(neigh, as.list)
    p1 <- unlist(unlist(rapply(p1, function(x) cent[x], how = 'replace'), recursive = F), recursive = F)
    p2 <- unlist(unlist(rapply(p2, function(x) cent[x], how = 'replace'), recursive = F), recursive = F)
    bearing12 <- unlist(Map(bearing, p2, p1)) + 180 - 22.5
    e <- which(sapply(pone, length)==0)
    if  (length(e)!=0) { 
      pone[e]  <- as.list(e)
      for (k in 1:length(e)) {
        if (e[k]==1) {
          bearing12 <- c(NA, bearing12)
        } else {
          bearing12 <- c(bearing12[1:length(unlist(pone[1:(e[k]-1)]))], NA,
                         bearing12[length(unlist(pone[1:(e[k])])):length(bearing12)])
        }
      }
    } 
    
    # Start sampling neighbors starting at leftmost polygon,
    # going in northeast direction until last polygon no longer has neighbors to sample
    # Algorithm is greedy in sampling multiple neighbors at one time
    left    <- st_bbox(spatial_df)[1]
    bb      <- lapply(st_geometry(spatial_df), st_bbox)
    which(sapply(sapply(bb, function(x) which(x == left)), length) != 0) -> polyleft

    # put bearing vector back to a list object per polygon i
    # so each element is the bearings of neighbors of polygon i
    bears <- split(bearing12, unlist(pone)) 
    universe <- 1:nrow(spatial_df)
    # panel is the list whose element is a vector of ordered sampled neighbors
    panel <- list()
    # chopped_map is the shapefile after neighbors have been sampled
    chopped_map <- spatial_df
    chopped_map$no <- 1:nrow(chopped_map)
    
    # Can set tolerance to non-zero, but that will leave some counties unsampled :(
    while (nrow(chopped_map)>(nrow(spatial_df)*tol)) {
      
      p <- polyleft
      sample <- polyleft
      m <- 1
      
      while(length(polyleft)!=0) {
        
        # find p's next neighbor (.r = row number, .d = bearing)
        nx.r <- neigh[[p]] 
        nx.d <- bears[[p]]  
        
        # sort bearings of neighbors in the clockwise order starting from north
        sortma  <- sort(nx.d)
        # when there are multiple neighbors ...>
        if (length(sortma)!=1) {
          #>... grab all neighbors -- this is what I meant by ``greedy'' --
          #     in the counter-clockwise direction
          #     so that at the end of sampling, we're back to the upperleftmost polygon
          #     but want to stop greedy neighbor extraction if sequence of neighbor isn't contiguous
          #     so need to check for 2nd degree neighbors too
          m <- rev(match(sortma, nx.d))
          nm <- lapply(m, function(x) neigh[[nx.r[x]]])     # 2nd degree neighbors of each neighbor of i
          checknm <- Map(match, as.list(lead(nx.r[m])), nm) # check for 2nd degree neighbors in the set of 2nd degree neighbors
          checknm <- unlist(checknm)[-length(checknm)]      # NA means in the set of neighbors of i = {A, B, C, D}, B isn't neighbor of A
          wnm <- which(is.na(checknm))
          if (length(wnm)!=0) {
            
            if (max(wnm) %in% c(1,length(checknm))) {
              m <- match(sortma, nx.d)[[1]]
            } else {
              m0 <- max(wnm)
              m  <- m[(m0+1):length(checknm)]
            }
            
          } 
          
          #>... otherwise just grab that neighbor
        } else {
          m <- which(nx.d == sortma)
        }
        
        # remove sampled polygons from neigh, bears, universe 
        # add sampled polygons to object sample
        if(length(m)!=0) {
          
          nxneigh <- nx.r[m]
          wn <- which(neigh[[p]]%in%nxneigh)
          universe[c(p, nxneigh)] <- NA
          ns <- lapply(neigh, function(x) which(x %in% c(p,nxneigh)))
          wns <- which(sapply(ns, length)!=0)
          neigh[wns] <- Map(function(x, y) replace(x, y, NA),
                            neigh[wns], ns[wns])
          bears[wns]  <- Map(function(x, y) replace(x, y, NA),
                             bears[wns], ns[wns])
          sample <- c(sample, nxneigh)
          
        } else {
          
          # when there is no more polygons to be sampled, break the while loop
          # and end this current sequence of neighbor sampling ...(*)
          polyleft <- integer(0)
          nxneigh  <- integer(0)
        }
        
        # define polygon for next iteration
        # nxneigh is all neighbors in current iteration
        # p is element of nxneigh for next iteration
        if (length(nxneigh)>1) {
          p <- nxneigh[length(nxneigh)]
        } else if (length(nxneigh)==1) {
          p <- nxneigh
        } else {
          p <- integer(0)
        }
        
      }
      
      # (*)... redefine next sequence of neighbors
      chopped_map <- chopped_map[-which(chopped_map$no %in% sample),]
      left    <- st_bbox(chopped_map)[1]
      which(sapply(sapply(bb, function(x) which(x == left)), length) != 0) -> polyleft
      if (length(polyleft)>1) polyleft <- polyleft[which(!polyleft%in%unlist(panel))][[1]]
      panel <- append(panel, list(spatial_df$rowno[sample]))
      
    }
    
    # Convert list of sequences of neighbors into dataframe
    # that contains sequence number and ordered polygons within the sequence
    panelgroup <- Map(rep, as.list(1:length(panel)), lapply(panel, length))
    paneldf <- data.frame(group = unlist(panelgroup),
                          rowno = unlist(panel),
                          stringsAsFactors = F)
    paneldf$angle <- dg
    # since ordering is paramount, need to have index of ordered rows
    paneldf$index <- 1:nrow(paneldf) 
    robSDF <- rbind(robSDF, paneldf)
    
  }
  
  # Produce transformed variables
  dvars <- setNames(as.data.frame(as.data.frame(SP)[robSDF$rowno,dependent_var]),
                    dependent_var)
  if (length(independent_vars)==1) {
    ivars <- setNames(as.data.frame(as.data.frame(SP)[robSDF$rowno,independent_vars]),
                      independent_vars)
  } else {
    ivars <- as.data.frame(SP)[robSDF$rowno,independent_vars]
  }
  robSDF <- cbind(robSDF, dvars, ivars)
  robSDF <- robSDF[,c('angle','group','rowno','index',dependent_var,independent_vars)]
  rowno  <- robSDF[,"rowno"]
  robSDF <- robSDF[,-grep('rowno',colnames(robSDF))]
  minuslag <- function(x,n) x - dplyr::lag(x, n)
  A <- robSDF %>%
    group_by(angle, group) %>%
    arrange(angle, group, index) %>%
    mutate_all(minuslag, 1)
  B <- A %>%
    mutate_all(minuslag, 1)
  A <- setNames(as.data.frame(A[,c(dependent_var,independent_vars)]),
                paste0('sfd.', c(dependent_var,independent_vars)))
  B <- setNames(as.data.frame(B[,c(dependent_var,independent_vars)]),
                paste0('sdd.', c(dependent_var,independent_vars)))
  robSDF$rowno <- rowno
  robSDF <- cbind(robSDF, A, B)
  rownames(robSDF) <- 1:nrow(robSDF)
  
  return(robSDF)
  
}

run_sfd_models <- function(city_name, dep, ind, model_name){
  
  sfd_save <- paste0("./cfho_analysis/sfd/", city_name, "_prepped_sfd.geojson")
  df_city <- st_read(sfd_save)
  
  # Check sensitivity of model to targeted &  reduced controls - do we arrive
  # at similar findings, if we limit to only immediate neighbors of treated sites 
  # for constructing comparisons?
  
  if (analysis_type == "sensitivity_test_neighbors"){
    
    df_city['index'] <- 1:dim(df_city)[1]
    treat_sites <- df_city[df_city$cfho_any == 1,]$index
    
    neigh <- st_intersects(df_city, df_city)
    neigh <- neigh[treat_sites]
    neigh <- unique(c(unlist(neigh), treat_sites))
    
    df_city <- df_city[df_city$index %in% neigh,]
      
  }
  
  sfd_city <- construct_sfd(spatial_df = df_city,
                               dependent_var = dep,
                               independent_vars = ind,
                               rotation = seq(0, 330, 30),
                               plot = F,
                               tol = 0)
  
  # Set up summary results table which holds overall estimates across models
  # and individual model results for a given angle:
  full_results <- list()
  full_results[['sfd_data']] <- sfd_city
  
  # For each angle, use calculated spatial first-differences within a OLS
  # regression using conley standard errors (i.e. account for clustered spatial
  # correlations)
  
  # Calculate centroid for each block, then add as lat/lon to the original dataset:
  df_conley <- st_centroid(df_city)
  
  # Note: We need to index rows in the original dataset to merge w/ the SFD
  # calculations:
  setDT(df_conley)
  
  df_conley[, lon := st_coordinates(geometry)[, 1]]
  df_conley[, lat := st_coordinates(geometry)[, 2]]
  
  df_conley[, rowno := .I]
  
  df_conley <- df_conley[, .(block_geoid, lon, lat, rowno)]
  
  # Add on location
  sfd_city <- join(sfd_city, df_conley, by = "rowno", type = "left")

  spec <- paste0(paste0("sfd.", dep), " ~ ", paste0("sfd.", ind, collapse = " + "))
  
  for (a in unique(sfd_city$angle)){
    
      df_a <- sfd_city[sfd_city$angle == a,]
      
      # Get dataframe of coefficient table
      lm_con <- conleyreg(spec, data = df_a, lat = 'lat', lon = 'lon', dist_cutoff = 200, 
                          gof = T, vcov = F)
      
      coef_res <- as.data.table(lm_con$coef[,])
      coef_res[, var := c("intercept", ind)]
      coef_res[, angle := a]
      coef_res[, city := city_name]
      
      sum_res <- data.table("nobs" = lm_con$nobs,
                            "rsquared" = lm_con$r.squared,
                            "rsquared_adj" = lm_con$adj.r.squared,
                            "model_name" = model_name)
      
      # Get dataframe of covariance-variance matrix
      lm_con <- conleyreg(spec, data = df_a, lat = 'lat', lon = 'lon', dist_cutoff = 200, 
                          gof = T, vcov = T)
      
      vc_res <- as.data.table(lm_con$vcov)
      vc_res[, var := c("intercept", ind)]
      vc_res[, angle := a]
      vc_res[, city := city_name]
      
      res_list = list("coef" = coef_res, "vcov" = vc_res, "sum" = sum_res)
      
      results_name <- paste0("angle=", a)
      full_results[[results_name]] <- res_list
      
  }
  
  return(full_results)
  
}

cfho_cities <- c("fremont", "hayward", "riverside", "san_diego_unincorporated")

if (prep_cfho_locs){
  lapply(cfho_cities, construct_analysis_df)
}

dep_var <- "evict_count"

adjusted_spec <- c("number_rental_units_100", "pop_black", "per_capita_income_10k",
                   "pop_asian", "pop_white", "pop_native_american", "pop_latin_hispanic")

sfd_results_any_unadj <- lapply(cfho_cities, run_sfd_models, 
                            dep = dep_var, 
                            ind = c("cfho_any", "number_rental_units_100"),
                            model_name = "cfho_any_unadj")

names(sfd_results_any_unadj) <- cfho_cities

setwd(model_dir)
saveRDS(sfd_results_any_unadj, "cfho_any_unadj.rds")
setwd(main_dir)

sfd_results_any_adj   <- lapply(cfho_cities, run_sfd_models, 
                            dep = dep_var, 
                            ind = c("cfho_any", adjusted_spec),
                            model_name = "cfho_any_adj")

names(sfd_results_any_adj) <- cfho_cities

setwd(model_dir)
saveRDS(sfd_results_any_adj, file = "cfho_any_adj.rds")
setwd(main_dir)

sfd_results_count_unadj <- lapply(cfho_cities, run_sfd_models, 
                            dep = dep_var,
                            ind = c("cfho_num", "number_rental_units_100"),
                            model_name = "cfho_count_unadj")

names(sfd_results_count_unadj) <- cfho_cities

setwd(model_dir)
saveRDS(sfd_results_count_unadj, "cfho_count_unadj.rds")
setwd(main_dir)

sfd_results_count_adj   <- lapply(cfho_cities, run_sfd_models, 
                            dep = dep_var,
                            ind = c("cfho_num", adjusted_spec),
                            model_name = "cfho_count_adj")

names(sfd_results_count_adj) <- cfho_cities

setwd(model_dir)
saveRDS(sfd_results_count_adj, "cfho_count_adj.rds")
setwd(main_dir)

if (save_tables){
  
  library(Hmisc)

  # Generate model tables and save separately. Also save "channels" by model to investigate gradient of first-difference trends.
  models <- list(sfd_results_any_unadj, sfd_results_any_adj, sfd_results_count_unadj, sfd_results_count_adj)
  
  # For each city and model, construct a model tables with y-axis for ordered covariates, and x-axis for
  # model angles. Each cell contains the mean estimate and se, along with a star for significance. Along 
  # the bottom row, add number of units, number of treated units, and R^2 adjusted:
  
  old_varnames <- c("cfho_num", "cfho_any", "number_rental_units_100", "per_capita_income_10k", "pop_asian", "pop_black", "pop_white",
                    "pop_native_american", "pop_latin_hispanic", "treat_n", "total_n", "r2_adj")
  
  new_varnames  <- c("Number of CFHP units", "One or more CFHP units", "Rental Units (in hundreds)", "Per Capita Income (in $10,000 dollars)", 
                     "Asian (Pop %)", "White (Pop %)", "Black (Pop %)", "American Indian/Alaskan Native (Pop %)", "Latin/Hispanic (Pop %)", 
                     "Treated units", "Total units", "Adjusted R-squared")
  
  old_modelnames <- c("cfho_any_unadj", "cfho_any_adj", "cfho_count_unadj", "cfho_count_adj")
  new_modelnames <- c("Unadjusted - binary treatement", "Adjusted - binary treatment", "Unadjusted - continuous treatment",
                      "Adjusted - continuous treatment")
  
  old_citynames <- cfho_cities
  new_citynames <- c("Fremont", "Hayward", "Riverside", "San Diego County")
  
  sig_figs = 2
  options(scipen = sig_figs)
  setwd(model_dir)

  # Save tables. 

  for (model in models){
    for (city in cfho_cities){
      
      res <- model[[city]]
      dd  <- setDT(res$sfd_data)
      
      collected_res <- list()
      
      treat_var <- names(dd)[names(dd) %like% 'cfho'][1]
      
      for (a in names(res)[names(res) %like% "angle"]){
        
        # Hold onto strings for the angle rotation and model name
        ang <- unique(res[[a]][['coef']]$angle)
        mod <- res[[a]]$sum$model_name
        
        betas <- data.table("variable" = res[[a]][['coef']]$var,
                            "beta_mean" = res[[a]][['coef']]$Estimate,
                            "beta_se" = res[[a]][['coef']]$`Std. Error`,
                            "p" = res[[a]][['coef']]$`Pr(>|t|)`,
                            "angle" = a) 
        
        # Keep results readable but without truncating a value to zero. Always
        # display same number of decimal places across units to make table
        # easier to read.
        betas[, `:=`(beta_mean = format(round(beta_mean, sig_figs), nsmall = sig_figs),
                     beta_se   = format(round(beta_se, sig_figs), nsmall = sig_figs))]
        
        # Convert numerics into a nicely formatted table cell, with stars for p-values
        # and separate parentheticals for se
        
        #ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p <= 0.1, "*", ""))),
        
        betas[, table_print := paste0(beta_mean, 
                                      ifelse(p < 0.05, "*", ""),
                                      " (", beta_se, ")")]
        
        betas <- betas[, .(variable, table_print, angle)]
        
        # Separately, get goodness of fit statistic and counts of total units.
        treat_n <- ifelse(treat_var == "cfho_any", 
                          dim(dd[angle == ang & cfho_any == 1,])[1],
                          dim(dd[angle == ang & cfho_num >= 1,])[1])
        
        gof <- data.table("treat_n" = treat_n,
                          "total_n" = dim(dd[angle == ang  & !is.na(sfd.evict_count),])[1],
                          "r2_adj" = res[[a]]$sum$rsquared_adj,
                          "angle" = a)
        
        gof <- setDT(melt(gof, id.vars = "angle", variable.name = "variable", value.name = "table_print"))
        
        #Make all numbers display
        gof[, table_print := formatC(round(table_print, sig_figs), drop0trailing = T)]
        
        # Combine coefficients and gof into single dataset
        interim_res <- setDT(rbind(betas, gof))
        collected_res[[a]] <- interim_res
        
      }
      
      # Collect all results by angle and cast wide
      df_res <- rbindlist(collected_res)
      df_res <- dcast(df_res, formula = variable ~ angle, value.var = "table_print")
      
      # Remove intercept from final output and add on nice names for variables:
      df_res <- df_res[variable != "intercept"]
      
      df_res[, variable := mapvalues(variable, old_varnames, new_varnames)]
      df_res[, variable := factor(variable, levels = new_varnames)]
      df_res <- plyr::arrange(df = df_res, variable)
  
      #Set order of columns so that angles are ordered smallest to largest:
      angle_order <- paste0("angle=", seq(0, 330, 30))
      df_res <- df_res[, c("variable", angle_order), with = F]
      
      setnames(df_res, names(df_res), c("", paste0(seq(0, 330, 30), "°")))
      
      write.csv(df_res, paste0(model_dir, "/", mod, "_", city, ".csv"), 
                row.names = F, quote = )
      
    }  
  }
}
