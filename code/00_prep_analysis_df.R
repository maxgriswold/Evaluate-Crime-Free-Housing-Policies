# Prep analysis dataset for DID and conduct DID analyses
# Max Griswold
# 3/21/22

rm(list = ls())

library(data.table)
library(plyr)
library(dplyr)
library(tidycensus)
library(Hmisc)

setwd("./data")
census_key <- ""

# Create plots?
make_plots <- F
check_ucr  <- F

# Load concatennated files from Jacob Kaplan. Too large to store in Repo, needs to be downloaded.

# https://www.openicpsr.org/openicpsr/project/100707/version/V17/view
# https://www.icpsr.umich.edu/web/ICPSR/studies/35158

ucr <- fread("./data/ucr/ucr_1960_2020_offenses_known.csv")

# Get treatment data
treat <- fread("./data/treated_locs.csv", fill = T)

# Load decennial census estimates.
# Obtained through IPUMS NHGIS
df_decennial <- fread('./data/decennial_processed.csv')


# UCR distinguishes between two broad categories of crime: those committed
# against persons & those committed against property. Following suggestions for aggregation in
# IPCSR UCR documentation:

keep_vars <- c("ori", "year", "state", "crosswalk_agency_name",
               "number_of_months_missing", "geoid", 
               names(ucr)[names(ucr) %like% "actual_"])

# Reconstruct geoid from fips codes. Only hold onto the first six characters
# for GEOID, which correspond to CDP (essentially, we're upscaling more
# granular geographies to be at the CDP-level). Make sure to include
# leading zeros so that fips place code is always six characters

ucr[, fips_place_code := formatC(fips_place_code, width = 5,
                                 format = "d", flag = "0")]

ucr[, geoid := paste0(fips_state_code, fips_place_code)]
ucr[, geoid := sub("^(\\d{6})", "\\1", geoid)]

ucr <- ucr[, keep_vars, with = F]

# Hold onto rows if data is in California:
#ucr <- ucr[state == "california",]

# Remove locations missing a geoid or agency_name
ucr <- ucr[(!((geoid %like% "NA")|crosswalk_agency_name == ""))]

ucr <- ucr[year %in% c(2000, 2009:2020), ]

# Melt data by crime_type, then calculate totals for each crime type by 
# census place. 
ucr <- melt(ucr, id.vars = c("ori", "year", "state", "number_of_months_missing",
                             "geoid", "crosswalk_agency_name"), 
            variable.name = "crime_type", value.name = "crime_count")

# Only hold onto specific crime types
ucr[, crime_type := gsub("actual_", "", crime_type)]

# Collapse counts by year, CDP, crime-type. Also aggregate number of LEAs informing
# counts, along with cumulative months missing:
ucr[, number_of_reporting_agencies := .N, by = c("crime_type", "geoid", "year")]
ucr[, number_of_months_missing := sum(.SD$number_of_months_missing), by = c("crime_type", "geoid", "year")]

ucr[, crime_count := sum(.SD$crime_count, na.rm = T), by = c("crime_type", "geoid", "year")]
ucr <- setDT(unique(ucr[, .(geoid, year, crime_type, crime_count, 
                            number_of_months_missing, number_of_reporting_agencies)]))

# Only hold onto cities within CA using TIGER files (so drop tribal areas and non-cbsa). 
# Simultaneously, extract pop counts from census to construct rate vars.

extract_census <- function(y, geo = "place"){
  
  # Grab tigris shapefiles w/ renter pop and total pop
  var_list <- data.table(variable_name = 
                           c("median_age", "median_age_male", "median_age_female",
                             "total_population", "pop_white", "pop_black", "pop_asian",
                             "pop_pacific_islander", "pop_native_american", "pop_other", "pop_immigrant",
                             "pop_male", "pop_female", "pop_latin_hispanic", "number_rental_units",
                             "number_housing_units", "number_vacant_units", "median_income", 
                             "agg_mortgage", "median_rent", "ratio_income_poverty_total",
                             "number_income_poverty_below_0.49", "number_income_poverty_0.5_to_0.99", 
                             "number_income_poverty_1_to_1.24", "number_income_poverty_1.25_to_1.49",
                             "pop_tenants", "pop_owners", "renter_household_total",
                             "renter_black", "renter_native_american", "renter_asian",
                             "renter_white", "renter_hispanic_latine", "renter_pacific_islander",
                             "renter_multiple_race", "renter_other", "renter_white_alone"),
                         variable_id   = 
                           c("B01002_001", "B01002_002", "B01002_003",
                             "B01003_001", "B02008_001", "B02009_001", "B02011_001",
                             "B02012_001", "B02010_001", "B02013_001", "B05012_003",
                             "B01001_002", "B01001_026", "B03001_003", "B25002_002",
                             "B25002_001", "B25002_003", "B19013_001",
                             "B25082_002", "B25064_001", "C17002_001",
                             "C17002_002", "C17002_003", 
                             "C17002_004", "C17002_005", 
                             "B25033_008", "B25003_002", "B25003_003",
                             "B25003B_003", "B25003C_003", "B25003D_003",
                             "B25003A_003", "B25003I_003", "B25003E_003",
                             "B25003G_003", "B25003F_003", "B25003H_003"))
  
  education_list <- data.table(variable_name = paste0("education=", rep(c(0, 2.5, 6.5, 7.5, 10, 11, 12, 12, 13, 14, 15, 15, 17, 19, 19, 21), 2)),
                               variable_id = sprintf("B15002_0%s", sprintf("%02d", c(seq(3, 18, 1), seq(20, 35, 1)))))
  
  var_list <- rbind(var_list, education_list)
  
  # Crosswalk 2020 census block groups to census tracts prior to 2020.
  # Based on NHGIS guidelines: https://www.nhgis.org/geographic-crosswalks
  if (y >= 2020 & geo == "tract"){
    census_vars <- get_acs(geography = "block group", year = y, variables = var_list$variable_id, 
                           state = "CA", geometry = F, key = census_key, survey = "acs5") %>%
      setDT %>%
      .[, .(GEOID, variable, estimate, moe)] %>%
      setnames(., names(.), c("geoid_new", "variable_id", "mean", "margin_of_error")) %>%
      .[, `:=`(mean = as.numeric(mean),
               margin_of_error = as.numeric(margin_of_error))] %>%
      join(., var_list, by = c("variable_id"), type = "left")
    
    tract_names <- get_acs(geography = "tract", year = 2019, variables = "B01002_001", 
                           state = "CA", geometry = F, key = census_key, survey = "acs5") %>%
      setDT %>%
      .[, .(GEOID, NAME)] %>%
      setnames(., names(.), c("geoid", "location"))
    
    cw <- fread("./data/nhgis_bg2020_tr2010.csv")
    cw <- cw[, .(bg2020ge, tr2010ge, wt_adult)]
    setnames(cw, names(cw), c("geoid_new", "geoid", "weight_adult_interpolated"))
    
    # NHGIS provides fips codes w/o a leading zero. So I'll set the census extraction to
    # a numeric to remove the leading zero in this dataset.
    cw[, geoid_new := as.numeric(geoid_new)]
    cw[, geoid := as.numeric(geoid)]
    
    census_vars[, geoid_new := as.numeric(geoid_new)]
    tract_names[, geoid := as.numeric(geoid)]
    
    census_vars <- join(census_vars, cw, by = "geoid_new", type = "left")
    census_vars <- join(census_vars, tract_names, by = "geoid", type = "left")
    
    census_vars[, mean := weighted.mean(.SD$mean, w = .SD$weight_adult_interpolated), by = c("geoid", "variable_name")]
    census_vars <- unique(census_vars[, .(geoid, location, variable_id, mean, variable_name)])
    
  }else{
    census_vars <- get_acs(geography = geo, year = y, variables = var_list$variable_id, 
                           state = "CA", geometry = F, key = census_key, survey = "acs5") %>%
      setDT %>%
      .[, .(GEOID, NAME, variable, estimate, moe)] %>%
      setnames(., names(.), c("geoid", "location", "variable_id", "mean", "margin_of_error")) %>%
      .[, `:=`(mean = as.numeric(mean),
               margin_of_error = as.numeric(margin_of_error))] %>%
      join(., var_list, by = c("variable_id"), type = "left")
  }
  
  #Process the education variables separately.
  education_vars <- census_vars[variable_name %like% "education", ]
  pop_vars       <- census_vars[variable_name %like% "pop",]
  renter_vars    <- census_vars[variable_name %like% "renter",]
  census_vars    <- census_vars[!(variable_name %like% "education" | variable_name %like% "pop" ), ]
  
  #Combine education variables into average years of schooling
  education_vars[, years := as.numeric(gsub(".*\\=", "", variable_name))]
  education_vars[, pop_years := sum(.SD$mean*.SD$years, na.rm = T), by = "location"]
  education_vars[, pop := sum(.SD$mean, na.rm = T), by = "location"]
  
  education_vars[, avg_edu := pop_years/pop]
  education_vars <- unique(education_vars[, .(location, avg_edu)])
  
  #Convert counts of renter housholds into percentages
  renter_vars <- dcast(renter_vars, location ~ variable_name, value.var = "mean")
  for (household in names(renter_vars)[names(renter_vars) %like% "renter" & !(names(renter_vars) %like% "total")]){
    
    renter_vars[, (household) := get(household)/renter_household_total]
    
  }
  
  #Convert counts of population into percentages
  pop_vars <- dcast(pop_vars, location ~ variable_name, value.var = "mean")
  for (populate in names(pop_vars)[names(pop_vars) %like% "pop" & !(names(pop_vars) %like% "total")]){
    
    pop_vars[, (populate) := get(populate)/total_population]
    
  }
  
  census_vars <- dcast(census_vars, geoid + location ~ variable_name, value.var = "mean")
  #census_vars[, unemployed_percentage := number_unemployed/number_labor_force]
  
  #Combine poverty percentage categories
  census_vars[, percent_poverty_150 := (number_income_poverty_below_0.49 + 
                                          number_income_poverty_0.5_to_0.99 + 
                                          number_income_poverty_1_to_1.24 + 
                                          number_income_poverty_1.25_to_1.49)/ratio_income_poverty_total]
  
  census_summary <- census_vars[, .(geoid, location, median_income, median_rent, 
                                    percent_poverty_150, number_rental_units)]
  census_summary <- join(census_summary, education_vars, by = "location", type = "left")
  census_summary <- join(census_summary, pop_vars, by = "location", type = "left")
  census_summary <- join(census_summary, renter_vars, by = "location", type = "left")
  
  
  census_summary[, year := y]
  census_summary[, geoid := as.numeric(geoid)]
  
  return(census_summary)
  
}

years <- 2009:2020
census_vars <- setDT(ldply(years, extract_census))
census_vars[, geoid := as.character(as.numeric(geoid))]

# Merge on rent CPI and calculate rent prices & income in 2020 dollars:
cpi_rent <- fread("./cpi_rent_shelter.csv")
setnames(cpi_rent, names(cpi_rent), c("year", "cpi_b82"))
cpi_rent[, year := as.numeric(gsub("1/1/", "", year))]

# Recalculate CPI so that base year is 2020 rather than 2017:
cpi_rent[, cpi_b20 := cpi_rent[year == 2020, cpi_b82]]
cpi_rent[, cpi_b20 := cpi_b82/cpi_b20]

census_vars  <- join(census_vars, cpi_rent, by = "year", type = "left")

census_vars[, median_rent_adj := median_rent*cpi_b20]
census_vars[, median_income_adj := median_income*cpi_b20]

pop_vars <- census_vars[, .(geoid, location, total_population, year)]

# Reshape UCR long so each crime category is a separate column
ucr <- dcast(ucr, geoid + year + number_of_months_missing + number_of_reporting_agencies ~ crime_type, value.var = "crime_count")

test <- setDT(join(pop_vars, ucr, by = c("geoid", "year"), type = "right"))
test <- test[!is.na(all_crimes)]

# For purposes of checking the merge to make sure it behaved appropriately, recast
# crime-type wide in the UCR dataset. Then, once I've confirmed the merge is
# working, I'll comment out this code:

# ucr <- dcast(ucr, geoid + year ~ crime_type, value.var = "crime_count")
# 
# test <- setDT(join(total_population, ucr, by = c("geoid", "year"), type = "left"))
# check <- test[is.na(all_crimes), unique(location)]
# 
# # From "View"ing check, I could see most of these locations are unincorporated places
# # (i.e. would not have a local government passing a CFHO and would be serviced
# # by the county sheriff). Let's look and see which locations are missing and
# # aren't CDP:
# 
# check[!(check %like% "CDP")]
# 
# # Some notable misses, about 66 cities. Some of these might be locations that
# # aren't reporting to UCR. However, I can see Burbank is in this set and
# # they should be reporting to UCR, based on googling.
# 
# # Confirmed burbank is in the UCR database. Accordingly, the merge seems to be
# # messing up if there's a leading zero in the fips_place_code (which is implied
# # by checking the burbank fips code in the UCR dataset, which is missing a leading
# # zero, based on census fips database).
# 
# # Based on this finding, I'll add leading zeros to all fips_place_codes so that
# # they have 5 characters. I'll also leave this investigation as comments for
# # posterity sake. After doing this, there were only 5 small cities left w/o
# # ucr counts

# There should be 480 cities in CA (459 cities + 21 towns). We have 505 units.
# Are the extras CDPs w/ LEAs?

length(unique(test$location)) - length(unique(test[location %like% "CDP"]$location))

#Yes, minus the 5 missing cities! We're in the clear

# Let's now reshape long, investigate extreme values, and make some plots:

ucr <- melt(ucr, id.vars = c("geoid", "location", "year", "total_population", 
                             "number_of_months_missing", "number_of_reporting_agencies"), 
            variable.name = "crime_type", value.name = "crime_count")

# Keep crime_type as a character (not a factor):
ucr[, crime_type := as.character(crime_type)]

# Remove "California" from location name
ucr[, location := gsub(", California", "", location)]

# Add on total population for locations in 2000:
df_pop_old <- df_decennial[, .(geoid, total_population_2000)]
df_pop_old[, year := 2000]

ucr <- join(ucr, df_pop_old, by = c('year', 'geoid'), type = 'left')
ucr[year == 2000, total_population := total_population_2000]

# Convert counts to rates, along with adding a standardized transformation
ucr[, crime_rate_10k := (crime_count/total_population)*10000]
ucr[, crime_z := (crime_rate_10k - mean(.SD$crime_rate_10k, na.rm = T))/sd(.SD$crime_rate_10k, na.rm = T),
    by = c("crime_type")]

# Look at decile tables & density distributions for each crime type across
# cities

percentiles <- seq(0, 1, 0.1)
percentile_names <- factor(paste0(percentiles*100, "th"), 
                           levels = paste0(percentiles*100, "th"))

summary_percentiles <- copy(ucr)
summary_percentiles <- summary_percentiles[, .(year, crime_type, crime_z)]

# Standardize rates by crime_type:
summary_percentiles <- summary_percentiles[, quantile(.SD$crime_z, 
                                                      probs = percentiles,
                                                      na.rm = T), 
                                           by = c("crime_type", "year")]

setnames(summary_percentiles, "V1", "crime_z")

summary_percentiles[, percentile := rep(percentile_names, 
                                        times = dim(summary_percentiles)[1]/length(percentile_names))]

# Cut continuous z-scores into reasonable bins for displaying on heatmap. Maintain
# actual values to fill in squares
summary_percentiles[, crime_binned := cut(crime_z, breaks = c(-1, -0.1, -0.5, -0, 0.5, 1, 5, 100),
                                          include.lowest = T, right = T)]

if (make_plots){
  percentile_heatmap <- function(y){
    
    df_plot <- summary_percentiles[year == y,]
    
    percentile_heatmap <- ggplot(df_plot, aes(x = percentile, y = crime_type, fill = crime_binned)) +
      geom_raster() +
      labs(title = y) +
      scale_fill_manual(values = c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000')) +
      theme_bw()
  }
  
  heatmaps <- lapply(unique(summary_percentiles$year), percentile_heatmap)
  heatmaps <- marrangeGrob(heatmaps, nrow = 1, ncol = 1)
  
  ggsave("plots/ucr_crime_type_heatmaps.pdf", heatmaps, width = 21, height = 29.7, units = "cm")
  
}

# Based on heatmaps above,
# Cut off density plots at the 90th percentile to make 
# distributions clearer (for example, a few cities
# w/ small pop counts, like Industry or Vernon, have
# rates exceeding 10k):

#Check "extreme" values
sort(table(ucr[crime_z >= 2, ]$location))
sort(table(ucr[crime_z >= 2, ]$crime_type))
#View(ucr[crime_z >= 2, ])

# Observations on above:
# City of Industry and Vernon are primary locations w/ large values
# Most locations with large values are in low-pop areas (some exceptions)
# Lots of extreme values are in the arson category (so small numbers again)

# Likely want to outlier Vernon & Industry + exclude arson from any analysis
# (which seems fine - I can't imagine arson is the primary endpoint LEAs
# are considering when they implement the policy)

if (make_plots){
  crime_type_density_by_year <- function(crime){
    
    df_plot <- ucr[crime_type == crime,]
    
    # Temporarily remove obs above 5 sd
    df_plot <- df_plot[crime_z < 5]
    
    plot <- ggplot(df_plot, aes(x = crime_rate_10k)) +
      geom_histogram(bins = 50) +
      facet_wrap(~year) +
      theme_bw() +
      labs(title = crime,
           x = "Crime rate per 10k",
           y = "# of obs") +
      theme(strip.background = element_blank())
    
    return(plot)
  }
  
  density_plots <- lapply(unique(ucr$crime_type), crime_type_density_by_year)
  density_plots <- marrangeGrob(density_plots, nrow = 1, ncol = 1)
  
  ggsave("./plots/ucr_histograms_crime_type.pdf", density_plots, width = 21, height = 29.7, units = "cm")
}

# Based on above densities, lots of crime types that should likely be excluded:
# e.g. remove arson crime-types, which aren't likely to be modified by policy 
# (and don't seem to be collected frequently)

# Take a peek at time-series by location to see if there's any
# notable jumps. 

# In addition, break out crime types into parent/child categories:
ucr[, crime_subtype := crime_type]
ucr[, crime_type := NA]

ucr[crime_subtype %like% c("burg"), crime_type := 'burglary_total']
ucr[crime_subtype %like% c("robbery"), crime_type := 'robbery_total']
ucr[crime_subtype %like% c("assault"), crime_type := 'assault_total']
ucr[crime_subtype %like% c("index_total"), crime_type := 'index_total']
ucr[crime_subtype %like% c("mtr"), crime_type := 'motor_theft_total']

ucr <- ucr[!is.na(crime_type)]
ucr[, `:=`(location = NULL, total_population_2000 = NULL)]

if (make_plots){
  crime_type_by_year <- function(type, city){
    
    if (type != "pop" & type != "months missing"){
      df_plot <- ucr[location == city & crime_type == type,]
    }
    
    if (type != "pop" & type != "months missing"){
      plot <- ggplot(df_plot, aes(x = year, y = crime_rate_10k, color = crime_subtype)) +
        geom_line(size = 1) +
        labs(title = paste0(type, " rates per 10k in ", city),
             x = "Year",
             y = "Rate per 10k") +
        theme_bw() +
        theme(strip.background = element_blank())
    } else if (type == "pop"){
      plot <- ggplot(df_plot, aes(x = year, y = total_population)) +
        geom_line(size = 1) +
        labs(title = paste0("Population estimates in ", city),
             x = "Year",
             y = "Population") +
        theme_bw() +
        theme(strip.background = element_blank())
    } else{
      plot <- ggplot(df_plot, aes(x = year, y = number_of_months_missing)) +
        geom_bar() +
        labs(title = paste0("Number of months not reported across ", 
                            unique(df_plot$number_of_reporting_agencies),
                            " agencies"),
             x = "Year",
             y = "Number of months") +
        theme_bw() +
        theme(strip.background = element_blank())
    }
    return(plot)
    
  }
  
  multiple_args <- expand.grid("type" = c(unique(ucr$crime_type), c("pop", "missingness")), "city" = unique(ucr$location))
  
  crimes_year <- mlply(multiple_args, crime_type_by_year)
  crimes_year <- marrangeGrob(crimes_year, nrow = 3, ncol = 2)
  
  ggsave("plots/crimes_by_city_year.pdf", crimes_year, height = 21, width = 29.7, units = "cm")
  
}

# Reviewing the above plots, consider the following:

# Outlier: Amador, Belvedere, Blue Lake, Bradbury, Calipatria, Colfax, Corte Madera,
# Del Rey Oaks, Dorris, Dunsmuir, Etna, Ferndale, Fort Jones, Half Moon Bay,
# Live Oak, Loyalton, Maricopa, Millbrae,
# Consider outlier: Clayton
# Analyze: Assault, motor vehicle thefts, burglary, robberies
# Consider: Take 3/5-year rolling average, given some random years w/ missing
# reports

# Other notes:
# Odd robbery results in Anderson; odd decreasing trends; Angels, Del Mar, Dos Palos,
# Fairfax, Fillmore, Fowler, Gridley, Healdsberg
# Atherton, Guadalupe are all over the place
# Odd burglary, robbery results in Avalon

# Restrict UCR data to aggregate categories, reshape wide, then add on census
# vars, treatment var, and FHEO/sundown town data

ucr_agg <- ucr[, .(geoid, year, number_of_months_missing, number_of_reporting_agencies,
                   crime_type, crime_subtype, crime_rate_10k)]

# Remove total categories and "simple" categories (reason for aggregating myself rather
# than using index - assault totals exclude assaults with a gun)
ucr_agg <- ucr_agg[!(crime_subtype %like% "assault_total"|crime_subtype %like% "burg_total"|crime_subtype %like% "robbery_total"|
                       crime_subtype %like% "mtr_veh_theft_total")]

ucr_agg[, crime_rate_10k := sum(.SD$crime_rate_10k), by = c("geoid", "year", "crime_type")]

ucr_agg <- setDT(unique(ucr_agg[, .(geoid, year, number_of_months_missing, 
                                    number_of_reporting_agencies, crime_type, crime_rate_10k)]))

ucr_agg <- dcast(ucr_agg, geoid + year + number_of_months_missing + number_of_reporting_agencies ~ crime_type,
                 value.var = "crime_rate_10k")

# Also calculate changes in crime rate by location between 2000 & 2010:
change_vars <- c("assault_total", "burglary_total", "index_total") 

calc_2000_2010_change <- function(var){
  
  var_2000 <- paste0(var, "_2000")
  var_2010 <- paste0(var, "_2010")
  
  # Take 5-year mean for both 2000 & 2010
  df_2000 <- ucr_agg[year == 2000, c("geoid", var), with = F]
  setnames(df_2000, var, var_2000)
  
  df_2010 <- ucr_agg[year == 2010, c("geoid", var), with = F]
  setnames(df_2010, var, var_2010)
  
  df_change <- join(df_2000, df_2010, type = 'inner')
  
  var_change <- paste0(var, "_2000_2010")
  df_change[, (var_change) := (get(var_2010) - get(var_2000))/get(var_2000)]
  df_change <- df_change[, c("geoid", var_change), with = F]
  
  hold <- join(ucr_agg, df_change, by = 'geoid', type = 'left')
  hold <- hold[, var_change, with = F]
  
  return(hold)
  
}

for (v in change_vars){
  ucr_agg <- cbind(ucr_agg, calc_2000_2010_change(v))
}

df_analysis <- setDT(join(ucr_agg, census_vars, by = c('geoid', "year")))

# Add on decennial changes:
loc_names <- unique(census_vars[, .(geoid, location)])
loc_names[, location := gsub(", California", "", location)]

df_analysis <- join(df_analysis, loc_names, by = c("geoid"), type = "left")
df_analysis <- join(df_analysis, df_decennial, by = c("geoid"), type = "left")

# Clean up census names
df_analysis[, location := gsub(" city| town", "", location)]
df_analysis[, location := tolower(location)]

# Look at overlap between treatment city names and census/ucr:
treat_cities <- unique(treat[treated == 1, ]$location)
cov_cities   <- unique(df_analysis$location)

treat_cities[!(treat_cities %in% cov_cities)]

# Only eastvale isn't matching up covariate data, which makes sense
# since eastvale was incorporated after 2010 census and relies on county
# to provide policing (so no direct UCR estimates in harmonized cities).
# It would be nice to include eastvale but not sure how I could do so,
# given lack of covariate data

df_analysis <- join(df_analysis, treat, by = c("geoid"), type = "left")
df_analysis <- df_analysis[!is.na(location) & !(location %like% "cdp"),]

# Set treatment status to "1" if implementation date is unknown; otherwise
# set treatment status to "0" prior to treatment, then set treatment status
# to 1 post treatment (set to 0 if repealed)

setorder(df_analysis, location, year)

df_analysis[is.na(treated), treated := 0]

df_analysis[year < year_implemented, treated := 0]
df_analysis[year == year_implemented & policy_type == "cfmhp", treated := 1]

#df_analysis[year >= repeal_date, treated:= 0]

# Also, make sure to turn off treatment for all locations
# where primary policy type is not cfmhp:

df_analysis[policy_type == "cfmhp"|policy_type == "good_neighbor", ever_treated := 1]
df_analysis[is.na(policy_type), ever_treated := 0]

# Make same judgments for nuisance ordinances
df_analysis[is.na(nuisance)|nuisance == 0, ever_nuisance := 0]
df_analysis[nuisance == 1, ever_nuisance := 1]

df_analysis[year < year_implemented & ever_nuisance == 1, nuisance := 0]

#Outlier vernon and industry since no one lives there and they create huge outliers
# in the outcome variables
df_analysis <- df_analysis[!(location %in% c("industry", "vernon"))]

# For covariates, set level equal to what's observed in 2009 (i.e. do not
# condition on post-treatment covariates)
covariates <- df_analysis[year == 2009, .(geoid, pop_female, pop_white, pop_immigrant,
                                          median_income, percent_poverty_150)]

new_names <- paste0(names(covariates)[names(covariates) != "geoid"], "_2009")
setnames(covariates, names(covariates), c("geoid", new_names))

df_analysis <- join(df_analysis, covariates, by = "geoid", type = 'left')
df_analysis <- df_analysis[year >= 2009,]

# Outlier sites missing estimates for the outcome variable (i.e. reported 0
# crimes over the time period)

df_analysis[, location := gsub(", california", "", location)]
outlier <- c("amador city", "blue lake", "calipatria", "half moon bay", "millbrae", "orange cove", "san anselmo", "san carlos",
             "colfax", "corte madera", "live oak", "loyalton", "maricopa", "plymouth", "point arena", 
             "portola", "portola valley", "san joaquin", "san juan bautista", "trinidad", "tehama", "wasco", "woodside")  

df_analysis <- df_analysis[!(location %in% outlier), ]

# Add on eviction records as well for summary table purposes:
df_evict <- fread("./data/evict_census_place.csv")
df_evict[, geoid := fips_code]

df_evict <- df_evict[, .(geoid, year, evict_count, evict_rate)]

df_analysis <- join(df_analysis, df_evict, by = c("geoid", "year"), type = 'left')

# Rename outcome variables to be more accurate to what they represent:
setnames(df_analysis, c("index_total", "assault_total", "burglary_total"),
         c("total_crime_rate_10k", "assault_rate_10k", "burglary_rate_10k"))

write.csv(df_analysis, "./data/df_did_analysis.csv", row.names = F)
