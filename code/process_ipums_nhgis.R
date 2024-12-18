# Process decennial census records:
# Max Griswold
# 1/9/23

library(data.table)
library(plyr)

# You'll need to change the below to the path of the github repostiory:
setwd("C:/path/to/your/repository")

df <- fread('./data/nhgis/nhgis_place_decennial.csv')

df <- df[(STATE == "California") & (DATAYEAR %in% c(1990, 2000, 2010)),]
df <- df[, .(GISJOIN, DATAYEAR, PLACE, CL8AA, 
             CM1AA, CM1AB, CY8AA, 
             CM1AD, CM1AC, CP4AB,
             CN1AB, CM7AA)]

setnames(df, names(df), c("geoid", "year", "location", "total_population",
                          "pop_white", "pop_black", "pop_white_nonhispanic",
                          "pop_asian", "pop_indigenous", "pop_hispanic_latino",
                          "rental_unit_count", "total_unit_count"))

poverty_ratios <- c(paste0("C20A", LETTERS[1:8]))

df_nominal <- fread('./data/nhgis/nhgis_nominal.csv')
df_nominal <- df_nominal[(STATE == "California") & (YEAR %in% c("1990", "2000", "2008-2012")),]
df_nominal <- df_nominal[, c("GISJOIN", "YEAR", "PLACE", "AT5AB", "BD5AA", poverty_ratios), with = F]

# Calculate joint poverty ratio <200% percentage:
below_200 <- rowSums(df_nominal[, c(poverty_ratios), with = F])

df_nominal[, pop_200_fpl := below_200]

setnames(df_nominal, c("GISJOIN", "YEAR", "PLACE", "AT5AB", "BD5AA"), c("geoid", "year", "location", "pop_non_citizen", "per_capita_income"))
df_nominal[year == '2008-2012', year := 2010]

df_nominal <- df_nominal[, .(geoid, year, location, pop_non_citizen, per_capita_income, pop_200_fpl)]

df <- join(df, df_nominal, by = c('geoid', 'year', 'location'), type = 'left')

# Convert population counts to percentages:
df[, percent_pop_white := pop_white/total_population]
df[, percent_pop_black := pop_black/total_population]
df[, percent_pop_asian := pop_asian/total_population]
df[, percent_pop_indigenous := pop_indigenous/total_population]

df[, percent_non_white_non_hispanic := 1 - (pop_white_nonhispanic/total_population)]
df[, percent_renter := rental_unit_count/total_unit_count]
df[, percent_non_citizen := pop_non_citizen/total_population]

# Remove CDP locations to match CFHO database:
df <- df[!(location %like% "CDP")]

# Calculate changes over time, 1990-2000 & 2000-2010:
df <- df[, .(geoid, location, year, total_population, pop_black, pop_asian, pop_indigenous,
             pop_hispanic_latino, pop_white,
             rental_unit_count, pop_non_citizen, per_capita_income,
             pop_200_fpl)]

change_vars <- names(df)[3:11]
df <- dcast(geoid + location ~ year, data = df, value.var = change_vars)

percent_change <- function(decade, var, dd = df){
  
  new <- decade
  old <- decade - 10
  
  old <- paste0(var, "_", old)
  new <- paste0(var, "_", new)
  
  old <- dd[, c(old), with = F]
  new <- dd[, c(new), with = F]
  
  return((new - old)/old)
  
}

# Leave out change in total population from
# percent change calculator:
change_vars <- change_vars[2:length(change_vars)]

for (v in change_vars){
  for (y in c(2010, 2000)){
    new_name <- paste0('percent_change_', v, "_", y - 10, "_", y)
    
    df[, (new_name) := percent_change(y, v)]
    
  }
}

keep <- c("geoid", "total_population_2000", "total_population_2010", names(df)[names(df) %like% "change"])
df <- df[, keep, with = F]

# Modify geoid to match other census series:
df[, geoid := gsub("G0", "", geoid)]
df[, geoid := gsub("60", "6", geoid)]

write.csv(df, './data/acs5/decennial_processed.csv', row.names = F)

# Prep block-group data in 2019:
df <- fread('./acs5/california_block_group_2019.csv')

df <- df[, .(GEOID, ALT0E001, ALT0E002, ALUCE002, ALUCE003, 
             ALUCE004, ALUCE005, ALULE003, 
             ALWVE008, ALX5E001,
             ALZ1E001, ALZ1E008,
             AL05E001, AL2OE005)]

new_names <- c("geoid", "total_population", "pop_male", "pop_white", "pop_black",
               "pop_native_american", "pop_asian", "pop_latin_hispanic",
               "pop_above_2_fpl", "per_capita_income",
               "number_housing_units", "number_rental_units",
               "median_rent", "pop_immigrant")

setnames(df, names(df), new_names)

# Convert population counts to percentages:
for (populate in names(df)[names(df) %like% "pop" & !(names(df) %like% "total")]){
  
  df[, (populate) := get(populate)/total_population]
  
}

df[, pop_below_2_fpl := 1 - pop_above_2_fpl]
df[, geoid := as.numeric(gsub("15000US0", "", geoid))]

write.csv(df, "./data/acs5/california_block_group_2019_processed.csv", row.names = F)
