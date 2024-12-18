# Construct summary tables, figures, and maps
# Max Griswold
# 3/17/23

library(data.table)
library(units)
library(ggplot2)
library(ggstance)
library(plyr)
library(dplyr)
library(parallel)
library(gridExtra)

library(extrafont)
font_import()
loadfonts(device = "win")

library(tidycensus)
library(sf)
library(scales)

options(scipen = 999)

census_key <- "bdb4891f65609f274f701e92911b94365992028a"

setwd("C:/users/griswold/Desktop/datasets/")

df_analysis <- fread("./cfho_analysis/did/df_did_analysis.csv")

formatted_results <- function(m, s, ci = F, eff_df = NA, sig = 0.025, t_dist = F, digs = 2){
  
  # If assuming normality, use z-score for critical value (t = F). 
  # If effective degrees-of-freedom are set and t_dist = T, use student t to 
  # determine critical value:
  if (ci == T){
    if (t_dist == T){
      cv <- qt(p = sig, df = eff_df, lower.tail = F)
    }else{
      cv <- qnorm(p = sig, lower.tail = F)
    }
    
    return(paste(round(m, digs), paste0("(", round(m - cv*s, digs), ", ", round(m + cv*s, digs), ")")))
  }else{
    return(paste(round(m, digs), paste0("(", round(s, digs), ")")))
  }
}

# Use Welch–Satterthwaite approximation to estimate degrees-of-freedom for
# t-dist as a means to derive critical value for estimating CI.

# This is relatively pendantic; critical values are typically in range of
# 2.03 to 1.98 rather than routine critical value from z of 1.96.
welch_satt <- function(sd0, sd1, n0, n1){
  
  d <- ((sd0/n0) + (sd1/n1))^2/(((sd0^2)/((n0^2)*(n0 - 1))) + ((sd1^2)/((n1^2)*(n1 - 1))))
  return(d)
  
}

########################################### 
# Changes in demographics prior to policy #
###########################################

# Investigate difference between to-be-treated locations & control sites,
# based on demographic changes between 2000 - 2010:

df_change    <- df_analysis[year == 2010 & (ever_treated == 0|
                                            (ever_treated == 1 & year_implemented >= 2010)),]

# Alternatively, see if removing nuisance ordinance sites changes interpretation:
# df_change <- df_change[!is.na(ever_treated),]

compare_vars <- c("index_total_2000_2010", "assault_total_2000_2010", "burglary_total_2000_2010",
                  "percent_change_pop_asian_2000_2010", "percent_change_pop_black_2000_2010", "percent_change_pop_hispanic_latino_2000_2010",
                  "percent_change_pop_indigenous_2000_2010", "percent_change_pop_non_citizen_2000_2010", "percent_change_pop_white_2000_2010",          
                  "percent_change_rental_unit_count_2000_2010", "percent_change_total_population_2000_2010", "inquiry_rate_1000_pop",                       
                  "inquiry_count_race", "inquiry_count_national_origin")

# Several sites went from having no people of color to having some, causing
# percent to be infinite. This only effects two sites, with small populations (<1k people
# in these cities). So for now, set any Inf to NA:
df_change[df_change == Inf] <- NA

# Convert inquiry counts into rates per 10k people:
df_change[, inquiry_rate_1000_pop := (inquiry_rate_1000_pop*10)]
df_change[, inquiry_count_race := (inquiry_count_race*10000)/total_population]
df_change[, inquiry_count_national_origin := (inquiry_count_national_origin*10000)/total_population]

# Try outliering Crescent city:
df_change <- df_change[location != "crescent city"]

# Calculate CI instead of using quantiles given reviewer comments
compare_table <- df_change[, sapply(.SD, function(x) list(mean = mean(x, na.rm = T), 
                                                          sd = sd(x, na.rm = T))),
                                                          by = "ever_treated", 
                                    .SDcols = compare_vars]

new_names <- paste0(levels(interaction(c("mean_", "sd_"), compare_vars, sep = "")))

names(compare_table) <- c("treated", new_names)

compare_table <- melt(compare_table, id.vars = "treated", 
                      measure  = patterns("mean_", "sd_"),
                      value.name = c("mean", "sd"))

new_names <- c("% change: total crime rate (per 10k people)",
               "% change: assault rate (per 10k people)", 
               "% change: burglary rate (per 10k people)",
               "% change: population asian", 
               "% change: population black",
               "% change: population hispanic/latine", 
               "% change: population indigenous",
               "% change: population immigrant", 
               "% change: population white",
               "% change: Number of rental units ", 
               "% change: Total population",
               "FHEO inquiry rate per 10,000 people (total), 2013-2020", 
               "FHEO inquiry rate per 10,000 people  (race discrimination), 2013-2020",
               "FHEO inquiry rate per 10,000 people  (national origin discrimination), 2013-2020")

setattr(compare_table$variable, "levels", new_names)

compare_table <- dcast(compare_table, variable ~ treated, value.var = c("mean", "sd"))

# Calculate relative mean difference assuming normally distributed in each pop.
# To calculate pooled variance, need to also incorporate sample size for each
# respective population
compare_table[, `:=`(n_0 = dim(df_change[ever_treated == 0,])[1],
                     n_1 = dim(df_change[ever_treated == 1,])[1])]

compare_table[, mean_diff := mean_1 - mean_0]

# Calculate pooled standard deviation, based off two-sample unpaired t-test
# (assuming here that sample is normally distributed, which is mostly true
# across these variables and that samples are independent, which is likely 
# unreasonable but not something easily corrected. Reviewer wanted this).
compare_table[, sd_diff := sqrt(((sd_0^2)/n_0) + ((sd_1^2)/n_1))]
compare_table[, df_diff := welch_satt(sd_0, sd_1, n_0, n_1)]

compare_table[, Control := formatted_results(mean_0, sd_0)]
compare_table[, `To-be-treated` := formatted_results(mean_1, sd_1)]
compare_table[, Difference := formatted_results(mean_diff, sd_diff, ci = T, eff_df = df_diff)]

compare_table <- compare_table[, c("variable", "Control", "To-be-treated", "Difference"), with = F]
setnames(compare_table, names(compare_table), c("", "Control", "To-be-treated", "Difference"))

compare_table <- rbind(compare_table,
                       data.table("N", dim(df_change[ever_treated == 0,])[1], dim(df_change[ever_treated == 1,])[1], ""),
                       use.names = F)

write.csv(compare_table, "./cfho_analysis/da/pre_policy_changes_ci_revised.csv", row.names = F)

###########################################
# Summary stats: DiD units in 2009 & 2019 #
###########################################

df_did <- copy(df_analysis[year %in% c(2009, 2019) & !is.na(ever_treated)])
df_did <- df_did[ever_treated == 0|(year_implemented >= 2010 & ever_treated == 1), ]

df_did <- df_did[, .(location, year, total_population,
                     ever_treated, total_crime_rate_10k, assault_rate_10k,
                     burglary_rate_10k, median_income_adj, median_rent_adj,
                     avg_edu, number_rental_units, percent_poverty_150,
                     pop_latin_hispanic, pop_asian, pop_black, 
                     pop_immigrant, pop_native_american, pop_white)]

# Convert number of rental units into a rate:
df_did[, number_rental_units := (number_rental_units*100)/total_population]

df_did <- melt(df_did, id.vars = c("location", "year", "ever_treated"))
df_did[, ever_treated := as.factor(ever_treated)]

df_did[, `:=`(mean = mean(.SD$value, na.rm = T),
              sd   = sd(.SD$value, na.rm = T),
              n    = .N), by = c("year", "ever_treated", "variable")]

df_did <- unique(df_did[, .(ever_treated, year, variable, mean, sd, n)])
df_did <- dcast(df_did, variable + year ~ ever_treated, value.var = c("mean", "sd", "n"))

df_did[, mean_diff := mean_1 - mean_0]

# Calculate pooled standard deviation, based off two-sample unpaired t-test
# (assuming here that sample is normally distributed, which is mostly true
# across these variables and that samples are independent, which is likely 
# unreasonable but not something easily corrected. Reviewer wanted this).
df_did[, sd_diff := sqrt(((sd_0^2)/n_0) + ((sd_1^2)/n_1))]
df_did[, df_diff := welch_satt(sd_0, sd_1, n_0, n_1)]

df_did[, Control := formatted_results(mean_0, sd_0)]
df_did[, Treated := formatted_results(mean_1, sd_1)]
df_did[, Difference := formatted_results(mean_diff, sd_diff, ci = T, eff_df = df_diff)]

df_did <- df_did[, c("variable", "year", "Control", "Treated", "Difference"), with = F]

# Recode some variable names for the summary plot
new_names <- c("Total population", 
               "Total crime rate (per 10k people)", "Assault rate (per 10k people)", "Burglary rate (per 10k people)",
               "Median income (in 2020 dollars)", "Median rent price (in 2020 dollars)", "Average years of education",
               "Number of rental units (per 100 people)", "Population: % below 150% of federal poverty line",
               "Population: % Hispanic", "Population: % Asian", "Population: % Black",
               "Population: % immigrant", "Population: % indigenous", "Population: % white")

df_did[, variable := factor(mapvalues(variable, unique(df_did$variable), new_names), levels = new_names)]
df_did[, year := factor(year, levels = c(2019, 2009))]

df_did <- dcast(df_did, variable ~ year, value.var = c("Control", "Treated", "Difference"))

# Reorder columns:
df_did <- df_did[, .(variable, Control_2009, Treated_2009, Difference_2009, Control_2019, Treated_2019, Difference_2019)]

write.csv(df_did, "./cfho_analysis/did/did_summary_stats.csv", row.names = F)

#####################################
# Summary stats: SFD blocks by city #
#####################################

sfd_files <- "./cfho_analysis/sfd/"
files <- list.files(sfd_files)[list.files(sfd_files) %like% ".geojson"]
files <- paste0(sfd_files, files)
  
cfho_cities <- c("fremont", "hayward", "riverside", "san_diego_unincorporated")

compare_vars <- c("evict_rate", "pop_black", "pop_asian", "pop_white",
                  "pop_native_american", "pop_latin_hispanic", 
                  "per_capita_income", "number_rental_units")

calc_sfd_summaries <- function(city){
  
  file <- files[files %like% city]
  df_city <- setDT(st_read(file))
  
  df_city[, evict_rate := (evict_count*100)/number_rental_units]
  df_city <-  melt(df_city, id.vars = c("block_geoid", "cfho_any"), measure.vars = compare_vars)

  df_city <- df_city[, `:=`(mean = mean(.SD$value, na.rm = T),
                             sd   = sd(.SD$value, na.rm = T),
                             n    = .N), 
                             by = c("cfho_any", "variable")]
  
  new_names <- c("Eviction rate per 100 rental units", "Population: % black",
                 "Population: % asian", "Population: % white", "Population: % indigenous",
                 "Population: % Hispanic/Latine", "Median income (2020 dollars)", "Number of rental units")
  
  df_city[, variable := factor(mapvalues(variable, unique(df_city$variable), new_names), levels = new_names)]
  df_city <- unique(df_city[, .(cfho_any, variable, mean, sd, n)])
  
  control_n <- unique(df_city[cfho_any == 0]$n)
  treat_n   <- unique(df_city[cfho_any == 1]$n)
  
  df_city <- dcast(df_city, variable ~ cfho_any, value.var = c("mean", "sd", "n"))
  
  df_city[, mean_diff := mean_1 - mean_0]
  
  df_city[, sd_diff := sqrt(((sd_0^2)/n_0) + ((sd_1^2)/n_1))]
  df_city[, df_diff := welch_satt(sd_0, sd_1, n_0, n_1)]
  
  df_city[, Control := formatted_results(mean_0, sd_0)]
  df_city[, Treated := formatted_results(mean_1, sd_1)]
  df_city[, Difference := formatted_results(mean_diff, sd_diff, ci = T, eff_df = df_diff)]
  
  df_city <- df_city[, c("variable", "Control", "Treated", "Difference"), with = F]
  
  df_city <- rbindlist(list(df_city, list("N", control_n, treat_n, "")), use.names = F)
  
  df_city[, City := city]
  
  return(df_city)
  
}

sfd_summary <- rbindlist(lapply(cfho_cities, calc_sfd_summaries))
sfd_summary <- sfd_summary[, .(City, variable, Control, Treated, Difference)]

sfd_summary[, City := mapvalues(City, unique(sfd_summary$City),
                                c("Fremont", "Hayward", "Riverside", "San Diego (unincorporated cities)"))]

write.csv(sfd_summary, "./cfho_analysis/da/city_summaries_2019.csv", row.names = F)

# Produce example of a SFD path through a map:
plot_map <- function(city){
  
  df_map <- setDT(st_read(paste0("./cfho_analysis/sfd/", city, "_prepped_sfd.geojson")))
  df_map[, evict_count := cut(evict_count, breaks = c(0, 1, 5, 10, 999), 
                              labels = c("0", "1 - 5", "6 - 10", "11+"),
                              right = F, include.lowest = T)]
  
  df_map <- st_as_sf(df_map)
  
  # Construct random points for each CFHO locations within respective polygons:
  sample_points <- function(geo, s){
    as.data.table(st_sample(geo, size = s, type = "random", exact = T))
  }
  
  df_cfho <- copy(df_map[df_map$cfho_num != 0,])
  
  df_cfho <- mcmapply(sample_points, 
                        geo = df_cfho$geometry, 
                        s = df_cfho$cfho_num,
                        SIMPLIFY = F)
  
  df_cfho <- setDT(rbindlist(df_cfho))
  df_cfho <- st_as_sf(df_cfho) %>%
             st_set_crs(., 4326)
  
  city <- ifelse(city == "fremont", "Fremont", ifelse(city == "hayward", "Hayward", ifelse(city == "riverside", "Riverside", "San Diego County")))
  
  map_plot <- ggplot(df_map) + 
              geom_sf(aes(fill = evict_count)) + 
              geom_sf(data = df_cfho, aes(color = cfho_num), size = 2, color = 'black', shape = 17) +
              theme_void() +
              scale_fill_manual(values = c("#ffffff", "#d9d9d9", "#969696", "#525252")) +
              labs(fill = "Number of Evictions",
                   title = "",
                   caption = "Note: The location of CFHP-covered rental units \nwas randomized within block-groups.") +
              guides(fill = guide_legend(nrow = 1, byrow = T),
                     color = guide_legend(nrow = 1)) +
              theme(legend.position = "bottom",
                    legend.direction = "vertical",
                    legend.title.align = 0.5,
                    plot.title = element_text(hjust = 0.5, family = 'sans', size = 16),
                    legend.text = element_text(family = 'sans', size = 12),
                    legend.title = element_text(family = 'sans', size = 14))
  
  return(map_plot)
  
}

calc_avg_units <- function(city){
  
  print(city)
  
  df_map <- st_read(paste0("./cfho_analysis/sfd/", city, "_prepped_sfd.geojson"))
  setDT(df_map)
  print(df_map[cfho_any == 1, mean(cfho_num, na.rm = T)])
  
}

cities     <- c("fremont", 'hayward', 'riverside', "san_diego_unincorporated")
nice_names <- c("Fremont", "Hayward", "Riverside", "San Diego County")

calc_avg_units <- l_ply(cities, calc_avg_units)
  
maps <- lapply(cities, plot_map)

plot_names <- paste0("Number of evictions and CFHP-covered rental units by block-group, in ", nice_names)

plots <- marrangeGrob(grobs = maps, nrow = 1, ncol = 1, 
                      top = quote(grid::textGrob(paste0("Exhibit 2-", g, ". ", plot_names[g]), 
                                                 gp = grid::gpar(fontsize = 24, font = 2, family = "sans"))))

ggsave("./cfho_analysis/sfd/plots/exhibit_2.pdf", plots, 
       device = "pdf", width = 11.69, heigh = 8.27, units = "in", scale = 1.4)

plot_path <- function(city){
  
  mod <- readRDS(sprintf("./cfho_analysis/sfd/main_analysis/cfho_any_unadj.Rds"))[[city]][['sfd_data']]
  df_map <- st_read(paste0("./cfho_analysis/sfd/main_analysis/", city, "_prepped_sfd.geojson"))
  
  setDT(mod)
  mod <- mod[angle == a, .(index, rowno, group)]
  
  setDT(map)
  map[, rowno := .I]
  
  map <- join(map, mod, by = "rowno", type = 'left')
  setorder(map, index)
  
  map[, cfho_any := mapvalues(cfho_any, c(0, 1), c("Control", "One or more CFHO properties"))]

  map <- st_as_sf(map) %>%
         st_transform(crs = 4269)
  
  map$X <- st_coordinates(st_centroid(map))[, "X"]
  map$Y <- st_coordinates(st_centroid(map))[, "Y"]
  
  map_colors <- c('#66c2a5', '#8da0cb')
  
  map$cfho_num[map$cfho_num == 0] <- NA
  map$evict_count[map$evict_count == 0] <- NA
  
  basemap <- ggplot(map) + 
    geom_sf(aes(fill = cfho_num)) + 
    # geom_point(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group, color = "black"),  size = 2.5) + 
    # geom_path(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group, color = "black"),
    #           linewidth = 1) +
    #scale_fill_manual(values = map_colors) +
    theme_void() +
    scale_fill_gradient(na.value = 'transparent', low = "lightgrey", high = "black", breaks = seq(1, 9, 2)) +
    labs(fill = "Number of CFHP-covered rental units") +
    guides(fill = guide_legend(nrow = 1, byrow = T),
           color = "none") +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.title.align = 0.5,
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(family = 'serif', size = 14))
  
  ggsave(basemap, filename = "./cfho_analysis/sfd/plots/riverside_treat.svg", device = "svg", 
         bg = "transparent", width = 11.7, height = 8.3, units = "in")
  
  basemap2 <- ggplot(map) + 
    geom_sf(aes(fill = evict_count)) + 
    # geom_point(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group, color = "black"),  size = 2.5) + 
    # geom_path(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group, color = "black"),
    #           linewidth = 1) +
    #scale_fill_manual(values = map_colors) +
    theme_void() +
    scale_fill_continuous(na.value = 'transparent', low = "lightgrey", high = "black", breaks = seq(1, 14, 2)) +
    labs(fill = "Number of evictions") +
    guides(fill = guide_legend(nrow = 1, byrow = T),
           color = "none") +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.title.align = 0.5,
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(family = 'serif', size = 14))
  
  ggsave(basemap2, filename = "./cfho_analysis/sfd/plots/riverside_out.svg", device = "svg", 
         bg = "transparent", width = 11.7, height = 8.3, units = "in")
  
  path_1 <- ggplot(map) + 
          geom_sf(aes(fill = cfho_num), alpha = 0.4) + 
          geom_point(data = map[map$group >= 1 & map$group <= 1, ], aes(x = X, y = Y, group = group), color = "black",  size = 2.5) +
          geom_path(data = map[map$group >= 1 & map$group <= 1, ], aes(x = X, y = Y, group = group), color = "black",
                    linewidth = 1) +
          scale_fill_gradient(na.value = 'transparent', low = "lightgrey", high = "black", breaks = seq(1, 9, 2)) +
          theme_void() +
          labs(fill = "Number of CFHP-covered rental units") +
          guides(fill = guide_legend(nrow = 1, byrow = T),
                 color = "none") +
          theme(legend.position = "bottom",
                legend.direction = "vertical",
                legend.title.align = 0.5,
                plot.title = element_text(hjust = 0.5),
                legend.text = element_text(family = 'serif', size = 14))
  
  ggsave(path_1, filename = "./cfho_analysis/sfd/plots/path1_90.svg", device = "svg", 
         bg = "transparent", width = 11.7, height = 8.3, units = "in")
  
  path_2 <- ggplot(map) + 
            geom_sf(aes(fill = cfho_num), alpha = 0.4) + 
            geom_point(data = map[map$group >= 1 & map$group <= 2, ], aes(x = X, y = Y, group = group), color = "black",  size = 2.5) +
            geom_path(data = map[map$group >= 1 & map$group <= 2, ], aes(x = X, y = Y, group = group), color = "black",
                      linewidth = 1) +
            scale_fill_gradient(na.value = 'transparent', low = "lightgrey", high = "black", breaks = seq(1, 9, 2)) +
            theme_void() +
            labs(fill = "Number of CFHP-covered rental units") +
            guides(fill = guide_legend(nrow = 1, byrow = T),
                   color = "none") +
            theme(legend.position = "bottom",
                  legend.direction = "vertical",
                  legend.title.align = 0.5,
                  plot.title = element_text(hjust = 0.5),
                  legend.text = element_text(family = 'serif', size = 14))
  
  ggsave(path_2, filename = "./cfho_analysis/sfd/plots/path2_90.svg", device = "svg", 
         bg = "transparent", width = 11.7, height = 8.3, units = "in")
  
  path_3 <- ggplot(map) + 
            geom_sf(aes(fill = cfho_num), alpha = 0.4) + 
            geom_point(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group), color = "black",  size = 2.5) +
            geom_path(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group), color = "black",
                      linewidth = 1) +
            scale_fill_gradient(na.value = 'transparent', low = "lightgrey", high = "black", breaks = seq(1, 9, 2)) +
            theme_void() +
            labs(fill = "Number of CFHP-covered rental units") +
            guides(fill = guide_legend(nrow = 1, byrow = T),
                   color = "none") +
            theme(legend.position = "bottom",
                  legend.direction = "vertical",
                  legend.title.align = 0.5,
                  plot.title = element_text(hjust = 0.5),
                  legend.text = element_text(family = 'serif', size = 14))
  
  ggsave(path_3, filename = "./cfho_analysis/sfd/plots/path3_90.svg", device = "svg", 
         bg = "transparent", width = 11.7, height = 8.3, units = "in")
  
  path <- ggplot(map) + 
            geom_sf(aes(fill = cfho_num), alpha = 0.4) + 
            geom_point(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group), color = "black",  size = 2.5) +
            geom_path(data = map[map$group >= 1 & map$group <= 3, ], aes(x = X, y = Y, group = group), color = "black",
                      linewidth = 1) +
            scale_fill_gradient(na.value = 'transparent', low = "lightgrey", high = "black", breaks = seq(1, 9, 2)) +
            theme_void() +
            labs(title = paste0("First three paths taken by SFD algorithm in Riverside, for angle = ", a),
                  fill = "Treatment group") +
            guides(fill = guide_legend(nrow = 2, byrow = T),
                   color = "none") +
            theme(legend.position = "bottom",
                  legend.direction = "vertical",
                  legend.title.align = 0.5,
                  plot.title = element_text(hjust = 0.5))
  
  return(path)
  
}

path_0  <- plot_path(df_map, model, 0)
path_90 <- plot_path(df_map, model, 90)

ggsave(path_0, filename = "./cfho_analysis/sfd/plots/sfd_path_riverside_0.svg", device = "svg", 
       bg = "transparent", width = 11.7, height = 8.3, units = "in")

ggsave(path_90, filename = "./cfho_analysis/sfd/plots/sfd_path_riverside_90.svg", device = "svg", 
       bg = "transparent", width = 11.7, height = 8.3, units = "in")


####################################
# Summary stats: All units in 2019 #
####################################

df_ever <- copy(df_analysis[year %in% c(2019) & !is.na(ever_treated)])

df_ever <- df_ever[, .(location, year, ever_treated, total_population, evict_rate,
                       total_crime_rate_10k, assault_rate_10k, burglary_rate_10k,
                       pop_black, pop_asian, pop_white, 
                       pop_latin_hispanic, pop_native_american,
                       percent_poverty_150, avg_edu, median_rent_adj, median_income_adj,
                       number_rental_units, pop_tenants, pop_immigrant, 
                       inquiry_rate_1000_pop, inquiry_count_race, inquiry_count_national_origin)]

# Convert number of rental units into a rate:
df_ever[, number_rental_units := (number_rental_units*100)/total_population]
df_ever[, inquiry_rate_1000_pop := (inquiry_rate_1000_pop*10)]
df_ever[, inquiry_count_race := (inquiry_count_race*10000)/total_population]
df_ever[, inquiry_count_national_origin := (inquiry_count_national_origin*10000)/total_population]

df_ever <- melt(df_ever, id.vars = c("location", "year", "ever_treated"))
df_ever[, ever_treated := as.factor(ever_treated)]

df_ever[, `:=`(mean = mean(.SD$value, na.rm = T),
              sd   = sd(.SD$value, na.rm = T),
              n    = .N), by = c("year", "ever_treated", "variable")]

df_ever <- unique(df_ever[, .(ever_treated, year, variable, mean, sd, n)])
df_ever <- dcast(df_ever, variable + year ~ ever_treated, value.var = c("mean", "sd", "n"))

df_ever[, mean_diff := mean_1 - mean_0]

# Calculate pooled standard deviation, based off two-sample unpaired t-test
# (assuming here that sample is normally distributed, which is mostly true
# across these variables and that samples are independent, which is likely 
# unreasonable but not something easily corrected. Reviewer wanted this).
df_ever[, sd_diff := sqrt(((sd_0^2)/n_0) + ((sd_1^2)/n_1))]
df_ever[, df_diff := welch_satt(sd_0, sd_1, n_0, n_1)]

df_ever[, Control := formatted_results(mean_0, sd_0)]
df_ever[, Treated := formatted_results(mean_1, sd_1)]
df_ever[, Difference := formatted_results(mean_diff, sd_diff, ci = T, eff_df = df_diff)]

df_ever <- df_ever[, c("variable", "Control", "Treated", "Difference"), with = F]

# Recode variable names for the summary plot
new_names <- c("Total population", "Evict rate (per 100 rental units)", 
               "Total crime rate (per 10k people)", "Assault rate (per 10k people)", "Burglary rate (per 10k people)",
               "Population proportion: Black", "Population proportion: Asian", "Population proportion: White",
               "Population proportion: Hispanic", "Population proportion: indigenous", 
               "Population proportion: <150% of federal poverty line",  
               "Average years of education", "Median rent price (in 2020 dollars)",
               "Median income (in 2020 dollars)", "Number of rental units (per 100 people)",
               "Population proportion: renters", "Population proportion: immigrant",
               "FHEO inquiry rate per 10,000 people (total), 2013-2020", 
               "FHEO inquiry rate per 10,000 people  (race discrimination), 2013-2020",
               "FHEO inquiry rate per 10,000 people  (national origin discrimination), 2013-2020")

df_ever[, variable := factor(mapvalues(variable, unique(df_ever$variable), new_names), levels = new_names)]

# Save as nicely formatted table as well:
df_ever <- rbindlist(list(df_ever, list("N", 
                                        dim(df_analysis[year == 2019 & ever_treated == 1,])[1], 
                                        dim(df_analysis[year == 2019 & ever_treated == 0,])[1], "")), 
                     use.names = F)

write.csv(df_ever, "./cfho_analysis/da/ever_treated_summary_stats_2019.csv", row.names = F)

# Look at subset of data and plot raw data + CI for each group:
df_z <- copy(df_analysis[year == 2019 & !is.na(ever_treated), .(location, ever_treated, pop_black, pop_asian,
                         pop_latin_hispanic, pop_white, total_population)])

df_z <- melt(df_z, 
             id.vars = c("location", "ever_treated", "total_population"), 
             value.name = "prop",
             variable.names = "var")

df_z[, ever_treated := factor(ever_treated)]
df_z[, variable := factor(variable, levels = rev(c("pop_black", "pop_asian", "pop_latin_hispanic", 
                                                   "pop_white")))]

#  Calculate notches to compare relative differences in medians. 
# (using guidelines in McGill, Tukey, and Larsen, 1978):
# One could also run a two-sided Wilcox rank test to determine if the medians are 
# significantly different (I checked: results are the same as the visual).

#e.g.
#wilcox.test(df_z[ever_treated == 1 & variable == "pop_black"]$prop,
#            df_z[ever_treated == 0 & variable == "pop_black"]$prop,
#            two.sided = T)

# Then calculate median and lower/upper notches:
df_z[, mid := mean(.SD$prop), by = c("variable", "ever_treated")]
df_z[, lower := t.test(.SD$prop)$conf.int[1], by = c("variable", "ever_treated")]
df_z[, upper := t.test(.SD$prop)$conf.int[2], by = c("variable", "ever_treated")]

plot_colors = c("#ba9c7d", "#807dba")

# Remap variables:
df_z[, ever_treated := ifelse(ever_treated == 1, "Yes", "No")]

new_names <- c("Black", "Asian", 
                   "Hispanic", "White")

df_z[, variable := factor(mapvalues(variable, unique(df_z$variable), new_names), levels = new_names)]

pop_prop <- ggplot(df_z, aes(x = prop, y = ever_treated, color = ever_treated)) +
              facet_wrap(~variable,  nrow = 4, strip.position = 'right') +
              geom_point(alpha = 0.3, position = position_jitterdodge()) +
              geom_crossbar(aes(xmin = lower, x = mid, xmax = upper)) +
              theme_bw() +
              labs(x = "Proportion of total population",
                   y = "CFHP \npolicy") +
              scale_color_manual(values = plot_colors) +
              scale_x_continuous(breaks = seq(0, 1, 0.2), labels = percent) +
              theme(strip.text.y.right = element_text(size = 14, family = 'serif', angle = 0),
                    strip.background = element_blank(),
                    legend.position = 'none',
                    axis.ticks = element_line(linewidth = 1),
                    axis.ticks.length = unit(5.6, "points"),
                    axis.text = element_text(family = "serif", size = 12),
                    axis.title = element_text(family = "serif", size = 18),
                    axis.title.y = element_text(angle = 0, vjust = 1),
                    panel.spacing = unit(0, "lines"))

ggsave(pop_prop, filename = "./cfho_analysis/figs_collected/figure_2_alt.svg", 
       device = "svg", bg = "transparent", width = 11.7, height = 8.3, units = "in")

plot_colors <- c("#4d9221", "#8856a7")

ever_plot <- ggplot(df_ever, aes(x = q50, y = ever_treated, color = ever_treated)) +
                geom_linerange(aes(xmin = q025, xmax = q975), lwd = 1, position = position_dodgev(height = 0.4)) +
                geom_point(size = 3, position = position_dodgev(height = 0.4)) +
                facet_wrap(~variable, scale = "free", ncol = 3) + 
                theme_bw() +
                labs(title = "Summary statistics for all cities, in 2019",
                     y = "Ever treated") +
                scale_color_manual(values = plot_colors) +
                theme(plot.title = element_text(hjust = 0.5),
                      strip.background = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      legend.position = "bottom")

ggsave(ever_dot, filename = "./cfho_analysis/da/figure_2.svg", device = "svg", bg = "transparent", width = 11.7, height = 8.3, units = "in")

###############################################
# Maps of control & treated sites by analysis #
###############################################

cal_state <- get_acs(geography = "state", year = 2020, variables = "B01003_001", 
                     state = "CA", geometry = T, key = census_key, survey = "acs5") %>%
                    setDT %>%
                    .[, .(geometry)] %>%
                    st_as_sf()

cal_county <- get_acs(geography = "county", year = 2020, variables = "B01003_001", 
                      state = "CA", geometry = T, key = census_key, survey = "acs5") %>%
                      setDT %>%
                      .[, .(geometry)] %>%
                      st_as_sf()

cal_place <- get_acs(geography = "place", year = 2020, variables = "B01003_001", 
                    state = "CA", geometry = T, key = census_key, survey = "acs5") %>%
                    setDT %>%
                    .[, .(GEOID, geometry)] %>%
                    setnames(., "GEOID", "geoid") %>%
                    .[, geoid := as.numeric(geoid)]

df_map <- copy(df_analysis)
df_map <- df_map[!is.na(ever_treated)]

sfd_sites <- c("fremont", "riverside", "hayward", "vista", 'san marcos', 'cardiff',
               'solana beach', 'del mar', 'encinitas', 'poway', 'san diego',
               'santee', 'lakeside', 'el cajon', 'la mesa', 'lemon grove',
               'fallbrook', 'ramona')

df_map[, treat_status := ifelse(ever_treated == 0, "Control site (DiD models)", "CFHO policy (excluded)")]
df_map[ever_treated == 1 & year_implemented >= 2011, treat_status := "CFHO policy (DiD models)"]
df_map[ever_treated == 1 & location %in% sfd_sites, treat_status := "CFHO policy (SFD models)"]

df_map <- join(df_map, cal_place, by = "geoid", type = "right")
df_map <- df_map[year == 2019, ]

df_map[, treat_status := factor(treat_status, levels = c("Control site (DiD models)", "CFHO policy (excluded)", 
                                                         "CFHO policy (DiD models)", "CFHO policy (SFD models)"))]

setorder(df_map, treat_status)

df_map <- st_as_sf(df_map)

# Standardize city sizes on the map by getting the centroid of each city, then
# drawing an equal-length circumference around each respective city:
df_map <- st_transform(df_map, 3488) %>%
          st_centroid(df_map) %>%
          st_buffer(dist = set_units(5, "km")) %>%
          st_transform(st_crs(cal_state))

treat_color <- c('#fec44f','#74c476','#2b8cbe','#807dba')

treat_by_model  <- ggplot() + 
                    geom_sf(data = cal_state, fill = "#f7f7f7") + 
                    geom_sf(data = cal_county, fill = "transparent", alpha = 0.6, color = "#bdbdbd") + 
                    geom_sf(data = df_map[df_map$treat_status == "Control site (DiD models)",], aes(fill = "Control site (DiD models)")) +
                    geom_sf(data = df_map[df_map$treat_status == "CFHO policy (excluded)",], aes(fill = "CFHO policy (excluded)")) +
                    geom_sf(data = df_map[df_map$treat_status == "CFHO policy (DiD models)",], aes(fill = "CFHO policy (DiD models)")) +
                    geom_sf(data = df_map[df_map$treat_status == "CFHO policy (DiD models & SFD models)",], aes(fill = "CFHO policy (DiD models & SFD models)")) +
                    geom_sf(data = df_map[df_map$treat_status == "CFHO policy (SFD models)",], aes(fill = "CFHO policy (SFD models)")) +
                    scale_color_manual(values = treat_color) +
                    scale_fill_manual(values = treat_color,
                                      breaks = unique(df_map$treat_status)) +
                    labs(fill = "Treatment status") +
                    theme_void() +
                    guides(fill = guide_legend(nrow = 2, byrow = T)) +
                    theme(legend.position = "bottom",
                          legend.direction = "vertical",
                          legend.title.align = 0.5)

ggsave(treat_by_model, filename = "./cfho_analysis/treated_locs_by_model.svg", device = "svg", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

treat_locs  <- ggplot() + 
                    geom_sf(data = cal_state, fill = "#bdbdbd") + 
                    geom_sf(data = cal_county, fill = "transparent", alpha = 0.6, color = "#bdbdbd") + 
                    geom_sf(data = df_map[df_map$treat_status != "Control site (DiD models)",], fill = "#807dba") +
                    scale_color_manual(values = treat_color) +
                    labs(fill = "Cities with CFHP policy") +
                    theme_void() +
                    guides(fill = guide_legend(nrow = 2, byrow = T)) +
                    theme(legend.position = "bottom",
                          legend.direction = "vertical",
                          legend.title.align = 0.5)

###############
# Time trends #
###############

# Add stacked dotplot histogram above map to display years of implementation 
# and a line graph to display population exposed by year:

df_dot  <- unique(df_analysis[ever_treated == 1, .(location, year_implemented)])

plot_histdot <- ggplot(df_dot, aes(x = year_implemented)) +
                geom_dotplot(binwidth = 1, color = "#807dba", fill = "#807dba", method = 'histodot') +
                annotate('text', x = 1998, y = 0.8,
                         label = paste("Implemented CFHP \nbefore 2010: ", 
                                       length(unique(df_dot[is.na(year_implemented), location])),
                                       "cities")) +
                scale_x_continuous(limits = c(1995, 2019), 
                                   breaks = c(seq(1995, 2015, 5), 2019)) +
                scale_y_continuous(NULL, breaks = NULL) +
                coord_fixed(ratio = 10) +
                labs(x = "Implementation year") +
                theme_classic() +
                theme(axis.line.y = element_blank(),
                      panel.grid.minor.y = element_blank()) 

# Calculate population of renters living in cities w/ CFHO

# For locations where we know policy is implemented before to 2010 (but no exact date)
# set to 2009 for dot plot purpose:
df_year <- unique(df_analysis[year <= 2019, 
                              .(year, ever_treated, location, year_implemented, pop_tenants, total_population,
                                renter_white, renter_asian, renter_black,
                                renter_hispanic_latine, renter_white_alone,
                                renter_native_american, renter_multiple_race,
                                renter_pacific_islander, renter_other)])

# Calculate percentage of total renter population by race. Note: since ACS only reports race by head of household
# for rental units, I'm assuming that renter population will stratify similarly to renter household
# population.

df_year[, renter_total := .SD$total_population*.SD$pop_tenants, by = "year"]
df_year[, state_total := sum(.SD$total_population*.SD$pop_tenants), by = "year"]

df_year <- df_year[ever_treated == 1,]
df_year <- df_year[(year >= year_implemented)|(is.na(year_implemented)),]

df_year[, cfhp_renter_total := sum(.SD$renter_total), by = "year"]

df_year[, renter_hispanic_latine_prop := sum(.SD$renter_hispanic_latine*.SD$renter_total)/state_total, by = "year"]
df_year[, renter_black_prop := sum(.SD$renter_black*.SD$renter_total)/state_total, by = "year"]
df_year[, renter_white_alone_prop := sum(.SD$renter_white_alone*.SD$renter_total)/state_total, by = "year"]
df_year[, renter_asian_prop := sum(.SD$renter_asian*.SD$renter_total)/state_total, by = "year"]
df_year[, renter_other_prop := sum((.SD$renter_pacific_islander + .SD$renter_native_american + .SD$renter_multiple_race + .SD$renter_other)*.SD$renter_total)/state_total, by = "year"]

df_year <- unique(df_year[, .(year, renter_total, renter_hispanic_latine_prop, renter_black_prop, renter_white_alone_prop, renter_asian_prop, renter_other_prop)])
df_year <- melt(df_year, id.vars = "year", variable.name = "category", value.name = "Population")

new_names <- c("Total", "Hispanic", "Black", "White (alone)", "Asian", "Other")
df_year[, category := mapvalues(category, unique(df_year$category), new_names) ]

df_year <- unique(df_year)

plot_colors <- c("#637079", "#1473af", "#ca7ba7", "#199e74", "#d1602f", "#e59e37")

plot_pop <- ggplot(df_year[category != "Total",], aes(x = year, y = Population, color = category)) +
                  geom_line(linewidth = 1.2) +
                  geom_point(size = 3) +
                  #scale_y_continuous(limits = c(0.2, 0.4), labels = scales::percent, breaks = seq(0.2, 0.4, 0.05)) +
                  scale_x_continuous(breaks = seq(2010, 2019, 2)) +
                  scale_color_manual(values = plot_colors) +
                  labs(title = "Proportion of total renters in a demographic group that live in CFHP cities",
                      y = "Proportion of \ntotal renters in \ndemographic \ngroup",
                      x = "Year",
                      color = "Demographic group") +
                  theme_bw() +
                  theme(legend.position = "bottom",
                        legend.text = element_text(family = 'serif', size = 12),
                        axis.ticks = element_line(linewidth = 1),
                        axis.ticks.length = unit(5.6, "points"),
                        axis.text = element_text(family = "serif", size = 12),
                        axis.title = element_text(family = "serif", size = 18),
                        axis.title.y = element_text(angle = 0, vjust = 0.9),
                        plot.title = element_text(size = 22, family = 'serif'))
                
plot_pop_simple <- ggplot(df_year, aes(x = year, y = Population)) +
                    geom_line(linewidth = 1.2, color = "#807dba") +
                    geom_point(size = 3, color = "#807dba") +
                    scale_y_continuous(labels = label_number(scale = 1e-6), 
                                       breaks = seq(2e6, 6e6, 1e6),
                                       limits = c(2e6, 5e6)) +
                    scale_x_continuous(breaks = seq(2010, 2019, 2)) +
                    labs(y = "Population \n(in millions)",
                         x = "Year") +
                    theme_bw() +
                    theme(legend.position = "none",
                          axis.ticks = element_line(linewidth = 1),
                          axis.ticks.length = unit(5.6, "points"),
                          axis.text = element_text(family = "serif", size = 12),
                          axis.title = element_text(family = "serif", size = 18),
                          axis.title.y = element_text(angle = 0, vjust = 0.9))

ggsave(treat_locs, filename = "./cfho_analysis/da/figure_1a.svg", device = "svg", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(plot_histdot, filename = "./cfho_analysis/da/figure_1b.svg", device = "svg", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(plot_pop, filename = "./cfho_analysis/da/figure_1c.svg", device = "svg", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(plot_pop_simple, filename = "./cfho_analysis/da/figure_1c_alt.svg", device = "svg", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(treat_locs, filename = "./cfho_analysis/da/figure_1a.png", device = "png", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(plot_histdot, filename = "./cfho_analysis/da/figure_1b.png", device = "png", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(plot_pop, filename = "./cfho_analysis/da/figure_1c.png", device = "png", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")

ggsave(plot_pop_simple, filename = "./cfho_analysis/da/figure_1c_alt.png", device = "png", 
       bg = "transparent", height = 11.7, width = 11.7, units = "in")
