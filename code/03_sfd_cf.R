# Summarize estimates  and simulate counterfactual policy results
# Max Griswold
# 3/10/23

library(data.table)
library(sf)
library(plyr)

library(ggplot2)
library(ggridges)
library(gridExtra)
library(scales)
library(ggstance)

library(stringr)
library(MASS)

# Get sfd analysis dataset and results:

sfd_path <- "./cfho_analysis/sfd/"

#analysis_type <- "main_analysis/"
#analysis_type <- "sensitivity_test_neighbors/"

analysis_type <- "main_analysis/"

setwd(sfd_path)

sfd_datasets <- list.files()[list.files() %like% ".geojson"]

sfd_models <- list.files(paste0(sfd_path, analysis_type))
sfd_models <- sfd_models[sfd_models %like% ".rds" & !(sfd_models %like% "no_renter") & !(sfd_models %like% "interact")]

# Read in cities and combine into single file:
combine_sfd_df <- function(file){
  
  city_name <- gsub("_prepped_sfd.geojson", "", file)
  
  d <- st_read(file)
  
  setDT(d)
  d[, city_name := city_name]
  
  return(d)
  
}

df_sfd <- setDT(ldply(sfd_datasets, combine_sfd_df))


df_sfd <- df_sfd[, .(city_name, cfho_num, cfho_any, evict_count, number_rental_units_100,
                     pop_white, pop_black, pop_asian, pop_native_american, pop_latin_hispanic,
                     per_capita_income_10k)]


# Read in a model for a given city. Using estimated betas and vcov,
# simulate 1000 parameter draws from each model and combine into one data.frame

sim_betas <- function(a, mod, city){
  
  submod <- mod[[a]]
  
  vars  <- submod[['coef']]$var
  
  k <- length(vars)
  
  ests  <- submod[['coef']]$Estimate
  vcovs <- submod[['vcov']][1:k, 1:k]
  
  sims <- as.data.table(mvrnorm(1000, ests, vcovs))
  
  # Second column is always the treatment variable
  treatname <- vars[2]
  
  if (length(vars) > 2){
    setnames(sims, names(sims), c(vars[1], "est_beta", vars[3:length(vars)]))
  }else{
    setnames(sims, names(sims), c(vars[1], "est_beta"))
    
  }
  sims[, treat_name := treatname]
  sims[, angle := a]
  
  return(sims)
  
}

combine_sims <- function(city, mod){
  
  mod <- mod[[city]]
  
  angles <- names(mod)[names(mod) %like% 'angle']
  
  sims <- rbindlist(lapply(angles, sim_betas, mod = mod, city = city))
  
  sims[, angle := gsub("angle=", "", angle)]
  sims[, angle := factor(angle, levels = unique(angle))]
  
  sims[, city_name := city]
  
  return(sims)
  
}

sim_wrapper <- function(model){
  
  mod <- readRDS(paste0(analysis_type, model))
  cities <- names(mod)
  
  sim_cities <- rbindlist(lapply(cities, combine_sims, mod = mod))
  
  # Make model names human-readable
  model_name <- gsub(".Rds|.rds", "", model)
  model_name <- gsub("_", " ", model_name)
  model_name <- ifelse(model_name %like% "any", 
                        gsub("cfho any", "CFHO any,", model_name),
                        gsub("cfho count", "CFHO count,", model_name))
  
  model_name <- ifelse(model_name %like% "adj", 
                        gsub("adj", "adjusted", model_name),
                        gsub("unadj", "unadjusted", model_name))
  
  sim_cities[, model_name := model_name]
  
  return(sim_cities)
  
}

sim_results <- rbindlist(lapply(sfd_models, sim_wrapper), fill = T)

plot_sim_betas <- function(cityname, modelname, df_sim){
  
  df_sim <- df_sim[model_name == modelname & city_name == cityname,]
  
  # Second column is always the treatment variable
  estname <- unique(df_sim$treat_name)
  
  df_overall <- data.table("angle" = "Overall \nEffect", 
                           "est_beta" = quantile(df_sim$est_beta, 0.5),
                           "est_lower" = quantile(df_sim$est_beta, 0.025),
                           "est_upper" = quantile(df_sim$est_beta, 0.975))

  estname_human <- ifelse(estname == "cfho_any",
                          "One or More CFHP units",
                          "Number of CFHP units")
  
  cityname_human <- ifelse(cityname == "san_diego_unincorporated",
                           "San Diego County",
                           str_to_title(cityname))
  
  if (estname == "cfho_any"){bks <- seq(-0.25, 2.5, 0.25)}else{bks <- seq(-0.25, 1.5, 0.25)}
  if (estname == "cfho_any"){lms <- c(-0.25, 2.5)}else{lms <- c(-0.25, 1.5)}

  # Make sure factor variables are ordered the same way in each plot df:
  df_sim <- df_sim[, .(angle, est_beta)]
  df_sim <- rbindlist(list(df_sim, df_overall), fill = T)
  
  df_sim[, angle := as.character(angle)]
  df_sim[, angle := factor(angle, levels = rev(unique(df_sim$angle)))]
  
  print(modelname)
  
  p <- ggplot(df_sim, aes(x = est_beta, y = angle)) +

        geom_density_ridges(scale = 1.7, alpha = 0.9, rel_min_height=0.01) +
        geom_point(data = df_sim[angle == "Overall \nEffect"], size = 2.5, position = position_nudge(y = 0.4)) +
        geom_errorbarh(data = df_sim[angle == "Overall \nEffect"], aes(xmin = est_lower, xmax = est_upper), 
                       size = 1, height = 0.3, position = position_nudge(y = 0.4)) +
        geom_vline(xintercept = 0, color = 'grey', linetype = 2, size = 1) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0),
                           breaks = bks,
                           limits = lms) +
        labs(x = ifelse(cityname_human == "San Diego County", paste0("\n", estname_human), ""),
             y = ifelse(cityname_human == "Fremont", "Map \nRotation", ""),
             title = paste0(cityname_human, "\n")) +
        theme_ridges() +
        theme(legend.position = "none",
              axis.title.y = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}

extract_gof <- function(model, cityname){
  
  mod <- readRDS(paste0(analysis_type, model))
  
  mod_city <- mod[[cityname]]
  angles <- names(mod_city)[names(mod_city) %like% 'angle']
  
  res <- list()
  for (a in angles){
    
    df_hold <- mod_city[[a]]$sum
    df_hold[, angle := a]
    
    res[[a]] <- df_hold
    
  }
  
  res <- rbindlist(res)
  res[, city_name := cityname]
  
  return(res)
  
}

gof_args <- expand.grid(model = as.character(sfd_models),
                        cityname = unique(sim_results$city_name))

gof <- rbindlist(mapply(extract_gof, 
                        model = gof_args$model,
                        cityname = unique(sim_results$city_name),
                        SIMPLIFY = F))

# Order arguments so that 2x2 facets display ordered models, with
# each facet displaying a city
args <- expand.grid(cityname = unique(sim_results$city_name),
                    modelname = unique(sim_results$model_name)[c(2, 1, 4, 3)])

sim_beta_density_plots <- mapply(plot_sim_betas, 
                                 cityname = args$cityname,
                                 modelname = args$modelname, 
                                 SIMPLIFY = F,
                                 MoreArgs = list(df_sim = sim_results))

# Set names of each plot:
plot_names <- paste0("Estimated ATT effect by map rotation, ", 
                     c("for unadjusted model using binary treatment",
                       "for adjusted model using binary treatment",
                       "for unadjusted model using continuous treatment",
                       "for adjusted model using continuous treatment"))

plots <- marrangeGrob(grobs = sim_beta_density_plots, nrow = 2, ncol = 2, 
                      top = quote(grid::textGrob(paste0("Exhibit B-", g, ". ", plot_names[g]), gp = grid::gpar(fontsize = 24, font = 2, family = "sans"))))

ggsave(paste0(analysis_type, "plots/sim_beta_density_plots.pdf"), plots,
       device = "pdf", width = 11.69, heigh = 8.27, units = "in", scale = 1.5)

# pdf(paste0(analysis_type, "plots/sim_beta_density_plots.pdf"), paper = 'A4r')
# par(mar = c(0, 0, 0, 0))
# for (p in sim_beta_density_plots){
#   print(p)
# }
# dev.off()

# Plot simulated coefficients as ATT plot, like DiD models:
sim_att <- sim_results[, .(est_beta, treat_name, city_name, model_name)]

# Calculate lower ui, mean, and upper ui estimates for each model (adjust/unadjust),
# treatment-parameter (any v. count), and city:

sim_att[, est_mean  := quantile(.SD$est_beta, probs = 0.5), by = c("model_name", "city_name")]
sim_att[, est_lower := quantile(.SD$est_beta, probs = 0.025), by = c("model_name", "city_name")]
sim_att[, est_upper := quantile(.SD$est_beta, probs = 0.975), by = c("model_name", "city_name")]

# Also calculate across all sites:
sim_att_overall <- sim_results[, .(est_beta, treat_name, city_name, model_name)]
sim_att_overall[, city_name := "Across Sites"]

sim_att_overall[, est_mean  := quantile(.SD$est_beta, probs = 0.5), by = c("model_name")]
sim_att_overall[, est_lower := quantile(.SD$est_beta, probs = 0.025), by = c("model_name")]
sim_att_overall[, est_upper := quantile(.SD$est_beta, probs = 0.975), by = c("model_name")]

sim_att <- unique(sim_att[, .(city_name, treat_name, model_name, est_mean, est_lower, est_upper)])
sim_att_overall <- unique(sim_att_overall[, .(city_name, treat_name, model_name, est_mean, est_lower, est_upper)])

sim_att <- rbind(sim_att, sim_att_overall)

# Rename variables to make plot more interpretable for readers:
sim_att[, city_name := ifelse(city_name == "san_diego_unincorporated",
                              "San Diego county",
                              str_to_title(city_name))]

# Set order of models to be more intuitive for the reader (treatment-parameter type first, then
# model type):

if (sim_att$model_name[1] %like% "no renter"){
  old_model_names <- c("CFHO any, unadjusted no renter", "CFHO any, adjusted no renter", "CFHO count, unadjusted no renter",
                       "CFHO count, adjusted no renter")
  
  ord_mods <- c("Unadjusted model, \nbinary treatment", "Adjusted model, \nbinary treatment", 
                "Unadjusted model, \ncontinuous treatment", "Adjusted model, \ncontinuous treatment")
}else{
  old_model_names <- c("CFHO any, unadjusted", "CFHO any, adjusted", "CFHO count, unadjusted",
                       "CFHO count, adjusted")
  
  ord_mods <- c("Unadjusted model, \nbinary treatment", "Adjusted model, \nbinary treatment", 
                "Unadjusted model, \ncontinuous treatment", "Adjusted model, \ncontinuous treatment")
}

sim_att[, model_name := mapvalues(model_name, old_model_names, ord_mods)]
sim_att[, model_name := factor(model_name, levels = ord_mods)]

plot_colors    <- c("#4d9221", "#b8e186", "#2166ac", "#92c5de")
plot_colors_bw <- c("#d9d9d9", "#969696", "#737373", "#252525")

plot_att <- ggplot(sim_att, aes(x = est_mean, color = model_name, y = city_name)) +
            geom_errorbarh(aes(xmin = est_lower, xmax = est_upper, color = factor(model_name, levels = rev(ord_mods))), size = 2, height = 0, position = position_dodgev(height = 0.4)) +
            geom_point(aes(color = factor(model_name, levels = rev(ord_mods))), size = 3, position = position_dodgev(height = 0.4)) +
            geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
            labs(title = "Exhibit D. Estimated ATT by model and location",
                 x = "Estimated Treatment Effect",
                 y = "Location",
                 color = "Model") +
            theme_bw() +
            scale_color_manual(values = plot_colors_bw, limits = ord_mods) +
            theme(plot.title = element_text(hjust = 0.5, family = 'sans', size = 16),
                  strip.background = element_blank(),
                  legend.position = "bottom",
                  legend.text = element_text(family = 'sans', size = 10),
                  axis.ticks = element_line(linewidth = 1),
                  axis.ticks.length = unit(5.6, "points"),
                  axis.title = element_text(size = 12, family = 'sans'),
                  axis.title.y = element_text(size = 12, family = 'sans', angle = 0),
                  axis.text = element_text(size = 10, family = 'sans'),
                  axis.text.x = element_text(size = 10, family = 'sans',
                                             margin = margin(t = 5, r = 0, b = 10, l = 0)),
                  legend.title = element_text(family = 'sans', size = 14))

# Add on formatted rows for table in SFD paper:
sim_att[, format_text := paste0(round(est_mean, 3), " (", round(est_lower, 3), ", ", round(est_upper, 3), ")")]
setorder(sim_att, city_name, treat_name, model_name)

write.csv(sim_att, "formatted_att_res.csv", row.names = F)

ggsave(filename = "./plots/exhibit_d_no_renter.svg", plot = plot_att, device = "svg", height = 8.3, width = 11.7, units = "in")
ggsave(filename = "plots/summary_att.png", plot = plot_att, device = "png", height = 8.3, width = 11.7, units = "in")

# Using simulated coefficients, calculate change in predicted evictions when
# CFHO policies are removed:

sim_pred <- function(cityname, modelname){
  
  obs  <- df_sfd[city_name == cityname,]
  sims <- sim_results[city_name == cityname & model_name == modelname]
  
  # Only hold onto treated locations:
  obs <- obs[cfho_any == 1,]
  
  truth <- sum(obs$evict_count)
  
  treat_var <- unique(sims$treat_name)
  
  # Remove extraneous columns from sims:
  sims <- Filter(function(x)!all(is.na(x)), sims)
  sims[, `:=`(treat_name = NULL,
              angle = NULL,
              city_name = NULL,
              model_name = NULL)]
  
  setnames(sims, "est_beta", treat_var)
  
  # Add intercept column to obs and only hold onto columns corresponding to
  # estimated betas:
  
  obs[, intercept := 1]
  obs <- obs[, names(sims), with = F]
  
  # Calculate estimate from simulated betas using observed values & a
  # counterfactual where no observed locations are treated w/ CFHO:
  obs_cf <- copy(obs)
  obs_cf[, (treat_var) := 0]
  
  pred    <- colSums(as.matrix(obs) %*% t(as.matrix(sims)))
  pred_cf <- colSums(as.matrix(obs_cf) %*% t(as.matrix(sims)))
  
  df_cf <- data.table("obs" = truth, "status_quo" = pred, "cf" = pred_cf)
  df_cf[, percent_change := round((status_quo - cf)/status_quo, 2)]
  
  df_cf[, city_name := cityname]
  df_cf[, model_name := modelname]
    
  return(df_cf)
  
}

# Alternatively, only use the beta coefficient for treatment and estimate
# relative change in observed evictions, if that effect held:

sim_cf <- function(modelname, cityname){
  
  obs  <- df_sfd[city_name == cityname, ]
  sims <- sim_results[city_name == cityname & model_name == modelname,]
  
  # Only hold onto treated locations:
  obs <- obs[cfho_any == 1,]
  
  treat_var <- unique(sims$treat_name)
  
  truth <- sum(obs$evict_count)
  obs <- obs[, (treat_var), with = F]
  
  # Remove extraneous columns from sims:
  sims <- sims[, .(est_beta)]
  
  # Calculate change from truth, if estimated beta held as causal effect:
  cf_change <- truth - colSums(as.matrix(obs) %*% t(as.matrix(sims)))
  
  df_cf <- data.table("obs" = truth, "cf" = cf_change)
  df_cf[, percent_change := (obs - cf)/obs]
  
  df_cf[, city_name := cityname]
  df_cf[, model_name := modelname]
  
  return(df_cf)
  
}

sim_cf_alt <- function(modelname, cityname){
  
  obs  <- df_sfd[city_name == cityname & cfho_any == 1, ]
  sims <- sim_results[city_name == cityname & model_name == modelname,]
  
  treat_var <- unique(sims$treat_name)
  
  truth <- mean(obs$evict_count)
  
  # Set treatment status to either equal 1 if
  # treat variable is "CFHO any" or mean across blocks,
  # if treat variable is "CFHO num"
  obs <- obs[, mean(get(treat_var))]

  # Remove extraneous columns from sims:
  sims <- sims[, .(est_beta)]
  
  # Calculate change from truth, if estimated beta held as causal effect:
  cf_change <- as.numeric(truth - obs*sims$est_beta)
  
  df_cf <- data.table("obs" = truth, "cf" = cf_change)
  df_cf[, percent_change := (obs - cf)/obs]
  
  df_cf[, city_name := cityname]
  df_cf[, model_name := modelname]
  
  return(df_cf)
  
}

sim_cf_overall <- function(modelname){
  
  evict_observed <- 0
  evict_cf       <- 0
  
  for (city in unique(df_sfd$city_name)){
    
    obs  <- df_sfd[city_name == city, ]
    sims <- sim_results[city_name == city & model_name == modelname,]
    
    # Only hold onto treated locations:
    obs <- obs[cfho_any == 1,]
    
    treat_var <- unique(sims$treat_name)
    
    truth <- sum(obs$evict_count)
    obs <- obs[, (treat_var), with = F]
    
    # Remove extraneous columns from sims:
    sims <- sims[, .(est_beta)]
    
    # Calculate change from truth, if estimated beta held as causal effect:
    cf_change <- truth - colSums(as.matrix(obs) %*% t(as.matrix(sims)))
    
    evict_observed <- evict_observed + truth
    evict_cf       <- evict_cf + cf_change
    
  }
  
  df_cf <- data.table("obs" = evict_observed, 
                      "cf" = evict_cf)
  
  df_cf[, percent_change := (obs - cf)/obs]
  
  df_cf[, city_name := "Across sites"]
  df_cf[, model_name := modelname]
  
  return(df_cf)
  
}

sim_beta_cf <- rbindlist(mapply(sim_cf_alt, 
                                modelname = args$modelname,
                                cityname = args$cityname, 
                                SIMPLIFY = F))

sim_beta_cf_overall <- rbindlist(mapply(sim_cf_overall, 
                                modelname = args$modelname,
                                SIMPLIFY = F))

sim_beta_cf <- rbindlist(list(sim_beta_cf, sim_beta_cf_overall))

ggplot(sim_beta_cf, aes(x = percent_change, y = model_name, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale = 1.7, alpha = 0.8) +
        scale_fill_viridis_c(name = "Percentage Change in Executed Evictions", option = "C") +
        facet_wrap(~city_name) +
        theme_ridges() +
        theme(strip.background = element_blank())

summary_beta_cf <- sim_beta_cf[, quantile(.SD$percent_change, probs = c(0.05, 0.5, 0.95)), by = c("city_name", "model_name")]

percentiles <- rep(c("q025", "q50", "q975"), dim(args)[1] + dim(args)[1]/4)
summary_beta_cf[, percentile := percentiles]
summary_beta_cf <- dcast(summary_beta_cf, formula = city_name + model_name ~ percentile, value.var = "V1")

summary_beta_cf[, city_name := ifelse(city_name == "san_diego_unincorporated",
                                       "San Diego County",
                                       str_to_title(city_name))]

# Set order of models to be more intuitive for the reader:
ord_mods <- c("CFHO any, unadjusted", "CFHO any, adjusted", 
              "CFHO count, unadjusted", "CFHO count, adjusted")

summary_beta_cf[, model_name := factor(model_name, levels = ord_mods)]
write.csv(summary_beta_cf, "sim_results.csv", row.names = F)

beta_cf_plot <- ggplot(summary_beta_cf, aes(x = q50, color = model_name, y = city_name)) +
                        geom_errorbarh(aes(xmin = q025, xmax = q975), size = 2, height = 0, position = position_dodgev(height = 0.4)) +
                        geom_point(size = 3, position = position_dodgev(height = 0.4)) +
                        labs(title = "Estimated effect of CFHO policy on total evictions, within treated cities",
                             x = "Percent change in evictions within city",
                             y = "",
                             color = "Model",
                             caption = "Bars are 95% UI.") +
                        theme_bw() +
                        scale_color_manual(values = c("#4d9221", "#b8e186", "#2166ac", "#92c5de")) +
                        scale_x_continuous(labels = percent, breaks = seq(0, 1, 0.1)) +
                        theme(plot.title = element_text(hjust = 0.5),
                              strip.background = element_blank(),
                              legend.position = "bottom",
                              axis.title = element_text(size = 12),
                              axis.text.y = element_text(size = 12),
                              axis.text.x = element_text(size = 12,
                                                         margin = margin(t = 5, r = 0, b = 10, l = 0)))

paper_copy <- data.table(city_name = unique(summary_beta_cf$city_name),
                         model_name = "CFHO any, unadjusted",
                         q025 = c(0.155, 0.236, 0.009, 0.051, 0.151),
                         q50 = c(0.41, 0.371, 0.171, 0.272, 0.249),
                         q975 = c(0.674, 0.507, 0.33, 0.492, 0.346))

df_plot <- summary_beta_cf[model_name == "CFHO any, unadjusted"]

exhibit_3 <- ggplot(df_plot, aes(x = q50, color = model_name, y = city_name)) +
            geom_errorbarh(aes(xmin = q025, xmax = q975), size = 2, height = 0) +
            geom_point(size = 3) +
            geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
            labs(title = "Exhibit 3. Estimated effect of CFHPs on evictions, within each location",
                 x = "Percentage Change in Executed Evictions",
                 y = "Location",
                 color = "Model") +
            theme_bw() +
            scale_color_manual(values = c("#969696")) +
            scale_x_continuous(labels = percent, breaks = seq(0, 1, 0.1)) +
            theme(plot.title = element_text(hjust = 0.5, family = 'sans', size = 16),
                  strip.background = element_blank(),
                  legend.position = "none",
                  legend.text = element_text(family = 'sans', size = 10),
                  axis.ticks = element_line(linewidth = 1),
                  axis.ticks.length = unit(5.6, "points"),
                  axis.title = element_text(size = 12, family = 'sans'),
                  axis.title.y = element_text(size = 12, family = 'sans', angle = 0),
                  axis.text = element_text(size = 10, family = 'sans'),
                  axis.text.x = element_text(size = 10, family = 'sans',
                                             margin = margin(t = 5, r = 0, b = 10, l = 0)),
                  legend.title = element_text(family = 'sans', size = 14))

ggsave(paste0(analysis_type, "plots/exhibit_3.pdf"), plot = exhibit_3, device = "pdf", height = 8.3*1.1, width = 11.7*1.1, units = "in")

sim_pred_cf <- rbindlist(mapply(sim_pred, 
                           cityname = args$cityname, 
                           modelname = args$modelname,
                           SIMPLIFY = F))

sim_pred_cf <- melt(sim_pred_cf, measure.vars = c("status_quo", "cf"), value.name = "est", variable.name = "simulation")

df_truth <- unique(sim_pred_cf[, .(city_name, model_name, obs)])

ggplot(sim_pred_cf, aes(x = est, y = simulation)) +
  geom_density_ridges(scale = 1.7, alpha = 0.9) +
  geom_vline(data = df_truth, mapping = aes(xintercept = obs), linetype = 2, lwd = 1) +
  scale_fill_viridis_c(name = "Estimated evictions", option = "C") +
  facet_wrap(city_name ~ model_name,  scales = 'free')
