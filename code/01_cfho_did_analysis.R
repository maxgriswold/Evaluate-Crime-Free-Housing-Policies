rm(list = ls())

library(data.table)
library(plyr)
library(dplyr)

library(ggplot2)
library(scales)
library(ggstance)
library(ggpubr)

library(bacondecomp)
library(fixest)
library(did)

main_dir <- "C:/users/griswold/Desktop/datasets/"
setwd(main_dir)

analysis_type <- 'main_analysis'

model_dir     <- paste0("./cfho_analysis/did/", analysis_type, "/models/")

df_analysis <- fread("./cfho_analysis/did/df_did_analysis.csv") 

# Prep event study dataset:
df_es <- df_analysis[year >= 2009 & year <= 2019,]

if (analysis_type == "sensitivity_analysis_nuisance"){
   
   # Set sites with nuisance ordinances to treated
   df_es[is.na(ever_treated), ever_treated := 1]
   df_es[ever_treated == 1, policy_type := "cfmhp"]
   
}

if (analysis_type == "reviewer_analysis"){
  
  df_es[, total_crime_rate_10k := asinh(total_crime_rate_10k)]
  df_es[, assault_rate_10k := asinh(assault_rate_10k)]
  df_es[, burglary_rate_10k := asinh(burglary_rate_10k)]
}

df_es <- df_es[ever_treated == 0 | (year_implemented >= 2010 & policy_type == "cfmhp"), ]
df_es[ever_treated == 1, ttt := year - year_implemented]

df_es <- df_es[, .(geoid, location, year, year_implemented, ttt, treated, ever_treated,
                   total_crime_rate_10k, assault_rate_10k, burglary_rate_10k, total_population, pop_latin_hispanic,
                   pop_black, pop_asian, pop_native_american, pop_white, pop_immigrant,
                   median_income_adj)]

perform_did_analysis <- function(outcome, dd, adjusters = NULL){
   
   #######################
   # Bacon Decomposition #
   #######################
   
   if (is.null(adjusters)){cov <- "unadjusted"}else{cov <- "adjusted"}
   
   model_name <- paste0(outcome, "_", cov)
   
   print(paste0("Running estimators for ", model_name))
   
   if (is.null(adjusters)){
      bacon_form <- formula(paste0(outcome, " ~ treated"))
   }else{
      bacon_form <- formula(paste0(outcome, " ~ treated + ", paste0(adjusters, collapse = " + ")))
   }
   
   # Construct bacon-goodman decomposition for our two outcomes
   bacon_res  <- bacon(formula = bacon_form, data = dd, 
                       id_var = "location", time_var = "year",
                       quietly = T)
   
   if (is.null(adjusters)){
      setDT(bacon_res)
   }else{
      bacon_res <- bacon_res$two_by_twos
      setDT(bacon_res)
   }
   
   # ggplot(bacon_res) +
   #       aes(x = weight, y = estimate, shape = factor(type)) +
   #       labs(x = "Weight", y = "Estimate", shape = "Type") +
   #       geom_point()
   
   if (is.null(adjusters)){
      bacon_res <- bacon_res[, `:=`(mean_est = round(weighted.mean(.SD$estimate, na.rm = T, w = .SD$weight), 3), 
                                    weight = round(sum(.SD$weight, na.rm = T), 3)), by = "type"]
   }else{
      bacon_res <- bacon_res[, `:=`(mean_est = round(weighted.mean(.SD$estimate, na.rm = T, w = .SD$weight), 3), 
                                    weight = round(unique(.SD$weight, na.rm = T), 3)), by = "type"]
   }

   
   bacon_res <- unique(bacon_res[, .(type, weight, mean_est)])
   write.csv(bacon_res, paste0("./cfho_analysis/did/", analysis_type, "/bacon/", model_name, "_bacon_decomposition.csv"), row.names = F)
   
   ########
   # TWFE #
   ########
   
   # Recode never-treated categories so that TWFE + C&SA estimate appropriately:
   dd[is.na(year_implemented), year_implemented := 0]
   dd[is.na(ttt), ttt := -1000]
   
   if (is.null(adjusters)){
      twfe_event_form  <- formula(paste0(outcome, " ~ i(ttt, ref = c(-1, -11, -1000)) | location + year"))
      twfe_att_form    <- formula(paste0(outcome, " ~ treated | location + year"))
   }else{
      twfe_event_form <- formula(paste0(outcome, " ~ i(ttt, ref = c(-1, -11, -1000)) + ", paste0(adjusters, collapse = " + "), "| location + year"))
      twfe_att_form    <- formula(paste0(outcome, " ~ treated + ", paste0(adjusters, collapse = " + "), "| location + year"))
   }
   
   twfe <- feols(twfe_event_form, dd, cluster = "location")
   
   # R will not serialize R objects unless you are in the respective folder. So quickly
   # change to preferred save location, then back to the main directory
   setwd(model_dir)
   saveRDS(twfe, paste0("twfe_", model_name, ".rds"))
   setwd(main_dir)
   
   twfe_hold <- data.table("estimate" = coef(twfe), 
                           "se" = se(twfe),
                           "eventTime" = as.numeric(gsub("ttt::", "", names(twfe$coefficients))),
                           "model" = "twfe")
   
   twfe_hold <- twfe_hold[!is.na(eventTime)]
   
   twfe_att <- feols(twfe_att_form, data = df_es, cluster = 'location')
   twfe_att <- data.table('model' = model_name,
                           'estimator' = 'twfe',
                           'estimate' = coef(twfe_att)[['treated']],
                           'se' = se(twfe_att)[['treated']])
   
   #######################
   # Sun & Abraham model #
   #######################
   
   if (is.null(adjusters)){
      sa_form <- formula(paste0(outcome, " ~ sunab(year_implemented, year) | location + year"))
   }else{
      sa_form <- formula(paste0(outcome, " ~ sunab(year_implemented, year, ref.p = c(-1, -11)) + ", paste0(adjusters, collapse = " + "), "| location + year"))
   }
   
   sa <- feols(sa_form, data = dd, cluster = "location")

   setwd(model_dir)
   saveRDS(sa, paste0("sa_", model_name, ".rds"))
   setwd(main_dir)
   
   sa_hold <- data.table("estimate" = coef(sa), 
                        "se" = se(sa),
                        "eventTime" = as.numeric(gsub("year::", "", rownames(summary(sa)[["coeftable"]]))),
                        "model" = "sa")
   
   sa_hold <- sa_hold[!is.na(eventTime)]
   
   sa_att <- aggregate(sa, "att")
   
   ########################
   # Callaway & Sant'Anna #
   ########################
   
   dd[, year_implemented := as.double(year_implemented)]
   dd[year_implemented == 0, year_implemented := Inf]
   
   csa <- att_gt(yname = outcome, 
                 tname = "year", 
                 idname = "geoid", 
                 gname = "year_implemented", 
                 control_group = "notyettreated", 
                 data = dd, 
                 xformla = NULL,
                 est_method = "dr", 
                 bstrap = T,
                 biters = 1000,
                 base_period	= "universal",
                 allow_unbalanced_panel = F)
   
   setwd(model_dir)
   saveRDS(csa, paste0("csv_", model_name, ".rds"))
   setwd(main_dir)
   
   csa_hold <- aggte(csa, type = 'dynamic')
   
   csa_hold <- data.table("estimate" = csa_hold$att.egt,
                          "se" = csa_hold$se.egt,
                          "eventTime" = csa_hold$egt,
                          "model" = "csa")
   
   csa_att <- aggte(csa, type = 'dynamic', balance_e = 3)
   csa_att <- data.table('model' = model_name,
                         'estimator' = 'csa',
                         'estimate' = csa_att$overall.att,
                         'se' = csa_att$overall.se)
   
   # Restructure ATT for both SA/CSA:
   sa_att <- data.table('model' = model_name,
                        'estimator' = 'sa',
                        'estimate' = sa_att[1],
                        'se' = sa_att[2])
   
   att <- rbindlist(list(twfe_att, sa_att, csa_att))
   
   # Construct and save event history plots: 
   dd_history <- rbindlist(list(twfe_hold, csa_hold, sa_hold))
   
   # Determine y-axis bounds
   dd_bound <- dd_history[eventTime <= 4 & eventTime >= -4,]
   dd_bound[, upper := estimate + 1.96*se]
   dd_bound[, lower := estimate - 1.96*se]
   
   lower <- floor(min(dd_bound$lower, na.rm = T))
   upper <- ceiling(max(dd_bound$upper, na.rm = T))
   
   new_outcome_name <- data.table(old = outcomes,
                                  new = c("Total crime rate (per 10k people)",
                                          "Assault rate (per 10k people)",
                                          "Burglary rate (per 10k people)"))
   
   # Create one outcome name for the title on single line; another for axis
   # with two lines
   new_name          <- new_outcome_name[old == outcome, ]$new
   new_name_carriage <- gsub("\\(", "\\\n\\(", new_name)
   
   # Remove rows for the reference period:
   dd_history <- dd_history[eventTime != -1]
   
   hplot <- ggplot(dd_history, aes(x = eventTime, y = estimate, color = model)) +
               geom_point(size = 3, position = position_dodge(width = 0.7)) +
               geom_errorbar(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), 
                             width = 0, size = 1.5, position = position_dodge(width = 0.7)) +
               geom_hline(aes(yintercept = 0), linetype = 'dashed', size = 1) +
               scale_x_continuous(limits = c(-4.5, 4.5), breaks = seq(-4, 4)) +
               scale_color_manual(values=c("#4d9221", "#8856a7", "#2166ac"), 
                                  breaks = c('csa', "sa", "twfe"), 
                                  labels = c("Callaway and Sant'Anna", "Sun and Abraham", 'Two-way Fixed Effects')) +
               labs(title = paste0("Event history estimates for outcome: ", new_name),
                    y = new_name_carriage, x = "Years from policy adoption", color = "Estimator") +
               lims(y = c(lower, upper)) +
               theme_bw() +
               theme(legend.position = 'bottom',
                     plot.title = element_text(hjust = 0.5),
                     axis.ticks = element_line(linewidth = 1),
                     axis.title = element_text(family = "serif", size = 18),
                     axis.title.y = element_text(angle = 0, hjust = 0.9),
                     axis.ticks.length = unit(5.6, "points"),
                     axis.text = element_text(family = "serif", size = 12),
                     legend.text = element_text(family = 'serif', size = 12))
   
   ggsave(hplot, filename = paste0("./cfho_analysis/did/", analysis_type, "/event_study_plots/event_", model_name, ".svg"),
          device = "svg", height = 8.3, width = 11.7, units = "in")
   
   ggsave(hplot, filename = paste0("./cfho_analysis/did/", analysis_type, "/event_study_plots/event_", model_name, ".png"),
          device = "png", height = 8.3, width = 11.7, units = "in")
   
   write.csv(dd_history, paste0("./cfho_analysis/did/", analysis_type, "/event_study_plots/df_event_", model_name, ".csv"), row.names = F)
   
   return(att)
   
}

outcomes <- c("total_crime_rate_10k", "assault_rate_10k", "burglary_rate_10k")
covs     <- c('pop_white', 'pop_black', 'pop_asian', 'pop_native_american', 
                'pop_latin_hispanic', 'median_income_adj')

att_unadj <- rbindlist(llply(outcomes, perform_did_analysis, dd = df_es))
att_adj   <- rbindlist(llply(outcomes, perform_did_analysis, adjusters = covs, dd = df_es))

# Combine ATT into single dataset:
att_res <- setDT(rbind(att_unadj, att_adj))
att_res[, outcome := gsub("_adjusted|_unadjusted", "", model)]
att_res[, adjusted := ifelse(model %like% "unadjusted", "Unadjusted\nmodels", "Adjusted\nmodels")]
att_res[, adjusted := factor(adjusted, levels = c( "Unadjusted\nmodels", "Adjusted\nmodels"))]

att_res[, c025 := estimate - 1.96*se]
att_res[, c050 := estimate - 1.645*se]
att_res[, c950 := estimate + 1.645*se]
att_res[, c975 := estimate + 1.96*se]

# Modify output for simulation purposes later:
att_sim <- copy(att_res)
att_sim[, model := paste(model, estimator, sep = "_")]
att_sim <- att_sim[, .(model, outcome, estimator, adjusted, estimate,  c025, c050, c950, c975)]

# Set order of factor variables for plot and rename:
estimator_names <- c("Two-way Fixed Effects", "Sun & Abraham", "Callaway & Sant'Anna")
att_res[, estimator := mapvalues(estimator,  unique(att_res$estimator), estimator_names)]
att_res[, estimator := factor(estimator, levels = rev(estimator_names))]

outcome_names <- c("Total crime rate (per 10k people)", "Assault rate (per 10k people)", "Burglary rate (per 10k people)")
att_res[, outcome := mapvalues(outcome,  unique(att_res$outcome), outcome_names)]
att_res[, outcome := factor(outcome, levels = rev(outcome_names))]

att_res <- att_res[, .(outcome, adjusted, estimator, estimate, se, c025, c050, c950, c975)]
write.csv(att_res, paste0("./cfho_analysis/did/", analysis_type,"/att_results.csv"), row.names = F)

# Plot ATT estimates across outcomes, models, and estimators:
plot_att_est <- ggplot(att_res, aes(x = estimate, color = estimator, y = outcome)) +
                  geom_errorbarh(aes(xmin = c025, xmax = c975), height = 0, position = position_dodgev(height = 0.4)) +
                  geom_point(size = 3, position = position_dodgev(height = 0.4)) +
                  geom_vline(xintercept = 0, linetype = 2, size = 1) +
                  labs(title = "Estimated ATT effect of CFHO policy, by estimator and model specification",
                       x = "Estimated Treatment Effect",
                       y = "",
                       color = "Estimator") +
                  facet_wrap(~adjusted, ncol = 1) +
                  theme_bw() +
                  lims(x = c(-0.3, 0.1)) +
                  scale_color_manual(values = c("#4d9221", "#8856a7", "#2166ac")) +
                  theme(plot.title = element_text(hjust = 0.5),
                        strip.background = element_blank(),
                        legend.position = "bottom",
                        axis.title = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12,
                                                   margin = margin(t = 5, r = 0, b = 10, l = 0)))

ggsave(filename = paste0("./cfho_analysis/did/", analysis_type,"/att_results.png"), 
       plot = plot_att_est, device = "png", height = 8.3, width = 11.7, units = "in")

ggsave(filename = paste0("./cfho_analysis/did/", analysis_type,"/att_results.svg"), 
       plot = plot_att_est, device = "svg", height = 8.3, width = 11.7, units = "in")

# Simulate counterfactual change in outcome, assuming estimated ATTs held:
df_sim <- copy(df_es[, .(location, year, ever_treated, total_crime_rate_10k, assault_rate_10k, burglary_rate_10k, total_population)])

sim_cf <- function(outcome_name, model_name, df_obs, df_att){
   
   obs  <- df_obs[year == 2019 & ever_treated == 1, c("location", "ever_treated", "total_population", outcome_name), with = F]
   sims <- df_att[model == model_name & outcome == outcome_name,]
   
   # Calculate total count for outcome across treated sites in 2019:
   truth <- sum(obs[, ceiling((get(outcome_name)*total_population)/10000)])
   
   # Calculate change from truth, if estimated beta held as causal effect (in count space)
   cf_change <- truth - colSums(((obs$ever_treated %*% as.matrix(sims[,5:9]))*obs$total_population)/10000)
   df_cf <- as.data.table(cf_change)

   df_cf[, truth := truth]
   df_cf[, percent_change := round((truth - cf_change)/truth, 3)]
   
   df_cf[, ci := c("q50", "q025", "q05", "q95", "q975")]
   
   df_cf[, outcome := outcome_name]
   df_cf <- dcast(df_cf, outcome ~ ci, value.var = "percent_change")
   
   df_cf[, estimator := sims$estimator]
   df_cf[, adjusted := sims$adjusted]
   
   return(df_cf)
   
}

# Rather than simulate possible ATT,  apply ATT to each site, aggregate, then
# summarize percent changes across ATT draws, instead apply ATT draws to
# average rate to calculate a CF, then summarize.

sim_cf_alt <- function(outcome_name, model_name, df_obs, df_att){
  
  obs  <- df_obs[year == 2019 & ever_treated == 1, c("location", "ever_treated", "total_population", outcome_name), with = F]
  sims <- df_att[model == model_name & outcome == outcome_name,]
  
  # Calculate average outcome across treated sites:
  truth <- obs[, mean(get(outcome_name), na.rm = T)]
  
  # Calculate change from truth, if estimated beta held as causal effect (in count space)
  cf_change <- truth - sims[,5:9]
  df_cf <- as.data.table(cf_change)
  
  df_cf[, truth := truth]
  df_cf <- melt(df_cf, id.vars = "truth", value.name = "cf",
                variable.name = "ci")
  
  df_cf[, percent_change := round((truth - cf)/truth, 3)]
  
  df_cf <- dcast(df_cf, . ~ ci, value.var = "percent_change")
  
  df_cf <- df_cf[, .(estimate, c025, c975)]
  setnames(df_cf, names(df_cf), c("q50", "q025", "q975"))
  
  df_cf[, outcome := outcome_name]

  df_cf[, estimator := sims$estimator]
  df_cf[, adjusted := sims$adjusted]
  
  return(df_cf)
  
}

if (analysis_type == "counterfactual_analysis"){
  df_sim_cf <- rbindlist(mapply(sim_cf_alt, model_name = att_sim$model, outcome = att_sim$outcome,
                                MoreArgs = list(df_obs = df_analysis, df_att = att_sim),
                                SIMPLIFY = F))
}else{
  df_sim_cf <- rbindlist(mapply(sim_cf, model_name = att_sim$model, outcome = att_sim$outcome,
                                MoreArgs = list(df_obs = df_sim, df_att = att_sim),
                                SIMPLIFY = F))
}

# Change names like previous plot:
df_sim_cf <- df_sim_cf[, estimator := mapvalues(estimator, unique(df_sim_cf$estimator), estimator_names)]
df_sim_cf[, estimator := factor(estimator, levels = rev(estimator_names))]

outcome_names <- c("Total crime count", "Assault count", "Burglary count")
df_sim_cf <- df_sim_cf[, outcome := mapvalues(outcome, unique(df_sim_cf$outcome), outcome_names)]
df_sim_cf[, outcome := factor(outcome, levels = rev(outcome_names))]

write.csv(df_sim_cf, paste0("./cfho_analysis/did/", analysis_type,"/sim_cf_results.csv"), row.names = F)

plot_sim_cf <- ggplot(df_sim_cf, aes(x = q50, color = estimator, y = outcome)) +
                  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, position = position_dodgev(height = 0.4)) +
                  geom_point(size = 4, position = position_dodgev(height = 0.4)) +
                  geom_vline(xintercept = 0, linetype = 2, size = 1) +
                  labs(title = "Simulated effect of CFHO policy, by estimator, and model specification",
                       x = "Percent change",
                       y = "Outcome",
                       color = "Estimator",
                       caption = "Bars are 95% CI.") +
                  facet_wrap(~adjusted, ncol = 1) +
                  theme_bw() +
                  scale_color_manual(values = c("#4d9221", "#8856a7", "#2166ac")) +
                  scale_x_continuous(labels = percent, limits = c(-0.3, 0.1)) +
                  theme(plot.title = element_text(hjust = 0.5),
                        strip.background = element_blank(),
                        legend.position = "bottom",
                        axis.title = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12,
                                                   margin = margin(t = 5, r = 0, b = 10, l = 0)))

ggsave(filename = paste0("./cfho_analysis/did/", analysis_type, "/sim_cf.svg"),
       plot = plot_sim_cf, device = "svg", height = 8.3, width = 11.7, units = "in")

df_capstone <- df_sim_cf[estimator == "Callaway & Sant'Anna" & adjusted == "Unadjusted\nmodels", .(outcome, q025, q50, q975)]

df_add  <- fread("./cfho_analysis/sfd/sim_results.csv")
df_add[, outcome := "Eviction count"]
df_add <- df_add[city_name == "Across All Sites" & model_name == "CFHO count, adjusted", .(outcome, q025, q50, q975)]

df_capstone <- rbind(df_capstone, df_add)

capstone <- ggplot(df_capstone, aes(x = q50, y = outcome)) +
               geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 0.5) +
               geom_errorbarh(aes(xmin = q025, xmax = q975), size = 1.5, height = 0, position = position_dodgev(height = 0.4),
                              color = "#8856a7") +
                 geom_point(size = 4, position = position_dodgev(height = 0.4), color = "#8856a7") +
               labs(x = "Percent change",
                    y = "Outcome",
                    color = "Outcome") +
               theme_bw() +
               scale_color_manual(values = c("#4d9221", "#8856a7", "#2166ac", '#d95f0e')) +
               scale_x_continuous(labels = percent, breaks = seq(-0.3, 0.5, 0.1)) +
               theme(plot.title = element_text(hjust = 0.5),
                     strip.background = element_blank(),
                     legend.position = "none",
                     axis.title = element_text(family = "serif", size = 18),
                     axis.title.y = element_text(angle = 0, vjust = 1),
                     axis.text = element_text(family = "serif", size = 12),
                     axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 10, l = 0)),
                     axis.ticks.length = unit(5.6, "points"),
                     axis.ticks = element_line(linewidth = 1))

ggsave(filename = paste0("./cfho_analysis/figs_collected/figure_3.svg"),
       plot = capstone, device = "svg", height = 8.3, width = 11.7, units = "in")

ggsave(filename = paste0("./cfho_analysis/figs_collected/figure_3.png"),
       plot = capstone, device = "png", height = 8.3, width = 11.7, units = "in")

# # Simulate varying number of control units compared to treated & calculate
# # percent early v. late comparison weight percentages in goodman-bacon decomp:
# 
# sim_scenarios <- function(nc, nt, d){
#   
#   # Hold onto all treated units and a select number of controls
#   d <- d[geoid %in% nt|geoid %in% nc,]
#   
#   model_sa_assault  <- feols(burglary_total ~ sunab(ever_treated, ttt) | location + year, d, cluster = "location")
#   model_ols_assault <- feols(burglary_total ~ treated  | location + year, d, cluster = "location")
#   
#   sa   <- coef(summary(model_sa_assault, agg = "ATT"))[['ATT']]
#   twfe <- coef(model_ols_assault)[['treated']]
#   
#   #Calculate percentage of estimate driven by "forbidden comparisons"
#   res <- data.table('n_control' = length(nc),
#                     'sf' = sa,
#                     'twfe' = twfe)
#   
#   return(res)
#   
# }
# 
# control_scenarios <- list()
# control_sites <- unique(df_es[ever_treated == 0, geoid])
# 
# for (i in seq(10, length(control_sites), 5)){
#   
#   if (i == 10){
#     
#     add <- sample(control_sites, 5)
#     control_scenarios[[i]] <- c(add)
#     
#     remaining_controls <- control_sites[control_sites != add]
#     
#   }else{
#     
#     add <- sample(remaining_controls, 5)
#     control_scenarios[[i]] <- c(control_scenarios[[i - 5]], add)
#     
#     remaining_controls <- remaining_controls[remaining_controls != add]
#     
#   }
# }
# 
# treatment_scenarios <- list(c(sample(unique(df_es[year_implemented >= 2011 & year_implemented <= 2015, geoid]), 6),
#                               sample(unique(df_es[year_implemented > 2015, geoid]), 3)),
#                             sample(unique(df_es[year_implemented >= 2011, geoid]), 9),
#                             c(sample(unique(df_es[year_implemented >= 2011 & year_implemented <= 2015, geoid]), 3),
#                                 sample(unique(df_es[year_implemented > 2015, geoid]), 6)))
# 
# sim1 <- setDT(ldply(control_scenarios, sim_scenarios, d = df_es, nt = treatment_scenarios[[1]]))
# sim1[, scenario := "Few early-treated"]
# 
# sim2 <- setDT(ldply(control_scenarios, sim_scenarios, d = df_es, nt = treatment_scenarios[[2]]))
# sim2[, scenario := "Balanced treatment"]
# 
# sim3 <- setDT(ldply(control_scenarios, sim_scenarios, d = df_es, nt = treatment_scenarios[[3]]))
# sim3[, scenario := "Few late-treated"]
# 
# sims <- rbind(sim1, sim2, sim3)
# 
# full_est <- coef(summary(feols(burglary_total ~ sunab(ever_treated, ttt) | location + year, df_es, cluster = "location"), agg = 'ATT'))
# 
# sims[, diff := 100*(sf - twfe)/twfe]
# 
# ggplot(sims, aes(x = n_control, y = diff, color = scenario)) + 
#   geom_line(size = 1) + 
#   labs(title = "Difference in effect size: \nSun & Abraham \ncompared to TWFE estimate", 
#        x = "Number of \ncontrol units", 
#        y = "Relative percent \ndifference in \nestimated \neffect size",
#        color = "Treatment timing \nscenario") + 
#   lims(y = c(-125, 125)) +
#   theme_bw() + 
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 11),
#         axis.title.x = element_text(size = 11),
#         plot.title = element_text(hjust = 0.5, size = 14),
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10),
#         legend.position = "bottom")
# 
# ggsave("bacon_ex.jpg", height = 8, width = 10, units = "in")

plot(df_analysis[location == "claremont",]$year, df_analysis[location == "claremont",]$assault_rate_10k, xlab = "Year", ylab = "Violent crime rate (per 10k people)")
lines(df_analysis[location == "claremont",]$year, df_analysis[location == "claremont",]$assault_rate_10k)
