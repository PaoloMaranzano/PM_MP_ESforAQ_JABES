###################################################
########## Spatio-temporal Event Studies ##########
###################################################

setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/PM_MP_ESforAQ_JABES")

library(readr)
library(tidyverse)
library(moments)
library(lubridate)
library(xts)
library(forecast)
library(mgcv)

out_path <- "Data/"
Data <- read_csv(paste0(out_path,"HDGM_output.csv"))

Estim_model <- T
rolling_window_analysis <- T

source('Aux_funs/FN - ES_Estimation.R', encoding = 'UTF-8')
source('Aux_funs/FN - ES_AbnRes_Stats.R', encoding = 'UTF-8')
source('Aux_funs/FN - ES_Stats_Extract.R', encoding = 'UTF-8')

estim_window_start <- dmy("01/01/2018")
estim_window_end <- dmy("31/01/2020")
event_window_start <- dmy("01/02/2020")
event_window_end <- dmy("31/05/2020")
y_name <- "NO2_log"
date_name <- "Date"
idx_name <- "Stz_Code"
X_names <- c("weekend", "Holidays","Temperature", "Rainfall", "Pressure",
             "WindU", "WindV", "WindU_max", "WindV_max", 
             "VegetationHigh", "VegetationLow",
             "GeopotHeight", "RelHumid",
             "Lag",
             "MetrArea","Mountain","UrbPlain",
             "Ind","Traf","Rural",
             "_B1",
             "S1","S2","C1","C2"
)
X_names <- colnames(Data)[grepl(pattern = paste(X_names, collapse = "|"), x = colnames(Data)) & 
                            !grepl(pattern = paste(c("Summer","Winter","Autumn","Spring"), collapse = "|"), x = colnames(Data))]
TV_X_names <- colnames(Data)[grepl(pattern = paste(c("Temperature", "Rainfall", "Pressure",
                                                     "WindU", "WindV",
                                                     "RelHumid","SolarRadiation"), collapse = "|"), x = colnames(Data))]

cols <- c("Background" = "#CC33FF",
          "Rural" = "orange",
          "Traffic" = "#3399FF",
          "Industrial" = "#339900",
          #####
          "Rural plain" = "#FF9900",
          "Urban plain" = "#666666",
          "Metrop. area" = "#3399FF",
          "Mountains" = "#339900",
          #####
          "HDGM" = "red",
          "GAMM with AR(1) dynamics" = "orange",
          "GAMM" = "#99CCFF",
          "GAM" = "green",
          #####
          "Extreme weather event" = "green",
          "Lockdown" = "red",
          "Train" = "blue",
          "8th March to 18th May" = "blue",
          #####
          "10%" = "#CCCCCC",
          "5%" = "#666666",
          "1%" = "#333333"
)

#######################################
########## Models estimation ##########
#######################################
# Full dataset
Data <- Data %>%
  group_by(Stz_Code) %>%
  mutate(Time = seq(from=1,to=n()),
         TimeTrend = seq(from=1,to=n())/n(),
         Stz_Code = as.factor(Stz_Code)) %>%
  ungroup()
# Training set
Data_estim <- Data %>%
  filter(Window == "Estim") %>%
  group_by(Stz_Code) %>%
  mutate(Time = seq(from=1,to=n()),
         TimeTrend = seq(from=1,to=761)/761,
         Stz_Code = as.factor(Stz_Code))
# Evaluation set
Data_event <- Data %>%
  filter(Window == "Event") %>%
  group_by(Stz_Code) %>%
  mutate(Time = seq(from=1,to=n()),
         TimeTrend = seq(from=1,to=120)/120,
         Stz_Code = as.factor(Stz_Code))

##### Formulae
form_gam <- as.formula(paste0("NO2_log ~ ",paste(X_names,collapse = "+"),
                              "+ s(Stz_Long, Stz_Lat, bs='gp', m = c(2))"))
form_gamm2 <- as.formula(paste0("NO2_log ~ ",paste(X_names,collapse = "+"),
                                "+ te(Stz_Long, Stz_Lat, bs='tp')"))

if (Estim_model == T) {
  ########## Model 1: GAM with non-linear smooth (GP) spatial trend
  m_gam <- gam(formula = form_gam,
               # subset = Data_estim$Stz_Code %in% paste0("stz_",seq(20,40)),
               data = Data_estim,
               method = "REML")
  performance::performance(m_gam)
  gam.check(m_gam)
  th <- seq(from = 10, to = 90, by = 10)
  win.graph()
  par(mfrow = c(3,3))
  for (i in 1:9) {
    vis.gam(m_gam, theta=th[i], main = "M1: GAM", view = c("Stz_Long", "Stz_Lat"))
  }
  par(mfrow = c(1,1))
  # Extract large scale range
  range_gam <- m_gam$smooth[[1]]$gp.defn[2]
  # Compute predictions
  Data$pred_log_gam <- predict.gam(object = m_gam, newdata = Data, type = "response")
  Data <- Data %>%
    mutate(NO2_hat_log_GAM = pred_log_gam,
           NO2_hat_GAM = exp(NO2_hat_log_GAM))
  
  
  ########## Model 2: GAMM with non-linear smooth (GP) spatial trend and 
  ###                     independent random effects across space/sites
  m_gamm <- gamm(formula = form_gam,
                 # subset = Data_estim$Stz_Code %in% paste0("stz_",seq(20,35)),
                 random=list(Stz_Code=~1),
                 data = Data_estim,
                 method = "REML")
  # Extract large scale range
  range_gamm <- m_gamm$gam$smooth[[1]]$gp.defn[2]
  # Extract random effects
  r_gamm <- ranef(m_gamm$lme, level = 2)
  rownames(r_gamm) <- substr(rownames(r_gamm),start=3,stop=8)
  r_gammdf <- r_gamm %>%
    rownames_to_column() %>%
    rename(Stz_Code = 'rowname', GAMMeff = '(Intercept)')
  # Compute predictions without REs
  Data$pred_log_gamm <- predict(object = m_gamm$gam, newdata = Data, type = "response")
  Data <- left_join(x = Data, y = r_gammdf, by = c("Stz_Code"))
  # Add REs to predictions
  Data <- Data %>%
    mutate(NO2_hat_log_GAMM = pred_log_gamm + GAMMeff,
           NO2_hat_GAMM = exp(NO2_hat_log_GAMM))
  
  
  ########## Model 3: GAMM with non-linear smooth (GP) spatial trend and 
  ###                     independent random effects across space/sites and
  ###                     within-group AR(1) residuals
  m_gammAR1 <- gamm(formula = form_gam,
                    # subset = Data_estim$Stz_Code %in% paste0("stz_",seq(20,35)),
                    random=list(Stz_Code=~1),
                    correlation = corARMA(form=~1|Stz_Code, p=1),
                    data=Data_estim,
                    method = "REML")
  # Extract large scale range
  range_gammAR1 <- m_gammAR1$gam$smooth[[1]]$gp.defn[2]
  # Extract random effects
  rAR1 <- ranef(m_gammAR1$lme, level = 2)
  rownames(rAR1) <- substr(rownames(rAR1),start=3,stop=8)
  rAR1df <- rAR1 %>%
    rownames_to_column() %>%
    rename(Stz_Code = 'rowname', AR1eff = '(Intercept)')
  # Compute predictions without REs
  Data$pred_log_gamm_ar1 <- predict(object = m_gammAR1$gam, newdata = Data, type = "response")
  Data <- left_join(x = Data, y = rAR1df, by = c("Stz_Code"))
  # Add REs to predictions and AR(1) effects to the grouped residuals
  Data$res_log_gamm_ar1 <- Data$NO2_log - Data$pred_log_gamm_ar1
  Data <- Data %>%
    group_by(Stz_Code) %>%
    mutate(L1_res_log_gamm_ar1 = lag(res_log_gamm_ar1,n=1)) %>%
    ungroup() %>%
    mutate(L1_res_log_gamm_ar1 = replace_na(L1_res_log_gamm_ar1,0),
           NO2_hat_log_GAMM_ar1 = pred_log_gamm_ar1 + AR1eff + L1_res_log_gamm_ar1*0.6476572,
           NO2_hat_GAMM_ar1 = exp(NO2_hat_log_GAMM_ar1))
  
  
  ########## Save output
  save(Data,Data_estim,Data_event,m_gam,m_gamm,m_gammAR1,
       date_name,idx_name,X_names,TV_X_names,
       estim_window_start,estim_window_end,event_window_start,event_window_end,
       file = paste0(out_path,"OutputGAMs.RData"))
} else if (Estim_model == F) {
  load(paste0(out_path,"OutputGAMs.RData"))
}


#########################################################
########## Abnormal concentrations computation ##########
#########################################################
if (Estim_model == T) {
  y_name <- "NO2_log"
  univHDGM_ES <- ES_estimation(Data = Data,
                               estim_window_start,estim_window_end,event_window_start,event_window_end,
                               y_name,y_hat_name = "NO2_hat_log_HDGM",date_name,idx_name,X_names,
                               TV_X_names,
                               model_type = NULL)
  
  regGAM_ES <- ES_estimation(Data = Data,
                             estim_window_start,estim_window_end,event_window_start,event_window_end,
                             y_name,y_hat_name = "NO2_hat_log_GAM",date_name,idx_name,X_names,
                             TV_X_names = TV_X_names,
                             model_type = NULL)
  
  regGAMM_ES <- ES_estimation(Data = Data,
                              estim_window_start,estim_window_end,event_window_start,event_window_end,
                              y_name,y_hat_name = "NO2_hat_log_GAMM",date_name,idx_name,X_names,
                              TV_X_names = TV_X_names,
                              model_type = NULL)
  
  regGAMMar1_ES <- ES_estimation(Data = Data,
                                 estim_window_start,estim_window_end,event_window_start,event_window_end,
                                 y_name,y_hat_name = "NO2_hat_log_GAMM_ar1",date_name,idx_name,X_names,
                                 TV_X_names = TV_X_names,
                                 model_type = NULL)
  ##### Save output
  save.image(paste0(out_path,"Output4models.RData"))
} else {
  load(paste0(out_path,"Output4models.RData"))
}



############################################################
########## Abnormal concentrations representation ##########
############################################################

# GAM
res_GAM <- data.frame(regGAM_ES$epsilon_it)
res_GAM <- res_GAM %>%
  mutate(Model = "GAM",
         Date = index(regGAM_ES$epsilon_it)) %>%
  select(Date,Model,everything())
# GAMM
res_GAMM <- data.frame(regGAMM_ES$epsilon_it)
res_GAMM <- res_GAMM %>%
  mutate(Model = "GAMM",
         Date = index(regGAMM_ES$epsilon_it)) %>%
  select(Date,Model,everything())
# GAMM with AR(1)
res_GAMMar1 <- data.frame(regGAMMar1_ES$epsilon_it)
res_GAMMar1 <- res_GAMMar1 %>%
  mutate(Model = "GAMM_ar1",
         Date = index(regGAMMar1_ES$epsilon_it)) %>%
  select(Date,Model,everything())
# HDGM
res_HDGM <- data.frame(univHDGM_ES$epsilon_it)
res_HDGM <- res_HDGM %>%
  mutate(Model = "HDGM",
         Date = index(univHDGM_ES$epsilon_it)) %>%
  select(Date,Model,everything())

# Bind and reshape
res <- bind_rows(res_HDGM,res_GAMM,res_GAMMar1,res_GAM)
res <- res %>%
  pivot_longer(cols = 3:last_col(), names_to = "Station", values_to = "Value") %>%
  mutate(Period = case_when(between(Date,estim_window_start,estim_window_end) ~ "Estim. window",
                            between(Date,event_window_start,event_window_end) ~ "Event window",
                            TRUE ~ NA_character_)) %>%
  group_by(Model,Period) %>%
  mutate(Period_MeanValue = mean(Value,na.rm = T)) %>%
  ungroup()
# Add station information
Stz_info <- Data %>%
  select(Stz_Code,Stz_Type_rec,Stz_ARPA_zone_rec) %>%
  unique()
res <- left_join(x = res, y = Stz_info, by = c("Station" = "Stz_Code"))
res <- res %>%
  group_by(Date,Stz_Type_rec,Model) %>%
  mutate(Type_MeanValue = mean(Value,na.rm = T)) %>%
  ungroup() %>%
  group_by(Date,Stz_ARPA_zone_rec,Model) %>%
  mutate(ARPAZone_MeanValue = mean(Value,na.rm = T))

# Plot whole time series
win.graph()
p_AC <- res %>%
  filter(Date >= "2020-01-01") %>%
  ggplot(mapping = aes(x = Date, y = Value)) + 
  geom_line(mapping = aes(col = Station)) + 
  geom_hline(yintercept = 0, col = "black", linewidth = 2) + 
  geom_smooth(method = "gam", col = "green", linewidth = 2) + 
  ylim(-4, 2) + 
  facet_wrap( ~ Model) + 
  labs(y = latex2exp::TeX("$\\mu$g/$m^3$"),
       x = "First lockdown wave in 2020",
       title = "Estimated abnormal log-concentrations by station") + 
  geom_rect(data = data.frame(Period  = factor(c("Train",
                                                 "Lockdown",
                                                 "Extreme weather event",
                                                 "Extreme weather event",
                                                 "Extreme weather event"),
                                               levels = c("Train",
                                                          "Lockdown",
                                                          "Extreme weather event")),
                              start = c(as.POSIXct("2020-01-01"),
                                        as.POSIXct("2020-03-08"),
                                        as.POSIXct("2020-02-01"),
                                        as.POSIXct("2020-02-12"),
                                        as.POSIXct("2020-02-27")),
                              end   = c(as.POSIXct("2020-01-31"),
                                        as.POSIXct("2020-05-18"),
                                        as.POSIXct("2020-02-06"),
                                        as.POSIXct("2020-02-15"),
                                        as.POSIXct("2020-03-02"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  theme_classic() + 
  theme(title = element_text(size = 12,face = "bold"),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Station type", values = cols) + 
  scale_fill_manual("Period", values = cols)
ggpubr::ggexport(p_AC,width = 1800, height = 1200, res = 150,
                 filename = "AC_by_model.png")

# Plot TS by station type
p_type <- res %>%
  filter(Date >= "2020-01-01") %>%
  ggplot(mapping = aes(x = Date, y = Type_MeanValue, col = Stz_Type_rec)) + 
  geom_line() + 
  geom_hline(yintercept = 0, col = "black", linewidth = 2) + 
  geom_smooth(method = "gam") + 
  facet_wrap( ~ Model) + 
  labs(y = latex2exp::TeX("$\\mu$g/$m^3$"),
       x = "First lockdown wave in 2020",
       title = "Estimated average abnormal concentrations by station type") + 
  geom_rect(data = data.frame(Period  = factor(c("Train",
                                                 "Lockdown",
                                                 "Extreme weather event",
                                                 "Extreme weather event",
                                                 "Extreme weather event"),
                                               levels = c("Train",
                                                          "Lockdown",
                                                          "Extreme weather event")),
                              start = c(as.POSIXct("2020-01-01"),
                                        as.POSIXct("2020-03-08"),
                                        as.POSIXct("2020-02-01"),
                                        as.POSIXct("2020-02-12"),
                                        as.POSIXct("2020-02-27")),
                              end   = c(as.POSIXct("2020-01-31"),
                                        as.POSIXct("2020-05-18"),
                                        as.POSIXct("2020-02-06"),
                                        as.POSIXct("2020-02-15"),
                                        as.POSIXct("2020-03-02"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  theme_classic() + 
  theme(title = element_text(size = 12,face = "bold"),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Station type", values = cols) + 
  scale_fill_manual("Period", values = cols)
ggpubr::ggexport(p_type,width = 1800, height = 1200, res = 150,
                 filename = "AC_by_stattype.png")

# Plot TS by zoning
p_zone <- res %>%
  filter(Date >= "2020-01-01") %>%
  mutate(Stz_ARPA_zone_rec = case_when(Stz_ARPA_zone_rec == "ARPA_MetrArea" ~ "Metrop. area",
                                       Stz_ARPA_zone_rec == "ARPA_Mountain" ~ "Mountains",
                                       Stz_ARPA_zone_rec == "ARPA_RurPlain" ~ "Rural plain",
                                       Stz_ARPA_zone_rec == "ARPA_UrbPlain" ~ "Urban plain")) %>%
  ggplot(mapping = aes(x = Date, y = ARPAZone_MeanValue, col = Stz_ARPA_zone_rec)) + 
  geom_line() + 
  geom_hline(yintercept = 0, col = "black", linewidth = 2) + 
  geom_smooth(method = "gam") + 
  facet_wrap( ~ Model) + 
  labs(y = latex2exp::TeX("$\\mu$g/$m^3$"),
       x = "First lockdown wave in 2020",
       title = "Estimated average abnormal concentrations by ARPA zoning") + 
  geom_rect(data = data.frame(Period  = factor(c("Train",
                                               "Lockdown",
                                               "Extreme weather event",
                                               "Extreme weather event",
                                               "Extreme weather event"),
                                             levels = c("Extreme weather event",
                                                        "Lockdown",
                                                        "Train")),
                              start = c(as.POSIXct("2020-01-01"),
                                        as.POSIXct("2020-03-08"),
                                        as.POSIXct("2020-02-01"),
                                        as.POSIXct("2020-02-12"),
                                        as.POSIXct("2020-02-27")),
                              end   = c(as.POSIXct("2020-01-31"),
                                        as.POSIXct("2020-05-18"),
                                        as.POSIXct("2020-02-06"),
                                        as.POSIXct("2020-02-15"),
                                        as.POSIXct("2020-03-02"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  theme_classic() + 
  theme(title = element_text(size = 12,face = "bold"),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Zone", values = cols) + 
  scale_fill_manual("Period", values = cols)
ggpubr::ggexport(p_zone,width = 1800, height = 1200, res = 150,
                 filename = "AC_by_zone.png")

# Time series of observed NO2 concentrations
p_TS <- Data %>%
  ggplot(mapping = aes(x = Date, y = NO2)) + 
  geom_line(mapping = aes(col = Stz_Code), show.legend = F) + 
  labs(y = latex2exp::TeX("$\\mu$g/$m^3$"),
       x = "",
       title = "Observed NO2 concentrations",
       subtitle = "Overlapping grey lines are the time series at station level") + 
  geom_rect(data = data.frame(Period  = factor(c("8th March to 18th May",
                                                 "8th March to 18th May",
                                                 "Lockdown"),
                                               levels = c("Lockdown",
                                                          "8th March to 18th May")),
                              start = c(as.POSIXct("2018-03-08"),
                                        as.POSIXct("2019-03-08"),
                                        as.POSIXct("2020-03-08")),
                              end   = c(as.POSIXct("2018-05-18"),
                                        as.POSIXct("2019-05-18"),
                                        as.POSIXct("2020-05-18"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  theme_classic() + 
  theme(title = element_text(size = 12,face = "bold"),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Zone", values = cols) + 
  scale_fill_manual("Period", values = cols)
ggpubr::ggexport(p_TS,width = 1800, height = 1200, res = 150, filename = "TS_lockdown.png")





##############################################################
########## Rolling-window for event date estimation ##########
##############################################################
if (rolling_window_analysis == T) {
  y_name <- "NO2_log"
  estim_window_start <- dmy("01/01/2018")
  estim_window_end <- dmy("31/01/2020")
  event_window_start <- dmy("01/02/2020")
  max_days <- 120
  Stats_TS <- matrix(data = NA, nrow = 4, ncol = 13)
  Stats_list <- vector(mode = "list", length = max_days)
  for (j in 1:max_days) {
    print(paste0("***** Day ",j," of ",max_days," started at", Sys.time()," *****"))
    event_window_end <- event_window_start + days(j)
    univHDGM_ES <- ES_estimation(Data = Data,
                                 estim_window_start,estim_window_end,event_window_start,event_window_end,
                                 y_name,y_hat_name = "NO2_hat_log_HDGM",date_name,idx_name,X_names,
                                 TV_X_names,
                                 model_type = NULL)
    regGAMM_ES <- ES_estimation(Data = Data,
                                estim_window_start,estim_window_end,event_window_start,event_window_end,
                                y_name,y_hat_name = "NO2_hat_log_GAMM",date_name,idx_name,X_names,
                                TV_X_names = TV_X_names,
                                model_type = NULL)
    regGAMMar1_ES <- ES_estimation(Data = Data,
                                   estim_window_start,estim_window_end,event_window_start,event_window_end,
                                   y_name,y_hat_name = "NO2_hat_log_GAMM_ar1",date_name,idx_name,X_names,
                                   TV_X_names = TV_X_names,
                                   model_type = NULL)
    regGAM_ES <- ES_estimation(Data = Data,
                               estim_window_start,estim_window_end,event_window_start,event_window_end,
                               y_name,y_hat_name = "NO2_hat_log_GAM",date_name,idx_name,X_names,
                               TV_X_names = TV_X_names,
                               model_type = NULL)
    Stats_TS[1,] <- c("HDGM",event_window_end,univHDGM_ES$P1,univHDGM_ES$P2,univHDGM_ES$Corrado_Tuckey,
                      univHDGM_ES$Z_patell_adj,univHDGM_ES$Z_BMP_adj,univHDGM_ES$T_grank,
                      univHDGM_ES$Z_grank_adj,univHDGM_ES$CumRank,univHDGM_ES$CumRank_mod,
                      univHDGM_ES$CumRank_T,univHDGM_ES$CumRank_Z_adj)
    Stats_TS[2,] <- c("GAMM",event_window_end,regGAMM_ES$P1,regGAMM_ES$P2,regGAMM_ES$Corrado_Tuckey,
                      regGAMM_ES$Z_patell_adj,regGAMM_ES$Z_BMP_adj,regGAMM_ES$T_grank,
                      regGAMM_ES$Z_grank_adj,regGAMM_ES$CumRank,regGAMM_ES$CumRank_mod,
                      regGAMM_ES$CumRank_T,regGAMM_ES$CumRank_Z_adj)
    Stats_TS[3,] <- c("GAMM with AR(1) errors",event_window_end,regGAMMar1_ES$P1,regGAMMar1_ES$P2,regGAMMar1_ES$Corrado_Tuckey,
                      regGAMMar1_ES$Z_patell_adj,regGAMMar1_ES$Z_BMP_adj,regGAMMar1_ES$T_grank,
                      regGAMMar1_ES$Z_grank_adj,regGAMMar1_ES$CumRank,regGAMMar1_ES$CumRank_mod,
                      regGAMMar1_ES$CumRank_T,regGAMMar1_ES$CumRank_Z_adj)
    Stats_TS[4,] <- c("GAM",event_window_end,regGAM_ES$P1,regGAM_ES$P2,regGAM_ES$Corrado_Tuckey,
                      regGAM_ES$Z_patell_adj,regGAM_ES$Z_BMP_adj,regGAM_ES$T_grank,
                      regGAM_ES$Z_grank_adj,regGAM_ES$CumRank,regGAM_ES$CumRank_mod,
                      regGAM_ES$CumRank_T,regGAM_ES$CumRank_Z_adj)
    Stats_TS <- data.frame(Stats_TS)
    colnames(Stats_TS) <- c("Model","Date","P1","P2","Corrado-Tuckey",
                            "Patell Z adj.","BMP Z adj.","Grank T",
                            "Grank Z adj.", "CumRank", "CumRank modified",
                            "CumRank T", "CumRank Z adj.")
    Stats_list[[j]] <- Stats_TS
    print(paste0("***** Day ",j," of ",max_days," completed at", Sys.time()," *****"))
  }
  
  ##### Save output
  save(Stats_list, file = paste0(out_path,"Output_rolling.RData"))
} else if (rolling_window_analysis == T) {
  load(paste0(out_path,"Output_rolling.RData"))
}

Stats_df <- bind_rows(Stats_list)

##### Plots of recursive windows
p_rolling <- Stats_df %>%
  mutate(across(c("Date","P1","P2","Corrado-Tuckey",
                  "Patell Z adj.","BMP Z adj.","Grank T",
                  "Grank Z adj.", "CumRank", "CumRank modified",
                  "CumRank T", "CumRank Z adj."), as.numeric),
         Date = as.Date.numeric(Date)) %>%
  pivot_longer(cols = c("P1","P2","Corrado-Tuckey",
                        "Patell Z adj.","BMP Z adj.","Grank T",
                        "Grank Z adj.", "CumRank", "CumRank modified",
                        "CumRank T", "CumRank Z adj."),names_to = "Stat",values_to = "Value") %>%
  mutate(Model = case_when(Model == "GAMM with AR(1) errors" ~ "GAMM with AR(1) dynamics",
                           TRUE ~ Model)) %>%
  ggplot(mapping = aes(x = Date, y = Value)) + 
  geom_line(mapping = aes(col = Model), linewidth = 1) + 
  facet_wrap(~ Stat, scales = "free") + 
  theme_classic() + 
  theme(title = element_text(size = 12,face = "bold"),
        plot.subtitle = element_text(size = 9),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Model", values = cols) + 
  scale_fill_manual("Period", values = cols) + 
  labs(y = "ES statistic",
       x = "First lockdown wave in 2020",
       title = "ES statistics cumulated across days",
       subtitle = "Statistics are computed using AC in log-scale") + 
  geom_rect(data = data.frame(Period  = factor(c("Lockdown",
                                               "Extreme weather event",
                                               "Extreme weather event",
                                               "Extreme weather event"),
                                             levels = c("Lockdown",
                                                        "Extreme weather event")),
                              start = c(as.Date("2020-03-08"),
                                        as.Date("2020-02-01"),
                                        as.Date("2020-02-12"),
                                        as.Date("2020-02-27")),
                              end   = c(as.Date("2020-05-18"),
                                        as.Date("2020-02-06"),
                                        as.Date("2020-02-15"),
                                        as.Date("2020-03-02"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  ggnewscale::new_scale_color() +
  geom_hline(data = data.frame(Crit = factor(x = c("10%","5%","1%"),
                                             levels = c("10%","5%","1%")),
                               Vals = c(-qnorm(p = 0.90),
                                        -qnorm(p = 0.95),
                                        -qnorm(p = 0.99))),
             mapping = aes(yintercept = Vals, col = Crit),
             linewidth = 1.2) + 
  scale_color_manual("Crits", values = cols[16:18])
ggpubr::ggexport(p_rolling,width = 1800, height = 1200, res = 150,
                 filename = "ES_rolling.png")



##################################################################
########## Abnormal concentrations in estimation window ##########
##################################################################
HDGM_AbnRes_Stats <- ES_AbnRes_Stats(ES_AbnRes = univHDGM_ES$epsilon_it,
                                     Window_start = estim_window_start, Window_end = estim_window_end,
                                     Plot_TS = T, Latex = T, Digits = 2)
regGAM_AbnRes_Stats <- ES_AbnRes_Stats(ES_AbnRes = regGAM_ES$epsilon_it,
                                       Window_start = estim_window_start, Window_end = estim_window_end,
                                       Plot_TS = T, Latex = T, Digits = 2)
regGAMM_AbnRes_Stats <- ES_AbnRes_Stats(ES_AbnRes = regGAMM_ES$epsilon_it,
                                        Window_start = estim_window_start, Window_end = estim_window_end,
                                        Plot_TS = T, Latex = T, Digits = 2)
regGAMMar1_AbnRes_Stats <- ES_AbnRes_Stats(ES_AbnRes = regGAMMar1_ES$epsilon_it,
                                           Window_start = estim_window_start, Window_end = estim_window_end,
                                           Plot_TS = T, Latex = T, Digits = 2)

AbnRes_Stats <- cbind(HDGM_AbnRes_Stats$mat$mean,
                      regGAM_AbnRes_Stats$mat$mean,
                      regGAMM_AbnRes_Stats$mat$mean,
                      regGAMMar1_AbnRes_Stats$mat$mean)
rownames(AbnRes_Stats) <- rownames(HDGM_AbnRes_Stats$mat)
colnames(AbnRes_Stats) <- c("HDGM","GAM","GAMM","GAMM-ar1")
tab_latex <- xtable::xtable(AbnRes_Stats,
                            caption = "c", 
                            align=c("l","c","c","c","c"))
print(tab_latex,file = paste0(out_path,"Avg_Stats_Estim.tex"))
tab_latex


######### BHL data from ARPA
PBL <- read_delim("zic_dd_dcast.csv", delim = ";", 
                  escape_double = FALSE, trim_ws = TRUE)
colnames(PBL) <- c("Date","Arconate","Caiolo","Cinisello PN", "Mantova - Lunetta 2", "Pieve S. Giacomo",
                   "Vertemate")
p_PBL <- PBL %>%
  mutate(Date = dmy(Date)) %>%
  filter(Date >= "2020-01-01" & Date <= "2020-05-31") %>%
  pivot_longer(cols = 2:6, names_to = "Station", values_to = "PBL") %>%
  ggplot(mapping = aes(x = Date, y = PBL)) + 
  geom_line() + 
  facet_wrap( ~ Station) + 
  labs(y = "",
       x = "",
       title = "Measured boundary layer height (BHL) at 5 specific monitoring site in Lombardy",
       subtitle = "Data provided by ARPA Lombardia") + 
  geom_rect(data = data.frame(Period  = factor(c("Train",
                                                 "Lockdown",
                                                 "Extreme weather event",
                                                 "Extreme weather event",
                                                 "Extreme weather event"),
                                               levels = c("Lockdown",
                                                          "Extreme weather event",
                                                          "Train")),
                              start = c(as.Date("2020-01-01"),
                                        as.Date("2020-03-08"),
                                        as.Date("2020-02-01"),
                                        as.Date("2020-02-12"),
                                        as.Date("2020-02-27")),
                              end   = c(as.Date("2020-01-31"),
                                        as.Date("2020-05-18"),
                                        as.Date("2020-02-06"),
                                        as.Date("2020-02-15"),
                                        as.Date("2020-03-02"))),
            inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = +Inf, fill = Period),
            alpha = 0.1) + 
  theme_classic() + 
  theme(title = element_text(size = 12,face = "bold"),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom")
ggpubr::ggexport(p_PBL,width = 1800, height = 1200, res = 150, filename = "TS_PBL.png")


########## PACF of the ACs within the estimation window
m_list <- list(regGAM_ES,regGAMMar1_ES,regGAMMar1_ES,univHDGM_ES)
m_names <- c("GAM","GAMM","GAMM with AR(1) dynamics","HDGM")
pacf_list <- vector(mode = "list", length = length(m_list))
for (m in 1:length(m_list)) {
  for (s in 1:dim(univHDGM_ES$epsilon_it)[2]) {
    pacf_s <- Pacf(m_list[[m]]$epsilon_it[1:761,s],lag.max = 30,na.action = na.pass)
    df_pacf <- data.frame(lag=pacf_s$lag,pacf=pacf_s$acf)
    df_pacf$IDStations <- s
    df_pacf$Model <- m_names[m]
    ifelse(s == 1,
           all_df_pacf <- df_pacf,
           all_df_pacf <- rbind(all_df_pacf,df_pacf)
    )
  }
  pacf_list[[m]] <- all_df_pacf
}
p_PACF <- ggplot(data = bind_rows(pacf_list), aes(x = as.factor(lag), y = pacf, col = Model)) +
  geom_boxplot()+
  ylim(-1,+1) + 
  facet_wrap(~ Model) + 
  geom_hline(yintercept = 0, linetype=2) +
  labs(x = "Temporal lag", y = "PACF", title = "Boxplot of PACF for ACs") + 
  theme_classic() + 
  theme(title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Model", values = cols)
ggpubr::ggexport(p_PACF,width = 1800, height = 1200, res = 150, filename = "Box_PACF.png")


############### ACF of the ACs within the estimation window
m_list <- list(regGAM_ES,regGAMMar1_ES,regGAMMar1_ES,univHDGM_ES)
m_names <- c("GAM","GAMM","GAMM with AR(1) dynamics","HDGM")
acf_list <- vector(mode = "list", length = length(m_list))
for (m in 1:length(m_list)) {
  for (s in 1:dim(univHDGM_ES$epsilon_it)[2]) {
    acf_s <- Acf(m_list[[m]]$epsilon_it[1:761,s],lag.max = 30,na.action = na.pass)
    df_acf <- data.frame(lag=acf_s$lag,acf=acf_s$acf)
    df_acf$IDStations <- s
    df_acf$Model <- m_names[m]
    ifelse(s == 1,
           all_df_acf <- df_pacf,
           all_df_acf <- rbind(all_df_acf,df_acf)
    )
  }
  acf_list[[m]] <- all_df_acf
}
p_ACF <- ggplot(data = bind_rows(acf_list), aes(x = as.factor(lag), y = acf, col = Model)) +
  geom_boxplot() +
  ylim(-1,+1) + 
  facet_wrap(~ Model) + 
  geom_hline(yintercept = 0, linetype=2) +
  labs(x = "Temporal lag", y = "ACF", title = "Boxplot of ACF for ACs") + 
  theme_classic() + 
  theme(title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Model", values = cols)
ggpubr::ggexport(p_ACF,width = 1800, height = 1200, res = 150, filename = "Box_ACF.png")


###### Linear correlation of the ACs within the estimation window
m_list <- list(regGAM_ES,regGAMMar1_ES,regGAMMar1_ES,univHDGM_ES)
m_names <- c("GAM","GAMM","GAMM with AR(1) dynamics","HDGM")
corr_list <- vector(mode = "list", length = length(m_list))
for (m in 1:length(m_list)) {
    corr_mat <- cor(m_list[[m]]$epsilon_it[1:761,], use = "na.or.complete")
    df_corr <- corr_mat %>%
      as_tibble() %>%
      pivot_longer(cols = 1:last_col(),names_to = "Station", values_to = "Corr") %>%
      filter(Corr < 1)
    df_corr$Model <- m_names[m]
    df_corr$IDStat <- as.numeric(substr(x = df_corr$Station,start = 5,stop = 6))
    corr_list[[m]] <- df_corr
}
p_lincor <- ggplot(data = bind_rows(corr_list), aes(x = as.factor(IDStat), y = Corr, col = Model)) +
  geom_boxplot() +
  ylim(-1,+1) + 
  facet_wrap(~ Model) + 
  geom_hline(yintercept = 0, linetype=2) +
  labs(x = "Station ID", y = "Pearson's linear correlation",
       title = "Boxplot of Pearson's linear correlation (pairwise) by station") + 
  theme_classic() + 
  theme(title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom") + 
  scale_color_manual("Model", values = cols)
ggpubr::ggexport(p_lincor,width = 1800, height = 1200, res = 150, filename = "Box_LinCorr.png")

library(ggpubr)
p_comb <- ggarrange(p_ACF,p_lincor,ncol = 2,common.legend = T,legend = "bottom")
p_comb <- annotate_figure(p = p_comb,
                          top = text_grob("Diagnostics for ACs within the estimation window",
                                          size = 12,face = "bold"))
ggexport(p_comb,width = 1800, height = 1200, res = 150, filename = "Box_Combined.png")


###### Linear correlation of the observed concentrations within the estimation window
NO2_long <- Data %>%
  select(Date,Stz_Code,NO2) %>%
  # pivot_wider(names_from = Stz_Code, values_from = NO2) %>%
  tsibble::as_tibble() %>%
  tsbox::ts_ts()
corr_mat <- cor(NO2_long[1:761,], use = "na.or.complete")
df_corr <- corr_mat %>%
  as_tibble() %>%
  pivot_longer(cols = 1:last_col(),names_to = "Station", values_to = "Corr") %>%
  filter(Corr < 1)
df_corr$IDStat <- as.numeric(substr(x = df_corr$Station,start = 5,stop = 6))
p_lincor_no2 <- ggplot(data = df_corr, aes(x = as.factor(IDStat), y = Corr)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype=2) +
  labs(x = "Station ID", y = "Pearson's linear correlation",
       title = "Boxplot of pairwise linear correlation by station") + 
  theme_classic() + 
  theme(title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom")

##### Linear correlation as function of distance
Stz_df <- Data %>%
  select(Stz_Code,Stz_Name,longitude = Stz_Long,latitude = Stz_Lat) %>%
  unique()
df_corr <- corr_mat %>%
  as_tibble() %>%
  pivot_longer(cols = 1:last_col(),names_to = "Station", values_to = "Corr")

Staz_geodist <- geodist::geodist(Stz_df)
colnames(Staz_geodist) <- rownames(Staz_geodist) <- colnames(corr_mat)
Staz_geodist_df <- Staz_geodist %>%
  as_tibble() %>%
  pivot_longer(cols = 1:last_col(),names_to = "Station", values_to = "Dist_m") %>%
  mutate(Dist_m = Dist_m/1000)

Staz_geodist_df <- bind_cols(Staz_geodist_df, df_corr[,-1])
Staz_geodist_df <- Staz_geodist_df %>%
  filter(Dist_m > 0)

p_dist_corr_no2 <- Staz_geodist_df %>%
  ggplot(mapping = aes(x = Dist_m, y = Corr)) +
  geom_point() + 
  geom_smooth(method = "gam") + 
  labs(x = "Distance (km)", y = "Pearson's linear correlation",
       title = "Linear correlation (pairwise) and distance") + 
  theme_classic() + 
  theme(title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom")

##### ACF
for (s in 1:dim(NO2_long)[2]) {
  acf_no2_s <- Acf(NO2_long[1:761,s],lag.max = 30,na.action = na.pass)
  df_acf_no2 <- data.frame(lag=acf_no2_s$lag,
                           acf=acf_no2_s$acf)
  df_acf_no2$IDStations <- s
  ifelse(s == 1,
         all_df_acf <- df_acf_no2,
         all_df_acf <- rbind(all_df_acf,df_acf_no2)
  )
}
p_ACF_no2 <- ggplot(data = all_df_acf, aes(x = as.factor(lag), y = acf)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype=2) +
  labs(x = "Temporal lag", y = "ACF", title = "Boxplot of ACF for ACs") + 
  theme_classic() + 
  theme(title = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.position = "bottom")



p_comb <- ggarrange(p_lincor_no2,p_dist_corr_no2,p_ACF_no2,ncol = 3,common.legend = T,legend = "bottom")
p_comb <- annotate_figure(p = p_comb,
                          top = text_grob("Spatio-temporal characteristics of observed NO2 concentrations",
                                          size = 12,face = "bold"))
ggexport(p_comb,width = 1800, height = 1200, res = 150, filename = "Box_Combined_no2.png")




########## Models performances
performance::performance(m_gammAR1)
performance::performance(m_gamm)
performance::performance(m_gam)
