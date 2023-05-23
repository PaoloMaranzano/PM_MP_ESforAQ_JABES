ES_estimation <- function(Data,estim_window_start,estim_window_end,event_window_start,
                          event_window_end,y_name,y_hat_name,date_name,idx_name,X_names,TV_X_names,
                          model_type,max.p = 7,max.q = 7, p = 0, d = 0, q = 0) {
  
  if (F) {
    setwd("H:/Il mio Drive/Spatiotemporal_ES/Code/ES_univHDGM")
    Data <- readr::read_csv("HDGM_output_1608.csv")
    estim_window_start <- lubridate::dmy("01/01/2018")
    estim_window_end <- lubridate::dmy("31/01/2020")
    event_window_start <- lubridate::dmy("01/02/2020")
    event_window_end <- lubridate::dmy("31/05/2020")
    y_name <- "NO2"
    y_hat_name <- "NO2_hat_HDGM"
    # y_hat_name <- NULL
    date_name <- "Date"
    idx_name <- "Stz_Code"
    X_names <- c("weekend", "Temperature", "Rainfall", "Pressure", "WindU", "WindV", 
                 "VegetationHigh", "VegetationLow", "GeopotHeight", "RelHumid",
                 "SolarRadiation",
                 "Humidex",
                 "Lag",
                 "MetrArea","Mountain","UrbPlain",
                 "Ind","Traf","Rural")
    X_names <- colnames(Data)[grepl(pattern = paste(X_names, collapse = "|"), x = colnames(Data))]
    TV_X_names <- colnames(Data)[grepl(pattern = paste(c("Temperature", "Rainfall", "Pressure",
                                                         "WindU", "WindV", "RelHumid","SolarRadiation",
                                                         "Humidex"), collapse = "|"), x = colnames(Data))]
    model_type <- NULL
  }
  
  
  
  #################################################
  ########## Data management and reshape ##########
  #################################################
  
  ##### Filtering dates
  Data <- Data %>%
    rename(Y = y_name,
           Y_hat = y_hat_name,
           Date = date_name,
           Idx = idx_name) %>%
    filter(Date >= as_date(estim_window_start) & Date <= as_date(event_window_end))
  
  ##### Dimensions
  Y_wide <- Data %>%
    select(Date,Idx,Y) %>%
    pivot_wider(names_from = Idx, values_from = Y)
  # m = number of cross-sectional units
  m <- dim(Y_wide)[2]-1
  dates_estim <- Y_wide %>%
    filter(Date >= as_date(estim_window_start) & Date <= as_date(estim_window_end)) %>%
    pull(Date)
  # n0 = estimation window length
  n0 <- length(dates_estim)
  dates_event <- Y_wide %>%
    filter(Date >= as_date(event_window_start) & Date <= as_date(event_window_end)) %>%
    pull(Date)
  # n1 = event window length
  n1 <- length(dates_event)
  # n = complete window length
  n <- n0 + n1
  # Check on the right number of stocks and market indices
  if (n != nrow(Y_wide)) {
      stop("Length of the complete window is different the dates sequence length")
  }
  # K = number of estimated regression coefficients (intercept excluded)
  if (is.null(X_names)) {
    K <- 0
  } else {
    K <- length(X_names)
  }
  
  ##### Dates vector
  dates <- Y_wide$Date
  
  
  
  #######################################################################
  ########## Modelling step and abnormal residuals computation ##########
  ##########     (cross-sectional units specific)              ##########
  #######################################################################
  
  ##### cross-sectional variances
  sigma_i <- numeric(length = m)
  
  ##### Regression design matrices X and Y
  X_it <- X_it_estim <- X_it_event <- Y_it <- Y_it_estim <- vector("list", length = m)
  
  ##### Abnormal residuals matrix
  epsilon_it <- Y_wide
  epsilon_it <- epsilon_it %>%
    column_to_rownames("Date")
  
  ##### Modelling output
  aux_reg_list <- vector("list", length = m)
  
  ##### Abornal residuals computation
  for (i in 1:m) {
    ### Extract data associated with the i-th cross-sectional unit
    Idx_cod <- colnames(epsilon_it)[i]
    # Complete window observations
    X_it[[i]] <- Data %>%
      filter(Idx == Idx_cod) %>%
      select(X_names)
    Y_it[[i]] <- Data %>%
      filter(Idx == Idx_cod) %>%
      select(Y)
    X_it[[i]] <- xts::xts(X_it[[i]],order.by = dates)
    Y_it[[i]] <- xts::xts(Y_it[[i]],order.by = dates)
    # Estimation window observations
    X_it_estim[[i]] <- Data %>%
      filter(Idx == Idx_cod,
             Date >= as_date(estim_window_start) & Date <= as_date(estim_window_end)) %>%
      select(X_names)
    X_it_estim[[i]] <- xts::xts(X_it_estim[[i]],order.by = dates_estim)
    Y_it_estim[[i]] <- window(Y_it[[i]],start = estim_window_start, end = estim_window_end)
    # Event window observations
    X_it_event[[i]] <- Data %>%
      filter(Idx == Idx_cod,
             Date >= as_date(event_window_start) & Date <= as_date(event_window_end)) %>%
      select(X_names)
    X_it_event[[i]] <- xts(X_it_event[[i]],order.by = dates_event)
    
    ### Modelling and abnormal residuals computation
    if (!is.null(y_hat_name)) {
      # Case 1: AR from external model (need Y and Y_hat columns)
      print(paste0("Computing abnormal residuals using external model on index ",i," of ", m))
      epsilon_it[,i] <- Data %>%
        filter(Idx == Idx_cod) %>%
        mutate(epsilon_t = Y - Y_hat) %>%
        pull(epsilon_t)
      epsilon_it <- xts(epsilon_it,order.by = dates)
      sigma_i[i] <- sqrt(var(epsilon_it[,i],na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- "External model used: epsilon_it from external model/source"
    } else if ((is.null(model_type) | model_type == "lm") & is.null(y_hat_name) & is.null(X_names)) {
      # Case 2: linear regression without covariates (demeaning-centering)
      print(paste0("Estimating linear regression without covariates on index ",i," of ", m))
      aux_reg <- forecast::Arima(y = Y_it_estim[[i]],
                                 order = c(0,0,0), method = "CSS-ML")
      epsilon_it[,i] <- forecast::Arima(model = aux_reg,
                                        y = Y_it[[i]])$residuals
      epsilon_it <- xts::xts(x = epsilon_it, order.by = dates)
      sigma_i[i] <- sqrt(var(aux_reg$residuals,na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- aux_reg
    } else if (model_type == "lm" & is.null(y_hat_name) & !is.null(X_names)) {
      # Case 3: linear regression with covariates
      print(paste0("Estimating linear regression on index ",i," of ", m))
      X_it_estim_timevary <- X_it_estim[[i]][,TV_X_names]
      X_it_timevary <- X_it[[i]][,TV_X_names]
      aux_reg <- forecast::Arima(y = Y_it_estim[[i]],
                                 xreg = X_it_estim_timevary,
                                 order = c(0,0,0), method = "CSS-ML")
      epsilon_it[,i] <- forecast::Arima(model = aux_reg,
                                        xreg = X_it_timevary,
                                        y = Y_it[[i]])$residuals
      sigma_i[i] <- sqrt(var(epsilon_it[,i],na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- aux_reg
    } else if (model_type == "auto.arima" & is.null(y_hat_name) & !is.null(X_names)) {
      # Case 4: linear regression with ARIMA errors (regARIMA)
      print(paste0("Estimating linear regression with ARIMA errors (regARIMA) on index ",i," of ", m))
      X_it_estim_timevary <- X_it_estim[[i]][,TV_X_names]
      X_it_timevary <- X_it[[i]][,TV_X_names]
      aux_reg <- forecast::auto.arima(y = Y_it_estim[[i]],
                                      xreg = X_it_estim_timevary,
                                      max.p = max.p,max.q = max.q, method = "ML", ic="aicc")
      epsilon_it[,i] <- forecast::Arima(model = aux_reg,
                                        xreg = X_it_timevary,
                                        y = Y_it[[i]])$residuals
      sigma_i[i] <- sqrt(var(epsilon_it[,i],na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- aux_reg
    } else if (model_type == "auto.arima" & is.null(y_hat_name) & is.null(X_names)) {
      # Case 5: ARIMA model with optimal (AICc) order
      print(paste0("Estimating ARIMA model on index ",i," of ", m))
      aux_reg <- forecast::auto.arima(y = Y_it_estim[[i]],
                                      max.p,max.q, method = "ML", ic="aicc")
      epsilon_it[,i] <- forecast::Arima(model = aux_reg,
                                        y = Y_it[[i]])$residuals
      sigma_i[i] <- sqrt(var(epsilon_it[,i],na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- aux_reg
    } else if (model_type == "arima" & is.null(y_hat_name) & !is.null(X_names)) {
      # Case 6: linear regression with covariates and specified ARIMA order
      print(paste0("Estimating linear regression with ARIMA(",p,",",d,",",q,") on index ",i," of ", m))
      X_it_estim_timevary <- X_it_estim[[i]][,TV_X_names]
      X_it_timevary <- X_it[[i]][,TV_X_names]
      aux_reg <- forecast::Arima(y = Y_it_estim[[i]],
                                 xreg = X_it_estim_timevary,
                                 order = c(p,d,q), method = "ML")
      epsilon_it[,i] <- forecast::Arima(model = aux_reg,
                                        xreg = X_it_timevary,
                                        y = Y_it[[i]])$residuals
      sigma_i[i] <- sqrt(var(epsilon_it[,i],na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- aux_reg
    } else if (model_type == "arima" & is.null(y_hat_name) & is.null(X_names)) {
      # Case 7: ARIMA model with specified order
      print(paste0("Estimating ARIMA(",p,",",d,",",q,") on index ",i," of ", m))
      aux_reg <- forecast::Arima(y = Y_it_estim[[i]],
                                 order = c(p,d,q), method = "ML")
      epsilon_it[,i] <- forecast::Arima(model = aux_reg,
                                        y = Y_it[[i]])$residuals
      sigma_i[i] <- sqrt(var(epsilon_it[,i],na.rm = T)*(n0-1)/(n0-2))
      aux_reg_list[[i]] <- aux_reg
    }
  }
  epsilon_it <- xts::xts(x = epsilon_it, order.by = dates)
  
  
  
  #########################################################
  ########## Event studies parametric statistics ##########
  #########################################################
  
  ##### Parametric t-test for single stock (not-adjuted for cross-correlation)
  # for testing H0: AR_it = 0
  s2_AR_i <- sigma_i^2
  t_AR_it <- sweep(epsilon_it,2,s2_AR_i,'/')
  # for testing H0: CAR_i = 0
  CAR_i <- colSums(window(epsilon_it,start = event_window_start, end = event_window_end),na.rm=T)
  s2_CAR <- n1 * s2_AR_i
  t_CAR_i <- CAR_i / sqrt(s2_CAR)
  
  
  ##### Cross-sectional t-statistic for multiple stocks (not-adjuted for cross-correlation)
  # for testing H0: AAR=0
  AAR_t <- xts::xts(rowMeans(epsilon_it), order.by = dates)
  s2_AAR_t <- 1 / (m - 1) * rowSums(sweep(epsilon_it,1,AAR_t,'-')^2,na.rm = T)
  t_AAR_t <- sqrt(m) * AAR_t / sqrt(s2_AAR_t)
  # for testing H0: CAAR=0
  CAR_i <- colSums(window(epsilon_it,start = event_window_start, end = event_window_end),na.rm=T)
  CAAR <- mean(CAR_i,na.rm=T)
  s2_CAAR <- 1 / (m - 1) * sum((CAR_i - CAAR)^2, na.rm=T)
  t_CAAR <- CAAR / sqrt(s2_CAAR / m)
  
  
  ##### Time-Series Standard Deviation or Crude Dependence Test (not-adjuted for cross-correlation)
  AAR_t <- xts::xts(rowMeans(epsilon_it,na.rm = T), order.by = dates)
  AAR_bar <- 1 / n0 * sum(window(x = AAR_t, start = estim_window_start, end = estim_window_end),na.rm=T)
  s2_AAR <- 1 / (n0 - 2) * sum((AAR_t - AAR_bar)^2, na.rm=T)
  # H0: AAR=0
  T_AAR_t <- sqrt(m) * AAR_t / sqrt(s2_AAR)
  # H0: CAAR=0
  T_CAAR <- CAAR / (sqrt(n1) * sqrt(s2_AAR))
  
  
  ##### Patell Z test (Patell, 1976 & 1979 or Kolari & Pynnonen, 2010 eq. 12)
  ### Not adjusted for cross-correlation
  # for testing H0: AAR = 0
  # !!! when n0 becomes large, the second and third terms tend to zero and the correction tends to 1 --> s2_AR_it = s2_AR_i
  SAR_it <- epsilon_it
  for (i in 1:m) {
    # Select on time-varying covariates (otherwise singular matrix)
    X_it_estim_timevary <- X_it_estim[[i]][,TV_X_names]
    X_it_event_timevary <- X_it_event[[i]][,TV_X_names]
    # X_it_event_timevary <- X_it_event[[i]][,apply(X_it_event[[i]], 2, var, na.rm=TRUE) != 0]
    # Correcting forecast errors (epsilon_it in event window) for their unbiased variance 
    # See Pelagatti & Maranzano (QJFA, 2021): sigma_i * ( 1 + Z'*inv(X'X)*Z)
    #    where Z is the design matrix in the event window
    #    where X is the design matrix in the estimation window
    FrcErr_corr_i <- c(1 + diag(X_it_estim_timevary %*%
                                  solve(t(X_it_estim_timevary) %*% X_it_estim_timevary) %*%
                                  t(X_it_estim_timevary)),
                       1 + 1/n0 + diag(X_it_event_timevary %*%
                                         solve(t(X_it_estim_timevary) %*% X_it_estim_timevary) %*%
                                         t(X_it_event_timevary)))
    s2_AR_it <- s2_AR_i[i] * FrcErr_corr_i
    SAR_it[,i] <- epsilon_it[,i] / sqrt(s2_AR_it)
  }
  ASAR_t <- rowSums(SAR_it, na.rm=T)
  s2_ASAR_t <- (n0 - K - 1) / (n0 - K - 3) * m    # when n0 becomes large, (n0-2)/(n0-4) tends to 1 --> s2_ASAR_t tends to m
  z_patell_t <- xts::xts(ASAR_t / sqrt(s2_ASAR_t), order.by = dates)
  # for testing H0: CAAR = 0
  CSAR_i <- colSums(window(SAR_it,start=event_window_start,end=event_window_end),na.rm = T)
  s2_CSAR_i <- n1 * (n0 - K - 1) / (n0 - K - 3)  # when n0 becomes large, (n0-2)/(n0-4) tends to 1 --> s2_CSAR_i tends to n1
  z_patell <- sum(CSAR_i / sqrt(m * s2_CSAR_i))
  ### Adjusted for cross-correlation Patell Z (Kolari & Pynnonen, 2010 eq. 13)
  # propose a modification to the Patell-test to account for cross-correlation of the abnormal returns.
  # It uses the standardized abnormal returns (SAR_it) and defining  r_bar the average of the sample cross-correlation
  # of the estimation period abnormal returns.
  # If the correlation r_bar is zero, the adjusted test statistic reduces to the original Patell test statistic.
  # The test statistic for testing H0: AAR = 0
  r_ij <- cor(window(epsilon_it,start = estim_window_start, end = estim_window_end), method = "pearson",
              use = "complete.obs")
  r_bar <- 2*sum(abs(r_ij[lower.tri(r_ij,diag = F)])) / (m*(m-1))
  z_adj_patell_t <- z_patell_t / sqrt(1 + (m-1)*r_bar)
  # Assuming the square-root rule holds for the standard deviation of different return periods, this test can be
  # used when considering Cumulated Abnormal Returns
  # The test statistic for testing H0: CAAR = 0
  z_adj_patell <- z_patell / sqrt(1 + (m-1)*r_bar)
  
  
  ##### Standardized Cross-Sectional test or BMP test
  # Ref: Boehmer, Musumeci and Poulsen, 1991
  # Ref: Kolari & Pynnonen, 2010 (eq. 6)
  ### Not adjusted for cross-correlation
  # They propose a standardized cross-sectional method which is robust to the variance induced by the event.
  # Test statistics on day t for testing H0: AAR = 0
  s2_BMP_ASAR_t <- 1 / (m -1) * rowSums((SAR_it - rowMeans(SAR_it,na.rm = T))^2, na.rm = T)
  s2_BMP_ASAR_t_adj <- s2_BMP_ASAR_t * (1 + (m-1)*r_bar) / (1 - r_bar)
  z_BMP_t <- xts::xts(ASAR_t / sqrt(m * s2_BMP_ASAR_t),order.by = dates)
  # for testing H0: CAAR = 0
  # The Mikkelson and Partch (1988) correction adjusts for each firm the test statistic for serial correlation in the returns.
  # The correction terms are
  # s2_CAR_i <- s2_AR_i * n1                  # Market Adjusted Model: sum of n1 independent variables
  # s2_CAR_i <- s2_AR_i * (n1 + n1^2 / n0)    # Comparison Period Mean Adjusted Model
  s2_CAR_i <- sigma_i
  for (i in 1:m) {
    # Select on time-varying covariates (otherwise singular matrix)
    X_it_estim_timevary <- X_it_estim[[i]][,apply(X_it_estim[[i]], 2, var, na.rm=TRUE) != 0]
    X_it_event_timevary <- X_it_event[[i]][,names(X_it_estim_timevary)]
    # X_it_event_timevary <- X_it_event[[i]][,apply(X_it_event[[i]], 2, var, na.rm=TRUE) != 0]
    MikkPatch_corr_i <- (n1 + n1^2/n0 + sum(diag(X_it_event_timevary %*% 
                                                   solve(t(X_it_estim_timevary) %*% X_it_estim_timevary) %*% 
                                                   t(X_it_event_timevary))))
    s2_CAR_i[i] <- s2_AR_i[i] * MikkPatch_corr_i
  }
  SCAR_i <- CAR_i / sqrt(s2_CAR_i)
  SCAR_bar  <- mean(SCAR_i, na.rm=T)
  s2_SCAR <- 1 / (m-1) * sum((SCAR_i - SCAR_bar)^2, na.rm=T)
  SCAR_i_prime <- SCAR_i / sqrt(s2_SCAR)
  GSAR_it <- rbind(as.matrix(window(SAR_it, start = estim_window_start, end = estim_window_end)),
                   SCAR_i_prime)
  s2_SCAR_bar <- s2_SCAR / m
  z_BMP <- SCAR_bar / sqrt(s2_SCAR_bar)
  ### Adjusted for cross-correlation BMP test (Kolari & Pynnonen, 2010 eq. 11)
  # Propose a modification to the BMP-test to account for cross-correlation of the abnormal returns.
  # Using the standardized abnormal returns (SAR_it) defined as in the previous section, and defining r_bar as
  # the average of the sample cross-correlation of the estimation period abnormal returns,
  # The test statistic for testing H0: AAR = 0
  z_adj_BMP_t <- ASAR_t / sqrt(m * s2_BMP_ASAR_t_adj)   # = z_BMP_t * sqrt((1-r_bar) / (1 + (m-1)*r_bar))
  # If the correlation r_bar is zero, the adjusted test statistic reduces to the original BMP test statistic.
  # Assuming the square-root rule holds for the standard deviation of different return periods, this test can
  # be used when considering Cumulated Abnormal Returns:
  # The test statistic for testing H0: CAAR = 0
  s2_SCAR_bar_adj <- s2_SCAR_bar * (1 + (m-1)*r_bar) / ( 1 - r_bar)
  z_adj_BMP <- SCAR_bar / sqrt(s2_SCAR_bar_adj)   # = z_BMP * sqrt((1-r_bar) / (1 + (m-1)*r_bar))
  
  
  ##### Skewness Corrected Test
  # The skewness-adjusted t-test, introduced by Hall 1992, corrects the cross-sectional t-test for skewed abnormal
  # return distribution. This test is applicable for averaged abnormal return (H0:AAR=0),
  # the cumulative averaged abnormal return (H0:CAAR=0), and the averaged buy-and-hold abnormal return (H0:ABHAR=0).
  # In the following, we are limited by the situation of cumulative averaged abnormal returns.
  # First, let's revisit the cross-sectional standard deviation (unbiased by sample size):
  s2_CAAR <- 1 / (m - 1) * sum((CAR_i - CAAR)^2, na.rm=T)
  # The skewness estimation (unbiased by sample size) is given by:
  gamma <- m / ((m-2)*(m-1)) * sum(((CAR_i - CAAR)^3)/sqrt(s2_CAAR)^3, na.rm=T)
  S <- CAAR / sqrt(s2_CAAR)
  # then the skewness adjusted test statistic for CAAR is given by
  t_skew <- sqrt(m) * (S + 1/3*gamma*S^2 + 1/27*gamma^2*S^3 + 1/(6*m)*gamma)
  # which is asymptotically standard normal distributed.
  # For a further discussion on skewness transformation we refer to Hall (1992) and for further discussion
  # on unbiased estimation of the second and third moment we refer to Cramer (1961) or Rimoldini (2013).
  
  
  
  #############################################################
  ########## Event studies non-parametric statistics ##########
  #############################################################
  
  ##### CumRank T and Z statistics --> refer to eq. (5)-(23) Luoma, 2010
  s2_SAR_t <- xts::xts(apply(SAR_it, 1, FUN = var,na.rm=T), order.by = dates)
  # SAR_it_prime <- rbind(window(SAR_it,start = estim_window_start, end = estim_window_end),
  #                     window(sweep(SAR_it,1,sqrt(s2_SAR_t),"/"),start = event_window_start, end = event_window_end))
  SAR_it_prime <- SAR_it
  R_it <- xts::xts(Rfast::colRanks(SAR_it_prime), order.by = zoo::index(SAR_it_prime))
  S_ievent <- window(R_it,start=event_window_start,end=event_window_end)    # (10)
  K_it <- R_it / (n+1)        # (14)
  U_it <- S_ievent / (n+1)    # (20)-(21)
  U_ievent <- colSums(window(K_it,start=event_window_start,end=event_window_end),na.rm=T)
  K_bar_t <- xts::xts(rowMeans(K_it,na.rm=T),order.by = zoo::index(SAR_it_prime))
  U_bar_event <- sum(window(K_bar_t,start=event_window_start,end=event_window_end),na.rm=T)   # U_bar_event = mean(U_ievent)
  s2_K_bar <- mean((K_bar_t - 1/2)^2,na.rm=T)
  ## Cumrank-Z with cross-sectional independence ~ N(0,1)
  #   !! References: Luoma, 2010 (30)-(32)
  #   !! This statistic uses the theoretical variance of VAR(U_bar_event) without correction for cross-correlation
  s2_U_bar_event <- (n1 * (n - n1)) / (12 * (n+1) * m)
  CumRank_Z <- (U_bar_event - n1/2) / sqrt(s2_U_bar_event)
  ## Cumrank-Z with cross-sectional dependence ~ N(0,1)
  #   !! References: Luoma, 2010 (30)-(32)
  #   !! This statistic uses the theoretical variance of VAR(U_bar_event) correcting for cross-correlation
  #   !! Uguale (ordine dei millesimi) alla statistica MC nello script di Pelagatti: il motivo Ã¨ che MC Ã¨ la versione
  #         asintotica di CumRank_Z, perchÃ¨ quando n -> Inf allora n0/(n+1) -> 1
  # rho_rank_ij <- cor(window(SAR_it_prime,start = estim_window_start, end = estim_window_end), method = "spearman")
  rho_rank_ij <- cor(window(K_it,start = estim_window_start, end = estim_window_end),
                     use = "complete.obs", method = "pearson")
  rho_rank_bar <- 2*sum(abs(rho_rank_ij[lower.tri(rho_rank_ij,diag = F)]),na.rm=T) / (m*(m-1))
  s2_U_bar_event_cross_corr <- (n1 * (n - n1)) / (12 * (n+1) * m) * (1 + (m-1)*rho_rank_bar)
  CumRank_Z_adj <- (U_bar_event - n1/2) / sqrt(s2_U_bar_event_cross_corr)
  ## Campbell & Wasley (1993) - with implicit cross-sectional dependence and biased variance ~ N(0,1)
  #   !! References: eq. (37)-(38) Kolari Pynnonnen, 2011
  #   !! References: eq. (35)-(36) Luoma, 2010
  #   !! s2_tilde_event is a biased estimator of the population variance VAR(U_bar_event) and implicitly accounts for the cross-corr
  #   !! This statistic is the multi-period version of Corrado (1989) Corrado-Zyvney test (1992) with biased variance
  #            if n1 = 1 --> t_rank_t
  #   !! Uguale (ordine dei millesimi) alla statistica MC nello script di Pelagatti
  s2_tilde_event <- n1 * s2_K_bar
  CumRank <- (U_bar_event - n1/2) / sqrt(s2_tilde_event)
  ## Modified Campbell & Wasley (1993) - with implicit cross-sectional dependence and unbiased variance ~ N(0,1)
  #   !! References: Luoma, 2010 (41-Z3) ... Kolari Pynnonen, 2011 (41)
  #   !! s2_hat_event is a unbiased estimator of the population variance VAR(U_bar_event) and implicitly accounts for the cross-corr
  #   !! This statistic is the multi-period version of Corrado-Zyvney test (1992) with unbiased variance if n1 = 1 --> t_rank_t
  #   !! Uguale (ordine dei millesimi) alla statistica MC nello script di Pelagatti
  s2_hat_event <- s2_tilde_event * (n - n1) / (n -1)
  CumRank_mod <- (U_bar_event - n1/2) / sqrt(s2_hat_event)
  # CumRank-T - Adjusted for cross-correlation ~ t_{n-2}
  #   !! References: Luoma, 2010 (42-Z4) ... Kolari Pynnonen, 2011 (42)
  #   !! This statsitic is more robust against cross-sectional correlation w.r.t to Cumrank_Z, because the variance estimator
  #   !! s2_hat_event implicitly accounts the possible cross-correlation.
  CumRank_T <- CumRank_mod * sqrt((n-2)/(n-1-CumRank_mod^2))
  
  
  ##### Generalized rank statistics --> eq. (3)-(23) Kolari Pynnonnen 2011
  CAR_ievent <- CAR_i
  s2_CAR_ievent <- s2_CAR_i
  SCAR_ievent <- CAR_i / sqrt(s2_CAR_ievent)
  s2_SCAR <- 1 / (m-1) * sum((SCAR_ievent - mean(SCAR_ievent))^2,na.rm=T)   # = s2_SCAR_bar
  SCAR_ievent_star <- SCAR_ievent / sqrt(s2_SCAR)
  # GSAR_it <- xts::xts(rbind(as.matrix(window(SAR_it, start = estim_window_start, end = estim_window_end)),
  #                           SCAR_ievent_star), order.by = dates[1:(n0+1)])
  GSAR_it <- xts::xts(rbind(as.matrix(window(SAR_it, start = estim_window_start, end = estim_window_end)),
                            SCAR_ievent),order.by = dates[1:(n0+1)])
  U_it <- K_it <- xts::xts(Rfast::colRanks(GSAR_it) / (n0+1+1) - 1/2,
                           order.by = zoo::index(GSAR_it))   # (9)
  U_bar_t <- rowMeans(U_it,na.rm=T)   # K_bar_t <- rowMeans(K_it)
  # grank-T: Kolari, 2011 (12)
  s2_U_bar <- mean(U_bar_t^2)   # s2_K_bar <- mean(K_bar_t^2)
  U_bar_zero <- U_bar_t[n0+1]   # K_bar_zero <- K_bar_t[n0+1]
  Z <- U_bar_t[n0+1] / sqrt(s2_U_bar)   # K_bar_t[n0+1] / sqrt(s2_K_bar)
  t_grank <- Z * sqrt((n0-1) / (n0-1+1 - Z^2))    # ~ t_(n0+1-2)
  # grank-Z with cross-sectional independence: Kolari, 2011 (20)
  s2_grank_z <- n0 / (12 * m * (n0+2))
  z_grank <- U_bar_zero / sqrt(s2_grank_z)   # ~ N(0,1)
  # grank-Z test with cross-sectional dependence: Kolari Pynnonnen, 2011 (21)-(23)
  rho_rank_ij <- cor(window(K_it,start = estim_window_start, end = estim_window_end),
                     method = "pearson", use = "complete.obs")
  rho_rank_bar <- 2*sum(abs(rho_rank_ij[lower.tri(rho_rank_ij,diag = F)]),na.rm=T) / (m*(m-1))
  z_grank_crosscorr <- U_bar_zero / sqrt(s2_grank_z*(1+(m-1)*rho_rank_bar))
  
  
  ##### Multivariate Corrado with Tuckey's correction for small samples and with cross-sectional dependence (2011)
  ## References: Corrado (2011) eq. (10)
  # 4.91 is the asymptotic variance of V, where V = K^0.14 - (1-K)^0.14 --> VAR(V) = 4.91
  # sum(rrcor) = m*(1 + (m-1)*rho_bar_K),where rho_bar_K = cor(V)
  R_it <- xts::xts(Rfast::colRanks(SAR_it), order.by = zoo::index(SAR_it))
  # R_it <- xts::xts(Rfast::colRanks(SAR_it_prime), order.by = zoo::index(SAR_it_prime))
  # mU <- (R_it - 1/2)/n
  K_it <- R_it / (n+1)
  rrcor <- cor(window((K_it^0.14 - (1-K_it)^0.14), start = estim_window_start, end = estim_window_end),
               use = "complete.obs",method = "pearson")
  s2_CT <- n1 * sum(abs(rrcor),na.rm=T)
  Corrado_Tuckey_adj <-  4.91 / sqrt(s2_CT) * sum(window((K_it^0.14 - (1-K_it)^0.14),
                                                         start = event_window_start,
                                                         end = event_window_end),na.rm=T)
  
  
  ##### Pelagatti 1
  e_t <- xts::xts(rowSums(SAR_it,na.rm = T), order.by = zoo::index(SAR_it_prime))
  # e_t <- xts::xts(rowSums(SAR_it_prime), order.by = zoo::index(SAR_it_prime))
  R_t <- xts::xts(Rfast::colRanks(e_t), order.by = zoo::index(SAR_it_prime))
  s2_P1 <- n0*n1/n
  P1 <- sum(qnorm((window(R_t, start = event_window_start, end = event_window_end) - 0.50) / n),na.rm=T) / sqrt(s2_P1)   # ~ N(0,1)
  # vU <- (rank(rowSums(e_t)) - 0.5)/n
  # P1 <- sum(qnorm(vU[(n0+1):(n0+n1)]))*sqrt(n/(n0*n1))
  
  
  ##### Pelagatti 2
  R_it <- xts::xts(Rfast::colRanks(SAR_it), order.by = zoo::index(SAR_it_prime))
  # R_it <- xts::xts(Rfast::colRanks(SAR_it_prime), order.by = zoo::index(SAR_it_prime))
  K_it <- R_it / (n+1)
  # rho_ij <- cor(qnorm(window(K_it, start = estim_window_start, end = estim_window_end)))
  # s2_P2 <- n1 * sum(abs(rho_ij))
  rho_ij <- cor(qnorm(window(K_it, start = estim_window_start, end = estim_window_end)),
                use = "complete.obs",method = "pearson")
  rho_ij <- cor(qnorm(window(K_it, start = event_window_start, end = event_window_end)),
                use = "complete.obs",method = "pearson")
  rho_rank_bar <- 2*sum(abs(rho_ij[lower.tri(rho_ij,diag = F)]),na.rm=T) / (m*(m-1))
  s2_P2 <- n1*m*(1+(m-1)*rho_rank_bar)
  P2 <- sum(rowSums(qnorm(window(K_it, start = event_window_start, end = event_window_end)))) / sqrt(s2_P2)  # ~ N(0,1)
  # mU <- (R_it - 1/2)/n
  # mPhi <- qnorm(mU)
  # phicov <- cov(mPhi[1:n0,])
  # s2_P2 <- n1 * sum(phicov)
  # P2 <- sum(window(mPhi, start = event_window_start, end = event_window_end)) / sqrt(s2_P2)
  
  return(list(epsilon_it = epsilon_it,
              SAR_it = SAR_it,
              aux_reg_list = aux_reg_list,
              P1 = P1, P2 = P2,
              T_test_t = t_AR_it,
              T_test = t_CAR_i,
              cross_T_test_t = t_AAR_t,
              cross_T_test = t_CAAR,
              crude_dep_T_test_t = T_AAR_t,
              crude_dep_T_test = T_CAAR,
              T_test_skew = t_skew,
              Z_patell_t = z_patell_t,
              Z_patell = z_patell,
              Z_adj_patell_t = z_adj_patell_t,
              Z_patell_adj = z_adj_patell,
              Z_BMP_t = z_BMP_t,
              Z_BMP = z_BMP,
              Z_adj_BMP_t = z_adj_BMP_t,
              Z_BMP_adj = z_adj_BMP,
              T_grank = t_grank,  Z_grank = z_grank, Z_grank_adj = z_grank_crosscorr,
              CumRank = CumRank, CumRank_mod = CumRank_mod,
              CumRank_T = CumRank_T, CumRank_Z = CumRank_Z, CumRank_Z_adj = CumRank_Z_adj,
              Corrado_Tuckey = Corrado_Tuckey_adj))
  
  
}
