ES_AbnRes_Stats <- function(ES_AbnRes,Window_start,Window_end,Latex = F, Digits = 2,
                            Plot_TS = T) {
  
  mean_cor <- function(cor_mat) {
    2*sum(abs(cor_mat[lower.tri(cor_mat,diag = F)]),na.rm=T) / (ncol(cor_mat)*(ncol(cor_mat)-1)) 
  }
  
  res <- ES_AbnRes
  
  # Filtering window of interest
  res <- window(res, start = Window_start, end = Window_end)
  dates <- zoo::index(res)

  ##### ACF
  phi <- matrix(NA,nrow = 84,ncol = 5)
  phi[,1] <- round(apply(res,2,FUN = function(x) acf(imputeTS::na_mean(x), lag.max = 21)$acf[1+1]),3)
  phi[,2] <- round(apply(res,2,FUN = function(x) acf(imputeTS::na_mean(x), lag.max = 21)$acf[3+1]),3)
  phi[,3] <- round(apply(res,2,FUN = function(x) acf(imputeTS::na_mean(x), lag.max = 21)$acf[7+1]),3)
  phi[,4] <- round(apply(res,2,FUN = function(x) acf(imputeTS::na_mean(x), lag.max = 21)$acf[14+1]),3)
  phi[,5] <- round(apply(res,2,FUN = function(x) acf(imputeTS::na_mean(x), lag.max = 21)$acf[21+1]),3)
  
  ##### Linear (Person) correlation
  rho <- cor(res, use = "na.or.complete")
  
  ##### Seasonality (Fs) and Trend (Ft) strength and outlier percentage
  Fs <- Ft <- out_perc <- numeric(dim(res)[2])
  for (s in 1:dim(res)[2]) {
    res1 <- tsbox::ts_ts(res[,s])
    out_perc[s] <- length(forecast::tsoutliers(res1)$index)
    # res1 <- imputeTS::na_interpolation(x = res1, option = "spline")
    res1 <- imputeTS::na_mean(x = res1)
    # res1 <- imputeTS::na_kalman(x = res1)
    ##### Estrarre con extracts di Matteo
    dec_add <- decompose(x = res1,type="additive")
    Fs[s] <- max(0,1-var(dec_add$random,na.rm = T)/var(dec_add$random+dec_add$seasonal,na.rm = T))
    Ft[s] <- max(0,1-var(dec_add$random,na.rm = T)/var(dec_add$random+dec_add$trend,na.rm = T))
  }
  
  ##### Skewness
  sk <- apply(res,2,FUN = moments::skewness, na.rm=T)
  
  ##### Kurtosis
  k <- apply(res,2,FUN = moments::kurtosis, na.rm=T)
  
  ##### Standard deviation
  sigma <- apply(res,2,FUN = sd, na.rm=T)
  
  ##### Average
  m <- apply(res,2,FUN = mean, na.rm=T)
  
  ##### Table of statistics
  mat <- matrix(NA,13,6)
  mat[1,] <- c(min(rho),quantile(rho,0.25),mean_cor(rho),quantile(rho,0.5),quantile(rho,0.75),max(rho[lower.tri(rho,diag = F)]))
  mat[2,] <- c(min(m),quantile(m,0.25),mean(m),quantile(m,0.5),quantile(m,0.75),max(m))
  mat[3,] <- c(min(sigma),quantile(sigma,0.25),mean(sigma),quantile(sigma,0.5),quantile(sigma,0.75),max(sigma))
  mat[4,] <- c(min(sk[sk>-19]),quantile(sk[sk>-19],0.25),mean(sk[sk>-19]),quantile(sk[sk>-19],0.5),quantile(sk[sk>-19],0.75),max(sk[sk>-19]))
  mat[5,] <- c(min(k),quantile(k,0.25),mean(k),quantile(k,0.5),quantile(k,0.75),max(k))
  mat[6:10,] <- c(apply(phi,2,min),
                  apply(phi,2,function(x) quantile(x,0.25)),
                  apply(phi,2,mean),
                  apply(phi,2,function(x) quantile(x,0.50)),
                  apply(phi,2,function(x) quantile(x,0.75)),
                  apply(phi,2,max))
  mat[11,] <- c(min(Fs),quantile(Fs,0.25),mean(Fs),quantile(Fs,0.5),quantile(Fs,0.75),max(Fs))
  mat[12,] <- c(min(Ft),quantile(Ft,0.25),mean(Ft),quantile(Ft,0.5),quantile(Ft,0.75),max(Ft))
  mat[13,] <- c(min(out_perc),quantile(out_perc,0.25),mean(out_perc),quantile(out_perc,0.50),quantile(out_perc,0.75),max(out_perc))
  colnames(mat) <- c("min","25p","mean","median","75p","max")
  rownames(mat) <- c("$rho$","$mu$","$sigma$","Skewness","Kurtosis","$phi_1$","$phi_3$","$phi_7$","$phi_{14}$","$phi_{21}$","Fs","Ft","Outliers")
  mat <- round(mat,Digits)
  
  mat <- data.frame(mat)
  
  ##### TS plot
  if (Plot_TS == T) {
    hlinemean <- mean(res, na.rm=T)
    hlinemean <- xts(x = rep(hlinemean,dim(window(res, start = Window_start, end = Window_end))[1]),
                     order.by = zoo::index(res))
    p <- plot(res, main = "Abnormal residuals (concentrations) from regAR(1) for NO2 in Lombardy")
    lines(hlinemean, col = "yellow", lwd=7)
    print(p)
  }
  
  ##### Latex table
  if (Latex == T) {
    tab_latex <- xtable::xtable(mat,
                                caption = "c", 
                                align=c("l","c","c","c","c","c","c"))
    tab_latex
  }
  
  ##### Output
  return(list(mat = mat, tab_latex = tab_latex))
}
