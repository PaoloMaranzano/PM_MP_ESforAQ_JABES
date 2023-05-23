ES_Stats_Extract <- function(ES_estim_out,Digits = 2) {
  
  if (F) {
    ES_estim_out <- reg_ES
  }
  
  ##### Output matrix
  mat <- matrix(NA, nrow = 18, ncol = 3)
  
  ##### Statitics extraction
  stats <- c(ES_estim_out$P1,
             ES_estim_out$P2,
             ES_estim_out$cross_T_test,
             ES_estim_out$crude_dep_T_test,
             ES_estim_out$T_test_skew,
             ES_estim_out$Z_patell,
             ES_estim_out$Z_patell_adj,
             ES_estim_out$Z_BMP,
             ES_estim_out$Z_BMP_adj,
             ES_estim_out$T_grank,
             ES_estim_out$Z_grank,
             ES_estim_out$Z_grank_adj,
             ES_estim_out$CumRank,
             ES_estim_out$CumRank_mod,
             ES_estim_out$CumRank_T,
             ES_estim_out$CumRank_Z,
             ES_estim_out$CumRank_Z_adj,
             ES_estim_out$Corrado_Tuckey)
  
  ##### Output
  mat[,1] <- c("P1","P2","Corrado_Tuckey",
               "cross_T_test","crude_dep_T_test","T_test_skew",
               "Z_patell","Z_patell_adj","Z_BMP","Z_BMP_adj",
               "T_grank","Z_grank","Z_grank_adj",
               "CumRank","CumRank_mod","CumRank_T","CumRank_Z","CumRank_Z_adj")
  mat[,2] <- c("Adj","Adj","Adj",
               "Unadj","Unadj","Unadj",
               "Unadj","Adj","Unadj","Adj",
               "Adj","Unadj","Adj",
               "Adj","Adj","Adj","Unadj","Adj")
  mat[,3] <- round(as.numeric(stats),Digits)
  colnames(mat) <- c("ES_stat","CD","Value")
  rownames(mat) <- mat[,1]
  mat <- as.data.frame(mat)
  
  ##### Output
  return(list(mat = mat))
}
