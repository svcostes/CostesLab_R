read_C_simulation_Ianik <- function(file_path) {
  # This funtion reads the RIFs.dat file in current directory file_path and extract
  # 1000 simulations num_DSb array made of DSB, DSB image, Ncomb, RIF, pmut and psurv
  # The two last values are derived from the gamma terms used for death and mutation
  # probability
  # 
  #
  # Author: Sylvain Costes, September 2019, NASA
  
  require("arrayhelpers")
  
  tmp = read.delim(file_path,sep='\t',header=TRUE)
  singlesim_df = data.frame(matrix(ncol = 5, nrow = 1000)) # initialize
  colnames(singlesim_df) = c('dsb','ncomb','rif','psurv','pmut')
  singlesim_df$dsb = tmp[,5]
  singlesim_df$ncomb = tmp[,7]
  singlesim_df$rif = tmp[,6]
  singlesim_df$psurv = tmp[,8]
  singlesim_df$pmut = tmp[,11]
  return(singlesim_df)
}
