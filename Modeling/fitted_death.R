fitted_death_V4 <- function(cell_line,LET) {
  # This funtions fittted_death <- function(cell_line,LET) will compute gam1 and gam2 for cell death
  # As a function of cell_type and LET
  #
  # Author: Sylvain Costes, Jul 2020
  
  if (cell_line == 'HF19') {
     gam1 = 0.020 + (0.179-0.020)*(1-exp(-9.83e-3*LET))
     gam2 = 2.28e-2*exp(-0.5*((LET-133)/51.7)^2)
  }
  
  if (cell_line == 'V79') {
    gam1 = 0.01 + (0.0930-0.01)*(1-exp(-2.25e-3*LET))
    gam2 = 0.005
  }
  
  return(c(gam1,gam2))
}