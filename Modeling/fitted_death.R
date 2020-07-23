fitted_death <- function(cell_line,LET) {
  # This funtions fittted_death <- function(cell_line,LET) will compute gam1 and gam2 for cell death
  # As a function of cell_type and LET
  #
  # Author: Sylvain Costes, Jul 2020
  
  if (cell_line == 'HF19') {
    if (LET<=160) {
      gam1 = 1.186e-3*LET}
    else { gam1 = 0.2}
    gam2 = 0.024
  }
  
  if (cell_line == 'V79') {
    gam1 = 0.017
    gam2 = 0.01
  }
  
  return(c(gam1,gam2))
}