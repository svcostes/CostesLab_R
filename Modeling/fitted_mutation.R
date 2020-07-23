fitted_mutation <- function(cell_line,LET) {
  # This funtions fittted_mutation <- function(cell_line,LET) will compute lg1 and gam2 for cell death
  # As a function of cell_type and LET
  #
  # Author: Sylvain Costes, Jul 2020
  
  if (cell_line == 'HF19') {
    if (LET>50) {
      lg2 = 6.8e-6}
    else { lg2 = 1.8e-7}
    lg1 = 5.141e-8*LET+3.743e-7
  }
  
  if (cell_line == 'V79') {
    if (LET>110) {
      lg1 = 2.7e-7}
    else { lg1 = 9.349e-10*LET+8.199e-8}
    lg2 = 2.240e-7*LET+2.837e-7
  }
  
  return(c(lg1,lg2))
}