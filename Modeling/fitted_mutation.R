fitted_mutation_V4 <- function(cell_line,LET) {
  # This funtions fittted_mutation <- function(cell_line,LET) will compute lg1 and gam2 for cell death
  # As a function of cell_type and LET
  #
  # Author: Sylvain Costes, Jul 2020
  
  if (cell_line == 'HF19') {
      if (LET<=90) {
        lg1 = 7.94e-7*exp(0.0350*LET)
      }
      else if (LET < 160) {lg1=-2.85e-7*(LET-90) + 1.99e-5}
      else if (LET < 200) {lg1=4.99e-8*(LET-160) + 1e-9}
        else {lg1 = 2e-6}

      if (LET<=90) {lg2=1e-9}
      else if (LET<110) {lg2 = 6.29e-8*(LET-90) + 1e-9}
      else if (LET<200) {lg2 = -1.26e-4 + 1.56e-6*LET - 3.66e-9*LET^2}
      else {lg2 = 3.98e-5}
  }
  
  if (cell_line == 'V79') {
    lg1 = 8.80e-8 + (5.11e-7-8.80e-8)*(1-exp(-2.73e-2*LET))
    lg2 = 6.31e-8 + (4.78e-6-6.31e-8)*(1-exp(-5.88e-3*LET))
  }
  
  return(c(lg1,lg2))
}