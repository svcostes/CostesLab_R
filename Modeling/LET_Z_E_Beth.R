LET_Z_E_Beth <- function(E,Z) {
  # This funtions LET = LET_Z_E_Beth(E,Z) will compute the LET for energy
  #  E (MeV/n) and charge Z
  # Unit of LET is in eV/um
  #
  #
  # Author: Sylvain Costes, Feb 2012, Lawrence Berkeley National Laboratory
  m_nucleon = 931.5# %MeV
  uen = 31e8# % Mu_en/rho in um^2/kg for 1.25 MeV photons 
  Iexc = 73e-6# % Mean excitation potential in water in MeV
  c = 29.9792# % Speed of light in 10^9 cm/sec
  me = 0.511# % mass of electron in MeV/c2
  ne = 3.34e29# % electron density of water in m-3
  re = 2.818e-9# % classical electronic radius in um
  let_par = 16.83# %4.pi.me.c^2.re^2.ne in eV/um
  w_dens = 1e3# % water density, in kg/m3#
  ev = 1.602189e-19# % eV in J
  gamma = E/m_nucleon + 1
  beta = sqrt(gamma^2-1)/gamma
  Zeff= Z*(1.-exp(-125*sqrt(beta^2)*Z^(-2./3.)))
  LET = let_par/beta^2*Zeff^2*(log(2*me*beta^2/Iexc/(1-beta^2))-beta^2)#  % in eV/um
  return(LET)
}
