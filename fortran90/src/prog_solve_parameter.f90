  integer, parameter :: n_init = 512
  integer, parameter :: betaopt_init = 2
  real(rk), parameter :: rho_init = 0.99_rk ! "_rk" is important

  integer, parameter :: n = 512
  integer, parameter :: nstep = 2
  integer, parameter :: betaopt = 2
  real(rk), parameter :: rho_term = 0.995_rk
  real(rk), parameter :: drho = 0.001_rk
