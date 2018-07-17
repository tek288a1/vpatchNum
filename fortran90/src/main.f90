program main
  use constants
  use myroutines
  implicit none
  integer, parameter :: n=128, betaopt=0
  character(len=*), parameter :: in_file='../data/vp_b0_n02_r0500.dat'
  character(len=*), parameter :: out_file='../data/vp_end.dat'
  integer :: nstep, n_old
  real(rk) :: rho_init, rho_term, beta_init, beta_term
  real(rk), dimension(n) :: a_init, a_term

  nstep = 9
  rho_term = 0.5_rk
  call timestamp()
  call readdata (n, in_file, n_old, rho_init, beta_init, a_init)
  call vsolver(n, rho_init, rho_term, nstep, &
       beta_term, betaopt, a_init, a_term)
  call writedata(out_file, n, rho_term, beta_term, a_term)
  write(*, '(a)') ''
  call timestamp()
  print *, linspaceh(0.0_rk, 1.0_rk, 0.1_rk)
  stop
end program main
