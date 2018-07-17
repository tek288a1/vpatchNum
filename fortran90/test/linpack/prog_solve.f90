program prog_solve
  use constants
  use myroutines
  implicit none

  ! declaration of parameters
  integer, parameter :: n_init = 4
  integer, parameter :: betaopt_init = 0
  real(rk), parameter :: rho_init = 0.05_rk ! "_rk" is important

  integer, parameter :: n = 128
  integer, parameter :: nstep = 2
  integer, parameter :: betaopt = 0
  real(rk), parameter :: rho_term = 0.3_rk
  real(rk), parameter :: drho = 0.05_rk
  !  include 'prog_solve_parameter.f90'
  ! n, betaopt, nstep, rho_init, rho_term, drho, file_init declared

  ! local variables
  integer :: i, n_old
  real(rk) :: rho, beta
  real(rk) :: a(n), a1(n)
  real(rk), dimension(:), allocatable :: rhov
  character(len=28) :: filename

  call timestamp()
  rhov = linspaceh(rho_init, rho_term, drho)
  filename = genfilename(n_init, rho_init, betaopt_init)
  call readdata(n, filename, n_old, rho, beta, a)
  do i = 2,size(rhov)
     call vsolver(n, rhov(i-1), rhov(i), nstep, &
          beta, betaopt, a, a1)
     filename = genfilename(n, rhov(i), betaopt)
     call writedata(filename, n, rhov(i), beta, a1)
     a = a1
  end do
  print *, ''
  call timestamp()
end program prog_solve
