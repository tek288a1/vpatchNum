program prog_solve
  use constants
  use myroutines
  use vproutines
  implicit none

  ! declaration of parameters
  integer, parameter :: n_init = 128
  integer, parameter :: betaopt_init = 1
  real(rk), parameter :: rho_init = 0.8_rk ! "_rk" is important

  integer, parameter :: n = 256
  integer, parameter :: nstep = 1
  integer, parameter :: betaopt = 1
  real(rk), parameter :: rho_term = 0.95_rk
  real(rk), parameter :: drho = 0.01_rk
  !  include 'prog_solve_parameter.f90'
  ! n, betaopt, nstep, rho_init, rho_term, drho, file_init declared

  ! local variables
  integer :: i, n_read
  real(rk) :: rho, beta
  real(rk) :: a_init(n_init), a(n), a1(n)
  real(rk), dimension(:), allocatable :: rhov
  character(len=50) :: filename

  call timestamp()
  rhov = linspaceh(rho_init, rho_term, drho)
  filename = genfilename(n_init, rho_init, betaopt_init)
  call readdata(filename, n_read, rho, beta, a_init)
  a = 0.0_rk
  a(1:n_read) = a_init
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
