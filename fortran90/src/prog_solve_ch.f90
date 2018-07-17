program prog_solve_ch
  use constants
  use myroutines
  use vproutines
  implicit none

  ! declaration of parameters
  integer, parameter :: n_init = 128
  integer, parameter :: betaopt_init = 0
  real(rk), parameter :: rho_init = 0.01_rk ! "_rk" is important
  integer, parameter :: n = 128
  integer, parameter :: nstep = 1
  integer, parameter :: betaopt = 0

  ! local variables
  integer :: i, n_read, ndp
  real(rk) :: rho, beta
  real(rk) :: a_init(n_init), a(n), a1(n)
  real(rk), dimension(:), allocatable :: rhov
  character(len=28) :: filename

  call timestamp()

  filename = genfilename(n_init, rho_init, betaopt_init)
  call readdata(filename, n_read, rho, beta, a_init)
  a = 0.0_rk
  a(1:n_read) = a_init
  open( unit=1, file='../data/rho_cheb.txt', status='unknown' )
  read(1,*) ndp
  allocate(rhov(0:ndp))
  rhov(0) = rho_init
  do i = 1,ndp
     read(1, *) rhov(i)
  end do
  close(1)

  do i = 1,ndp
     call vsolver(n, rhov(i-1), rhov(i), nstep, &
          beta, betaopt, a, a1)
     filename = genfilename(n, rhov(i), betaopt)
     call writedata(filename, n, rhov(i), beta, a1)
     a = a1
  end do
  print *, ''
  call timestamp()
end program prog_solve_ch
