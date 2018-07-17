program transmat
  ! This program is meant to test the routines calculating transfer
  ! matrix M relating positive and negative modes
  use constants
  use myroutines
  use vproutines
  use fftw
  implicit none
  ! parameters
  integer, parameter :: n = 128, nm1 = n-1, n4 = 4*n
  integer, parameter :: betaopt = 0
  real(rk), parameter :: rho = 0.4_rk
  ! local variables
  character(len=28) :: filename
  integer :: n_read, i
  real(rk) :: beta, rho_read, a(n), an(nm1)
  complex(rk), dimension(n4) :: eta, etah, zeta, z, zd, dthdnu
  ! fft
  integer(8) :: plan_forward
  complex(rk) :: in(n4), out(n4)
  real(rk) :: rout(n4)
  ! matrices
  real(rk) :: M(n, n)

  call timestamp()
  filename = genfilename(n, rho, betaopt)
  call readdata(filename, n_read, rho_read, beta, a)
  call ptval(n4, rho, beta, a, eta, etah, dthdnu, zeta, z, zd)
  call trmatval(n, rho, beta, M)

  open( unit=1, file='../tmp/adata.dat', status='replace')
  do i = 1,n
     write(1, '(es30.16E3)') a(i)
  end do
  close(1)

  open( unit=1, file='../tmp/mdata.dat', status='replace')
  do i = 1,n
     write(1, '(128es30.16E3)') M(i, :)
  end do
  close(1)

  write(*, *) beta

  ! fft check: replace z by za = z - i(zeta-rho)/(zeta+rho) and
  ! calculate its Fourier coefficients by fft; compare it against a.
  in = z - ii*(zeta-rho)/(zeta+rho)
  call dfftw_plan_dft_1d_ ( plan_forward, n4, in, out, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_forward )
  rout = aimag(out/real(n4, rk))

  ! positive coefficients: checked
  write(*, '(/a/, es12.5)') &
       'norm difference of positive modes: a data vs. fft data', &
       vnorm(rout(2:n)-a(2:n), 'inf')

  ! nonpositive coefficients
  an = -matmul( M(1:nm1, 1:nm1), a(2:n) )
  write(*, '(/a/, es12.5/)') &
       'norm difference of nonpositive modes: fft vs transfer matrix', &
       vnorm(an-rout((/1, (i, i=n4,n4-n+3,-1) /)), 'inf')
  call timestamp()
  stop
end program transmat
