program fft2d
  use constants
  use myroutines
  use vproutines
  implicit none

  ! parameters
  integer, parameter :: n = 128, n4 = 4*n
  integer, parameter :: betaopt = 0
  real(rk), parameter :: rho = 0.5_rk

  ! local variables
  character(len=28) :: filename
  integer :: ntmp
  real(rk) :: a(n), beta, rhotmp
  complex(rk), dimension(4*n) :: eta, etah, zeta, z, zd, dthdnu

  ! a0val and c0val
  real(rk) :: u, c0
  complex(rk) :: a0(4*n), q0(4*n)

  call timestamp()

  filename = genfilename(n, rho, betaopt)
  call readdata(filename, ntmp, rhotmp, beta, a)
  u = a(1)
  print *, n4
  call ptval(n4, rho, beta, a, eta, etah, dthdnu, zeta, z, zd)
  call a0val(n4, u, z, zd, a0, q0)
  call c0val(a0, zd, c0)

  ! ! check
  ! do i = 1,n
  !    write(*, '(2es22.10)') a0(2*i)
  ! end do
  ! print *, size(a0)
  ! print *, c0

  call timestamp()
  stop
end program fft2d
