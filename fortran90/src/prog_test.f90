program test
  use constants
  use myroutines
  use vproutines
  implicit none

  ! parameters
  integer, parameter :: n_input = 256
  integer, parameter :: n = n_input
  integer, parameter :: n4 = 4*n
  integer, parameter :: n2 = 2*n
  integer, parameter :: nm1 = n-1
  integer, parameter :: n2m1 = n2-1
  integer, parameter :: nj =200
  integer, parameter :: betaopt = 1
  real(rk), parameter :: rho = 0.85_rk

  ! local variables
  character(len=50) :: filename
  integer :: n_read
  integer :: i, j, k, m
  integer :: ind(n2)
  real(rk) :: a(n), beta, rho_read
  complex(rk), dimension(n4) :: eta, etah, zeta, z, zd
  complex(rk), dimension(n4) :: jt, dthdnu, omega0
  ! integer :: plan_backward
  ! integer, dimension(n) :: io = (/ (i,i=1,n2,2) /)

  ! calculated quantities
  real(rk), dimension(n,n) :: trm
  complex(rk) :: c0
  complex(rk), dimension(n4) :: a0, q0, lu
  complex(rk), dimension(n2,n2) :: j1, j2, j3, h4
  complex(rk), dimension(n2,n2) :: a1, a2, a3, b4
  ! complex(rk), dimension(-n:nm1, -n:nm1) :: t0c
  ! complex(rk), dimension(-n2:n2m1,-n2:n2m1) :: a1c, a2c, a3c, b4c

  ! allocatable arrays
  ! complex(rk), allocatable :: tmp(:)
  complex(rk), allocatable :: t_hat(:,:), t1(:,:)
  complex(rk), allocatable :: a1jm(:,:), ajm(:,:)
  ! Note: passing an array with non-standard indexing, i.e., starting
  ! from an integer other than 1, requires allocatability.
  allocate(t_hat(-n:nm1,-n:nm1))
  allocate(t1(-n:nm1,-n:nm1))
  allocate(a1jm(0:nm1,1:n))
  allocate(ajm(1:nm1,1:nm1))

  call timestamp()
  filename = genfilename(n_input, rho, betaopt)
  call readdata(filename, n_read, rho_read, beta, a)
  call ptval(n4, rho, beta, a, eta, etah, dthdnu, zeta, z, zd)
  call a0val(n4, a(1), z, zd, a0, q0)
  call c0val(a0, zd, c0)
  call trmval(n, rho, beta, trm)
  call jtval(nj, rho, zeta, jt)
  omega0 = q0*dthdnu*(c0 - a(1)*jt) ! u = a(1)
  lu = q0*(jt*dthdnu + omega0)
  call j1val_sh(eta, z, zd, omega0, j1, a1)
  call j2val_sh(eta, z, zd, omega0, j2, a2)
  call j3val_sh(eta, z, zd, omega0, j3, a3)
  call h4val_sh(z, zd, h4, b4)
  call thval(a1, a2, a3, t_hat)
  call t1val(a1, a2, a3, b4, dthdnu, t1)
  call a1jmval(q0, t1, trm, a1jm)
  call ajmval(q0, t1, trm, lu, ajm)

  open(unit=1, file='../tmp/Ajm.dat', status='unknown')
  write(1, '(2es26.16)') ajm
  close(1)

  deallocate(t_hat, t1, a1jm, ajm)
  write(*, '(a)') ''
  call timestamp()

  stop
end program test
