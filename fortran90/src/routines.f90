!-----------------------------------------------------------------------------80
! reading data from a file
!-----------------------------------------------------------------------------80
subroutine readdata(filename, n, rho, beta, a)
  use constants, only: rk
  implicit none
  ! arguments declaration
  character(len=*) :: filename
  integer :: n
  real(rk) :: rho
  real(rk) :: beta
  real(rk) :: a(n)
  ! local variables
  integer :: i, ist             ! ist: I/O status
  real(rk) :: tmp

  ! reading in data: rho, beta, n_old, and a = (u, a_1, a_2, ..., a_nm1)
  open( unit=1, file=filename, status='unknown' )
  read(1,*) n
  read(1,*) rho
  read(1,*) beta
  do i = 1,n          ! truncate if n < n_old; pad with 0 if n > n_old
     read(1, *, iostat=ist) tmp
     if (ist == 0) then
        a(i) = tmp
     else                       ! ist /= 0 : no more lines in the file
        a(i) = 0.0_rk
     end if
  end do
  close(1)

end subroutine readdata
!-----------------------------------------------------------------------------80
! writing data to a file
!-----------------------------------------------------------------------------80
subroutine writedata(filename, n, rho, beta, a)
  use constants, only: rk
  ! argument declaration
  character(len=*), intent(in) :: filename
  integer, intent(in) :: n
  real(rk), intent(in) :: rho
  real(rk), intent(in) :: beta
  real(rk), intent(in) :: a(n)
  ! local variables
  integer :: i
  open( unit=1, file=filename, status='unknown' )
  write(1, '(i24)') n
  write(1, '(f24.4)') rho
  write(1, '(f24.4)') beta
  do i=1,n
     write(1, '(es24.16)') a(i) ! es: format in scientific notation
  end do
  close(1)
end subroutine writedata
!-----------------------------------------------------------------------------80
! point evaluation subroutine
!-----------------------------------------------------------------------------80
subroutine ptval(npt, la, rho, beta, a, eta, etah, zeta, z, zd)
  use constants
  use myroutines
  implicit none
  ! input: npt, la, rho, beta, a
  ! output: eta, zeta, z, zd
  ! aux var: nu, eta, zeta, etah, dth_dnu, fm, fmd, al, ah, ald, ahd
  integer, intent(in) :: npt ! number of uniformly spaced-out points on unit circle
  integer, intent(in) :: la  ! length of a = 1 (U_0) + ncoef
  real(rk), intent(inout) :: a(la)
  complex(rk), dimension(npt), intent(inout) :: eta, etah, zeta, z, zd
  integer :: j, ncoef
  real(rk) :: rho, beta, tmp
  complex(rk), dimension(npt) :: nu, dth_dnu, fm, fmd
  complex(rk), dimension(npt) :: al, ald, ah, ahd

  ncoef = la - 1
  ! preliminary calculations
  nu = linspace1(0.0_rk, 2*pi, npt)
  eta = exp(ii*nu)
  zeta = (eta - beta)/(1.0_rk - beta*eta)
  etah = (rho**2 + beta*zeta)/(zeta + beta*rho**2)
  dth_dnu =( zeta*(rho**2)*(beta**2 - 1.0_rk) )/ &
       ( (rho**2 + beta*zeta)*(zeta + beta*rho**2) )
  dth_dnu = dth_dnu*eta*(1.0_rk - beta**2)/ &
       ( (1.0_rk - beta*eta)*(eta - beta))

  ! mobius part
  fm = (zeta - rho)/(zeta + rho)
  fmd = 2.0_rk*rho*zeta/(zeta + rho)**2 * &
       (1.0_rk - beta**2)*eta/( zeta*(1.0_rk - beta*eta)**2 )

  ! analytical part, alpha(eta) and alpha_hat(eta)
  al = a(ncoef+sh)           ! this sets all elements to be the constant
  ah = a(ncoef+sh)
  ald = ncoef*al
  ahd = ncoef*ah
  tmp = a(0+sh)
  a(0+sh) = 0.0_rk
  do j = ncoef-1,0,-1
     al = al*eta + a(j+sh)
     ah = ah*etah + a(j+sh)
     ald = ald*eta + j*a(j+sh)
     ahd = ahd*etah + j*a(j+sh)
  end do
  a(0+sh) = tmp
  z = ii*(fm + al - ah)
  zd = -(fmd + ald - ahd*dth_dnu)
end subroutine ptval

!-----------------------------------------------------------------------------80
! velocity calculation (quadrature routine): deprecated
!-----------------------------------------------------------------------------80
!
! A_0(\nu) is calculated for \nu between 0 and \pi (strictly upper
! half eta-plane) and is meant to be used in residual
! calculation. This routine is now deprecated; use =a0val= instead.
!
subroutine a0val_old(n, u, z, zd, a0, q0)
  use constants
  use myroutines
  implicit none
  ! input: u, z, zd, n
  ! output: a0, q0
  integer :: n, i
  real(rk) :: u
  complex(rk) :: z(4*n), zd(4*n), a0(n), q0(n)
  complex(rk) :: zc(n), zq(2*n), zdq(2*n)

  zc = z((/ (i, i=2,2*n,2) /))
  zq = z((/ (i, i=1,4*n,2) /))
  zdq = zd((/ (i, i=1,4*n,2) /))
  do i = 1,n
     a0(i) = sum( (abs(zc(i))**2 - abs(zq)**2)/(zc(i)**2 - zq**2)*zdq )
     a0(i) = a0(i)/real(2*n, rk)
     q0(i) = 1.0_rk/(u + a0(i))
  end do
end subroutine a0val_old
!-----------------------------------------------------------------------------80
! velocity calculation (quadrature routine): current, updated
!-----------------------------------------------------------------------------80
!
! This routine calculated A_0(nu) and q_0(nu) using at points
! corresponding to the input z values using alternating trapezoidal
! quadrature method.
!
!   A_0(\nu) = 1/(2\pi) \int_{0}^{2\pi}
! 	( |z(\nu)|^2 - |z(\nu')|^2 )/( z(\nu)^2 - z(\nu')^2 ) z'(\nu) d\nu
!
!   q_0(\nu) = 1/( U_0 + A_0(\nu))
!
subroutine a0val(n2, u, z, zd, a0, q0)
  use constants
  use myroutines
  implicit none
  ! argument declaration
  integer :: n2
  real(rk) :: u
  complex(rk) :: z(n2), zd(n2), a0(n2), q0(n2)
  ! local variables
  integer :: n, i, j
  complex(rk) :: zc, zq(n2/2), zdq(n2/2)

  n = int(n2/2)                 ! n2 must be even
  do i = 1,n2
     zc = z(i)
     if (mod(i,2)==1) then      ! odd colloc / even quad
        zq = z((/ (j, j=2,n2,2) /))
        zdq = zd((/ (j, j=2,n2,2) /))
     else                       ! even colloc / odd quad
        zq = z((/ (j, j=1,n2,2) /))
        zdq = zd((/ (j, j=1,n2,2) /))
     end if
     a0(i) = sum( (abs(zc)**2 - abs(zq)**2)/(zc**2 - zq**2)*zdq )
     a0(i) = a0(i)/real(n, rk)
  end do
  q0 = 1/(u + a0)
end subroutine a0val
!-----------------------------------------------------------------------------80
! residual calculation
!-----------------------------------------------------------------------------80
function resval(n, rho, beta, a)
  use constants
  implicit none
  ! input: rho, beta, a, n
  ! output: resval
  integer :: n, i, iv(n)
  real(rk) :: rho, beta, a(n), resval(n), u
  complex(rk), dimension(4*n) :: a0, q0
  complex(rk), dimension(4*n) :: eta, etah, zeta, z, zd
  call ptval(4*n, n, rho, beta, a, eta, etah, zeta, z, zd)
  u = a(1)
  iv = (/ (i, i=2,2*n,2) /)
  call a0val(4*n, u, z, zd, a0, q0)
  resval = aimag( (a0(iv)+u)*zd(iv) )
end function resval
!-----------------------------------------------------------------------------80
! vpatch solver
!-----------------------------------------------------------------------------80
subroutine vsolver &
     (n, rho_init, rho_term, nstep, beta_term, betaopt, a_init, a_term)
  use constants
  use myroutines
  implicit none
  ! input: rho_init, rho_term, nstep, betaopt, a_init, n
  ! output: beta_term, a_term

  ! parameters
  real(rk), parameter :: tol_res = 1e-10   ! residual
  real(rk), parameter :: tol_dsol = 1e-10  ! change in solution
  real(rk), parameter :: h = 1e-8     ! forward difference step size
  integer, parameter :: max_iter = 10 ! maximum number of iteration

  ! variable declaration
  integer :: nstep, betaopt, iter, j, k, n
  real(rk) :: rho, rho_init, rho_term, beta, beta_term
  real(rk) :: resmax, dsolmax
  real(rk), dimension(0:nstep+1) :: rho_vec
  real(rk), dimension(n) :: a_init, a_term, a, anew, res, respert, resnew
  real(rk) :: id(n,n), jac(n,n)

  ! for lapack routines
  integer :: lda, lwork, ldb
  integer :: ipiv(n), iwork(n), info
  real(rk) :: anorm, rcond, work(4*n)

  ! ! for linpack routines
  ! integer :: ipvt(n)            ! pivot indices
  ! integer :: job = 0
  ! real(rk) :: rcond
  ! real(rk) :: z(n)

  interface
     function resval(n, rho, beta, a)
       use constants
       integer :: n
       real(rk) :: rho, beta, a(n), resval(n)
     end function resval
  end interface

  rho_vec = (/ linspace(rho_init, rho_term, nstep+1), rho_term /)
  a = a_init
  do k = 1,nstep+1           ! continuation loop, note k starts from 2
     rho = rho_vec(k)
     if (betaopt==0) then
        beta = 0.0_rk
     else if (betaopt==1) then
        beta = (1-sqrt(1-rho**2))/rho
     else if (betaopt==2) then
        beta = (1-2*sqrt(1-rho**2))/rho
     end if

     ! initial iterate
     iter = 0
     res = resval(n, rho, beta, a)

     ! print to terminal
     write(*, '(/ a8, f6.4, a12, f6.4, a9, i4)') & ! / for a new line
          'rho = ', rho, 'beta = ', beta, 'n = ', n
     write(*, '(2x, a43)') line
     write(*, '(a6, a13, a13, a13)') &
          'iter', '|res|', '|dsol|', 'rcond'

     ! preallocation
     id = eye(n)
     jac = zeros(n,n)

     ! Newton iteration: this is how to write an infinite loop
     do
        do j = 1,n ! Jacobian calculation
           respert = resval(n, rho, beta, a+h*id(:,j))
           jac(:, j) = (respert - res)/h
        end do

        ! ! linear algebra routines (LINPACK) ----------------------------
        ! call DGECO(jac, n, n, ipvt, rcond, z)
        ! ! LU-factors JAC and estimates RCOND; IPVT is the pivot
        ! ! indices; z is a work vector.
        ! call DGESL(jac, n, n, ipvt, res, job)
        ! ! Solves JAC*X = RES; on return, RES is the solution; JOB = 0
        ! ! for non-transposed problem

        ! linear algebra routines (LAPACK) -------------------------------
        lda = n; lwork = 4*n; ldb = n; anorm = mnorm(n, n, jac)
        ! calculating infinity norm of matrix JAC
        call DGETRF(n, n, jac, lda, ipiv, info)
        call DGECON('O', n, jac, lda, anorm, rcond, work, iwork, info)
        call DGETRS('N', n, 1, jac, lda, ipiv, res, ldb, info)

        ! iterative formula
        iter = iter + 1
        anew = a - res
        resnew = resval(n, rho, beta, anew)
        resmax = vnorm(resnew, 'inf')
        dsolmax = vnorm(a-anew, 'inf')

        ! print to terminal
        write(*, '(i6, es13.4, es13.4, es13.4)') &
             iter, resmax, dsolmax, rcond

        ! stopping criteria
        if ((resmax <= tol_res).and.(dsolmax <= tol_dsol)) then
           a = anew
           write(*, '(a45, /)') &
                '  ******** Solution has been found. *********'
           exit
        else if (iter == max_iter) then
           write(*, '(a45, /)') &
                '  ********** No solution is found. **********'
           stop
        end if

        ! update iterates for next iteration
        a = anew
        res = resnew
     end do
  end do
  beta_term = beta
  a_term = a
end subroutine vsolver
