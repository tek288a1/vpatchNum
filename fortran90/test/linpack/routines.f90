!
! reading data from a file
!
subroutine readdata(n_new, filename, n_old, rho, beta, a)
  use constants, only: rk
  implicit none
  ! arguments declaration
  integer, intent(in) :: n_new
  character(len=*), intent(in) :: filename
  integer, intent(out) :: n_old
  real(rk), intent(out) :: rho
  real(rk), intent(out) :: beta
  real(rk), intent(out) :: a(n_new)
  ! local variables
  integer :: i, ist             ! ist: I/O status
  real(rk) :: tmp

  ! reading in data: rho, beta, n_old, and a = (u, a_1, a_2, ..., a_nm1)
  open( unit=1, file=filename, status='unknown' )
  read(1,*) n_old
  read(1,*) rho
  read(1,*) beta
  do i = 1,n_new          ! truncate if n < n_old; pad with 0 if n > n_old
     read(1, *, iostat=ist) tmp
     if (ist == 0) then
        a(i) = tmp
     else                       ! ist /= 0 : no more lines in the file
        a(i) = 0.0_rk
     end if
  end do
  close(1)
end subroutine readdata
!
! writing data to a file
!
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
!
! constructing data file name
!

!
! point evaluation subroutine
!
subroutine ptval(n, rho, beta, a, eta, zeta, z, zd)
  use constants
  use myroutines
  implicit none
  ! input: rho, beta, a, n
  ! output: eta, zeta, z, zd
  ! aux var: nu, eta, zeta, etah, dth_dnu, fm, fmd, al, ah, ald, ahd
  integer, intent(in) :: n
  real(rk), dimension(n), intent(inout) :: a
  complex(rk), dimension(4*n), intent(out) :: eta, zeta, z, zd
  integer :: j, n4, nm1
  real(rk) :: rho, beta, tmp
  complex(rk), dimension(4*n) :: nu, etah, dth_dnu, fm, fmd
  complex(rk), dimension(4*n) :: al, ald, ah, ahd

  n4 = 4*n
  nm1 = n-1
  ! preliminary calculations
  nu = linspace1(0.0_rk, 2*pi, n4)
  eta = exp(ii*nu)
  zeta = (eta - beta)/(1.0_rk - beta*eta)
  etah = (rho**2 + beta*zeta)/(zeta + beta*rho**2)
  dth_dnu =( zeta*(rho**2)*(beta**2 - 1.0_rk) )/ &
       ( (rho**2 + beta*zeta)*(zeta + beta*rho**2) )
  dth_dnu = dth_dnu*eta*(1.0_rk - beta**2)/ &
       ( (1.0_rk - beta*eta)*(eta - beta))

  ! mobius part
  fm = (zeta - rho)/(zeta + rho)
  fmd = 2*rho*zeta/(zeta + rho)**2 * &
       (1.0_rk - beta**2)*eta/( zeta*(1.0_rk - beta*eta)**2 )

  ! analytical part, alpha(eta) and alpha_hat(eta)
  al = a(nm1+sh)           ! this sets all elements to be the constant
  ah = a(nm1+sh)
  ald = nm1*al
  ahd = nm1*ah
  tmp = a(0+sh)
  a(0+sh) = 0.0_rk
  do j = nm1-1,0,-1
     al = al*eta + a(j+sh)
     ah = ah*etah + a(j+sh)
     ald = ald*eta + j*a(j+sh)
     ahd = ahd*etah + j*a(j+sh)
  end do
  a(0+sh) = tmp
  z = ii*(fm + al - ah)
  zd = -(fmd + ald - ahd*dth_dnu)
end subroutine ptval
!
! Velocity calculation (quadrature routine)
!
subroutine a0val(n, u, z, zd, a0, q0)
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
     q0(i) = 1/(u + a0(i))
  end do
end subroutine a0val
!
! residual calculation
!
function resval(n, rho, beta, a)
  use constants
  implicit none
  ! input: rho, beta, a, n
  ! output: resval
  integer :: n, i
  real(rk) :: rho, beta, a(n), resval(n), u
  complex(rk), dimension(n) :: a0, q0
  complex(rk), dimension(4*n) :: eta, zeta, z, zd
  call ptval(n, rho, beta, a, eta, zeta, z, zd)
  u = a(1)
  call a0val(n, u, z, zd, a0, q0)
  resval = aimag( (a0+u)*zd((/ (i, i=2,2*n,2) /)) )
end function resval
!
! vpatch solver
!
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

  ! ! for lapack routines
  ! integer :: lda, lwork, ldb
  ! integer :: ipiv(n), iwork(n), info
  ! real(rk) :: anorm, rcond, work(4*n)

  ! for linpack routines
  integer :: ipvt(n)            ! pivot indices
  integer :: job = 0
  real(rk) :: rcond
  real(rk) :: z(n)

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

        ! linear algebra routines (LINPACK) ----------------------------
        call DGECO(jac, n, n, ipvt, rcond, z)
        ! LU-factors JAC and estimates RCOND; IPVT is the pivot
        ! indices; z is a work vector.
        call DGESL(jac, n, n, ipvt, res, job)
        ! Solves JAC*X = RES; on return, RES is the solution; JOB = 0
        ! for non-transposed problem

        ! ! linear algebra routines (LAPACK) -------------------------------
        ! lda = n; lwork = 4*n; ldb = n; anorm = mnorm(n, n, jac)
        ! ! calculating infinity norm of matrix JAC
        ! call DGETRF(n, n, jac, lda, ipiv, info)
        ! call DGECON('O', n, jac, lda, anorm, rcond, work, iwork, info)
        ! call DGETRS('N', n, 1, jac, lda, ipiv, res, ldb, info)

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
