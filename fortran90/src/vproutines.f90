module vproutines
  use constants
  use myroutines
  use fftw
contains
  !---------------------------------------------------------------------------80
  ! reading data from a file
  !---------------------------------------------------------------------------80
  subroutine readdata(filename, n, rho, beta, a)
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
  !---------------------------------------------------------------------------80
  ! writing data to a file
  !---------------------------------------------------------------------------80
  subroutine writedata(filename, n, rho, beta, a)
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
    write(1, '(f24.16)') beta
    do i=1,n
       write(1, '(es24.16)') a(i) ! es: format in scientific notation
    end do
    close(1)
  end subroutine writedata
  !---------------------------------------------------------------------------80
  ! point evaluation subroutine: deprecated
  ! ---------------------------------------------------------------------------80
  ! Since this routine is part of a module file, we can utilized
  ! "assumed shape" technique to eliminate one of the arguments (la)
  subroutine ptval_old(npt, la, rho, beta, a, eta, etah, zeta, z, zd)
    implicit none
    ! input: npt, la, rho, beta, a
    ! output: eta, zeta, z, zd
    ! aux var: nu, eta, zeta, etah, dth_dnu, fm, fmd, al, ah, ald, ahd
    integer, intent(in) :: npt ! number of uniform pts on unit circle
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
  end subroutine ptval_old

  !---------------------------------------------------------------------------80
  ! point evaluation subroutine - updated
  !---------------------------------------------------------------------------80
  !
  ! This routine calculates point values of various quantities some of
  ! which are related to quasi-solution such as
  !
  ! 	eta = uniformly spaced-out points on unit circle
  !
  ! 	etah = pre-image of rho^2/zeta under eta-Mobius transformation defined by
  ! 		( rho^2/zeta + beta )/( 1 + beta*rho^2/zeta )
  !
  ! 	dthdnu = d(theta)/d(nu) where theta and nu are angular
  !	 	parameters of the unit zeta and eta circle respectively:
  ! 		eta ( 1 - beta^2 )/( (1 - beta*eta) (eta - beta) )
  !
  !	zeta = ( eta - beta )/( 1 - beta*eta )
  !
  !	z = f_0(nu)
  !
  ! 	zd = f_0'(nu)
  !
  subroutine ptval(npt, rho, beta, a, eta, etah, dthdnu, zeta, z, zd)
    implicit none
    ! input: npt, la, rho, beta, a
    ! output: eta, zeta, z, zd
    ! aux var: nu, eta, zeta, etah, dthdnu, fm, fmd, al, ah, ald, ahd
    integer, intent(in) :: npt ! number of uniform pts on unit circle
    real(rk), intent(inout) :: a(:)
    complex(rk), dimension(npt), intent(inout) :: eta, etah, zeta, z, zd
    integer :: j, ncoef
    real(rk) :: rho, beta, tmp
    complex(rk), dimension(npt) :: nu, fm, fmd, dnuhdnu, dthdnu
    complex(rk), dimension(npt) :: al, ald, ah, ahd

    ncoef = size(a) - 1
    ! preliminary calculations
    nu = linspace1(0.0_rk, 2*pi, npt)
    eta = exp(ii*nu)
    zeta = (eta - beta)/(1.0_rk - beta*eta)
    etah = (rho**2 + beta*zeta)/(zeta + beta*rho**2)
    ! the following is d\theta/d\nu
    dthdnu = eta*(1.0_rk - beta**2)/( (1.0_rk - beta*eta)*(eta - beta))
    ! we are multiplying the above by d\hat{\nu}/d\theta to form d\hat{\nu}/d\nu
    dnuhdnu = dthdnu*( zeta*(rho**2)*(beta**2 - 1.0_rk) )/ &
         ( (rho**2 + beta*zeta)*(zeta + beta*rho**2) )

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
    zd = -(fmd + ald - ahd*dnuhdnu)
  end subroutine ptval

  !---------------------------------------------------------------------------80
  ! velocity calculation (quadrature routine): deprecated
  !---------------------------------------------------------------------------80
  !
  ! A_0(\nu) is calculated for \nu between 0 and \pi (strictly upper
  ! half eta-plane) and is meant to be used in residual
  ! calculation. This routine is now deprecated; use =a0val= instead.
  !
  subroutine a0val_old(n, u, z, zd, a0, q0)
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
  !---------------------------------------------------------------------------80
  ! velocity calculation (quadrature routine): current, updated
  !---------------------------------------------------------------------------80
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
  !---------------------------------------------------------------------------80
  ! residual calculation
  !---------------------------------------------------------------------------80
  function resval(n, rho, beta, a)
    implicit none
    ! input: rho, beta, a, n
    ! output: resval
    integer :: n, i, iv(n)
    real(rk) :: rho, beta, a(n), resval(n), u
    complex(rk), dimension(4*n) :: a0, q0
    complex(rk), dimension(4*n) :: eta, etah, zeta, z, zd, dthdnu
    ! call ptval_old(4*n, n, rho, beta, a, eta, etah, zeta, z, zd)
    call ptval(4*n, rho, beta, a, eta, etah, dthdnu, zeta, z, zd)
    u = a(1)
    iv = (/ (i, i=2,2*n,2) /)
    call a0val(4*n, u, z, zd, a0, q0)
    resval = aimag( (a0(iv)+u)*zd(iv) )
  end function resval
  !---------------------------------------------------------------------------80
  ! vpatch solver
  !---------------------------------------------------------------------------80
  subroutine vsolver &
       (n, rho_init, rho_term, nstep, beta_term, betaopt, a_init, a_term)
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
  !---------------------------------------------------------------------------80
  ! c_0 calculation
  !---------------------------------------------------------------------------80
  !
  ! The real constant c_0 defined by
  !
  ! 	c_0 = 1/(2\pi) \int_{0}^{2\pi} A_0(\nu) f_0'(\nu) \; d\nu
  !
  ! is calculated using a simple trapezoidal quadrature.
  !
  ! TODO: write it as a function?
  subroutine c0val(a0, zd, c0)
    ! input: a0, zd (same size)
    ! output: c0
    implicit none
    ! argument declaration
    complex(rk) :: a0(:), zd(:), c0
    ! local variables
    integer :: n
    n = size(a0)
    c0 = 1/real(n, rk)*sum( a0*zd ) ! this is supposed to be real; check
  end subroutine c0val
  !---------------------------------------------------------------------------80
  ! transfer matrix calculation: vectorized version; less efficient
  !---------------------------------------------------------------------------80
  !
  ! The "transfer matrix" M relates the positive Fourier coefficients
  ! of a function in S^2 to its negative modes. The matrix entries are
  ! calculated by the formula
  !
  !   M(k,j) = 1/(2\pi) \int_{0}^{2\pi} \hat{\eta}^j \eta^{k-1} d\nu
  !
  ! for 1 \le j,k \le N. Since eta = exp(i\nu), we note that M(k,j) =
  ! F[\hat{\eta}^j](1-k). So we can use FFT to calculate the entries
  ! efficiently.
  !
  ! TODO: Write it as a function?
  !
  ! NOTE: For a large n, M(k,j) for large j becomes very small,
  ! smaller than the smallest number representable in floating point
  ! system, which needs to be remembered.
  subroutine trmval_old(n, rho, beta, M)
    implicit none
    ! argument declaration
    integer :: n
    real(rk) :: rho, beta, M(n,n)
    ! local variables
    integer :: k, j
    complex(rk) :: eta(n), zeta(n), etah(n), integ(n,n), etahM(n,n)
    eta = exp(ii*linspace1(0.0_rk, 2*pi, n))
    zeta = (eta-beta)/(1.0_rk-beta*eta)
    etah = (rho**2/zeta+beta)/(1.0_rk+beta*(rho**2/zeta))
    M = zeros(n, n)
    integ = cmplx(spread(eta, 1, n) ** spread((/ (k-1, k=1,n) /), 2, n), kind=rk)
    etahM = spread(etah, 1, n)
    do j = 1,n
       integ = cmplx(integ*etahM, kind=rk)
       M(:,j) = 1/real(n, rk)*real(sum(integ, 2), rk)
    end do
    ! print *, 1/real(n, rk)
  end subroutine trmval_old
  !---------------------------------------------------------------------------80
  ! transfer matrix calculation: fft version; more efficient
  !---------------------------------------------------------------------------80
  !
  ! The "transfer matrix" M relates the positive Fourier coefficients
  ! of a function in our solution space to its negative modes. The
  ! matrix entries are calculated by the formula
  !
  !   M(k,j) = 1/(2\pi) \int_{0}^{2\pi} \hat{\eta}^j \eta^{k-1} d\nu
  !
  ! for 1 \le j,k \le N. Since eta = exp(i\nu), we note that M(k,j) =
  ! F[\hat{\eta}^j](1-k). So we can use FFT to calculate the entries
  ! efficiently.
  !
  ! TODO: Write it as a function?
  !
  ! NOTE: For a large n, M(k,j) for large j becomes very small,
  ! smaller than the smallest number representable in floating point
  ! system, which needs to be remembered.
  subroutine trmval(n, rho, beta, M)
    implicit none
    ! argument declaration
    integer :: n
    real(rk) :: rho, beta
    real(rk) :: M(n,n)
    ! local variables
    integer :: j, n2, plan_forward
    complex(rk), dimension(:), allocatable :: eta, zeta, etah, out
    complex(rk), dimension(:,:), allocatable :: P
    n2 = 2*n
    allocate(eta(n2), zeta(n2), etah(n2), out(n2))
    allocate(P(n2,n))
    eta = exp(ii*linspace1(0.0_rk, 2*pi, n2))
    zeta = (eta-beta)/(1.0_rk-beta*eta)
    etah = (rho**2/zeta+beta)/(1.0_rk+beta*(rho**2/zeta))
    P = spread(etah, 2, n) ** spread( (/(j, j=1,n)/), 1, n2)
    do j = 1,n
       call dfftw_plan_dft_1d_(plan_forward, n2, P(:,j), out, &
            FFTW_FORWARD, FFTW_ESTIMATE)
       call dfftw_execute_(plan_forward)
       P(:,j) = out/real(n2, rk)
       call dfftw_destroy_plan_( plan_forward )
    end do
    ! NOTE: be careful to match dimensions. P(n2,n), M(n,n)
    M(1, :) = real(P(1, :), rk)
    M(2:n, :) = real(P(n2:n+2:-1, :), rk)
  end subroutine trmval
  !---------------------------------------------------------------------------80
  ! j(nu) evaluation: related to Jacobi-Theta function
  !---------------------------------------------------------------------------80
  ! TODO: Write it using fast-convergent series representation. See
  ! Jacobi-Theta function.
  subroutine jtval(nterm, rho, zeta, jt)
    implicit none
    ! argument declaration
    integer :: nterm
    real(rk) :: rho
    complex(rk) :: zeta(:)
    complex(rk) :: jt(:)
    ! local variables
    integer :: i
    real(rk) :: rhotmp
    rhotmp = rho
    jt = cmplx(0, kind=rk)
    do i = 1,nterm
       jt = jt + rhotmp*(1/(zeta+rhotmp)**2 + 1/(1+rhotmp*zeta)**2)
       rhotmp = rhotmp*rho**2
    end do
    jt = 2*zeta*jt
  end subroutine jtval

  !---------------------------------------------------------------------------80
  ! j1(nu, theta) evaluation: original version with fft on shifted data
  !---------------------------------------------------------------------------80
  !
  ! This routine evaluates j1(nu, theta) at N x N (nu, theta)-points where
  !
  ! 	nu_j    = 2 pi (j-1) / N   for  j = 1..N ( normal data points)
  ! 	theta_k = 2 pi (k-1/2) / N for  k = 1..N (shifted data points)
  !
  ! and
  !
  ! 	j1(nu, theta) = -omega0(nu)*f0(nu)*conj( f0'(theta) )/( 2i*exp(i*nu) )
  !	   * ( exp(i*nu) - exp(i*theta) )/( f0(nu)^2 - f0(theta)^2 )
  !
  ! and calculate its Fourier coefficients
  !
  ! 	j1(nu, theta) = \sum_{j,k} A_{j,k}^{(1)} exp(i*(j*nu+k*theta))
  !
  subroutine j1val_sh(eta, z, zd, omega, j1, a1)
    implicit none
    ! argument declaration
    complex(rk) :: eta(:), z(:), zd(:), omega(:) ! inputs
    complex(rk) :: j1(:,:), a1(:,:)              ! output
    ! local variables
    integer :: n2, n, i
    integer, allocatable :: io(:), ie(:)
    real(rk), allocatable :: v(:)
    complex(rk), dimension(:,:), allocatable :: &
         znu, zth, zdnu, zdth, ph, omnu, einu, eith
    ! fft
    integer(8) :: plan_forward
    complex(rk), dimension(:,:), allocatable :: in2, out2

    n2 = size(z)
    n = n2/2
    allocate(io(n), ie(n), v(n), ph(n,n))
    allocate(znu(n,n), zth(n,n), zdnu(n,n), zdth(n,n))
    allocate(omnu(n,n), einu(n,n), eith(n,n))
    allocate(in2(n,n), out2(n,n))
    io = (/ (i, i=1,n2,2) /)    ! odd indices: nu
    ie = (/ (i, i=2,n2,2) /)    ! even indices: theta
    znu = spread( z(io), 2, n )
    zth = spread( z(ie), 1, n )
    zdth = spread( zd(ie), 1, n )
    omnu = spread( omega(io), 2, n )
    einu = spread( eta(io), 2, n )
    eith = spread( eta(ie), 1, n )

    ! point values j1 on staggered grid points
    j1 = -omnu*znu*conjg(zdth)/(2*ii*einu)
    j1 = j1*( einu-eith )/( znu**2 - zth**2 )

    ! 2d-fft with phase-factor adjustment
    in2 = j1
    call dfftw_plan_dft_2d_ ( plan_forward, n, n, in2, out2, FFTW_FORWARD, &
         FFTW_ESTIMATE)
    call dfftw_execute_ ( plan_forward )
    a1 = out2/real(n**2, rk)


    v = (/ (i,i=0,n/2), (i,i=-n/2+1,-1) /)/real(n, rk)
    ph = spread( exp(-ii*pi*v), 1, n )
    a1 = a1 * ph
    a1(n/2+1, :) = 0.0_rk
    a1(:, n/2+1) = 0.0_rk
    call dfftw_destroy_plan_ ( plan_forward )
    deallocate(io, ie, v, ph)
    deallocate(znu, zth, zdnu, zdth)
    deallocate(omnu, einu, eith)
    deallocate(in2, out2)
  end subroutine j1val_sh

  !---------------------------------------------------------------------------80
  ! j2(nu, theta) evaluation: original version with fft on shifted data
  !---------------------------------------------------------------------------80
  !
  ! This routine evaluates j2(nu, theta) at N x N (nu, theta)-points where
  !
  ! 	nu_j    = 2 pi (j-1) / N   for  j = 1..N ( normal data points)
  ! 	theta_k = 2 pi (k-1/2) / N for  k = 1..N (shifted data points)
  !
  ! and
  !
  ! 	j2(nu, theta) = omega0(nu)*f0(theta)*f0'(theta)/( -2i*exp(-i*nu) )
  !	   * ( exp(-i*nu) - exp(-i*theta) )/( f0(nu)^2 - f0(theta)^2 )
  !
  ! and calculate its Fourier coefficients
  !
  ! 	j2(nu, theta) = \sum_{j,k} A_{j,k}^{(2)} exp(i*(j*nu+k*theta))
  !
  subroutine j2val_sh(eta, z, zd, omega, j2, a2)
    implicit none
    ! argument declaration
    complex(rk) :: eta(:), z(:), zd(:), omega(:) ! inputs
    complex(rk) :: j2(:,:), a2(:,:)              ! output
    ! local variables
    integer :: n2, n, i
    integer, allocatable :: io(:), ie(:)
    real(rk), allocatable :: v(:)
    complex(rk), dimension(:,:), allocatable :: &
         znu, zth, zdnu, zdth, ph, omnu, einu, eith
    ! fft
    integer(8) :: plan_forward
    complex(rk), dimension(:,:), allocatable :: in2, out2

    n2 = size(z)
    n = n2/2
    allocate(io(n), ie(n), v(n), ph(n,n))
    allocate(znu(n,n), zth(n,n), zdnu(n,n), zdth(n,n))
    allocate(omnu(n,n), einu(n,n), eith(n,n))
    allocate(in2(n,n), out2(n,n))
    io = (/ (i, i=1,n2,2) /)    ! odd indices: nu
    ie = (/ (i, i=2,n2,2) /)    ! even indices: theta
    znu = spread( z(io), 2, n )
    zth = spread( z(ie), 1, n )
    zdth = spread( zd(ie), 1, n )
    omnu = spread( omega(io), 2, n )
    einu = spread( eta(io), 2, n )
    eith = spread( eta(ie), 1, n )

    ! point values j2 on staggered grid points
    j2 = omnu*zth*zdth/(-2*ii*conjg(einu))
    j2 = j2*( conjg(einu-eith) )/( znu**2 - zth**2 )

    ! 2d-fft with phase-factor adjustment
    in2 = j2
    call dfftw_plan_dft_2d_ ( plan_forward, n, n, in2, out2, FFTW_FORWARD, &
         FFTW_ESTIMATE)
    call dfftw_execute_ ( plan_forward )
    a2 = out2/real(n**2, rk)


    v = (/ (i,i=0,n/2), (i,i=-n/2+1,-1) /)/real(n, rk)
    ph = spread( exp(-ii*pi*v), 1, n )
    a2 = a2 * ph
    a2(n/2+1, :) = 0.0_rk
    a2(:, n/2+1) = 0.0_rk
    call dfftw_destroy_plan_ ( plan_forward )
    deallocate(io, ie, v, ph)
    deallocate(znu, zth, zdnu, zdth)
    deallocate(omnu, einu, eith)
    deallocate(in2, out2)
  end subroutine j2val_sh

  !---------------------------------------------------------------------------80
  ! j3(nu, theta) evaluation: original version with fft on shifted data
  !---------------------------------------------------------------------------80
  !
  ! This routine evaluates j3(nu, theta) at N x N (nu, theta)-points where
  !
  ! 	nu_j    = 2 pi (j-1) / N   for  j = 1..N ( normal data points)
  ! 	theta_k = 2 pi (k-1/2) / N for  k = 1..N (shifted data points)
  !
  ! and
  !
  ! 	j3(nu, theta) = i * omega0(nu)*conj( f0'(theta) )/( f0(nu)+f0(theta) )
  !
  ! and calculate its Fourier coefficients
  !
  ! 	j3(nu, theta) = \sum_{j,k} A_{j,k}^{(3)} exp(i*(j*nu+k*theta))
  !
  subroutine j3val_sh(eta, z, zd, omega, j3, a3)
    implicit none
    ! argument declaration
    complex(rk) :: eta(:), z(:), zd(:), omega(:) ! inputs
    complex(rk) :: j3(:,:), a3(:,:)              ! output
    ! local variables
    integer :: n2, n, i
    integer, allocatable :: io(:), ie(:)
    real(rk), allocatable :: v(:)
    complex(rk), dimension(:,:), allocatable :: &
         znu, zth, zdnu, zdth, ph, omnu, einu, eith
    ! fft
    integer(8) :: plan_forward
    complex(rk), dimension(:,:), allocatable :: in2, out2

    n2 = size(z)
    n = n2/2
    allocate(io(n), ie(n), v(n), ph(n,n))
    allocate(znu(n,n), zth(n,n), zdnu(n,n), zdth(n,n))
    allocate(omnu(n,n), einu(n,n), eith(n,n))
    allocate(in2(n,n), out2(n,n))
    io = (/ (i, i=1,n2,2) /)    ! odd indices: nu
    ie = (/ (i, i=2,n2,2) /)    ! even indices: theta
    znu = spread( z(io), 2, n )
    zth = spread( z(ie), 1, n )
    zdth = spread( zd(ie), 1, n )
    omnu = spread( omega(io), 2, n )
    einu = spread( eta(io), 2, n )
    eith = spread( eta(ie), 1, n )

    ! point values j3 on staggered grid points
    j3 = ii*omnu*conjg(zdth)/( znu + zth )

    ! 2d-fft with phase-factor adjustment
    in2 = j3
    call dfftw_plan_dft_2d_ ( plan_forward, n, n, in2, out2, FFTW_FORWARD, &
         FFTW_ESTIMATE)
    call dfftw_execute_ ( plan_forward )
    a3 = out2/real(n**2, rk)
    v = (/ (i,i=0,n/2), (i,i=-n/2+1,-1) /)/real(n, rk)
    ph = spread( exp(-ii*pi*v), 1, n )
    a3 = a3 * ph
    a3(n/2+1, :) = 0.0_rk
    a3(:, n/2+1) = 0.0_rk
    call dfftw_destroy_plan_ ( plan_forward )
    deallocate(io, ie, v, ph)
    deallocate(znu, zth, zdnu, zdth)
    deallocate(omnu, einu, eith)
    deallocate(in2, out2)
  end subroutine j3val_sh

  !---------------------------------------------------------------------------80
  ! h4(nu, theta) evaluation: new, updated, non-shifted data
  !---------------------------------------------------------------------------80
  !
  ! This routine evaluates h4(nu, theta) at N x N (nu, theta)-points where
  !
  ! 	nu_j    = 2 pi (j-1) / N for  j = 1..N
  ! 	theta_k = 2 pi (k-1) / N for  k = 1..N
  !
  ! and
  !
  ! 	h4(nu, theta) = f_0'(theta)
  !	   * ( |f0(nu)|^2 - |f0(theta)|^2 ) / ( f0(nu)^2 - f0(theta)^2 )
  !
  ! In case j = k, the integrand has a removable singularity; we
  ! replace h4(nu, theta) by the corresponding limit.
  !
  subroutine h4val(z, zd, h4, b4)
    implicit none
    ! argument declaration
    complex(rk) :: z(:), zd(:)      ! inputs
    complex(rk) :: h4(:,:), b4(:,:) ! output
    ! local variables
    integer :: n, i
    complex(rk), dimension(:,:), allocatable :: &
         znu, zth, zdth
    ! fft
    integer(8) :: plan_forward
    complex(rk), dimension(:,:), allocatable :: in2, out2

    n = size(z)
    allocate(znu(n,n), zth(n,n), zdth(n,n))
    allocate(in2(n,n), out2(n,n))
    znu = spread( z, 2, n )
    zth = spread( z, 1, n )
    zdth = spread( zd, 1, n )
    h4 = zdth*( abs(znu)**2 - abs(zth)**2 )/( znu**2 - zth**2 )
    ! The above line causes warning for floating point error, which
    ! gets fixed in the lines below.
    do i = 1,n
       h4(i,i) = zd(i)/2*( conjg(zd(i))/zd(i) + conjg(z(i))/z(i) )
    end do

    ! 2d-fft
    in2 = h4
    call dfftw_plan_dft_2d_ ( plan_forward, n, n, in2, out2, FFTW_FORWARD, &
         FFTW_ESTIMATE)
    call dfftw_execute_ ( plan_forward )
    b4 = out2/real(n**2, rk)
    b4(n/2+1, :) = 0.0_rk
    call dfftw_destroy_plan_ ( plan_forward )
    deallocate(znu, zth, zdth, in2, out2)
  end subroutine h4val

  !---------------------------------------------------------------------------80
  ! h4(nu, theta) evaluation: original version with fft on shifted data
  !---------------------------------------------------------------------------80
  !
  ! This routine evaluates h4(nu, theta) at N x N (nu, theta)-points where
  !
  ! 	nu_j    = 2 pi (j-1) / N   for  j = 1..N
  ! 	theta_k = 2 pi (k-1/2) / N for  k = 1..N
  !
  ! and
  !
  ! 	h4(nu, theta) = f_0'(theta)
  !	   * ( |f0(nu)|^2 - |f0(theta)|^2 ) / ( f0(nu)^2 - f0(theta)^2 )
  !
  subroutine h4val_sh(z, zd, h4, b4)
    implicit none
    ! argument declaration
    complex(rk) :: z(:), zd(:)      ! inputs
    complex(rk) :: h4(:,:), b4(:,:) ! output
    ! local variables
    integer :: n2, n, i
    integer, allocatable :: io(:), ie(:)
    real(rk), allocatable :: v(:)
    complex(rk), dimension(:,:), allocatable :: &
         znu, zth, zdnu, zdth, ph
    ! fft
    integer(8) :: plan_forward
    complex(rk), dimension(:), allocatable :: in, out
    complex(rk), dimension(:,:), allocatable :: in2, out2

    n2 = size(z)
    n = n2/2
    allocate(io(n), ie(n))
    allocate(znu(n,n), zth(n,n), zdnu(n,n), zdth(n,n))
    allocate(in(n), out(n))
    allocate(in2(n,n), out2(n,n))
    io = (/ (i, i=1,n2,2) /)
    ie = (/ (i, i=2,n2,2) /)
    znu = spread( z(io), 2, n )
    zth = spread( z(ie), 1, n )
    zdth = spread( zd(ie), 1, n )
    h4 = zdth*( abs(znu)**2 - abs(zth)**2 )/( znu**2 - zth**2 )
    b4 = cmplx(0.0, kind=rk)

    ! 1d-fft along columns then along rows
    ! do j = 1,n
    !    in = h4(:, j)
    !    call dfftw_plan_dft_1d_ &
    !         ( plan_forward, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    !    call dfftw_execute_ ( plan_forward )
    !    b4(:, j) = out/real(n, rk)
    ! end do
    ! do i = 1,n
    !    in = b4(i, :)
    !    call dfftw_plan_dft_1d_ &
    !         ( plan_forward, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    !    call dfftw_execute_ ( plan_forward )
    !    b4(i, :) = out/real(n, rk)
    ! end do
    ! call dfftw_destroy_plan_ ( plan_forward )

    ! 2d-fft
    in2 = h4
    call dfftw_plan_dft_2d_ ( plan_forward, n, n, in2, out2, FFTW_FORWARD, &
         FFTW_ESTIMATE)
    call dfftw_execute_ ( plan_forward )
    b4 = out2/real(n**2, rk)
    allocate(v(n), ph(n,n))
    v = (/ (i,i=0,n/2), (i,i=-n/2+1,-1) /)/real(n, rk)
    ph = spread( exp(-ii*pi*v), 1, n )
    b4 = b4 * ph
    b4(n/2+1, :) = 0.0_rk
    call dfftw_destroy_plan_ ( plan_forward )
    deallocate(znu, zth, zdnu, zdth, io, ie, in, out)
  end subroutine h4val_sh

  !---------------------------------------------------------------------------80
  ! construction of \hat{t}_{j,m}
  !---------------------------------------------------------------------------80
  subroutine thval(a1, a2, a3, t)
    implicit none
    ! argument declaration
    complex(rk), dimension(:,:) :: a1, a2, a3
    complex(rk), dimension(:,:), allocatable :: t

    ! local variables
    integer :: n, n2, i, j, k, m, nm1, n2m1
    integer, allocatable :: ind(:)
    complex(rk), allocatable :: tmp(:)
    complex(rk), dimension(:,:), allocatable :: a1c, a2c, a3c

    n2 = size(a1, 1)
    n = int(n2/2)
    nm1 = n-1
    n2m1 = n2-1

    ! allocate(t(-n:nm1,-n:nm1))
    allocate(ind(n2))
    allocate(a1c(-n2:n2m1,-n2:n2m1))
    allocate(a2c(-n2:n2m1,-n2:n2m1))
    allocate(a3c(-n2:n2m1,-n2:n2m1))

    ! fftshift: shifting indices so that the zero-th modes appear in the
    ! middle of spectrum
    a1c = cmplx(0.0, kind=rk)
    a2c = cmplx(0.0, kind=rk)
    a3c = cmplx(0.0, kind=rk)
    ind = (/ (i, i=n+1,n2), (i, i=1,n) /)
    a1c(-n:nm1,-n:nm1) = a1(ind,ind)
    a2c(-n:nm1,-n:nm1) = a2(ind,ind)
    a3c(-n:nm1,-n:nm1) = a3(ind,ind)

    ! Forming \tilde{t}_{j,m}
    do j = -n,nm1
       do m = -n,nm1
          allocate(tmp(abs(m)+1))
          ! tmp = (/(a1c(j-k,k-m) + a2c(j-k,k+m), k=0,m,sgn0(m))/)
          tmp = (/(a1c(j-k,k-m) + a2c(j+k,-k+m), k=0,m,sgn0(m))/)
          ! t(j+n+1,m+n+1) = &
          t(j,m) = &
               - a1c(j-m,0) + a1c(j,-m) &
               - a2c(j+m,0) + a2c(j,m) &
               - a3c(j,-m) &
               -sgn(m)*sum( tmp(1:abs(m)) + tmp(2:abs(m)+1) )
          deallocate(tmp)
       end do
    end do
    deallocate(ind, a1c, a2c, a3c)
  end subroutine thval

  !---------------------------------------------------------------------------80
  ! construction of \hat{t}_{j,m}: check version
  !---------------------------------------------------------------------------80
  subroutine thval_check(n, a1, a2, a3, t)
    implicit none
    ! argument declaration
    integer :: n
    complex(rk), dimension(2*n,2*n) :: a1, a2, a3
    complex(rk), dimension(-n:n-1,-n:n-1) :: t

    ! local variables
    integer :: n2, i, j, k, m, nm1, n2m1
    integer, allocatable :: ind(:)
    complex(rk), allocatable :: tmp(:)
    complex(rk), dimension(:,:), allocatable :: a1c, a2c, a3c

    n2 = 2*n
    nm1 = n-1
    n2m1 = n2-1

    allocate(ind(n2))
    allocate(a1c(-n2:n2m1,-n2:n2m1))
    allocate(a2c(-n2:n2m1,-n2:n2m1))
    allocate(a3c(-n2:n2m1,-n2:n2m1))

    ! fftshift: shifting indices so that the zero-th modes appear in the
    ! middle of spectrum
    a1c = cmplx(0.0, kind=rk)
    a2c = cmplx(0.0, kind=rk)
    a3c = cmplx(0.0, kind=rk)
    ind = (/ (i, i=n+1,n2), (i, i=1,n) /)
    a1c(-n:nm1,-n:nm1) = a1(ind,ind)
    a2c(-n:nm1,-n:nm1) = a2(ind,ind)
    a3c(-n:nm1,-n:nm1) = a3(ind,ind)

    ! print *, lbound(t), ubound(t)
    ! Forming \tilde{t}_{j,m}
    do j = -n,nm1
       do m = -n,nm1
          allocate(tmp(abs(m)+1))
          tmp = (/(a1c(j-k,k-m) + a2c(j+k,-k+m), k=0,m,sgn0(m))/)
          ! t(j+n+1,m+n+1) = &
          t(j,m) = &
               - a1c(j-m,0) + a1c(j,-m) &
               - a2c(j+m,0) + a2c(j,m) &
               - a3c(j,-m) &
               - sgn(m)*sum( tmp(1:abs(m)) + tmp(2:abs(m)+1) )
          deallocate(tmp)
       end do
    end do
    ! print *, lbound(t), ubound(t)
    deallocate(ind, a1c, a2c, a3c)
  end subroutine thval_check

  !---------------------------------------------------------------------------80
  ! t1val: construction of t^{(1)}_{j,m}
  !---------------------------------------------------------------------------80
  subroutine t1val(a1, a2, a3, b4, dthdnu, t)
    implicit none
    ! argument declaration
    complex(rk), dimension(:) :: dthdnu
    complex(rk), dimension(:,:) :: a1, a2, a3, b4
    complex(rk), dimension(:,:), allocatable :: t

    ! local variables
    integer :: n, n2, n4, i, j, k, m, nm1, n2m1, ltmp
    integer, allocatable :: ind(:)
    complex(rk), allocatable :: tmp(:)
    complex(rk), dimension(:,:), allocatable :: a1c, a2c, a3c, b4c

    ! fft
    integer(8) :: plan_forward  ! integer kind must be 8.
    complex(rk), dimension(:), allocatable :: v

    n4 = size(dthdnu, 1)
    n2 = size(a1, 1)
    n = int(n2/2)
    nm1 = n-1
    n2m1 = n2-1

    allocate(v(-n2:n2m1))
    allocate(a1c(-n2:n2m1,-n2:n2m1))
    allocate(a2c(-n2:n2m1,-n2:n2m1))
    allocate(a3c(-n2:n2m1,-n2:n2m1))
    allocate(b4c(-n2:n2m1,-n2:n2m1))

    ! fft on dthdnu to obtain v_j.
    ! - dthdnu is a vector of point values at 4N equi-spaced points.
    allocate(tmp(n4))
    call dfftw_plan_dft_1d_ ( plan_forward, n4, dthdnu, tmp, &
         FFTW_FORWARD, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_forward )
    tmp = tmp/real(n4, rk)
    call dfftw_destroy_plan_ ( plan_forward )

    ! fftshift: shifting indices so that the zero-th modes appear in the
    ! middle of spectrum
    allocate(ind(n4))
    ind = (/ (i, i=n2+1,n4), (i, i=1,n2) /)
    v = tmp(ind)
    deallocate(ind, tmp)
    a1c = cmplx(0.0, kind=rk)
    a2c = cmplx(0.0, kind=rk)
    a3c = cmplx(0.0, kind=rk)
    b4c = cmplx(0.0, kind=rk)
    allocate(ind(n2))
    ind = (/ (i, i=n+1,n2), (i, i=1,n) /)
    a1c(-n:nm1,-n:nm1) = a1(ind,ind)
    a2c(-n:nm1,-n:nm1) = a2(ind,ind)
    a3c(-n:nm1,-n:nm1) = a3(ind,ind)
    b4c(-n:nm1,-n:nm1) = b4(ind,ind)
    deallocate(ind)

    ! Forming \hat{t}_{j,m}
    do j = -n,nm1
       do m = -n,nm1
          ltmp = abs(m) + 1
          allocate(tmp(ltmp))
          tmp = (/(a1c(j-k,k-m) + a2c(j+k,-k+m), k=0,m,sgn0(m))/)
          t(j,m) = &
               - a1c(j-m,0) + a1c(j,-m) &
               - a2c(j+m,0) + a2c(j,m) &
               - a3c(j,-m) &
               - sgn(m)*sum( tmp(1:ltmp-1) + tmp(2:ltmp) )
          deallocate(tmp)
       end do
    end do

    ! Forming t^{(1)}_{j,m} for -N <= j,m < N
    allocate(tmp(-n:nm1))
    tmp = (/ (t(0,i), i=-n,nm1) /)
    do j = -n,nm1
       do m = -n,nm1
          t(j,m) = t(j,m) - ( tmp(m)-m*b4c(-m,0) ) * v(j)
       end do
    end do
    deallocate(v, a1c, a2c, a3c, b4c, tmp)
  end subroutine t1val

  !---------------------------------------------------------------------------80
  ! construction of a^{(1)}_{j,m}
  !---------------------------------------------------------------------------80
  subroutine a1jmval(q0, t1, trm, a1jm)
    implicit none
    ! declare arguments
    complex(rk) :: q0(:)
    complex(rk), allocatable :: t1(:,:)
    real(rk) :: trm(:,:)
    complex(rk), allocatable :: a1jm(:,:)

    ! print *, size(q0)
    ! print *, shape(t1), lbound(t1), ubound(t1)
    ! print *, shape(trm), lbound(trm), ubound(trm)
    ! print *, shape(a1jm), lbound(a1jm), ubound(a1jm)

    ! local variables
    integer :: n, n2, n4, nm1, n2m1, i, j, m
    integer, allocatable :: kk(:), ll(:)
    integer, allocatable :: ind(:)
    complex(rk), allocatable :: tmp(:)
    complex(rk), allocatable :: t1c(:,:)
    complex(rk), allocatable :: q(:)

    ! fft variables
    integer(8) :: plan_forward

    n4 = size(q0)
    n2 = int(n4/2)
    n = int(n2/2)
    n2m1 = n2-1
    nm1 = n-1
    allocate(t1c(-n2:n2m1,-n2:n2m1))
    allocate(q(-n2:n2m1))

    ! fft on q0 followed by fftshift
    allocate(tmp(n4), ind(n4))
    call dfftw_plan_dft_1d_ ( plan_forward, n4, q0, tmp, &
         FFTW_FORWARD, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_forward )
    tmp = tmp/real(n4, rk)
    call dfftw_destroy_plan_ ( plan_forward )
    ind = (/ (i, i=n2+1,n4), (i, i=1,n2) /)
    q = tmp(ind)
    deallocate(ind, tmp)

    ! padding t1 with zeros; recall that t1 already has zeroth modes
    ! in the middle of spectrum.
    t1c = cmplx(0.0, kind=rk)
    t1c(-n:nm1,-n:nm1) = t1

    ! forming a1jm
    allocate(tmp(n2), kk(n2), ll(n))
    kk = (/ (i, i=-n,nm1) /)
    ll = (/ (i, i=1,n) /)
    do j = 0,nm1
       do m = 1,n
          tmp = q(kk) * ( t1c(j-kk,m) - matmul(t1c(j-kk,1-ll), trm(ll,m)) )
          a1jm(j,m) = sum( tmp )
       end do
    end do
    deallocate(tmp, kk, ll)
  end subroutine a1jmval


  !---------------------------------------------------------------------------80
  ! construction of a_{j,m}
  !---------------------------------------------------------------------------80
  subroutine ajmval(q0, t1, trm, lu, ajm)
    implicit none
    ! declare arguments
    complex(rk) :: q0(:)
    complex(rk), allocatable :: t1(:,:)
    real(rk) :: trm(:,:)
    complex(rk) :: lu(:)
    complex(rk), allocatable :: ajm(:,:)

    ! print *, size(q0)
    ! print *, shape(t1), lbound(t1), ubound(t1)
    ! print *, shape(trm), lbound(trm), ubound(trm)
    ! print *, shape(ajm), lbound(ajm), ubound(ajm)

    ! local variables
    integer :: n, n2, n4, nm1, n2m1, i, j, m
    integer, allocatable :: kk(:), ll(:)
    integer, allocatable :: ind(:)
    complex(rk), allocatable :: tmp(:)
    complex(rk), allocatable :: t1c(:,:), a1jm(:,:)
    complex(rk), allocatable :: q(:), lt(:)

    ! fft variables
    integer(8) :: plan_forward

    n4 = size(q0)
    n2 = int(n4/2)
    n = int(n2/2)
    n2m1 = n2-1
    nm1 = n-1
    allocate(t1c(-n2:n2m1,-n2:n2m1))
    allocate(a1jm(0:nm1,1:n))
    allocate(q(-n2:n2m1), lt(0:n4-1))


    ! fft on q0 followed by fftshift
    allocate(tmp(n4), ind(n4))
    call dfftw_plan_dft_1d_ ( plan_forward, n4, q0, tmp, &
         FFTW_FORWARD, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_forward )
    tmp = tmp/real(n4, rk)
    call dfftw_destroy_plan_ ( plan_forward )
    ind = (/ (i, i=n2+1,n4), (i, i=1,n2) /)
    q = tmp(ind)

    ! fft on lu followed by division by zeroth mode to form \tilde{l}_j
    call dfftw_plan_dft_1d_ ( plan_forward, n4, lu, tmp, &
         FFTW_FORWARD, FFTW_ESTIMATE )
    call dfftw_execute_ ( plan_forward )
    tmp = tmp/real(n4, rk)
    lt = -tmp/tmp(1)
    deallocate(tmp, ind)

    ! padding t1 with zeros; recall that t1 already has zeroth modes
    ! in the middle of spectrum.
    t1c = cmplx(0.0, kind=rk)
    t1c(-n:nm1,-n:nm1) = t1

    ! forming a1jm
    allocate(tmp(n2), kk(n2), ll(n))
    kk = (/ (i, i=-n,nm1) /)
    ll = (/ (i, i=1,n) /)
    do j = 0,nm1
       do m = 1,n
          tmp = q(kk) * ( t1c(j-kk,m) - matmul(t1c(j-kk,1-ll), trm(ll,m)) )
          a1jm(j,m) = sum( tmp )
       end do
    end do

    ! forming ajm
    do j = 1,nm1
       do m = 1,nm1
          ajm(j,m) = a1jm(j,m) + a1jm(0,m)*lt(j)
       end do
    end do
    deallocate(a1jm, lt, tmp, t1c, kk, ll, q)
  end subroutine ajmval

  !---------------------------------------------------------------------------80
  ! signum function
  !---------------------------------------------------------------------------80
  !
  ! sgn(n) = 1	if n > 0
  ! 	   = 0	if n = 0
  ! 	   = -1 if n < 0
  !
  function sgn(n)
    implicit none
    integer :: n, sgn
    if (n > 0) then
       sgn = 1
    else if (n == 0) then
       sgn = 0
    else
       sgn = -1
    end if
  end function sgn

  !---------------------------------------------------------------------------80
  ! signum0 function
  !---------------------------------------------------------------------------80
  !
  ! sgn0(n) = 1   if n >=0
  ! 	    = -1  if n < 0
  !
  function sgn0(n)
    implicit none
    integer :: n, sgn0
    if (n >= 0) then
       sgn0 = 1
    else
       sgn0 = -1
    end if
  end function sgn0
end module vproutines
