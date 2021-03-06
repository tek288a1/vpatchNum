module constants
  implicit none
  integer, parameter :: rk = 8 ! rk = 8 (double) ; rk = 16 (quadruple)
  integer, parameter :: sh = 1 ! shift
  real(rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
  real(rk), parameter :: eps = epsilon(1.0_rk)
  complex(rk), parameter :: ii = cmplx(0.0_rk, 1.0_rk)
  character(len=80), parameter ::  line = &
       '--------------------------------------------------------------------------------'
end module constants

module myroutines
  use constants
  implicit none
contains
  !
  ! subroutine printing complex data on terminal
  !
  subroutine write_complex_data(z)
    integer :: i, size_z
    complex(rk), dimension(:), intent(in) :: z

    size_z = size(z)
    do i = 1,size_z
       write(*,'(2es24.15)') real(z(i)), imag(z(i))
    end do
  end subroutine write_complex_data
  !
  ! linspace function: n equispaced points between a and b
  !
  function linspace(a, b, n)
    integer :: i
    integer, intent(in) :: n
    real(rk), intent(in) :: a, b
    real(rk) :: linspace(n)
    linspace = a + (b-a)/real(n-1, rk)*[(i, i=0,n-1)]
  end function linspace
  !
  ! linspace1 function: first n points in the set of n+1 equispaced
  !	 points between a and b
  function linspace1(a, b, n)
    integer :: i
    integer, intent(in) :: n
    real(rk) :: tmp(n+1)
    real(rk), intent(in) :: a, b
    real(rk) :: linspace1(n)
    tmp = linspace(a, b, n+1)
    linspace1 = tmp([(i, i=1,n)])
  end function linspace1
  !
  ! linspaceh function: points spaced-out by h between a and b
  !
  function linspaceh(a, b, h)
    real(8), intent(in) :: a, b, h
    real(8), dimension(:), allocatable :: linspaceh
    integer :: i, n
    n = int((b-a)/h)
    ! if (abs((b-a)/h-n)<eps) then  -----------> it doesn't work
    ! if ( abs( (b-a)-n*h ) < eps ) then ------> this works
    if ( abs(mod(b-a, h)) < eps ) then ! note we compare against the
                                       ! machine epsilon
       allocate(linspaceh(n+1))
       linspaceh = (/ (a + h*i, i=0,n) /)
    else
       allocate(linspaceh(n+2))
       linspaceh = (/ (a+h*i, i=0,n), b /)
    end if
  end function linspaceh
  !
  ! (square) identity matrix constructor
  !
  function eye(n)
    integer :: i, j, n
    real(rk) :: eye(n,n)
    do i = 1,n
       do j = 1,n
          if (i==j) then
             eye(i,j) = 1.0_rk
          else
             eye(i,j) = 0.0_rk
          end if
       end do
    end do
  end function eye
  !
  ! zero array constructor
  !
  function zeros(m, n)
    integer :: i, j, m, n
    real(rk) :: zeros(m,n)
    do i = 1,m
       do j = 1,n
          zeros(i,j) = 0.0_rk
       end do
    end do
  end function zeros
  !
  ! vector-norm evaluator
  !
  function vnorm(v, which)
    integer :: i, n
    real(rk), dimension(:) :: v
    real(rk) :: vnorm
    character(len=*) :: which

    n = size(v)
    if (which=='inf') then
       vnorm = maxval(abs(v))
    end if
  end function vnorm
  !
  ! genfilename: creating a string for data file name
  !
  function genfilename(n, rho, betaopt)
    use constants
    implicit none
    ! declaring arguments
    integer, intent(in) :: n
    real(rk), intent(in) :: rho
    integer, intent(in) :: betaopt
    ! output
    character(len=28) :: genfilename
    ! local variables
    character(len=*), parameter :: fmt = trim('(a, i1, a, i0.2, a, i0.4, a)')
    integer :: log2n
    log2n = int(log(real(n))/log(real(2)))
    write( genfilename, fmt ) &
         '../data/vp_b', betaopt, '_n', log2n, '_r', int(1.0d4*rho), '.dat'
  end function genfilename
  !
  ! time stamp
  !
  subroutine timestamp ( )
    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none

    character ( len = 8 ) ampm
    integer ( kind = 4 ) d
    integer ( kind = 4 ) h
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4 ) n
    integer ( kind = 4 ) s
    integer ( kind = 4 ) values(8)
    integer ( kind = 4 ) y

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
       ampm = 'AM'
    else if ( h == 12 ) then
       if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
       else
          ampm = 'PM'
       end if
    else
       h = h - 12
       if ( h < 12 ) then
          ampm = 'PM'
       else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
             ampm = 'Midnight'
          else
             ampm = 'AM'
          end if
       end if
    end if

    write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
  end subroutine timestamp
  !
  ! matrix infinity-norm calculator
  !
  function mnorm ( m, n, a )

    !*****************************************************************************80
    !
    !! R8MAT_NORM_LI returns the matrix L-infinity norm of an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is a two dimensional matrix of double precision real values.
    !
    !    The matrix L-infinity norm is defined as:
    !
    !      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
    !
    !    The matrix L-infinity norm is derived from the vector L-infinity norm,
    !    and satisifies:
    !
    !      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer M, the number of rows in A.
    !
    !    Input, integer N, the number of columns in A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose L-infinity
    !    norm is desired.
    !
    !    Output, real ( kind = 8 ) R8MAT_NORM_LI, the L-infinity norm of A.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real   ( kind = 8 ) a(m,n)
    integer ( kind = 4 ) i
    real ( kind = 8 ) mnorm

    mnorm = 0.0D+00

    do i = 1, m
       mnorm = max ( mnorm, sum ( abs ( a(i,1:n) ) ) )
    end do

    return
  end function mnorm
end module myroutines
