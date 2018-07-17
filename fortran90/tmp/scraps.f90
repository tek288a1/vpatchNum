  ! ! visual check
  ! write(*, '(/a)') '(1,1)-block'
  ! do i = 1,n2
  !    write(*, '(*(es9.2))') real(aa(-n2-1+i, -n2:-1))
  ! end do

  ! write(*, '(/a)') '(1,2)-block'
  ! do i = 1,n2
  !    write(*, '(*(es9.2))') real(aa(-n2-1+i, 0:n2m1))
  ! end do

  ! write(*, '(/a)') '(2,1)-block'
  ! do i = 1,n2
  !    write(*, '(*(es9.2))') real(aa(i-1, -n2:-1))
  ! end do

  ! write(*, '(/a)') '(2,2)-block'
  ! do i = 1,n2
  !    write(*, '(*(es9.2))') real(aa(i-1, 0:n2m1))
  ! end do

  ! do i = 0,4
  !    do j = 0,4
  !       write(*, '(2x,i4,2x,i4,2x,2es26.16)') i, j, aa(i,j)-a1(i+sh,j+sh)
  !    end do
  ! end do
  ! print *, ''
  ! do i = 1,n2
  !    write(*, '(*(es10.2))') real(a1(i,:))
  ! end do
  ! print *, 'n = ', n
  ! print *, 'shape of aa', shape(aa)
  ! print *, 'shape of a1', shape(a1)
  ! write(*, '(a, 2x, es15.5)') &
  !      'norm of imaginary part', maxval( reshape(aimag(aa), (/size(aa)/)) )

  ! call dfftw_plan_dft_1d_ (plan_backward, n, b4(:,1), aa(:,1), &
  !      FFTW_BACKWARD, FFTW_ESTIMATE)
  ! call dfftw_execute_ (plan_backward)
  ! ! do i = 1,n
  ! !    write(*, '(2x, i4, 2es16.6)') i, aa(i, 1)-a0(2*i-1)
  ! ! end do
  ! write(*, '(/a, es26.16)') 'difference between a0 ptval and a0 from fft', &
  !      maxval(abs(aa(:,1)-a0((/(i, i=1,n2,2)/))))
  ! call dfftw_destroy_plan_ (plan_backward)

  ! call h4val(z, zd, h4, b4)
  ! call h4valnew(z(io), zd(io), h41, b41)
  ! write(*, '(/,2x,a,2x,es12.5)') 'infinity-norm of r_0*dth/dnu:', &
  !      vnorm( abs(zd-q0*dthdnu*(c0-u*jt)), 'inf' )
  ! open( unit=1, file='../tmp/h4.dat', status='replace')
  ! do i = 1,n
  !    do j = 1,n
  !       write(1, '(2es26.16)') h4(i,j)
  !    end do
  ! end do
  ! close(1)
  ! open( unit=1, file='../tmp/b4.dat', status='replace')
  ! do i = 1,n
  !    do j = 1,n
  !       write(1, '(2es26.16)') b4(i,j)
  !    end do
  ! end do
  ! close(1)
  ! open( unit=1, file='../tmp/b41.dat', status='replace')
  ! do i = 1,n
  !    do j = 1,n
  !       write(1, '(2es26.16)') b41(i,j)
  !    end do
  ! end do
  ! close(1)
  ! print *, 'difference between b4 and b41: ', mnorm(n,n,abs(b4-b41))
  ! print *, 'difference between b4 and b41: ', maxval(reshape(abs(b4-b41), (/size(b4)/)))

  ! open( unit=1, file='../tmp/jval.dat', status='replace')
  ! open( unit=2, file='../tmp/zeta.dat', status='replace')
  ! write(1, '(2es25.16E3)') jt
  ! write(2, '(2es25.16E3)') zeta
  ! close(1)
  ! close(2)

  ! open( unit=1, file='../tmp/jval.dat', status='replace' )
  ! open( unit=2, file='../tmp/zeta.dat', status='replace' )
  ! write(1, '(2es25.16)') jt
  ! write(2, '(2es25.16)') zeta
  ! close(1)
  ! write(*, *) shape(jt)


  !
  ! some checking
  !
  ! open( unit=1, file='../tmp/zdata.dat', status='unknown' )
  ! write(1, '(2es24.16)') z
  ! close(1)
  ! open( unit=1, file='../tmp/zddata.dat', status='unknown')
  ! write(1, '(2es24.16)') zd
  ! close(1)
  ! open( unit=1, file='../tmp/zetadata.dat', status='unknown')
  ! write(1, '(2es24.16)') zeta
  ! close(1)
  ! in = z(io)

  ! ! comparison
  ! open( unit=1, file='../tmp/adata.mat', status='unknown' )
  ! do i = 1,n
  !    read(1, *) a1(i)
  !    ! write(*, *)  a1(i)-a(i)
  ! end do
  ! print *, ''
  ! print *, vnorm(a-a1, 'inf')
  ! close(1)

  ! call ptval(npt, n, rho, beta, a1, eta, zeta, z1, zd)
  ! do i = 1,n
  !    write(*, '(2es24.16)') z(i)-z1(i)
  ! end do
  ! print *, ''
  ! print *, vnorm(abs(z-z1), 'inf')
