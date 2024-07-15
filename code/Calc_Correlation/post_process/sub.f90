! 
! calculating average and error
! 
subroutine calc_av_er(N,data_,av,er) 
  implicit none 
  integer , intent(in)  :: N 
  real(8) , intent(in)  :: data_(N)
  real(8) , intent(out) :: av, er
  integer :: i

  av = 0d0
  er = 0d0
  do i = 1, N
    av = av + data_(i)
  end do ! i
  av = av / dble(N)

  do i = 1, N
    er = er + (data_(i)-av)**2
  end do ! i 
  er = er / dble(N) / dble(N-1)
  er = dsqrt(er)

  return 
end subroutine 
! 
!
!
subroutine bond_map(L,N,i0,i1,i2) 
  implicit none 
  integer , intent(in)  :: L, N 
  integer , intent(in)  :: i0
  integer , intent(in)  :: i1
  integer , intent(out) :: i2

  integer :: ix, iy 
  ! 
  ! i0: 1-2 => NN  (J1) 
  ! i0: 3-4 => NNN (J2) 
  ! 
  ix = (i1-1) / L 
  iy = (i1-1) - ix * L  
  if( ix < 0  ) stop 'ix<0' 
  if( ix >= L ) stop 'ix>=L' 
  if( iy < 0  ) stop 'iy<0' 
  if( iy >= L ) stop 'iy>=L' 

  if( (i1-1) /= ix*L+iy ) stop 'i1 wrong'

  if( i0 == 1 ) then 
    ix = ix + 1    
    if( ix >= L ) ix = ix - L 
  else if( i0 == 2 ) then 
    iy = iy + 1    
    if( iy >= L ) iy = iy - L 
  else if( i0 == 3 ) then 
    ix = ix + 1    
    if( ix >= L ) ix = ix - L 
    iy = iy + 1    
    if( iy >= L ) iy = iy - L 
  else 
    if( i0 /= 4 ) stop 'i0/=4' 
    ix = ix + 1    
    if( ix >= L ) ix = ix - L 
    iy = iy - 1    
    if( iy < 0  ) iy = iy + L 
  end if

  i2 = ix*L + iy + 1 
  if( i2 < 1 ) stop 'i2<1' 
  if( i2 > N ) stop 'i2>N' 

  return 
end subroutine 
