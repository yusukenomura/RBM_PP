!
! spin configuration to EleIdx, Ele_Cfg 
! this subroutine is needed to fit the mVMC code
!
subroutine x_to_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg) 
  use mod_main , only : N
  use mod_PP   , only : Npair
  implicit none
  integer , intent(in)  :: x(N) 
  integer , intent(out) :: Ele_Idx(N) 
  integer , intent(out) :: Ele_Cfg(N) 
  integer :: i 
  integer :: i1, i2 

  i1 = 0; i2 = Npair
  do i = 1, N 
    if( x(i) == 1 ) then 
      i1 = i1 + 1
      Ele_Idx(i1) = i
      Ele_Cfg(i) = i1
    else if( x(i) == -1 ) then 
      i2 = i2 + 1
      Ele_Idx(i2) = i
      Ele_Cfg(i) = i2
    else 
      stop 'x is wrong (x_to_Ele_Idx_Ele_Cfg)'
    end if 
  end do ! i 
  if( i1 /= Npair ) stop 'i1/=Npair (x_to_Ele_Idx_Ele_Cfg)'
  if( i2 /= N     ) stop 'i2/=N (x_to_Ele_Idx_Ele_Cfg)'
  
  return 
end subroutine
!
! check Ele_Idx and Ele_Cfg
!
subroutine check_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg)
  use mod_main , only : N
  use mod_PP   , only : Npair
  implicit none
  integer , intent(in) :: x(N) 
  integer , intent(in) :: Ele_Idx(N) 
  integer , intent(in) :: Ele_Cfg(N) 
  integer :: i 
  logical :: chk(N) 

  chk(:) = .false.
  do i = 1, Npair
    if( Ele_Cfg(Ele_Idx(i)) /= i ) stop 'error1 (check_Ele_Idx)'
    if( x(Ele_Idx(i)) /= 1 ) stop 'error2 (check_Ele_Idx)'
    chk(Ele_Idx(i)) = .true. 
  end do ! i  

  do i = Npair+1, N
    if( Ele_Cfg(Ele_Idx(i)) /= i ) stop 'error3 (check_Ele_Idx)'
    if( x(Ele_Idx(i)) /= -1 ) stop 'error4 (check_Ele_Idx)'
    chk(Ele_Idx(i)) = .true. 
  end do ! i  

  do i = 1, N 
    if(.not.chk(i)) stop 'error5 (check_Ele_Idx)'
  end do ! i

  return 
end subroutine
!
! Depending on the order of fermionic operator, sign changes
! In spin model, we use the order: c^dagger_{1 s1} c^dagger_{2 s2} ...
! This subroutine gives additional sign coming from order exchange  
!
subroutine sign_PP(Ele_Idx,sgn)
  use mod_main , only : N, N_MP_Trans, subB_or_not
  use mod_PP   , only : Npair, MP_Trans_Idx
  implicit none 
  integer , intent(in)  :: Ele_Idx(N) 
  integer , intent(out) :: sgn(N_MP_Trans) 
  integer :: i1, i2, mpidx
  integer :: Tmp_Ele_Idx(N)
  integer :: sgn_tmp(N_MP_Trans)

!$omp parallel do default(shared) private(mpidx,i1,i2,Tmp_Ele_Idx)
  do mpidx = 1, N_MP_Trans
    do i1 = 1, N 
      i2 = Ele_Idx(i1)
      Tmp_Ele_Idx(i1) = MP_Trans_Idx(i2,mpidx)
    end do ! i1
    sgn(mpidx) = 1 
    sgn_tmp(mpidx) = 1 
    do i1 = 2, N 
      if( i1 > Npair ) then 
        if( subB_or_not(Tmp_Ele_Idx(i1)) ) sgn(mpidx) = -sgn(mpidx)
        if( subB_or_not(Tmp_Ele_Idx(i1)) ) sgn_tmp(mpidx) = -sgn_tmp(mpidx)
      end if 
      do i2 = 1, i1-1 
        if( Tmp_Ele_Idx(i2) > Tmp_Ele_Idx(i1) ) sgn(mpidx) = -sgn(mpidx)
      end do ! i2
    end do ! i1
  end do ! mpidx
!$omp end parallel do

  do mpidx = 2, N_MP_Trans
    if( sgn_tmp(mpidx) /= sgn_tmp(1) ) stop 'sgn_tmp is wrong'
  end do ! mpidx

  return 
end subroutine
! 
! define sublattice 
! 
subroutine define_sublattice()
  use mod_main, only : L, N, subB_or_not
  implicit none
  integer :: i1, i2
  integer :: i

  i = 0
  do i2 = 0, L-1
  do i1 = 0, L-1
    i = i + 1
    if( mod(i1+i2,2) == 1 ) then 
      subB_or_not(i) = .true.
    else 
      subB_or_not(i) = .false.
    end if
  end do ! i1
  end do ! i2
  if( i /= N ) stop 'i/=N (map_2D)'

  i1 = 0 
  do i = 1, N 
    if( subB_or_not(i) ) i1 = i1 + 1
  end do ! i 
  if( i1 /= N/2 ) stop 'i1 /= N/2 (map_2D)'

  return
end subroutine

