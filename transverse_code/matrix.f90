!
! calculate pfaffian from scratch
! we also calculate inverse matrix when calc_inv = .true.
! algorithm taken from mVMC package
!
subroutine calculate_MAll(calc_inv,Tmp_Ele_Idx,pfaff,invM)
  use mod_main , only : N, N_QP_Full, N_MP_Trans
  use mod_PP   , only : Npair, N_SP_Gauss_Leg, Slater_Elm, lwork_pfa 
  implicit none 
  logical , intent(in)  :: calc_inv
  integer , intent(in)  :: Tmp_Ele_Idx(N) 
  real(8) , intent(out) :: pfaff(N_QP_Full)
  real(8) , intent(out) :: invM(N,N,N_QP_Full)
  integer :: mpidx, spidx, qpidx
  integer :: msi, msj
  integer :: rsi, rsj
  real(8) :: bufM1(N,N)

  if( N /= 2*Npair ) stop 'N/=2*Npair, not implemented'
  ! 
  ! set lwork_pfa
  ! 
  if( lwork_pfa == 0 ) then 
    do msj = 1, N 
    do msi = 1, N 
      rsj = Tmp_Ele_Idx(msj) + ((msj-1)/Npair)*N 
      rsi = Tmp_Ele_Idx(msi) + ((msi-1)/Npair)*N 
      bufM1(msi,msj) = -Slater_Elm(rsi,rsj,1) !!! = Slter_Elm(rsj,rsi,1)
    end do ! msi 
    end do ! msj 
    call get_lwork_pfa(N,lwork_pfa,bufM1,pfaff(1))
    if( lwork_pfa < 0 ) stop 'lwork_pfa < 0' 
  else
    if( lwork_pfa < 0 ) stop 'lwork_pfa < 0' 
  end if

!$omp parallel do default(shared) private(qpidx,mpidx,spidx)
  do mpidx = 1, N_MP_Trans
    qpidx = (mpidx-1)*N_SP_Gauss_Leg 
    do spidx = 1, N_SP_Gauss_Leg
      qpidx = qpidx + 1 
      call calculate_MAll_child(qpidx,calc_inv,Tmp_Ele_Idx,pfaff(qpidx),invM(:,:,qpidx))
    end do ! spidx
    if( qpidx /= mpidx*N_SP_Gauss_Leg ) stop 'qpidx /= mpidx*N_SP_Gauss_Leg'
  end do ! mpidx
!$omp end parallel do

  return 
end subroutine  
!
! calculate pfaffian from scratch
! we also calculate inverse matrix when calc_inv = .true.
! algorithm taken from mVMC package
!
subroutine calculate_MAll_child(qpidx,calc_inv,Tmp_Ele_Idx,pfaff,invM)
  use mod_main , only : N
  use mod_PP   , only : Npair, Slater_Elm, lwork_pfa 
  implicit none 
  integer , intent(in)  :: qpidx
  logical , intent(in)  :: calc_inv
  integer , intent(in)  :: Tmp_Ele_Idx(N) 
  real(8) , intent(out) :: pfaff
  real(8) , intent(out) :: invM(N,N)
  integer :: msi, msj
  integer :: rsi, rsj
  real(8) :: bufM1(N,N)
  real(8) :: bufM2(N,N)
  real(8) :: r1, r2, r3
  ! 
  ! preparation of matrix
  ! 
  do msj = 1, N 
  do msi = 1, N 
    rsj = Tmp_Ele_Idx(msj) + ((msj-1)/Npair)*N 
    rsi = Tmp_Ele_Idx(msi) + ((msi-1)/Npair)*N 
    bufM1(msi,msj) = -Slater_Elm(rsi,rsj,qpidx) !!! = Slter_Elm(rsj,rsi,qpidx)
  end do ! msi 
  end do ! msj 

  bufM2(:,:) = bufM1(:,:) 

  ! 
  ! calculate pfaffian
  ! 
  call calc_pfaffian(N,lwork_pfa,bufM2,pfaff)

  if(.not.calc_inv) return
  if(.not.calc_inv) stop 'invM is calculated although calc_inv = .false.'
  ! 
  ! calculate inverse matrix 
  !
  call inv(N,bufM1)
  ! 
  ! check and symmetrize
  ! 
  r1 = maxval(abs(bufM1))
  do msj = 1, N 
  do msi = 1, N 
    r2 = 0.5d0*(bufM1(msj,msi)-bufM1(msi,msj))
    r3 = 0.5d0*(bufM1(msj,msi)+bufM1(msi,msj))
    if( abs(r3/r1) > 1.0d-4 ) then 
      write(6,'(a,3E13.4)') 'Warning: InvM is not skew-symmetric (calculate_MAll)', r1, r2, r3
      if( abs(r3/r1) > 2.0d-3 ) stop 'InvM is not skew-symmetric (calculate_MAll)'
    end if 
    invM(msi,msj) = r2
  end do ! msi 
  end do ! msj 

  return 
end subroutine  
! 
! get optimal lwork for pfaffian package  
!
subroutine get_lwork_pfa(N,lwork,mat,pfaff)
  implicit none 
  integer , intent(in)    :: N 
  integer , intent(out)   :: lwork
  real(8) , intent(inout) :: mat(N,N)
  real(8) , intent(out)   :: pfaff
  integer :: iwork(N)
  real(8) :: work(1)
  integer :: info

  lwork = -1
  call dskpfa('U','P',N,mat,N,pfaff,iwork,work,lwork,info)
  lwork = idnint(work(1))  

  return 
end subroutine
! 
! calculate pfaffian using pfapack
!
subroutine calc_pfaffian(N,lwork,mat,pfaff)
  implicit none 
  integer , intent(in)    :: N 
  integer , intent(in)    :: lwork
  real(8) , intent(inout) :: mat(N,N)
  real(8) , intent(out)   :: pfaff
  integer :: iwork(N)
  real(8) :: work(lwork)
  integer :: info

  call dskpfa('U','P',N,mat,N,pfaff,iwork,work,lwork,info)

  if( info /= 0 ) then 
    write(6,*) 'info (subroutine calc_pfaffian):' , info
    stop 
  end if 

  return 
end subroutine 
!
! calculate inverse of matrix
!
subroutine inv(nm,mat)
  implicit none 
  integer , intent(in) :: nm
  real(8) , intent(inout) :: mat(nm,nm)
  integer :: ipiv(nm)
  integer :: Lwork 
  real(8) , allocatable :: work(:)
  integer :: info 

  Lwork = 10*nm
  allocate (work(Lwork))
  info = 0
  call dgetrf(nm,nm,mat,nm,ipiv,info)
  call dgetri(nm,mat,nm,ipiv,work,Lwork,info)


  if(info /= 0) then
    write(6,*) 'info (subrouitine inv):' , info
    stop
  end if 
  deallocate(work)
  return 
end subroutine
