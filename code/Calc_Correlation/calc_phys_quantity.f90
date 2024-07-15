!
! calculate physical quantities 
!
!
! <H> = \sum_x p_x Eloc_x
!
!   p_x = | psi(x) |**2
!   Eloc_x = \sum_x'  <x|H|x'> * ( psi(x')/psi(x) )
!
!  
! gvec(k) = 2 <H*O_k> - 2 <H> <O_k>   
!   (for details, see e.g., Tahara and Imada, J. Phys. Soc. Jpn. 77, 114701 (2008)) 
! 
!   <H*O_k> = \sum_x p_x Eloc_x Ovec_loc_x(k) 
!   Ovec(k) = <O_k> = \sum_x p_x Ovec_loc_x(k) 
!   
!
! Smat(k1,k2) = <O_k1*O_k2> - <O_k1><O_k2>
!
subroutine calc_phys_quantity_opt(myrank,Etot,Etot2,gvec,Ovec,Smat,norm,wf_av)
  use mod_main   , only : N, Nv, N_QP_Full, N_MP_Trans, shift
  use mod_RBM    , only : Mc, Mr, alphac, alphar, Ncopy, W_map
  use mod_PP     , only : N_Slater, N_SP_Gauss_Leg
  use mod_sample 
  implicit none 
  integer    , intent(in)   :: myrank
  real(8)    , intent(out)  :: Etot, Etot2 
  real(8)    , intent(out)  :: gvec(Nv)
  real(8)    , intent(out)  :: Ovec(Nv)
  real(8)    , intent(out)  :: Smat(Nv,Nv)
  real(8)    , intent(out)  :: norm
  real(8)    , intent(out)  :: wf_av

  complex(8) :: psi_x(0:N_MP_Trans*Ncopy*2) 
  complex(8) :: psi_x_wo_pf(N_QP_Full)
  complex(8) :: thetac(Mc)
  real(8)    :: thetar(Mr)
  integer    :: Ele_Idx(N)
  integer    :: Ele_Cfg(N)
  integer    :: sgn_PP(N_MP_Trans)
  real(8)    :: Pf_M(N_QP_Full)
  real(8)    :: Inv_M(N,N,N_QP_Full)

  integer :: offset
  integer :: isample
  integer :: i, j, k, f, iw, icopy, itmp 
  integer :: k1, k2
  integer :: mpidx, spidx, qpidx, orbidx
  integer :: t1
  real(8) :: t(11)
  real(8) :: reweight_factor 
  real(8) :: deriv_PfM(N_Slater,N_QP_Full)

  complex(8) :: c1, c2, c3
  complex(8) :: Eloc_x
  complex(8) :: SdotS(N,4)
  complex(8) :: Ovec_loc_x(Nv)  
  complex(8) :: Ovec_Store(Nv,Nsample)  
  complex(8) :: OO(Nv,Nv)
   
  complex(8) :: ci = (0d0,1d0)

  Etot  = 0d0 
  Etot2 = 0d0
  norm  = 0d0
  wf_av = 0d0
  gvec(:) = 0d0
  Ovec(:) = 0d0
  SdotS(:,:) = 0d0
  !Smat(:,:) = 0d0
  t(:) = 0d0
  do isample = 1, Nsample
    ! 
    ! recalculate wave function 
    ! 
    call system_clock(t1)
    call calc_theta(x_sample(:,isample),thetac,thetar)
    call measure_time(t1,t(1))


    call system_clock(t1)
    call x_to_Ele_Idx_Ele_Cfg(x_sample(:,isample),Ele_Idx,Ele_Cfg)
    call check_Ele_Idx_Ele_Cfg(x_sample(:,isample),Ele_Idx,Ele_Cfg)
    call measure_time(t1,t(2))


    call system_clock(t1)
    call calculate_MAll(.true.,Ele_Idx,Pf_M,Inv_M)
    call measure_time(t1,t(3))


    call system_clock(t1)
    call sign_PP(Ele_Idx,sgn_PP)
    call measure_time(t1,t(4))


    call system_clock(t1)
    call calc_amplitude_RBM_PP(thetac,thetar,Pf_M,sgn_PP,psi_x_wo_pf,psi_x)
    call measure_time(t1,t(5))
    ! 
    ! calculate Eloc(x) and take sum over x for Etot
    ! H = J ( SzSz - SxSx - SySy ) after the gauge transformation
    !
    call system_clock(t1)
    reweight_factor = 1d0/(abs(psi_x(0))+shift)
    call calc_local_energy_correlation(x_sample(:,isample),psi_x(0),thetac,thetar,& 
                                       Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M,Eloc_x,SdotS)
    Eloc_x = Eloc_x * reweight_factor
    Etot   = Etot  + reweight_factor * dble(conjg(psi_x(0))*Eloc_x)  
    Etot2  = Etot2 + dble(conjg(Eloc_x)*Eloc_x)  
    norm   = norm + (abs(psi_x(0))*reweight_factor)**2
    wf_av  = wf_av + abs(psi_x(0))
    call measure_time(t1,t(6))
    !
    ! calculate Ovec_loc_x 
    !
    call system_clock(t1)
    Ovec_loc_x(:) = 0d0
    ! 
    !$omp parallel default(shared) & 
    !$omp private(f,k1,itmp,mpidx,spidx,icopy,j,offset,iw,i,c1,c2,c3)
    !$omp do 
    do f = 1, alphac
      k1 = 2*f-1
      itmp = 0 
      do mpidx = 1, N_MP_Trans 
      do spidx = 1, -1, -2
      do icopy = 1, Ncopy
        j = itmp * alphac + f
        itmp = itmp + 1 
        c2 = dcmplx(tanh(dble(thetac(j))),tan(dimag(thetac(j))))
        c3 = dcmplx(1d0,tanh(dble(thetac(j)))*tan(dimag(thetac(j))))
        Ovec_loc_x(k1) = Ovec_loc_x(k1) + (c2/c3)*psi_x(itmp)  
      end do ! icopy
      end do ! spidx
      end do ! mpidx
      if( itmp /= N_MP_Trans*Ncopy*2 ) stop 'itmp /= N_MP_Trans*Ncopy*2'  
      Ovec_loc_x(k1+1) = Ovec_loc_x(k1)*ci
    end do ! f 
    !$omp end do 

    !$omp do 
    do f = 1, alphac
      offset = 2*alphac + (f-1)*N*2
      itmp = 0 
      do mpidx = 1, N_MP_Trans 
      do spidx = 1, -1, -2
      do icopy = 1, Ncopy
        j = itmp * alphac + f
        itmp = itmp + 1 
        c2 = dcmplx(tanh(dble(thetac(j))),tan(dimag(thetac(j))))
        c3 = dcmplx(1d0,tanh(dble(thetac(j)))*tan(dimag(thetac(j))))
        c1 = (c2/c3)*psi_x(itmp)
        do i = 1, N  
          iw = W_map(i,icopy,mpidx)
          k1 = offset + 2*iw - 1   !!! please carefully check  
          Ovec_loc_x(k1) = Ovec_loc_x(k1) + c1 * dble(x_sample(i,isample)) * dble(spidx)
          if( mpidx == N_MP_Trans .and. icopy == Ncopy ) Ovec_loc_x(k1+1) = Ovec_loc_x(k1) * ci
        end do ! i
      end do ! icopy
      end do ! spidx 
      end do ! mpidx 
      if( itmp /= N_MP_Trans*Ncopy*2 ) stop 'itmp /= N_MP_Trans*Ncopy*2'  
    end do ! f 
    !$omp end do 

    !$omp do 
    do f = 1, alphar
      k1 = 2*alphac*(N+1) + f
      itmp = 0 
      do mpidx = 1, N_MP_Trans 
      do spidx = 1, -1, -2
      do icopy = 1, Ncopy
        j = itmp * alphar + f
        itmp = itmp + 1 
        Ovec_loc_x(k1) = Ovec_loc_x(k1) + tanh(thetar(j))*psi_x(itmp)
      end do ! icopy
      end do ! spidx
      end do ! mpidx
      if( itmp /= N_MP_Trans*Ncopy*2 ) stop 'itmp /= N_MP_Trans*Ncopy*2'  
    end do ! f 
    !$omp end do 

    !$omp do 
    do f = 1, alphar
      offset = 2*alphac*(N+1) + alphar + (f-1)*N
      itmp = 0 
      do mpidx = 1, N_MP_Trans 
      do spidx = 1, -1, -2
      do icopy = 1, Ncopy
        j = itmp * alphar + f
        itmp = itmp + 1 
        c1 = tanh(thetar(j))*psi_x(itmp)
        do i = 1, N  
          iw = W_map(i,icopy,mpidx)
          k1 = offset + iw 
          Ovec_loc_x(k1) = Ovec_loc_x(k1) + c1 * dble(x_sample(i,isample)) * dble(spidx)
        end do ! i
      end do ! icopy
      end do ! spidx
      end do ! mpidx 
      if( itmp /= N_MP_Trans*Ncopy*2 ) stop 'itmp /= N_MP_Trans*Ncopy*2'  
    end do ! f 
    !$omp end do 
    !$omp end parallel 
    call measure_time(t1,t(7))


    call system_clock(t1)
    offset = 2*alphac*(N+1) + alphar*(N+1)
    call calculate_derivative_Pfaffian(deriv_PfM,Ele_Idx,Pf_M,Inv_M)
    do orbidx = 1, N_Slater 
      qpidx = 0
      k1 = offset + orbidx
      do mpidx = 1, N_MP_Trans
      do spidx = 1, N_SP_Gauss_Leg
        qpidx = qpidx + 1 
        Ovec_loc_x(k1) = Ovec_loc_x(k1) + deriv_PfM(orbidx,qpidx) * psi_x_wo_pf(qpidx) 
      end do ! spidx 
      end do ! mpidx
      if( qpidx /= N_QP_Full ) stop 'qpidx /= N_QP_Full'
    end do ! orbidx 
    if( offset + N_Slater /= Nv ) stop 'offset is wrong 3'
    Ovec_loc_x(:) = Ovec_loc_x(:) * reweight_factor
    call measure_time(t1,t(8))
    !
    ! summation over x for gvec and Ovec and Ovec
    !
    call system_clock(t1)
    !$omp parallel do default(shared) private(k)
    do k = 1, Nv
      gvec(k) = gvec(k) + 2d0 * conjg(Eloc_x) * Ovec_loc_x(k)  
      Ovec(k) = Ovec(k) + reweight_factor * dble(conjg(psi_x(0))*Ovec_loc_x(k))   
      Ovec_Store(k,isample) = Ovec_loc_x(k)  
    end do ! k 
    !$omp end parallel do 
    call measure_time(t1,t(9))

  end do ! isample

  call system_clock(t1)
  call zgemm('N','C',Nv,Nv,Nsample,(1d0,0d0),Ovec_Store,Nv,Ovec_Store,Nv,(0d0,0d0),OO,Nv)
  call measure_time(t1,t(10))


  call system_clock(t1)
  !$omp parallel do default(shared) private(k1,k2)
  do k2 = 1, Nv
  do k1 = 1, Nv
    Smat(k1,k2) = dble(OO(k1,k2))
    if( k1 == k2 .and. dimag(OO(k1,k2)) > 1d-5 ) stop 'OO is wrong'  
  end do ! k1 
  end do ! k2 
  !$omp end parallel do
  call measure_time(t1,t(11))

  if( mod(myrank,100) == 0 ) write(50000+myrank,'(11F10.6)') t(:)

  return 
end subroutine 
!
!
!
subroutine calc_phys_quantity_after_opt(Etot,Etot2,Chi_sp,Chi_vbs,norm,wf_av,wf_sgn)
  use mod_main   , only : N, Nchi, N_QP_Full, N_MP_Trans, shift
  use mod_RBM    , only : Mc, Mr, Ncopy
  use mod_sample 
  implicit none 
  real(8)    , intent(out)  :: Etot, Etot2 
  real(8)    , intent(out)  :: Chi_sp (2,Nchi)
  real(8)    , intent(out)  :: Chi_vbs(Nchi,4)
  real(8)    , intent(out)  :: norm
  real(8)    , intent(out)  :: wf_av
  complex(8) , intent(out)  :: wf_sgn
  
  integer    :: isample
  integer    :: i0, i1, i2, i3 
  complex(8) :: Eloc_x
  complex(8) :: SdotS(N,4)
  complex(8) :: Chi_loc_x(2,Nchi)
  real(8)    :: reweight_factor 

  complex(8) :: psi_x(0:N_MP_Trans*Ncopy*2) 
  complex(8) :: psi_x_wo_pf(N_QP_Full) 
  complex(8) :: thetac(Mc)
  real(8)    :: thetar(Mr)
  integer    :: Ele_Idx(N)
  integer    :: Ele_Cfg(N)
  integer    :: sgn_PP(N_MP_Trans)
  real(8)    :: Pf_M(N_QP_Full)
  real(8)    :: Inv_M(N,N,N_QP_Full)

  Etot  = 0d0
  Etot2 = 0d0
  SdotS  (:,:) = 0d0
  Chi_sp (:,:) = 0d0
  Chi_vbs(:,:) = 0d0
  norm  = 0d0
  wf_av = 0d0
  wf_sgn = 0d0
  do isample = 1, Nsample
    ! 
    ! recalculate wave function 
    ! 
    call calc_theta(x_sample(:,isample),thetac,thetar)
    call x_to_Ele_Idx_Ele_Cfg(x_sample(:,isample),Ele_Idx,Ele_Cfg)
    call check_Ele_Idx_Ele_Cfg(x_sample(:,isample),Ele_Idx,Ele_Cfg)
    call calculate_MAll(.true.,Ele_Idx,Pf_M,Inv_M)
    call sign_PP(Ele_Idx,sgn_PP)
    call calc_amplitude_RBM_PP(thetac,thetar,Pf_M,sgn_PP,psi_x_wo_pf,psi_x)
    ! 
    ! calculate Eloc(x) and take sum over x for Etot
    ! H = J ( SzSz - SxSx - SySy ) after the gauge transformation
    !
    reweight_factor = 1d0/(abs(psi_x(0))+shift)
    call calc_local_energy_correlation(x_sample(:,isample),psi_x(0),thetac,thetar,& 
                                       Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M,Eloc_x,SdotS)
    Eloc_x = Eloc_x * reweight_factor
    Etot   = Etot  + reweight_factor * dble(conjg(psi_x(0))*Eloc_x)  
    Etot2  = Etot2 + dble(conjg(Eloc_x)*Eloc_x)  
    norm   = norm + (abs(psi_x(0))*reweight_factor)**2
    wf_av  = wf_av + abs(psi_x(0))
    wf_sgn = wf_sgn + psi_x(0)/abs(psi_x(0))
    SdotS(:,:) = SdotS(:,:) * reweight_factor
    do i3 = 1, 4  
      i0 = 0 
      do i2 = 1, N
      do i1 = 1, i2
        i0 = i0 + 1 
        Chi_vbs(i0,i3) = Chi_vbs(i0,i3) + dble(conjg(SdotS(i1,i3))*SdotS(i2,i3)) 
      end do ! i2 
      end do ! i1 
      if( i0 /= Nchi ) stop 'i0/=Nchi'
    end do ! i3 
    ! 
    ! calclate correlation function
    ! 
    call calc_local_correlation(x_sample(:,isample),psi_x(0),thetac,thetar,& 
                                Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M,Chi_loc_x)
    Chi_sp(:,:) = Chi_sp(:,:) + reweight_factor * reweight_factor * dble(conjg(psi_x(0))*Chi_loc_x(:,:))  
  end do ! isample

  return
end subroutine 
! 
! 
! 
subroutine calc_local_energy_correlation(x,psi_x,thetac,thetar,Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M,Eloc_x,SdotS)
  use mod_main , only : N, N_MP_Trans, N_QP_Full, J1, J2, subB_or_not 
  use mod_RBM  , only : Mc, Mr, Ncopy
  implicit none
  integer    , intent(in)  :: x(N)
  complex(8) , intent(in)  :: psi_x 
  complex(8) , intent(in)  :: thetac(Mc)
  real(8)    , intent(in)  :: thetar(Mr)
  integer    , intent(in)  :: Ele_Idx(N)
  integer    , intent(in)  :: Ele_Cfg(N)
  integer    , intent(in)  :: sgn_PP(N_MP_Trans)
  real(8)    , intent(in)  :: Pf_M(N_QP_Full)
  real(8)    , intent(in)  :: Inv_M(N,N,N_QP_Full)
  complex(8) , intent(out) :: Eloc_x
  complex(8) , intent(out) :: SdotS(N,4)

  integer    :: i0, i1, i2, i3, i4, jtmp 
  complex(8) :: psi_xp(0:N_MP_Trans*Ncopy*2)
  complex(8) :: psi_xp_wo_pf(N_QP_Full)
  complex(8) :: thetapc(Mc)
  real(8)    :: thetapr(Mr)
  integer    :: Tmp_Ele_Idx(N)
  real(8)    :: Tmp_Pf_M(N_QP_Full)
  integer    :: tmp_sgn_PP(N_MP_Trans)
  real(8)    :: J, sgn

  Eloc_x = 0d0
  SdotS(:,:) = 0d0
  do i0 = 1, 4
    if( i0 <= 2 ) then 
      J = J1; sgn = -1d0
    else
      J = J2; sgn =  1d0
    end if
    do i1 = 1, N
      call bond_map(i0,i1,i2)
      !
      ! SzSz contribution 
      ! 
      Eloc_x = Eloc_x + J*x(i1)*x(i2)*psi_x 
      SdotS(i1,i0) = SdotS(i1,i0) + x(i1)*x(i2)*psi_x 
      !
      ! ( SxSx + SySy ) = 2 * ( S+S- + S-S+ ) contribution 
      ! 
      if( x(i1) /= x(i2) ) then 

        !$omp parallel do default(shared) private(jtmp) 
        do jtmp = 1, Mc 
          thetapc(jtmp) = thetac(jtmp)
        end do ! jtmp 
        !$omp end parallel do 
        thetapr(:) = thetar(:)
        call update_theta(x,i1,i2,thetapc,thetapr)

        i3 = Ele_Cfg(i1)
        i4 = Ele_Cfg(i2)
        if( Ele_Idx(i3) /= i1 ) stop 'Ele_Idx(i3) /= i1'
        if( Ele_Idx(i4) /= i2 ) stop 'Ele_Idx(i4) /= i2'
        Tmp_Ele_Idx(:) = Ele_Idx(:)
        Tmp_Ele_Idx(i3) = i2
        Tmp_Ele_Idx(i4) = i1
        call calculate_new_PfM_two(i3,i4,Tmp_Ele_Idx,Inv_M,Pf_M,Tmp_Pf_M)

        tmp_sgn_PP(:) = -sgn_PP(:)
        if( subB_or_not(i1) ) tmp_sgn_PP(:) = -tmp_sgn_PP(:) 
        if( subB_or_not(i2) ) tmp_sgn_PP(:) = -tmp_sgn_PP(:) 
        call calc_amplitude_RBM_PP(thetapc,thetapr,Tmp_Pf_M,tmp_sgn_PP,psi_xp_wo_pf,psi_xp)

        Eloc_x = Eloc_x + 2d0*J*sgn*psi_xp(0)
        SdotS(i1,i0) = SdotS(i1,i0) + 2d0*sgn*psi_xp(0)

      end if
    end do ! i1      
  end do ! i0  

  return 
end subroutine 
!
!
!
subroutine bond_map(i0,i1,i2) 
  use mod_main , only : L, N
  implicit none 
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
!
! calculate local correlation function
!
subroutine calc_local_correlation(x,psi_x,thetac,thetar,Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M,Chi_loc_x)
  use mod_main , only : N, Nchi, N_MP_Trans, N_QP_Full, subB_or_not 
  use mod_RBM  , only : Mc, Mr, Ncopy
  implicit none
  integer    , intent(in)  :: x(N)
  complex(8) , intent(in)  :: psi_x 
  complex(8) , intent(in)  :: thetac(Mc)
  real(8)    , intent(in)  :: thetar(Mr)
  integer    , intent(in)  :: Ele_Idx(N)
  integer    , intent(in)  :: Ele_Cfg(N)
  integer    , intent(in)  :: sgn_PP(N_MP_Trans)
  real(8)    , intent(in)  :: Pf_M(N_QP_Full)
  real(8)    , intent(in)  :: Inv_M(N,N,N_QP_Full)
  complex(8) , intent(out) :: Chi_loc_x(2,Nchi)

  integer    :: i0, i1, i2, i3, i4, jtmp 
  complex(8) :: psi_xp(0:N_MP_Trans*Ncopy*2)
  complex(8) :: psi_xp_wo_pf(N_QP_Full)
  complex(8) :: thetapc(Mc)
  real(8)    :: thetapr(Mr)
  integer    :: Tmp_Ele_Idx(N)
  real(8)    :: Tmp_Pf_M(N_QP_Full)
  integer    :: tmp_sgn_PP(N_MP_Trans)
  real(8)    :: sgn

  i0 = 0 
  do i2 = 1, N
  do i1 = 1, i2
    i0 = i0 + 1
    !
    ! SzSz contribution 
    ! 
    Chi_loc_x(1,i0) = x(i1)*x(i2)*psi_x 
    !
    ! ( SxSx + SySy ) = 2 * ( S+S- + S-S+ ) contribution 
    ! 
    Chi_loc_x(2,i0) = 0d0
    if( i1 == i2 ) Chi_loc_x(2,i0) = psi_x
    if( x(i1) /= x(i2) ) then 

      !$omp parallel do default(shared) private(jtmp) 
      do jtmp = 1, Mc 
        thetapc(jtmp) = thetac(jtmp)
      end do ! jtmp 
      !$omp end parallel do 
      thetapr(:) = thetar(:)
      call update_theta(x,i1,i2,thetapc,thetapr)

      i3 = Ele_Cfg(i1)
      i4 = Ele_Cfg(i2)
      if( Ele_Idx(i3) /= i1 ) stop 'Ele_Idx(i3) /= i1'
      if( Ele_Idx(i4) /= i2 ) stop 'Ele_Idx(i4) /= i2'
      Tmp_Ele_Idx(:) = Ele_Idx(:)
      Tmp_Ele_Idx(i3) = i2
      Tmp_Ele_Idx(i4) = i1
      call calculate_new_PfM_two(i3,i4,Tmp_Ele_Idx,Inv_M,Pf_M,Tmp_Pf_M)

      tmp_sgn_PP(:) = -sgn_PP(:)
      if( subB_or_not(i1) ) tmp_sgn_PP(:) = -tmp_sgn_PP(:) 
      if( subB_or_not(i2) ) tmp_sgn_PP(:) = -tmp_sgn_PP(:) 
      call calc_amplitude_RBM_PP(thetapc,thetapr,Tmp_Pf_M,tmp_sgn_PP,psi_xp_wo_pf,psi_xp)

      sgn = 1 
      if( subB_or_not(i1) ) sgn = -sgn
      if( subB_or_not(i2) ) sgn = -sgn 
      Chi_loc_x(2,i0) = sgn*psi_xp(0)

    end if
  end do ! i1      
  end do ! i2  
  if( i0 /= Nchi ) stop 'i0/=Nchi'

  return 
end subroutine 
