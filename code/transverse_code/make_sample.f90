!
! make sample  
!
subroutine make_sample(myrank,first_iteration,naccept)
  use mod_main   , only : N, N_QP_Full, N_MP_Trans, Nupdate, subB_or_not, shift 
  use mod_RBM    , only : Mc, Mr, Ncopy
  use mod_PP     , only : Nfastup
  use mod_sample 
  implicit none 
  integer    , intent(in)    :: myrank 
  logical    , intent(in)    :: first_iteration
  integer    , intent(out)   :: naccept

  integer    :: x(N) 
  complex(8) :: psi_x(0:N_MP_Trans*Ncopy*2) 
  complex(8) :: psi_x_wo_pf(N_QP_Full) 
  complex(8) :: thetac(Mc)
  real(8)    :: thetar(Mr)
  integer    :: Ele_Idx(N)
  integer    :: Ele_Cfg(N)
  integer    :: sgn_PP(N_MP_Trans)
  real(8)    :: Pf_M(N_QP_Full)
  real(8)    :: Inv_M(N,N,N_QP_Full)

  complex(8) :: psi_xp(0:N_MP_Trans*Ncopy*2)
  complex(8) :: psi_xp_wo_pf(N_QP_Full)
  complex(8) :: thetapc(Mc)
  real(8)    :: thetapr(Mr)
  integer    :: Tmp_Ele_Idx(N)
  integer    :: Tmp_Ele_Cfg(N)
  integer    :: tmp_sgn_PP(N_MP_Trans)
  real(8)    :: Tmp_Pf_M(N_QP_Full)
  real(8)    :: Tmp_Inv_M(N,N,N_QP_Full)

  integer :: Nwarmup
  integer :: i0, i1, i2, i3, i4, j
  integer :: mpidx, qpidx  
  integer :: iupdate
  integer :: nacc_tmp
  integer :: Ntmp_accept
  integer :: t1
  real(8) :: r1, r2, r3
  real(8) :: t(10)
  real(8) :: diff_InvM, diff_PfM
  real(8) , external :: grnd 

  t(:) = 0d0

  call system_clock(t1)
  if( first_iteration ) then 
    call make_initial_x(x,psi_x,thetac,thetar,Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M) 
    if( mod(myrank,100) == 0 ) write(6,'(a,I12,2E15.5)') 'myrank, initial psi_x =', myrank, psi_x(0)
  else 
    x(:) = x_sample(:,Nsample)
    call calc_theta(x,thetac,thetar)
    call x_to_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg)
    call check_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg)
    call calculate_MAll(.true.,Ele_Idx,Pf_M,Inv_M)
    call sign_PP(Ele_Idx,sgn_PP)
    call calc_amplitude_RBM_PP(thetac,thetar,Pf_M,sgn_PP,psi_x_wo_pf,psi_x)
    if( grnd() < 0.05d0 ) call make_initial_x(x,psi_x,thetac,thetar,Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M) 
  end if 
  call measure_time(t1,t(1))


  Nwarmup = max(30,Nsample/100)

  naccept = 0
  nacc_tmp = 0
  Ntmp_accept = 0
  do i0 = -Nwarmup, Nsample

    call system_clock(t1)
    if( nacc_tmp > Nfastup ) then 
      diff_InvM = -1.0d0 
      diff_PfM  = -1.0d0 
  
      thetapc(:) = thetac(:)        
      thetapr(:) = thetar(:)        
      tmp_sgn_PP(:) = sgn_PP(:)
      call calc_theta(x,thetac,thetar)             
      call check_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg)
      call sign_PP(Ele_Idx,sgn_PP)
      do j = 1, Mc 
        if( abs(thetac(j)-thetapc(j)) > 1.0d-5 ) stop 'thetac is wrong'
      end do ! j  
      do j = 1, Mr
        if( abs(thetar(j)-thetapr(j)) > 1.0d-5 ) stop 'thetar is wrong'
      end do ! j
      do mpidx = 1, N_MP_Trans
        if( tmp_sgn_PP(mpidx) /= sgn_PP(mpidx) ) stop 'sgn_PP is wrong'
      end do ! mpidx 

      Tmp_Inv_M(:,:,:) = Inv_M(:,:,:)
      Tmp_Pf_M(:) = Pf_M(:)
      call calculate_MAll(.true.,Ele_Idx,Pf_M,Inv_M)
      call calc_amplitude_RBM_PP(thetac,thetar,Pf_M,sgn_PP,psi_x_wo_pf,psi_x) 
      r2 = 0d0 
      do qpidx = 1, N_QP_Full
        r2 = r2 + abs(Pf_M(qpidx))
      end do ! qpidx
      r2 = r2 / dble(N_QP_Full)
      do qpidx = 1, N_QP_Full
        diff_PfM = max( diff_PfM, abs(Tmp_Pf_M(qpidx)-Pf_M(qpidx))/r2 ) 
      end do ! qpidx
      do qpidx = 1, N_QP_Full
        r3 = maxval(abs(Inv_M(:,:,qpidx)))
        do i2 = 1, N
        do i1 = 1, N
          diff_InvM = max( diff_InvM, abs(Tmp_Inv_M(i1,i2,qpidx)-Inv_M(i1,i2,qpidx))/r3 ) 
        end do ! i1
        end do ! i2
      end do ! qpidx

      nacc_tmp = 0
      if( grnd() < 1.0d-5 ) & 
        write(6,'(a,2I8,3E12.4)') 'myrank, Nfastup, diff_InvM, diff_PfM, av_PfM', &  
                                   myrank, Nfastup, diff_InvM, diff_PfM, r2

      if( diff_InvM > 1.0d-5 .or. diff_PfM > 1.0d-5 ) then 
        Nfastup = max(1,idnint(dble(Nfastup)/1.2d0)) 
        if( diff_InvM > 10d0 .and. diff_PfM > 10d0 .and. grnd() < 1d-3 ) &  
            write(6,'(a,2I8,2E12.4)') 'diff is large, myrank, Nfastup, diff_InvM, diff_PfM:', & 
                                                      myrank, Nfastup, diff_InvM, diff_PfM
      else 
        Nfastup = min(100,idnint(dble(Nfastup)*1.1d0)) 
      end if
    end if
    call measure_time(t1,t(2))


    do iupdate = 1, Nupdate 

      call system_clock(t1)
      if( grnd() < 3.0d-4 ) then 
        ! global spin flip 
        x(:) = -x(:) 
        call calc_theta(x,thetapc,thetapr)              
        call x_to_Ele_Idx_Ele_Cfg(x,Tmp_Ele_Idx,Tmp_Ele_Cfg)
        call check_Ele_Idx_Ele_Cfg(x,Tmp_Ele_Idx,Tmp_Ele_Cfg)
        call calculate_MAll(.true.,Tmp_Ele_Idx,Tmp_Pf_M,Tmp_Inv_M)
        call sign_PP(Tmp_Ele_Idx,tmp_sgn_PP)
        call calc_amplitude_RBM_PP(thetapc,thetapr,Tmp_Pf_M,tmp_sgn_PP,psi_xp_wo_pf,psi_xp) 
        r1 = (abs(psi_xp(0))+shift)/(abs(psi_x(0))+shift)
        if( grnd() < r1**2 ) then 
          if( grnd() < 1d-6 ) write(6,'(a)') 'Global spin flip done'
          thetac(:)    = thetapc(:)
          thetar(:)    = thetapr(:)
          Ele_Idx(:)   = Tmp_Ele_Idx(:)
          Ele_Cfg(:)   = Tmp_Ele_Cfg(:)
          Pf_M(:)      = Tmp_Pf_M(:)
          Inv_M(:,:,:) = Tmp_Inv_M(:,:,:) 
          sgn_PP(:)    = tmp_sgn_PP(:)
          psi_x(:)     = psi_xp(:)
        else 
          x(:) = -x(:) 
        end if
      end if
      call measure_time(t1,t(3))


      call system_clock(t1)
      call make_spin_flip_candidate(x,i1,i2)
      call measure_time(t1,t(4))


      call system_clock(t1)
      !$omp parallel do default(shared) private(j) 
      do j = 1, Mc 
        thetapc(j) = thetac(j)
      end do ! j 
      !$omp end parallel do 
      thetapr(:) = thetar(:)
      call measure_time(t1,t(5))

 
      call system_clock(t1)
      call update_theta(x,i1,i2,thetapc,thetapr)
      call measure_time(t1,t(6))


      call system_clock(t1)
      i3 = Ele_Cfg(i1)
      i4 = Ele_Cfg(i2)
      if( Ele_Idx(i3) /= i1 ) stop 'Ele_Idx(i3) /= i1'
      if( Ele_Idx(i4) /= i2 ) stop 'Ele_Idx(i4) /= i2'
      Tmp_Ele_Idx(:) = Ele_Idx(:)
      Tmp_Ele_Idx(i3) = i2
      Tmp_Ele_Idx(i4) = i1
      call calculate_new_PfM_two(i3,i4,Tmp_Ele_Idx,Inv_M,Pf_M,Tmp_Pf_M)
      call measure_time(t1,t(7))


      call system_clock(t1)
      tmp_sgn_PP(:) = -sgn_PP(:)
      if( subB_or_not(i1) ) tmp_sgn_PP(:) = -tmp_sgn_PP(:) 
      if( subB_or_not(i2) ) tmp_sgn_PP(:) = -tmp_sgn_PP(:) 
      call calc_amplitude_RBM_PP(thetapc,thetapr,Tmp_Pf_M,tmp_sgn_PP,psi_xp_wo_pf,psi_xp)
      call measure_time(t1,t(8))


      call system_clock(t1)
      r1 = (abs(psi_xp(0))+shift)/(abs(psi_x(0))+shift)
      if( grnd() < r1**2 ) then 
        if( i0 > 0 ) naccept = naccept + 1
        nacc_tmp = nacc_tmp + 1
        Ntmp_accept = Ntmp_accept + 1 
        x(i1) = -x(i1)
        x(i2) = -x(i2)
        thetac(:) = thetapc(:)
        thetar(:) = thetapr(:)
        Ele_Cfg(i1) = i4
        Ele_Cfg(i2) = i3
        Ele_Idx(i3) = i2
        Ele_Idx(i4) = i1
        Pf_M(:) = Tmp_Pf_M(:)
        call calculate_new_InvM_two(i3,i4,Ele_Idx,Inv_M)
        sgn_PP(:)  = tmp_sgn_PP(:)
        psi_x(:)   = psi_xp(:)
      end if
      call measure_time(t1,t(9))


    end do ! iupdate
 

    call system_clock(t1)
    if( i0 > 0 ) then 
      x_sample(:,i0) = x(:)
      PfM_sample(:,i0) = Pf_M(:)
      psi_x_sample(:,i0) = psi_x(:)
    end if
    call measure_time(t1,t(10))

  end do ! i0

  if( mod(myrank,100) == 0 ) write(40000+myrank,'(I10,10F10.6)') Ntmp_accept, t(:)

  return 
end subroutine 
!
! x => xp (i1,i2 denotes the position for )  
!
subroutine make_spin_flip_candidate(x,i1,i2)
  use mod_main , only : L, N
  implicit none 
  integer , intent(in)  :: x(N) 
  integer , intent(out) :: i1, i2
  real(8) , external :: grnd 
  integer :: i
  integer :: ix, iy, dist_max
  logical :: lsuccess
  integer :: nran
  nran(i) = mod(int(dble(i)*grnd()),i) + 1  

  lsuccess = .false. 
  dist_max = nran(L-1) + 1 
  if( dist_max <= 1 ) stop 'dist_max <= 1' 
  do i = 1, 1000000   
    i1 = nran(N)

    ix = (i1-1) / L 
    iy = (i1-1) - ix * L  
    if( ix < 0  ) stop 'ix<0' 
    if( ix >= L ) stop 'ix>=L' 
    if( iy < 0  ) stop 'iy<0' 
    if( iy >= L ) stop 'iy>=L' 
    if( (i1-1) /= ix*L+iy ) stop 'i1 wrong'

    ix = ix + nran(dist_max) - 1  
    iy = iy + nran(dist_max) - 1 
    if( ix >= L ) ix = ix - L 
    if( iy >= L ) iy = iy - L 

    if( ix < 0  ) stop 'ix<0' 
    if( ix >= L ) stop 'ix>=L' 
    if( iy < 0  ) stop 'iy<0' 
    if( iy >= L ) stop 'iy>=L' 

    i2 = ix*L + iy + 1 
    if( i2 < 1 ) stop 'i2<1' 
    if( i2 > N ) stop 'i2>N' 

    if( x(i1) /= x(i2) ) then 
      lsuccess = .true.
      exit 
    end if
  end do ! i 

  if( .not. lsuccess ) stop 'make_spin_flip_candidate fails 1'
  if( i1 == i2 ) stop 'make_spin_flip_candidate fails 2'

  return 
end subroutine 
!
! make initial configuration x 
!
subroutine make_initial_x(x,psi_x,thetac,thetar,Ele_Idx,Ele_Cfg,sgn_PP,Pf_M,Inv_M) 
  use mod_main , only : N, N_QP_Full, N_MP_Trans 
  use mod_RBM  , only : Mc, Mr, Ncopy
  implicit none 
  integer    , intent(out) :: x(N)
  complex(8) , intent(out) :: psi_x(0:N_MP_Trans*Ncopy*2)
  complex(8) , intent(out) :: thetac(Mc)
  real(8)    , intent(out) :: thetar(Mr)
  integer    , intent(out) :: Ele_Idx(N)
  integer    , intent(out) :: Ele_Cfg(N)
  integer    , intent(out) :: sgn_PP(N_MP_Trans)
  real(8)    , intent(out) :: Pf_M(N_QP_Full)
  real(8)    , intent(out) :: Inv_M(N,N,N_QP_Full)

  complex(8) :: psi_x_wo_pf(N_QP_Full)
  logical :: lsuccess1, lsuccess2
  integer :: i, i1, i2, i3
  real(8) , external :: grnd
  integer :: nran
  nran(i) = mod(int(dble(i)*grnd()),i) + 1  

  lsuccess1 = .false. 

  do i1 = 1, 1000000 

    i3 = 0 
    x(:) = -1
    lsuccess2 = .false. 
    do i2 = 1, 1000000
      i = nran(N)
      if( x(i) == -1 ) then 
        i3 = i3 + 1 
        x(i) = 1 
      end if  
      if( i3 == N/2 ) then 
        lsuccess2 = .true.
        exit 
      end if 
    end do ! i2 
    if( .not. lsuccess2 ) stop 'make_initial_x fails 1'

    i2 = 0 
    do i = 1, N 
      i2 = i2 + x(i)
    end do ! i
    if( i2 /= 0 ) stop 'total magnetization /= 0'

    call calc_theta(x,thetac,thetar)
    call x_to_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg)
    call check_Ele_Idx_Ele_Cfg(x,Ele_Idx,Ele_Cfg)
    call calculate_MAll(.true.,Ele_Idx,Pf_M,Inv_M)
    call sign_PP(Ele_Idx,sgn_PP)
    call calc_amplitude_RBM_PP(thetac,thetar,Pf_M,sgn_PP,psi_x_wo_pf,psi_x)

    if( abs(psi_x(0)) > 1d-300 ) then 
      lsuccess1 = .true.
      exit
    end if

  end do ! i1 

  if( .not. lsuccess1 ) stop 'make_initial_x fails 2'

  return 
end subroutine 
