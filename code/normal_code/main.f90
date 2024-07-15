!
!  Solving 1D antiferromagnetic Heisenberg model by RBM 
!    using stochastic reconfiguration (SR) method for optimization  
!
!  Copyright (C) 2018 Yusuke Nomura
!
!
!  We assume translational symmetry in variational parameters 
!  Bias terms a_i b_j are set to be zero 
!  We only optimize interaction parameters W_ij
!
program RBM_solver 
  use mpi 
  use mod_main
  use mod_RBM
  use mod_PP
  use mod_sample
  use mod_pdposv
  implicit none 

  integer :: Nstep     ! total number of optimization steps
  real(8) :: delta_tau ! controls the amount of parameter update 

  integer :: f         ! index for independent neurons 
  integer :: i         ! index for visible units
  integer :: j         ! index for theta angles 
  integer :: k         ! index for variational parameters 
  integer :: iw        ! index for W interaction 
  integer :: mpidx     ! index for momentum projections
  integer :: orbidx    ! index for slater element
  integer :: iteration ! index for optimization iterations

  real(8) :: norm      
  real(8) :: wf_av_list(10)      
  real(8) :: Etot      ! total energy <H>
  real(8) :: Etot2     ! <H^2> 
  real(8) :: delW_max  ! maximum change in W at each optimization step

  real(8) , allocatable :: Smat(:,:)     ! metric  
  real(8) , allocatable :: gvec(:)       ! derivative of energy with respect to variational parameters 
  real(8) , allocatable :: Ovec(:)        
  real(8) , allocatable :: gamma_save(:,:) 
  logical , allocatable :: chk(:)

  real(8) , allocatable :: Etot_data(:), Etot2_data(:)
  real(8) :: Etot_er, Etot2_er
  integer :: naccept
  integer :: Nav, Nav_chk 
  real(8) :: accept_ratio
  real(8) :: max_change = 0.04d0
  real(8) :: rsmall_variance = 0.001d0
  integer :: Nsave = 1

  real(8) :: r1
  integer :: iseed
  integer :: isave, icopy, i1, i2, i3
  integer :: k1, k2, k_nonzero
  integer :: t1
  logical :: lsmall_variance
  real(8) :: t(6)
  real(8) , external :: grnd
  real(8) :: s_factor0, s_factor, s_factor_scale, s_cut 
  character(len=25) :: ch1
!---------- for MPI
  integer :: Nmpi, myrank, ierr
  real(8) , allocatable :: vec_tmp(:)
  real(8) , allocatable :: mat_tmp(:,:)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nmpi,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  ! 
  ! read inputs
  ! 
  open(unit=1,file='RBM.input',status='old')
  read(1,*)
  read(1,*) L, Nsample, spin_parity
  read(1,*)
  read(1,*) alphac, alphar, Nstep, Nav, delta_tau
  read(1,*)
  read(1,*) s_factor0, s_factor_scale, s_cut
  read(1,*)
  read(1,*) J1, J2
  close(1)
  N  = L * L                                     ! number of visible units 
  ! 
  ! read symmetry information
  ! 
  open(unit=1,file='sym_info.txt',status='old')
  read(1,*)
  read(1,*) N_MP_Trans, i1
  read(1,*)
  if( i1 /= Ncopy ) stop 'i1/=Ncopy'
  allocate( W_map(N,Ncopy,N_MP_Trans) ); W_map = -10
  allocate( para_mp_trans(N_MP_Trans) ); para_mp_trans = 0d0
  allocate( chk(N) ); chk = .false.
  do mpidx = 1, N_MP_Trans 
    read(1,*) i1, para_mp_trans(mpidx) 
    if( i1 /= mpidx ) stop 'i1/=mpidx' 
  end do ! mpidx
  read(1,*)
  do mpidx = 1, N_MP_Trans 
  do icopy = 1, Ncopy 

    do i = 1, N 
      read(1,*) i1, i2, i3, W_map(i,icopy,mpidx) 
      if( i1 /= mpidx ) stop 'i1/=mpidx (sym_info.txt)'
      if( i2 /= icopy ) stop 'i2/=icopy (sym_info.txt)'
      if( i3 /= i     ) stop 'i3/=i     (sym_info.txt)'
      if( W_map(i,icopy,mpidx) < 1 .or. W_map(i,icopy,mpidx) > N ) stop 'W_map wrong 1 (sym_info.txt)'
    end do ! i

    chk = .false.
    do i = 1, N 
      chk(W_map(i,icopy,mpidx)) = .true.
    end do ! i 
    do i = 1, N 
      if( .not. chk(i) ) stop 'W_map wrong 2 (sym_info.txt)'
    end do ! i 

  end do ! icopy
  end do ! mpidx
  deallocate(chk)
  close(1) 

  open(unit=1,file='orbitalidx.def',status='old')
  read(1,*)
  read(1,*) ch1, N_Slater 
  close(1)
  ! 
  ! constants
  ! 
  N_SP_Gauss_Leg = 2 
  Mc = alphac * N_MP_Trans * Ncopy * 2           ! number of independent theta angles = number of hidden unit 
  Mr = alphar * N_MP_Trans * Ncopy * 2           ! number of independent theta angles = number of hidden unit 
  Nv = alphac*2*(N+1) + alphar*(N+1) + N_Slater  ! number of variational parameters
  Nv_RBM = Nv - N_Slater                         ! number of RBM variational parameters
  Nsample = Nsample / Nmpi
  Nupdate = N !* 2
  N_QP_Full = N_SP_Gauss_Leg * N_MP_Trans
  Npair = N/2
  ! 
  ! write calculation condition
  ! 
  if( myrank == 0 ) then 
do mpidx = 1, N_MP_Trans
  write(6,*) para_mp_trans(mpidx)
end do ! mpidx
    write(6,*) 'Nmpi', Nmpi
    write(6,*) 'Ncopy', Ncopy
    write(6,*)
    write(6,'(a)')         '=================== calculation condition ===================='
    write(6,'(a,2I8)')     '  L (system size), N                 =', L, N
    write(6,'(a,I8)')      '  number of momentum projections     =', N_MP_Trans
    write(6,'(a,I8)')      '  number of spin projections         =', N_SP_Gauss_Leg
    write(6,'(a,I8)')      '  toptal number of MC sampling       =', Nsample
    write(6,'(a,I8)')      '  spin parity (+1/-1 for S even/odd) =', spin_parity
    write(6,'(a,I8)')      '  # of updates between measurement   =', Nupdate
    write(6,'(a,I8)')      '  alphac (# of independent neurons)  =', alphac
    write(6,'(a,I8)')      '  alphar (# of independent neurons)  =', alphar
    write(6,'(a,I8)')      '  total number of optimization steps =', Nstep
    write(6,'(a,F14.10)')  '  delta tau                          =', delta_tau 
    write(6,'(a,2F10.5)')  '  J1, J2                             =', J1, J2
    write(6,'(a)')         '=============================================================='
    write(6,*)
  end if 
  ! 
  ! allocate arrays 
  !  
  allocate( para_spin_proj(N_SP_Gauss_Leg) ); para_spin_proj = 0d0

  allocate( Wcirr(N,alphac)       ); Wcirr = 0d0
  allocate( Wrirr(N,alphar)       ); Wrirr = 0d0
  allocate( bcirr(alphac)         ); bcirr = 0d0
  allocate( brirr(alphar)         ); brirr = 0d0
  allocate( gvec(Nv)              ); gvec = 0d0
  allocate( Ovec(Nv)              ); Ovec = 0d0
  allocate( Smat(Nv,Nv)           ); Smat = 0d0

  allocate( vec_tmp(Nv)        ); vec_tmp = 0d0
  allocate( mat_tmp(Nv,Nv)     ); mat_tmp = 0d0

  allocate( x_sample(N,Nsample)                  ); x_sample = 0 
!! DELETE delete Delete
allocate( psi_x_sample(0:N_MP_Trans*Ncopy*2,Nsample) ); psi_x_sample = 0d0 
!! DELETE delete Delete
  allocate( PfM_sample(N_QP_Full,Nsample) ); PfM_sample = 0d0

  allocate( Wcirr_av(N,alphac) ); Wcirr_av = 0d0
  allocate( Wrirr_av(N,alphar) ); Wrirr_av = 0d0
  allocate( bcirr_av(alphac)   ); bcirr_av = 0d0
  allocate( brirr_av(alphar)   ); brirr_av = 0d0

  allocate( Orbital_Idx(N,N) ); Orbital_Idx = 0  
  allocate( Orbital_Sgn(N,N) ); Orbital_Sgn = 1  
  allocate( MP_Trans_Idx(N,N_MP_Trans) ); MP_Trans_Idx = 0 
  allocate( MP_Trans_Sgn(N,N_MP_Trans) ); MP_Trans_Sgn = 1 
  allocate( Slater(N_Slater) ); Slater = 0d0
  allocate( Slater_av(N_Slater) ); Slater_av = 0d0

  allocate( Slater_Elm(2*N,2*N,N_QP_Full) ); Slater_Elm = 0d0
  allocate( subB_or_not(N)                ); subB_or_not = .false.

  if( myrank == 0 ) then  
    allocate( gamma_save(Nv,Nsave) ); gamma_save = 0d0
  end if
  ! 
  ! preparation for pdposv
  ! 
  allocate( numar(0:Nmpi-1) ); numar = 0
  allocate( numac(0:Nmpi-1) ); numac = 0
  allocate( numbr(0:Nmpi-1) ); numbr = 0
  allocate( numbc(0:Nmpi-1) ); numbc = 0

  include "pre_define_pdposv.f90"

  allocate( amat(maxnumar,Nv)    ); amat   = 0d0
  allocate( bmat(maxnumbr,N_rhs) ); bmat   = 0d0
  allocate( s_save(maxnumar,Nv)  ); s_save = 0d0
  ! 
  ! creating different seed 
  !
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call system_clock(count=j)
  iseed = modulo(j,10000000) + 1
  if( myrank == 0 ) write(6,'(a,I12)') 'iseed (initial) =', iseed
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  iseed = iseed + myrank
  if( mod(myrank,100) == 0 ) write(6,'(a,2I12)') 'myrank, iseed =', myrank, iseed
  call sgrnd(iseed)
  do i = 0, Nsample*myrank
    r1 = grnd()
  end do ! i 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! 
  ! define global parameter for spin quantum projection for pair-product part
  ! 
  call InitSPWeight(myrank)
  ! 
  ! define sublattice 
  ! 
  call define_sublattice()
  ! 
  ! initialize variational parameters   
  !
  call initialize_parameters(myrank)
  call write_initial_parameters(myrank)
 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if( alphac > 0 ) then  
    call MPI_BCAST(bcirr,alphac,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Wcirr,alphac*N,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
  end if
  if( alphar > 0 ) then  
    call MPI_BCAST(brirr,alphar,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Wrirr,alphar*N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
  call MPI_BCAST(Slater,N_Slater,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  if( myrank == 0 ) then 

    k1 = 0    
    do f = 1, alphac     
      k1 = k1 + 1
      gamma_save(k1,Nsave) = dble (bcirr(f))
      k1 = k1 + 1
      gamma_save(k1,Nsave) = dimag(bcirr(f))
    end do ! f  
    do f  = 1, alphac     
    do iw = 1, N 
      k1 = k1 + 1
      gamma_save(k1,Nsave) = dble (Wcirr(iw,f))
      k1 = k1 + 1
      gamma_save(k1,Nsave) = dimag(Wcirr(iw,f))
    end do ! iw 
    end do ! f  

    do f = 1, alphar     
      k1 = k1 + 1
      gamma_save(k1,Nsave) = brirr(f)
    end do ! f  
    do f  = 1, alphar
    do iw = 1, N
      k1 = k1 + 1
      gamma_save(k1,Nsave) = Wrirr(iw,f)
    end do ! iw 
    end do ! f  

    do orbidx = 1, N_Slater 
      k1 = k1 + 1 
      gamma_save(k1,Nsave) = Slater(orbidx)
    end do ! orbidx 
    if( k1 /= Nv ) stop 'k1/=Nv'

  end if
  ! 
  ! initial wave function
  ! 
  call update_slater_elements()
  ! 
  ! start optimization  
  ! 
  lsmall_variance = .false.
  shift = 0d0
  shift_av = 0d0
  wf_av_list = 0d0
  Nav_chk = 0
  Nfastup = min(100,N*2) 
  if( myrank == 0 ) open(unit=10,file='Energy_vs_Iteration.txt',status='unknown')
  do iteration = 1, Nstep 
    if( myrank == 0 ) write(6,'(a,I6)')     'Iteration :', iteration
    t(:) = 0d0
    !
    ! making samples 
    !
    call system_clock(t1)
    call update_slater_elements()
    call measure_time(t1,t(1))


    call system_clock(t1)
    if( iteration == 1 ) then 
      call make_sample(myrank,.true.,naccept)
    else 
      call make_sample(myrank,.false.,naccept)
    end if
    call measure_time(t1,t(2))
    ! 
    ! calculate physical quantities 
    ! 
    call system_clock(t1)
    call calc_phys_quantity_opt(myrank,Etot,Etot2,gvec,Ovec,Smat,norm,wf_av)
    call measure_time(t1,t(3))
    ! 
    ! Dividing by Nsample*Nmpi, we get  
    !  gvec(k) = 2 <H*O_k> 
    !  Ovec(k) = <O_k>  
    !  Etot = <H>  
    !  Smat(k1,k2) = <O_k1*O_k2>
    !  
    call system_clock(t1)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
    call measure_time(t1,t(4))

 
    call system_clock(t1)
    call MPI_REDUCE(naccept,i1,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if( myrank == 0 ) accept_ratio = dble(i1)/dble(Nsample*Nmpi)/dble(Nupdate)

    call MPI_ALLREDUCE(wf_av,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    wf_av = r1 / dble(Nsample*Nmpi) 

    call MPI_ALLREDUCE(norm,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    norm = r1 

    call MPI_ALLREDUCE(Etot,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    Etot = r1 / norm!dble(Nsample*Nmpi)

    call MPI_ALLREDUCE(Etot2,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    Etot2 = r1 / norm!dble(Nsample*Nmpi)

    call MPI_ALLREDUCE(gvec,vec_tmp,Nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    gvec(:) = vec_tmp(:) / norm!dble(Nsample*Nmpi)
 
    call MPI_ALLREDUCE(Ovec,vec_tmp,Nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    Ovec(:) = vec_tmp(:) / norm!dble(Nsample*Nmpi)

    call MPI_ALLREDUCE(Smat,mat_tmp,Nv*Nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call measure_time(t1,t(5))
    ! 
    ! update variational parameters 
    ! 
    call system_clock(t1)
    include "SR_pdposv.f90"
    call measure_time(t1,t(6))


    if( myrank == 0 ) then 
      ! 
      ! write information 
      ! 
      write(6,'(a,F20.10)')     '   Total Energy         :', Etot
      write(6,'(a,F20.10)')     '   Energy Variance      :', (Etot2-Etot*Etot)/(Etot*Etot)
      write(6,'(a,F20.10)')     '   Acceptance Ratio     :', accept_ratio
      write(6,'(a,2F20.10)')    '   Stabilization Factor :', s_factor, s_cut
      write(6,'(a,2E20.5)')     '   wf_av, shift         :', wf_av, shift
      write(6,'(a,2E20.5)')     '   max fij parameter    :', maxval(abs(Slater(:)))
      write(6,'(a,F20.10)')     '   norm                 :', norm/dble(Nsample*Nmpi)
      write(6,'(a,2F20.10,I8)') '   Maximum change in W  :', min(delW_max,max_change), delW_max/max_change, kmax 
      write(6,'(a,6F16.8)')     '   Computational time   :', t(:)
      if( mod(iteration,200) == 0 ) then  
        write(6,'(a,I15)')           '   lwork_pfa            :', lwork_pfa
        write(6,'(a,2E15.5)')        '   psi_x (snapshot)     :', psi_x_sample(0,1)
        write(6,'(a,100(2F9.5,2x))') '   bias term            :', bcirr(:)
        write(6,'(a,100(2F9.5,2x))') '   bias term            :', brirr(:)
        write(6,'(a)') '   W parameters:'
        do f = 1, alphac 
          write(6,'(5x,200F9.5)')  (Wcirr(iw,f), iw = 1, N)
        end do ! f
        do f = 1, alphar 
          write(6,'(5x,200F9.5)')  (Wrirr(iw,f), iw = 1, N)
        end do ! f
      end if

      write(10,*) iteration, Etot, (Etot2-Etot*Etot)/(Etot*Etot)
      ! 
      ! take average 
      ! 
      if( iteration > Nstep-Nav ) then 
        Nav_chk = Nav_chk + 1
        shift_av = shift_av + shift
        do f = 1, alphac 
          Wcirr_av(:,f) = Wcirr_av(:,f) + Wcirr(:,f)
          bcirr_av(f)   = bcirr_av(f)   + bcirr(f)
        end do ! f 
        do f = 1, alphar 
          Wrirr_av(:,f) = Wrirr_av(:,f) + Wrirr(:,f)
          brirr_av(f)   = brirr_av(f)   + brirr(f)
        end do ! f 
        Slater_av(:) = Slater_av(:) + Slater(:)
      end if

      !if( accept_ratio < 0.5d0 ) Nupdate = N * 2
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if( alphac > 0 ) then 
      call MPI_BCAST(bcirr,alphac,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Wcirr,alphac*N,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    end if
    if( alphar > 0 ) then 
      call MPI_BCAST(Wrirr,alphar*N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(brirr,alphar,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    end if
    call MPI_BCAST(Slater,N_Slater,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nupdate,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    s_factor0 = s_factor0 * s_factor_scale

    do i1 = 2, 10
      wf_av_list(i1-1) = wf_av_list(i1)
    end do ! i1
    wf_av_list(10) = wf_av

    wf_av = 0d0 
    do i1 = 1, 10 
      wf_av = wf_av + wf_av_list(i1)
    end do ! i1
    wf_av = wf_av / 10d0
    shift = wf_av * 0.001d0 / dble(N) / dble(N)
    !shift = 0d0
    call MPI_BCAST(shift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  end do ! iteration

  if( myrank == 0 ) then 
    if( Nav_chk /= Nav ) stop 'Nav_chk /= Nav'
    shift_av      = shift_av      / dble(Nav)
    Wcirr_av(:,:) = Wcirr_av(:,:) / dble(Nav)
    Wrirr_av(:,:) = Wrirr_av(:,:) / dble(Nav)
    bcirr_av(:)   = bcirr_av(:)   / dble(Nav)
    brirr_av(:)   = brirr_av(:)   / dble(Nav)
    Slater_av(:)  = Slater_av(:)  / dble(Nav)
    close(10)
  end if
  ! 
  ! distribute avaraged parameters 
  !
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if( alphac > 0 ) then 
    call MPI_BCAST(bcirr_av,alphac,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Wcirr_av,alphac*N,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
  end if
  if( alphar > 0 ) then 
    call MPI_BCAST(brirr_av,alphar,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Wrirr_av,alphar*N,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
  call MPI_BCAST(Slater_av,N_Slater,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(shift_av,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  Wcirr(:,:) = Wcirr_av(:,:)
  Wrirr(:,:) = Wrirr_av(:,:)
  bcirr(:)   = bcirr_av(:)
  brirr(:)   = brirr_av(:)
  Slater(:)  = Slater_av(:)
  shift = shift_av
  ! 
  ! Final Energy  
  ! 
  call write_optimized_parameters(myrank)
  if( myrank == 0 ) then 
    open(unit=10,file='energy_after_opt.txt',status='unknown')
    allocate( Etot_data (Nav) ); Etot_data  = 0d0
    allocate( Etot2_data(Nav) ); Etot2_data = 0d0
  end if
  ! 
  ! Warm up
  ! 
  call make_sample(myrank,.true.,naccept)
  do iteration = 1, 2
    call make_sample(myrank,.false.,naccept)
  end do ! iteration 
  ! 
  ! Measure Energy 
  ! 
  do iteration = 1, Nav 
    t(:) = 0d0
    !
    ! making samples 
    !
    call system_clock(t1)
    call make_sample(myrank,.false.,naccept)
    call measure_time(t1,t(2))
    ! 
    ! calculate physical quantities 
    ! 
    call system_clock(t1)
    call calc_phys_quantity_after_opt(Etot,Etot2,norm,wf_av)
    call measure_time(t1,t(3))
    ! 
    ! gather data
    ! 
    call system_clock(t1)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
    call measure_time(t1,t(4))

 
    call system_clock(t1)
    call MPI_REDUCE(naccept,i1,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if( myrank == 0 ) accept_ratio = dble(i1)/dble(Nsample*Nmpi)/dble(Nupdate)

    call MPI_REDUCE(wf_av,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if( myrank == 0 ) wf_av = r1 / dble(Nsample*Nmpi) 

    call MPI_REDUCE(norm,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if( myrank == 0 ) norm = r1 

    call MPI_REDUCE(Etot,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if( myrank == 0 ) Etot = r1 / norm!dble(Nsample*Nmpi)

    call MPI_REDUCE(Etot2,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if( myrank == 0 ) Etot2 = r1 / norm!dble(Nsample*Nmpi)
    call measure_time(t1,t(5))
    ! 
    ! write info
    ! 
    if( myrank == 0 ) then 
      write(6,'(a,I6)')         'Iteration (Energy Measuement) :', iteration
      write(6,'(a,F20.10)')     '   total energy         :', Etot
      write(6,'(a,F20.10)')     '   energy variance      :', (Etot2-Etot*Etot)/(Etot*Etot)
      write(6,'(a,F20.10)')     '   acceptance ratio     :', accept_ratio
      write(6,'(a,2E20.5)')     '   wf_av, shift         :', wf_av, shift
      write(6,'(a,F20.10)')     '   norm                 :', norm/dble(Nsample*Nmpi)
      write(6,'(a,6F16.8)')     '   computational time   :', t(:)
      Etot_data (iteration) = Etot
      Etot2_data(iteration) = Etot2-Etot*Etot
      write(10,*) iteration, Etot, (Etot2-Etot*Etot)/(Etot*Etot)
    end if 

  end do ! iteration
 
  if( myrank == 0 ) then 
    close(10)
    call calc_av_er(Nav,Etot_data,Etot,Etot_er) 
    call calc_av_er(Nav,Etot2_data,Etot2,Etot2_er) 
    Etot2    = Etot2    / (Etot*Etot)
    Etot2_er = Etot2_er / (Etot*Etot)
    open(unit=10,file='Average_Energy.txt',status='unknown')
    write(10,'(3F25.15)') Etot,  Etot_er /dsqrt(dble(Nav-1)), Etot_er
    write(10,'(3E25.15)') Etot2, Etot2_er/dsqrt(dble(Nav-1)), Etot2_er
    write(10,'(4F20.12)') Etot,  Etot_er /dsqrt(dble(Nav-1)), & 
                          Etot2, Etot2_er/dsqrt(dble(Nav-1))
    close(10)
  end if
  ! 
  ! write optimized variational parameters
  ! 
  !call write_optimized_parameters(myrank)

  deallocate(W_map,para_mp_trans,para_spin_proj)
  deallocate(SP_GL_Cos,SP_GL_Sin,SP_GL_CosCos,SP_GL_CosSin,SP_GL_SinSin)
  deallocate(Wcirr,Wrirr,bcirr,brirr)
  deallocate(gvec,Ovec,Smat)
  deallocate(vec_tmp,mat_tmp)
  deallocate(Orbital_Idx,Orbital_Sgn,MP_Trans_Idx,MP_Trans_Sgn)
  deallocate(Slater,Slater_Elm,subB_or_not)
  deallocate(Wcirr_av,Wrirr_av,bcirr_av,brirr_av,Slater_av)
  deallocate(x_sample,PfM_sample,psi_x_sample)
  deallocate(numar,numac,numbr,numbc,amat,bmat,s_save)
  if( myrank == 0 ) deallocate(Etot_data,Etot2_data)
  if( myrank == 0 ) deallocate(gamma_save)
  call MPI_FINALIZE(ierr)
end program 
