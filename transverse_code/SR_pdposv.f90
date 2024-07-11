  ! 
  ! gvec(k) = 2 <H*O_k> - 2 <H> <O_k>  
  !  
  do k = 1, Nv
    gvec(k) = gvec(k) - 2d0 * Etot * Ovec(k)
  end do ! k
  ! 
  ! Smat(k1,k2) = <O_k1*O_k2> - <O_k1><O_k2>
  !  
  !$omp parallel do default(shared) private(k1,k2) 
  do k2 = 1, Nv
  do k1 = 1, Nv
    Smat(k1,k2) = mat_tmp(k1,k2)/norm - Ovec(k1)*Ovec(k2)
  end do ! k1
  end do ! k2
  !$omp end parallel do 

  Sdiag_av_RBM = 0d0
  Sdiag_av_PP  = 0d0
  i1 = 0; i2 = 0 
  do k1 = 1, Nv
    if( Smat(k1,k1) > 0d0 ) then 
      if( k1 <= Nv_RBM ) then 
        i1 = i1 + 1 
        Sdiag_av_RBM = Sdiag_av_RBM + Smat(k1,k1) 
      else 
        i2 = i2 + 1 
        Sdiag_av_PP  = Sdiag_av_PP  + Smat(k1,k1) 
      end if
    end if  
  end do ! k1
  if( i1 /= Nv_RBM ) stop 'i1/=Nv_RBM'
  if( i2 /= N_Slater - N_Slater/N ) stop 'i2 /= N_Slater - N_Slater/N'
  if( Nv_RBM > 0 ) then 
    Sdiag_av_RBM = Sdiag_av_RBM / dble(i1)
  else 
    Sdiag_av_RBM = 1d0
  end if
  Sdiag_av_PP = Sdiag_av_PP  / dble(i2)
  Smat_scaling_factor = dsqrt(Sdiag_av_RBM/Sdiag_av_PP)

  if( myrank == 0 .and. mod(Iteration,200) == 0 ) then 
    write(6,'(a,3E15.5)') '   Sdiag_av_RBM, Sdiag_av_PP, Smat_scaling_factor', &  
                              Sdiag_av_RBM, Sdiag_av_PP, Smat_scaling_factor
    rewind(901)
    do k1 = 1, Nv
      write(901,*) Smat(k1,k1), gvec(k1)
    end do ! k1
    close(901)
  end if
  ! 
  ! scaling Smat and gvec 
  ! 
  !$omp parallel default(shared) private(k1,k2)
  !$omp do 
  do k1 = Nv_RBM+1, Nv
    gvec(k1) = gvec(k1) * Smat_scaling_factor 
  end do ! k1
  !$omp end do 
  !$omp do 
  do k2 = Nv_RBM+1, Nv
    do k1 = 1, Nv 
      Smat(k1,k2) = Smat(k1,k2) * Smat_scaling_factor
    end do ! k1 
  end do ! k2  
  !$omp end do 
  !$omp do 
  do k1 = Nv_RBM+1, Nv
    do k2 = 1, Nv 
      Smat(k1,k2) = Smat(k1,k2) * Smat_scaling_factor
    end do ! k2
  end do ! k1
  !$omp end do 
  !$omp end parallel 
  ! 
  ! check Smat
  ! 
  do k1 = 1, Nv
    if( Smat(k1,k1) < 0d0 ) then 
      stop 'diagonal element of S matrix is not positive'
    else if( Smat(k1,k1) == 0d0 ) then 
      i2 = 0 
      do i1 = 1, N 
        if( k1 == Orbital_Idx(i1,i1)+Nv-N_Slater ) i2 = i2 + 1  
      end do ! i1 
      if( i2 == 0 ) stop 'diagonal element of S matrix is zero for other than f_ii parameters'
      if( abs(Ovec(k1)) > 1d-300 ) stop 'Ovec related with f_ii parameters is wrong'
      if( abs(gvec(k1)) > 1d-300 ) stop 'gvec related with f_ii parameters is wrong'
      do k2 = 1, Nv
        if( abs(Smat(k2,k1)) > 1d-300 ) stop 'Smat related with f_ii parameters is wrong'  
        if( abs(Smat(k1,k2)) > 1d-300 ) stop 'Smat related with f_ii parameters is wrong'  
      end do ! k2 
    else  
      i2 = 0 
      do i1 = 1, N 
        if( k1 == Orbital_Idx(i1,i1)+Nv-N_Slater ) i2 = i2 + 1  
      end do ! i1 
      if( i2 /= 0 ) stop 'diagonal element of S matrix is positive for f_ii parameters'
    end if
  end do ! k1 
  ! 
  ! calculate Sdiag_av and nSmat 
  ! 
  Sdiag_av = 0d0
  Sdiag_av_RBM = 0d0
  Sdiag_av_PP  = 0d0
  Smax = 0d0; Smin = exp(50d0)
  kmax = -100; kmin = -100
  i1 = 0; i2 = 0 
  do k1 = 1, Nv 
    if( Smat(k1,k1) < 0d0 ) stop 'e1 (SR)'
    if( Smat(k1,k1) > 0d0 ) then
      if( Smat(k1,k1) > Smax ) then 
        Smax = Smat(k1,k1) 
        kmax = k1 
      end if
      if( Smat(k1,k1) < Smin ) then 
        Smin = Smat(k1,k1) 
        kmin = k1 
      end if
      Sdiag_av = Sdiag_av + Smat(k1,k1)
      if( k1 <= Nv_RBM ) then 
        i1 = i1 + 1 
        Sdiag_av_RBM = Sdiag_av_RBM + Smat(k1,k1) 
      else 
        i2 = i2 + 1 
        Sdiag_av_PP  = Sdiag_av_PP  + Smat(k1,k1) 
      end if
    end if
  end do ! k1
  if( i1 /= Nv_RBM ) stop 'i1/=Nv_RBM'
  if( i2 /= N_Slater - N_Slater/N ) stop 'i2 /= N_Slater - N_Slater/N'
  Sdiag_av = Sdiag_av / dble(i1+i2)
  if( Nv_RBM > 0 ) then 
    Sdiag_av_RBM = Sdiag_av_RBM / dble(i1)
  else 
    Sdiag_av_RBM = 1d0
  end if
  Sdiag_av_PP = Sdiag_av_PP / dble(i2)
  if( abs(Sdiag_av_RBM/Sdiag_av_PP-1d0) > 1d-5 ) stop 'Sdiag_av_RBM/=Sdiag_av_PP' 


  rcut = Sdiag_av * s_cut 
  nSmat = 0 
  do k1 = 1, Nv
    if( Smat(k1,k1) > rcut ) nSmat = nSmat + 1  
  end do ! i 

  if( myrank == 0 .and. mod(Iteration,200) == 0 ) then 
    write(6,'(a,2E15.5)') '   Sdiag_av_RBM, Sdiag_av_PP', &  
                              Sdiag_av_RBM, Sdiag_av_PP
    rewind(902)
    do k1 = 1, Nv
      write(902,*) Smat(k1,k1), gvec(k1)
    end do ! k1
    close(902)
  end if
  ! 
  ! TODO cut treatment  !!!!!
  ! 
  if( myrank == 0 ) write(6,'(a,3I6,4E11.4)') '   nSmat,kmax,kmin,Smax,Smin,Sdiag_av,rcut',& 
                                                  nSmat,kmax,kmin,Smax,Smin,Sdiag_av,rcut

  ! 
  ! stabilization factor
  ! 
  s_factor = max(1d-4,s_factor0)
  do k1 = 1, Nv
    Smat(k1,k1) = Smat(k1,k1) + Sdiag_av*s_factor
  end do ! i     

  if( grnd() < 0.005d0 .and. myrank == 0 ) then 
    allocate( s_save_gomi(Nv,Nv) ); s_save_gomi(:,:) = Smat(:,:)
  end if
  ! 
  ! preparation of amat and bmat
  !
  offset = 0
  do i = 0, myrank-1 
    offset = offset + numar(i)
  end do ! i 

  do k2 = 1, Nv 
  do k1 = 1, numar(myrank)
    amat(k1,k2) = Smat(offset+k1,k2)
  end do ! k1 
  end do ! k2 
  s_save(:,:) = amat(:,:)

  do k1 = 1, numar(myrank)
    bmat(k1,1) = gvec(offset+k1)
  end do ! k1 

  if( myrank == 0 ) then 
    allocate( gchk(Nv) ); gchk = gvec
  end if
  ! 
  ! solve linear equation using pdposv 
  ! 
  call sl_init(ictxt,nprow,npcol)
  call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)
  if( myrow == -1 ) stop 'myrow == -1 (pdposv)'
  if( mycol == -1 ) stop 'mycol == -1 (pdposv)'

  mynumar = numroc(Nv,Nb,myrow,rsrc,nprow);
  mynumac = numroc(Nv,Nb,mycol,csrc,npcol);
  if( mynumar /= numar(myrank) ) stop 'mynumar /= numar(myrank) (pdposv)'
  if( mynumac /= numac(myrank) ) stop 'mynumac /= numac(myrank) (pdposv)'

  mynumbr = mynumar
  mynumbc = numroc(N_rhs,Nb_rhs,mycol,csrc,npcol);
  if( mynumbr /= numbr(myrank) ) stop 'mynumbr /= numbr(myrank) (pdposv)'
  if( mynumbc /= numbc(myrank) ) stop 'mynumbc /= numbc(myrank) (pdposv)'

  mylda = maxnumar
  info = 0
  call descinit(desca,Nv,Nv,Nb,Nb,rsrc,csrc,ictxt,mylda,info)
  if(info /= 0) then
    write(6,*) 'info (descinit for amat):', info
    stop
  end if 

  myldb = maxnumbr
  info = 0
  call descinit(descb,Nv,N_rhs,Nb,Nb_rhs,rsrc,csrc,ictxt,myldb,info)
  if(info /= 0) then
    write(6,*) 'info (descinit for bmat):', info
    stop
  end if 

  info = 0
  call pdposv('U',Nv,N_rhs,amat,ia,ja,desca,bmat,ib,jb,descb,info)
  if(info /= 0) then
    write(6,*) 'info (pdposv):', info
    stop
  end if 

  call blacs_gridexit(ictxt)
  ! 
  ! gather solution
  ! 
  gvec(:) = 0d0
  do k1 = 1, numbr(myrank)
    gvec(offset+k1) = bmat(k1,1)
  end do ! k1

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  call MPI_ALLREDUCE(gvec,vec_tmp,Nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  gvec(:) = vec_tmp(:) 

  lambda = 0d0
  lambda_RBM_RBM = 0d0
  lambda_PP_PP   = 0d0
  lambda_RBM_PP  = 0d0
  do k2 = 1, Nv 
  do k1 = 1, numar(myrank)
    r1 = gvec(offset+k1)*s_save(k1,k2)*gvec(k2)
    lambda = lambda + r1 
    if( k2 <= Nv_RBM .and. offset+k1 <= Nv_RBM ) then 
      lambda_RBM_RBM = lambda_RBM_RBM + r1 
    else if ( k2 > Nv_RBM .and. offset+k1 > Nv_RBM ) then  
      lambda_PP_PP = lambda_PP_PP + r1 
    else 
      lambda_RBM_PP = lambda_RBM_PP + r1 
    end if
  end do ! k1
  end do ! k2

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  call MPI_REDUCE(lambda,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if( myrank == 0 ) lambda = sqrt(r1)
  call MPI_REDUCE(lambda_RBM_RBM,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if( myrank == 0 ) lambda_RBM_RBM = r1
  call MPI_REDUCE(lambda_PP_PP,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if( myrank == 0 ) lambda_PP_PP = r1
  call MPI_REDUCE(lambda_RBM_PP,r1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if( myrank == 0 ) lambda_RBM_PP = r1
  !
  ! check
  !
  if( myrank == 0 .and. grnd() < 0.005d0 ) then
    call solve_linear_equation(Nv,1,Smat,gchk) 
    do k1 = 1, Nv 
      if( abs(gchk(k1)) > 1d-12 ) then 
        if( abs(gchk(k1)-gvec(k1))/abs(gchk(k1)) > 1d-5 .and. abs(gchk(k1)-gvec(k1)) > 1d-6 ) stop 'gchk/=gvec'
      end if
    end do ! k1
    write(6,*) 'SR check done'
  end if
  if( myrank == 0 ) deallocate(gchk)

  if( myrank == 0 ) then 

    if( allocated(s_save_gomi) ) then 
      lam_chk = 0d0
      do k2 = 1, Nv
      do k1 = 1, Nv
        lam_chk = lam_chk + gvec(k1)*s_save_gomi(k1,k2)*gvec(k2)  
      end do ! k1 
      end do ! k2
      lam_chk = sqrt(lam_chk) 
      if( abs(lam_chk-lambda) > 1d-5 ) stop 'lambda wrong'
      write(6,*) 'lambda check done'
    end if
    ! 
    ! scale delta_tau 
    ! 
    if( abs(lambda**2-lambda_RBM_RBM-lambda_PP_PP-lambda_RBM_PP) > 1d-5 ) stop 'lambda division is wrong' 
    lambda = lambda * 5d0 / dble(L)
    write(6,'(a,4E15.6)')  '   lambda(tot,RBM,PP,cross-RBM-PP)', lambda, lambda_RBM_RBM, lambda_PP_PP, lambda_RBM_PP  
    if( lambda > 1d0 ) gvec(:) = gvec(:) / lambda 

    if( grnd() < 0.1d0 .and. lambda > 1d0 .and. allocated(s_save_gomi) ) then 
      write(6,'(a)') '  norm of gSg is checked'
      lambda = 0d0
      do k2 = 1, Nv
      do k1 = 1, Nv
        lambda = lambda + gvec(k1)*s_save_gomi(k1,k2)*gvec(k2)  
      end do ! k1 
      end do ! k2
      lambda = sqrt(lambda) 
      if( abs(lambda-(dble(L)/5d0)) > 1.0d-6 ) stop 'error in g'
    end if 

    do k1 = Nv_RBM+1, Nv
      gvec(k1) = gvec(k1) * Smat_scaling_factor 
    end do ! k1

    delW_max = 0d0 
    kmax = -100
    k_nonzero = 0 
    do k = 1, Nv
      if( abs(delta_tau*gvec(k)) > delW_max ) then 
        delW_max = abs(delta_tau*gvec(k))
        kmax = k 
      end if
      if( abs(gvec(k)) > 1d-12 ) then 
        k_nonzero = k_nonzero + 1
      end if
    end do ! k 
    write(6,'(a,I10)')  '   k nonzero            :', k_nonzero

    do isave = 2, Nsave  
      gamma_save(:,isave-1) = gamma_save(:,isave)
    end do ! isave 

    do k = 1, Nv
      if( delW_max > max_change ) then  
        ! avoid large change in variational parameters to stabilize optimization
        gamma_save(k,Nsave) = gamma_save(k,Nsave) - delta_tau*gvec(k) * (max_change/delW_max) 
      else 
        gamma_save(k,Nsave) = gamma_save(k,Nsave) - delta_tau*gvec(k) 
      end if
    end do ! k

    offset = 2*alphac*(N+1) + alphar*(N+1)
    k1 = 0; r1 = 0d0 
    do k = offset+1, Nv
      k1 = k1 + 1 
      r1 = r1 + abs(gamma_save(k,Nsave)) 
    end do ! k 
    if( k1 /= N_Slater ) stop 'k1/=N_Slater'
    r1 = r1 / dble(N_Slater)
    do k = offset+1, Nv
      gamma_save(k,Nsave) = gamma_save(k,Nsave) * slater_abs_av / r1
    end do ! k 


    if( (Etot2-Etot*Etot)/(Etot*Etot) < rsmall_variance ) lsmall_variance = .true.
    i1 = max(Nsave*5,Nav)
    if( lsmall_variance .and. Iteration > Nsave .and. mod(Iteration,i1) == 1 ) then 
      write(6,'(a,7x,a)')       '   parameter averaging  :', 'done'
      do k = 1, Nv
        r1 = 0d0
        do isave = 1, Nsave
          r1 = r1 + gamma_save(k,isave)
        end do ! isave
        r1 = r1 / dble(Nsave)
        gamma_save(k,Nsave) = r1
      end do ! k 
    else 
      write(6,'(a,7x,a)')       '   parameter averaging  :', 'not performed'
    end if
 
    k1 = 0 
    do f = 1, alphac     
      k1 = k1 + 2
      bcirr(f) = dcmplx(gamma_save(k1-1,Nsave),gamma_save(k1,Nsave))
    end do ! f
    do f  = 1, alphac     
    do iw = 1, N 
      k1 = k1 + 2
      Wcirr(iw,f) = dcmplx(gamma_save(k1-1,Nsave),gamma_save(k1,Nsave))
    end do ! iw
    end do ! f
    do f = 1, alphar     
      k1 = k1 + 1
      brirr(f) = gamma_save(k1,Nsave)
    end do ! f
    do f  = 1, alphar     
    do iw = 1, N
      k1 = k1 + 1
      Wrirr(iw,f) = gamma_save(k1,Nsave)
    end do ! iw
    end do ! f
    r1 = 0d0 
    do orbidx = 1, N_Slater 
      k1 = k1 + 1 
      Slater(orbidx) = gamma_save(k1,Nsave)
      r1 = r1 + abs(Slater(orbidx))
    end do ! orbidx 
    r1 = r1 / dble(N_Slater)
    if( k1 /= Nv ) stop 'k1/=Nv'
    if( abs(r1-slater_abs_av) > 1d-6 ) stop 'r1/=slater_abs_av'

    rewind(900)
    do k1 = 1, Nv
      write(900,*) k1, gamma_save(k1,Nsave)
    end do ! k1
    close(900)

    if( mod(iteration,100) == 0 ) then
      rewind(800+iteration/100)
      do k1 = 1, Nv
        write(800+iteration/100,*) k1, gamma_save(k1,Nsave)
      end do ! k1
      close(800+iteration/100)
    end if 

  end if

  if( allocated(s_save_gomi) ) deallocate(s_save_gomi)
