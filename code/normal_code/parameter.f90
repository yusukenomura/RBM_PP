!
! read ***.def file
! definition follows that of mVMC 
!
subroutine initialize_parameters(myrank)
  use mod_main , only : N, N_MP_Trans, spin_parity
  use mod_RBM  , only : alphac, alphar, bcirr, brirr, Wcirr, Wrirr, W_map  
  use mod_PP   , only : l_AP, N_Slater, Orbital_Idx, Orbital_Sgn, &
                        MP_Trans_Idx, Slater, slater_abs_av 
  implicit none 
  integer , intent(in) :: myrank 
  character(len=25) :: ch1, ch2 
  real(8) :: r1, r2, r3
  real(8) , external :: grnd
  integer :: f, iw, i, i1, i2, i3, i4, mpidx, orbidx
  logical :: file_exist
  ! 
  ! initialize W (variational parameters)   
  ! put random numbers between -0.01 and 0.01
  !
  if( myrank == 0 ) then 
    do f = 1, alphac     
      r1 = 0.1d0*(grnd()-0.5d0)
      r2 = 0.1d0*(grnd()-0.5d0)
      bcirr(f) = dcmplx(r1,r2)
      if( spin_parity == -1 ) bcirr(f) = bcirr(f)*10d0
      do iw = 1, N 
        r1 = 0.1d0*(grnd()-0.5d0)
        r2 = 0.1d0*(grnd()-0.5d0)
        Wcirr(iw,f) = dcmplx(r1,r2)
      end do ! iw 
    end do ! f 
    do f = 1, alphar     
      r1 = 0.1d0*(grnd()-0.5d0)
      brirr(f) = r1
      if( spin_parity == -1 ) brirr(f) = brirr(f)*10d0
      do iw = 1, N
        r1 = 0.1d0*(grnd()-0.5d0)
        Wrirr(iw,f) = r1
      end do ! iw 
    end do ! f 
    file_exist = .false.
    INQUIRE(FILE="W.input", EXIST=file_exist) 
    if( file_exist ) then  
      write(6,'(a)') 'initial variational parameters: read from W.input'
      open(unit=1,file='W.input',status='old')
      read(1,*)
      do f = 1, alphac     
        read(1,*) ch1, ch2, i1, r1, r2
        bcirr(f) = dcmplx(r1,r2)
        if( i1 /= f ) stop 'i1/=f (W.input)'  
        do iw = 1, N 
          read(1,*) i1, i2, r1, r2 
          Wcirr(iw,f) = dcmplx(r1,r2)
          if( i1 /= f ) stop 'i1/=f (W.input)'  
          if( i2 /= iw ) stop 'i1/=f (W.input)'  
        end do ! iw 
        read(1,*)
      end do ! f 
      do f = 1, alphar     
        read(1,*) ch1, ch2, i1, r1, r2
        brirr(f) = r1
        if( abs(r2) > 1d-15 ) stop 'r2 /= 0 (W.input)'
        if( i1 /= f ) stop 'i1/=f (W.input)'  
        do iw = 1, N 
          read(1,*) i1, i2, r1, r2 
          Wrirr(iw,f) = r1
          if( abs(r2) > 1d-15 ) stop 'r2 /= 0 (W.input)'
          if( i1 /= f ) stop 'i1/=f (W.input)'  
          if( i2 /= iw ) stop 'i1/=f (W.input)'  
        end do ! iw 
        read(1,*)
      end do ! f 
      close(1)
    end if
  end if
  ! 
  ! orbital index 
  ! 
  l_AP = .false.
  open(unit=1,file='orbitalidx.def',status='old')
  read(1,*)
  read(1,*) ch1, i1 
  if( i1 /= N_Slater ) stop 'error1, orbitalidx' 
  read(1,*)
  read(1,*)
  read(1,*)
  do i1 = 1, N
  do i2 = 1, N
    if( l_AP ) then  
      read(1,*) i3, i4, Orbital_Idx(i1,i2), Orbital_Sgn(i1,i2)
    else 
      read(1,*) i3, i4, Orbital_Idx(i1,i2)
    end if
    Orbital_Idx(i1,i2) = Orbital_Idx(i1,i2) + 1 
    if( i1 /= i3+1 ) stop 'error2, orbitalidx' 
    if( i2 /= i4+1 ) stop 'error3, orbitalidx' 
    if( Orbital_Idx(i1,i2) < 1 .or. Orbital_Idx(i1,i2) > N_Slater ) stop 'error4, orbitalidx'
  end do ! i2
  end do ! i1
  close(1)
  ! 
  ! quantum projection for translation
  ! 
  do mpidx = 1, N_MP_Trans
    do i = 1, N 
      MP_Trans_Idx(i,mpidx) = W_map(i,1,mpidx)
    end do ! i 
  end do ! mpidx 
  close(1)
  ! 
  !  slater part 
  ! 
  open(unit=1,file='zqp_orbital_initial.dat',status='old')
  read(1,*)
  read(1,*) ch1, i1
  if( i1 /= N_Slater ) stop 'error in zqp_orbital_initial 1'
  read(1,*)
  read(1,*)
  read(1,*)
  r3 = 0d0 
  do orbidx = 1, N_Slater
    read(1,*) i1, r1, r2
    Slater(orbidx) = r1
    r3 = r3 + abs(Slater(orbidx))
    if( abs(r2) > 1.0d-10 ) stop 'imaginary zqp_orbital_initial, not implemented'
    if( orbidx /= i1+1 ) stop 'error in zqp_orbital_initial 2'   
  end do ! orbidx    
  r3 = r3 / dble(N_Slater)

  do orbidx = 1, N_Slater
    Slater(orbidx) = Slater(orbidx) * slater_abs_av / r3 
  end do ! orbidx    

  r3 = 0d0 
  do orbidx = 1, N_Slater
    r3 = r3 + abs(Slater(orbidx))
  end do ! orbidx    
  r3 = r3 / dble(N_Slater)
  if( abs(r3-slater_abs_av) > 1d-6 ) stop 'r3/=slater_abs_av'
  close(1)

  return 
end subroutine 
!
! write initial parameters
!
subroutine write_initial_parameters(myrank)
  use mod_main , only : N, N_MP_Trans, para_mp_trans
  use mod_RBM  , only : alphac, alphar, bcirr, brirr, Wcirr, Wrirr
  use mod_PP   , only : N_Slater, Orbital_Idx, Orbital_Sgn,  &
                        MP_Trans_Idx, MP_Trans_Sgn, Slater 
  implicit none  
  integer , intent(in) :: myrank 
  integer :: f, iw, i, i1, i2, mpidx, orbidx   

  if( myrank == 0 ) then 

    open(unit=1,file='initial_W.txt',status='unknown')
    write(1,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
    do f = 1, alphac     
      write(1,'(a,I6,2F22.15)') '# bias', f, bcirr(f)
      do iw = 1, N 
        write(1,'(2I6,2F22.15)') f, iw, Wcirr(iw,f)
      end do ! iw 
      write(1,*)
    end do ! f 
    do f = 1, alphar     
      write(1,'(a,I6,2F22.15)') '# bias', f, brirr(f), 0d0
      do iw = 1, N 
        write(1,'(2I6,2F22.15)') f, iw, Wrirr(iw,f), 0d0
      end do ! iw 
      write(1,*)
    end do ! f 
    close(1)  

    open(unit=1,file='chk_orbitalidx.def',status='unknown')
    write(1,*) N_Slater
    do i1 = 1, N 
    do i2 = 1, N 
      write(1,'(I5,2I7)') i1-1, i2-1, Orbital_Idx(i1,i2)-1 
      if( abs(Orbital_Sgn(i1,i2)) /= 1 ) stop 'Orbital_Sgn is wrong'
    end do ! i2
    end do ! i1
    close(1)

    open(unit=1,file='chk_qptransidx.def',status='unknown') 
    write(1,*) N_MP_Trans 
    do mpidx = 1, N_MP_Trans
      write(1,'(I3,2F20.12)') mpidx, para_mp_trans(mpidx)
    end do ! i 
    do mpidx = 1, N_MP_Trans
    do i = 1, N
      write(1,'(I5,2I7)') mpidx, i, MP_Trans_Idx(i,mpidx)
      if( abs(MP_Trans_Sgn(i,mpidx)) /= 1 ) stop 'MP_Trans_Sgn is wrong'
    end do ! i 
    end do ! mpidx 
    close(1) 


    open(unit=1,file='chk_zqp_orbital_opt.dat',status='unknown') 
    write(1,'(a,I10)') '# N_Slater', N_Slater
    do orbidx = 1, N_Slater
      write(1,'(I4,2e25.15)') orbidx-1, Slater(orbidx)
    end do ! orbidx 
    close(1)

  end if


  if( myrank == 4 ) then 

    open(unit=11,file='sub_initial_W.txt',status='unknown')
    write(11,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
    do f = 1, alphac     
      write(11,'(a,I6,2F22.15)') '# bias', f, bcirr(f)
      do iw = 1, N 
        write(11,'(2I6,2F22.15)') f, iw, Wcirr(iw,f)
      end do ! iw 
      write(11,*)
    end do ! f 
    do f = 1, alphar     
      write(11,'(a,I6,2F22.15)') '# bias', f, brirr(f), 0d0
      do iw = 1, N 
        write(11,'(2I6,2F22.15)') f, iw, Wrirr(iw,f), 0d0
      end do ! iw 
      write(11,*)
    end do ! f 
    close(11)  

  end if

  return 
end subroutine 
!
! write optimized parameters
!
subroutine write_optimized_parameters(myrank)
  use mod_main , only : N
  use mod_RBM  , only : alphac, alphar, bcirr, brirr, Wcirr, Wrirr
  use mod_PP   , only : N_Slater, Slater
  implicit none  
  integer , intent(in) :: myrank 
  integer :: f, iw, orbidx   

  if( myrank == 0 ) then 

    open(unit=1,file='optimized_W.txt',status='unknown')
    write(1,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
    do f = 1, alphac     
      write(1,'(a,I6,2F22.15)') '# bias', f, bcirr(f)
      do iw = 1, N 
        write(1,'(2I6,2F22.15)') f, iw, Wcirr(iw,f)
      end do ! iw 
      write(1,*)
    end do ! f 
    do f = 1, alphar     
      write(1,'(a,I6,2F22.15)') '# bias', f, brirr(f), 0d0
      do iw = 1, N 
        write(1,'(2I6,2F22.15)') f, iw, Wrirr(iw,f), 0d0
      end do ! iw 
      write(1,*)
    end do ! f 
    close(1)  

    open(unit=1,file='zqp_orbital_opt.dat',status='unknown') 
    write(1,*)
    write(1,'(a,I10)') ' Slater', N_Slater
    write(1,*)
    write(1,*)
    write(1,*)
    do orbidx = 1, N_Slater
      write(1,'(I5,2e25.15)') orbidx-1, Slater(orbidx), 0d0
    end do ! orbidx 
    close(1)

  end if


  if( myrank == 3 ) then 

    open(unit=11,file='sub_optimized_W.txt',status='unknown')
    write(11,'(a)') '# index for independent neurons (f), index for interaction (iw), irreducible interaction Wirr(iw,f)'
    do f = 1, alphac     
      write(11,'(a,I6,2F22.15)') '# bias', f, bcirr(f)
      do iw = 1, N 
        write(11,'(2I6,2F22.15)') f, iw, Wcirr(iw,f)
      end do ! iw 
      write(11,*)
    end do ! f 
    do f = 1, alphar     
      write(11,'(a,I6,2F22.15)') '# bias', f, brirr(f), 0d0
      do iw = 1, N 
        write(11,'(2I6,2F22.15)') f, iw, Wrirr(iw,f), 0d0
      end do ! iw 
      write(11,*)
    end do ! f 
    close(11)  

    open(unit=11,file='sub_zqp_orbital_opt.dat',status='unknown') 
    write(11,*)
    write(11,'(a,I10)') ' Slater', N_Slater
    write(11,*)
    write(11,*)
    write(11,*)
    do orbidx = 1, N_Slater
      write(11,'(I5,2e25.15)') orbidx-1, Slater(orbidx), 0d0
    end do ! orbidx 
    close(11)  

  end if

  return 
end subroutine 
