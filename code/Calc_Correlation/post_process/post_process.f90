program post_process
  integer :: L 
  integer :: N 
  integer :: Nav
  integer :: Nchi
  integer , allocatable :: map1(:,:) 
  integer , allocatable :: map2(:,:) 
  integer :: i0, i1, i2, i3, i4
  integer :: i  
  integer :: x, x1, x2, y, y1, y2
  integer :: kx, ky 
  real(8) :: r1, r2, r3
  real(8) :: r4, r5, r6
  real(8) :: r7, r8
  real(8) , allocatable :: Chi_R_SzSz(:,:,:)
  real(8) , allocatable :: Chi_R_SpSm(:,:,:)
  real(8) , allocatable :: Chi_q_SzSz(:,:,:)
  real(8) , allocatable :: Chi_q_SpSm(:,:,:)
  real(8) , allocatable :: Chi_SzSz_10(:,:,:)
  real(8) , allocatable :: Chi_SzSz_01(:,:,:)
  real(8) , allocatable :: Chi_SzSz_11(:,:,:)
  real(8) , allocatable :: Chi_SpSm_10(:,:,:)
  real(8) , allocatable :: Chi_SpSm_01(:,:,:)
  real(8) , allocatable :: Chi_SpSm_11(:,:,:)
  real(8) , allocatable :: Chi_R_VBS(:,:,:,:)
  real(8) , allocatable :: Chi_q_VBS(:,:,:,:)
  real(8) , allocatable :: Chi_tmp1(:,:)
  real(8) , allocatable :: Chi_tmp2(:)
  real(8) , allocatable :: Chi_tmp3(:)
  real(8) , allocatable :: Chi_tmp4(:)
  real(8) , allocatable :: Chi_tmp5(:)
  integer , allocatable :: chk(:,:)
  real(8) , allocatable :: Chi_R_SzSz_av(:,:)
  real(8) , allocatable :: Chi_R_SpSm_av(:,:)
  real(8) :: tpi = 2d0*dacos(-1d0)
  real(8) :: Etot, J1, J2

  open(unit=1,file='RBM.input',status='old')
  read(1,*)
  read(1,*) L
  read(1,*)
  read(1,*) i1, i2, i3, Nav
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) J1, J2
  close(1)
  N = L*L
  Nchi = N*(N+1)/2

  allocate( map1(2,N) )
  allocate( map2(0:L-1,0:L-1) )
  i = 0
  do i1 = 0, L-1
  do i2 = 0, L-1
    i = i + 1
    map1(1,i) = i1  
    map1(2,i) = i2  
    map2(i1,i2) = i  
    x = (i-1) / L 
    y = (i-1) - x * L  
    if( i1 /= x ) stop 'lattice definition not consistent'
    if( i2 /= y ) stop 'lattice definition not consistent'
  end do ! i2
  end do ! i1
  if( i /= N ) stop 'map wrong' 

  allocate( Chi_R_SzSz(Nav,0:L-1,0:L-1) ); Chi_R_SzSz = 0d0
  allocate( Chi_R_SpSm(Nav,0:L-1,0:L-1) ); Chi_R_SpSm = 0d0

  allocate( Chi_q_SzSz(Nav,0:L,0:L) ); Chi_q_SzSz = 0d0
  allocate( Chi_q_SpSm(Nav,0:L,0:L) ); Chi_q_SpSm = 0d0

  allocate( Chi_R_VBS(Nav,0:L-1,0:L-1,4) ); Chi_R_VBS = 0d0
  allocate( Chi_q_VBS(Nav,0:L  ,0:L  ,4) ); Chi_q_VBS = 0d0

  allocate( Chi_SzSz_10(Nav,0:L-1,N) ); Chi_SzSz_10 = 0d0
  allocate( Chi_SzSz_01(Nav,0:L-1,N) ); Chi_SzSz_01 = 0d0
  allocate( Chi_SzSz_11(Nav,0:L-1,N) ); Chi_SzSz_11 = 0d0

  allocate( Chi_SpSm_10(Nav,0:L-1,N) ); Chi_SpSm_10 = 0d0
  allocate( Chi_SpSm_01(Nav,0:L-1,N) ); Chi_SpSm_01 = 0d0
  allocate( Chi_SpSm_11(Nav,0:L-1,N) ); Chi_SpSm_11 = 0d0

  allocate( Chi_R_SzSz_av(0:L-1,0:L-1) ); Chi_R_SzSz_av = 0d0
  allocate( Chi_R_SpSm_av(0:L-1,0:L-1) ); Chi_R_SpSm_av = 0d0

  allocate( Chi_tmp1(N,N) ); Chi_tmp1 = 0d0
  allocate( Chi_tmp2(Nav) ); Chi_tmp2 = 0d0
  allocate( Chi_tmp3(Nav) ); Chi_tmp3 = 0d0
  allocate( Chi_tmp4(Nav) ); Chi_tmp4 = 0d0
  allocate( Chi_tmp5(Nav) ); Chi_tmp5 = 0d0

  allocate(chk(0:L-1,0:L-1) ); chk = 0

  open(300,file='m2_SS_vs_iteration.txt',status='unknown')
  open(301,file='Chi_VBS_vs_iteration.txt',status='unknown')
  open(302,file='sub_Chi_VBS_vs_iteration.txt',status='unknown')
  do i = 1, Nav
    !  
    ! read data
    ! 
    i0 = 0 
    chk = 0
    read(60000+i,*)
    read(70000+i,*)
    do i2 = 1, N
    do i1 = 1, i2
      i0 = i0 + 1
      read(60000+i,*) i3, r1, r2
      read(70000+i,*) i3, r3, r4, r5, r6
      Chi_tmp1(i1,i2) = 0.25d0 * (r1+2d0*r2)
      Chi_tmp1(i2,i1) = Chi_tmp1(i1,i2)
      x1 = map1(1,i1)
      y1 = map1(2,i1)
      x2 = map1(1,i2)
      y2 = map1(2,i2)
      x = x2-x1 
      y = y2-y1 
      if( x < 0 ) x = x + L
      if( y < 0 ) y = y + L
      if( x < 0 .or. x > L-1 ) stop 'x is wrong 1' 
      if( y < 0 .or. y > L-1 ) stop 'y is wrong 1'
      Chi_R_SzSz(i,x,y)  = Chi_R_SzSz(i,x,y)  + r1  
      Chi_R_SpSm(i,x,y)  = Chi_R_SpSm(i,x,y)  + r2  
      Chi_R_VBS(i,x,y,1) = Chi_R_VBS(i,x,y,1) + r3  
      Chi_R_VBS(i,x,y,2) = Chi_R_VBS(i,x,y,2) + r4  
      Chi_R_VBS(i,x,y,3) = Chi_R_VBS(i,x,y,3) + r5  
      Chi_R_VBS(i,x,y,4) = Chi_R_VBS(i,x,y,4) + r6  
      chk(x,y) = chk(x,y) + 1 
      if( y == 0 .and. x /= 0 ) then 
        Chi_SzSz_10(i,x,chk(x,y)) = r1
        Chi_SpSm_10(i,x,chk(x,y)) = r2
      end if 
      if( x == 0 .and. y /= 0 ) then 
        Chi_SzSz_01(i,y,chk(x,y)) = r1
        Chi_SpSm_01(i,y,chk(x,y)) = r2
      end if 
      if( x == y .and. x /= 0 ) then 
        Chi_SzSz_11(i,x,chk(x,y)) = r1
        Chi_SpSm_11(i,x,chk(x,y)) = r2
      end if 

      x = -x 
      y = -y 
      if( x < 0 ) x = x + L
      if( y < 0 ) y = y + L
      if( x < 0 .or. x > L-1 ) stop 'x is wrong 2'
      if( y < 0 .or. y > L-1 ) stop 'y is wrong 2'
      Chi_R_SzSz(i,x,y)  = Chi_R_SzSz(i,x,y)  + r1  
      Chi_R_SpSm(i,x,y)  = Chi_R_SpSm(i,x,y)  + r2  
      Chi_R_VBS(i,x,y,1) = Chi_R_VBS(i,x,y,1) + r3  
      Chi_R_VBS(i,x,y,2) = Chi_R_VBS(i,x,y,2) + r4  
      Chi_R_VBS(i,x,y,3) = Chi_R_VBS(i,x,y,3) + r5  
      Chi_R_VBS(i,x,y,4) = Chi_R_VBS(i,x,y,4) + r6  
      chk(x,y) = chk(x,y) + 1 
      if( y == 0 .and. x /= 0 ) then 
        Chi_SzSz_10(i,x,chk(x,y)) = r1
        Chi_SpSm_10(i,x,chk(x,y)) = r2
      end if 
      if( x == 0 .and. y /= 0 ) then 
        Chi_SzSz_01(i,y,chk(x,y)) = r1
        Chi_SpSm_01(i,y,chk(x,y)) = r2
      end if 
      if( x == y .and. x /= 0 ) then 
        Chi_SzSz_11(i,x,chk(x,y)) = r1
        Chi_SpSm_11(i,x,chk(x,y)) = r2
      end if 

    end do ! i1 
    end do ! i2
    if( i0 /= Nchi ) stop 'i0 /= Nchi' 
    close(60000+i)
    ! 
    ! average using translational symmetry 
    ! 
    do y = 0, L-1 
    do x = 0, L-1 
      if( x == 0 .and. y == 0 ) then 
        if( chk(x,y) /= 2*N ) stop 'chk(x,y) /= 2*N'
        Chi_R_SzSz(i,x,y)  = Chi_R_SzSz(i,x,y)  * 0.25d0 / dble(2*N)
        Chi_R_SpSm(i,x,y)  = Chi_R_SpSm(i,x,y)  * 0.25d0 / dble(2*N)
        Chi_R_VBS(i,x,y,1) = Chi_R_VBS(i,x,y,1) * 0.25d0 * 0.25d0 / dble(2*N)
        Chi_R_VBS(i,x,y,2) = Chi_R_VBS(i,x,y,2) * 0.25d0 * 0.25d0 / dble(2*N)
        Chi_R_VBS(i,x,y,3) = Chi_R_VBS(i,x,y,3) * 0.25d0 * 0.25d0 / dble(2*N)
        Chi_R_VBS(i,x,y,4) = Chi_R_VBS(i,x,y,4) * 0.25d0 * 0.25d0 / dble(2*N)
        if( abs(Chi_R_SzSz(i,x,y)-0.25d0) > 1d-5 ) stop 'onsite Chi_R wrong 1'
        if( abs(Chi_R_SpSm(i,x,y)-0.25d0) > 1d-5 ) stop 'onsite Chi_R wrong 2'
        Chi_R_SzSz(i,x,y) = 0.25d0
        Chi_R_SpSm(i,x,y) = 0.25d0
      else 
        if( chk(x,y) /= N ) stop 'chk(x,y) /= N'
        Chi_R_SzSz(i,x,y)  = Chi_R_SzSz(i,x,y)  * 0.25d0 / dble(N)
        Chi_R_SpSm(i,x,y)  = Chi_R_SpSm(i,x,y)  * 0.25d0 / dble(N)
        Chi_R_VBS(i,x,y,1) = Chi_R_VBS(i,x,y,1) * 0.25d0 * 0.25d0 / dble(N)
        Chi_R_VBS(i,x,y,2) = Chi_R_VBS(i,x,y,2) * 0.25d0 * 0.25d0 / dble(N)
        Chi_R_VBS(i,x,y,3) = Chi_R_VBS(i,x,y,3) * 0.25d0 * 0.25d0 / dble(N)
        Chi_R_VBS(i,x,y,4) = Chi_R_VBS(i,x,y,4) * 0.25d0 * 0.25d0 / dble(N)
      end if 
    end do ! x
    end do ! y 
    ! 
    ! Correction for Chi_R_VBS
    ! <BiBj> - <Bi><Bj>
    ! 
    do i0 = 1, 4 
      do y = 0, L-1 
      do x = 0, L-1 
        do i1 = 1, N  
          x1 = map1(1,i1)
          y1 = map1(2,i1)
          x2 = x1+x
          y2 = y1+y
          if( x2 > L-1 ) x2 = x2 - L
          if( y2 > L-1 ) y2 = y2 - L
          if( x2 < 0 .or. x2 > L-1 ) stop 'x2 is wrong' 
          if( y2 < 0 .or. y2 > L-1 ) stop 'y2 is wrong'
          i2 = map2(x2,y2)
          call bond_map(L,N,i0,i1,i3) 
          call bond_map(L,N,i0,i2,i4) 
          Chi_R_VBS(i,x,y,i0) = Chi_R_VBS(i,x,y,i0) - Chi_tmp1(i1,i3) * Chi_tmp1(i2,i4) / dble(N) 
        end do ! i1 
      end do ! x
      end do ! y 
    end do ! i0
    ! 
    ! Fourier transform correlation function 
    ! 
    do ky = 0, L 
    do kx = 0, L 
      do y = 0, L-1 
      do x = 0, L-1 
        r1 = dble(kx*x+ky*y) * tpi / dble(L)
        Chi_q_SzSz(i,kx,ky)  =  Chi_q_SzSz(i,kx,ky)  + Chi_R_SzSz(i,x,y)  * cos(r1) 
        Chi_q_SpSm(i,kx,ky)  =  Chi_q_SpSm(i,kx,ky)  + Chi_R_SpSm(i,x,y)  * cos(r1) 
        Chi_q_VBS(i,kx,ky,1) =  Chi_q_VBS(i,kx,ky,1) + Chi_R_VBS(i,x,y,1) * cos(r1) 
        Chi_q_VBS(i,kx,ky,2) =  Chi_q_VBS(i,kx,ky,2) + Chi_R_VBS(i,x,y,2) * cos(r1) 
        Chi_q_VBS(i,kx,ky,3) =  Chi_q_VBS(i,kx,ky,3) + Chi_R_VBS(i,x,y,3) * cos(r1) 
        Chi_q_VBS(i,kx,ky,4) =  Chi_q_VBS(i,kx,ky,4) + Chi_R_VBS(i,x,y,4) * cos(r1) 
      end do ! x
      end do ! y 
    end do ! kx 
    end do ! ky 
    !  
    ! write correlatin function at high-symmetry points
    ! 
    r1 = Chi_q_SzSz(i,L/2,L/2)+2d0*Chi_q_SpSm(i,L/2,L/2)
    r2 = Chi_q_SzSz(i,L/2,0)+2d0*Chi_q_SpSm(i,L/2,0)
    r3 = Chi_q_SzSz(i,0,L/2)+2d0*Chi_q_SpSm(i,0,L/2)
    write(300,'(I5,3F20.10)') i, r1/dble(N), r2/dble(N), r3/dble(N) 
    write(301,'(I5,6F20.10)') i, Chi_q_VBS(i,L/2,0  ,1), Chi_q_VBS(i,L/2,0  ,2), & 
                                 Chi_q_VBS(i,  0,L/2,1), Chi_q_VBS(i,  0,L/2,2), & 
                                 Chi_q_VBS(i,L/2,L/2,1), Chi_q_VBS(i,L/2,L/2,2)  
    write(302,'(I5,6F20.10)') i, Chi_q_VBS(i,L/2,0  ,3), Chi_q_VBS(i,L/2,0  ,4), & 
                                 Chi_q_VBS(i,  0,L/2,3), Chi_q_VBS(i,  0,L/2,4), & 
                                 Chi_q_VBS(i,L/2,L/2,3), Chi_q_VBS(i,L/2,L/2,4)  

  end do ! i
  close(300)
  close(301)
  close(302)
  ! 
  ! write spin correlation function
  ! 
  open(100,file='Chi_SS_R.txt',status='unknown')
  do y = 0, L-1 
    do x = 0, L-1 
      call calc_av_er(Nav,Chi_R_SzSz(:,x,y),r1,r2)
      call calc_av_er(Nav,Chi_R_SpSm(:,x,y),r3,r4)
      Chi_tmp2(:) = Chi_R_SzSz(:,x,y)+2d0*Chi_R_SpSm(:,x,y)
      call calc_av_er(Nav,Chi_tmp2,r5,r6)
      write(100,'(2I5,6F20.10)') x, y, r5, r6, r1, r2, r3, r4
      Chi_R_SzSz_av(x,y) = r1
      Chi_R_SpSm_av(x,y) = r3
    end do ! x 
    write(100,*)
  end do ! y 
  close(100)

  open(101,file='Chi_SS_q.txt',status='unknown')
  open(200,file='total_S.txt',status='unknown')
  open(201,file='m2_SS_pipi.txt',status='unknown')
  open(202,file='m2_SS_pi0.txt',status='unknown')
  open(203,file='m2_SS_0pi.txt',status='unknown')
  write(200,'(a)') '#  total  SzSz  SpSm'
  do ky = 0, L 
    do kx = 0, L 
      call calc_av_er(Nav,Chi_q_SzSz(:,kx,ky),r1,r2)
      call calc_av_er(Nav,Chi_q_SpSm(:,kx,ky),r3,r4)
      Chi_tmp2(:) = Chi_q_SzSz(:,kx,ky)+2d0*Chi_q_SpSm(:,kx,ky)
      call calc_av_er(Nav,Chi_tmp2,r5,r6)
      write(101,'(2I5,6F20.10)') kx, ky, r5, r6, r1, r2, r3, r4
      if( kx == 0 .and. ky == 0 ) then 
        write(200,'(6F20.10)') r5*dble(N), r6*dble(N), & 
                               r1*dble(N), r2*dble(N), &
                               r3*dble(N), r4*dble(N)
      else if( kx == L/2 .and. ky == L/2 ) then 
        write(201,'(6F20.10)') r5/dble(N), r6/dble(N), & 
                               r1/dble(N), r2/dble(N), &
                               r3/dble(N), r4/dble(N)
      else if( kx == L/2 .and. ky == 0 ) then 
        write(202,'(6F20.10)') r5/dble(N), r6/dble(N), & 
                               r1/dble(N), r2/dble(N), &
                               r3/dble(N), r4/dble(N)
      else if( kx == 0 .and. ky == L/2 ) then 
        write(203,'(6F20.10)') r5/dble(N), r6/dble(N), & 
                               r1/dble(N), r2/dble(N), &
                               r3/dble(N), r4/dble(N)
      end if
    end do ! kx 
    write(101,*)
  end do ! ky 
  close(101)
  close(200)
  close(201)
  close(202)
  close(203)
  ! 
  ! write VBS correlation function
  ! 
  open(400,file='Chi_VBS_R.txt',status='unknown')
  do y = 0, L-1 
    do x = 0, L-1 
      call calc_av_er(Nav,Chi_R_VBS(:,x,y,1),r1,r2)
      call calc_av_er(Nav,Chi_R_VBS(:,x,y,2),r3,r4)
      write(400,'(2I5,4F20.10)') x, y, r1, r2, r3, r4
    end do ! x 
    write(400,*)
  end do ! y 
  close(400)

  open(400,file='sub_Chi_VBS_R.txt',status='unknown')
  do y = 0, L-1 
    do x = 0, L-1 
      call calc_av_er(Nav,Chi_R_VBS(:,x,y,3),r1,r2)
      call calc_av_er(Nav,Chi_R_VBS(:,x,y,4),r3,r4)
      write(400,'(2I5,4F20.10)') x, y, r1, r2, r3, r4
    end do ! x 
    write(400,*)
  end do ! y 
  close(400)

  open(401,file='Chi_VBS_q.txt',status='unknown')
  open(501,file='m2_VBS_pipi.txt',status='unknown')
  open(502,file='m2_VBS_pi0.txt',status='unknown')
  open(503,file='m2_VBS_0pi.txt',status='unknown')
  do ky = 0, L 
    do kx = 0, L 
      call calc_av_er(Nav,Chi_q_VBS(:,kx,ky,1),r1,r2)
      call calc_av_er(Nav,Chi_q_VBS(:,kx,ky,2),r3,r4)
      write(401,'(2I5,4F20.10)') kx, ky, r1, r2, r3, r4
      r1 = r1 / dble(N)
      r2 = r2 / dble(N)
      r3 = r3 / dble(N)
      r4 = r4 / dble(N)
      if( kx == L/2 .and. ky == L/2 ) write(501,'(4F20.10)') r1, r2, r3, r4
      if( kx == L/2 .and. ky == 0   ) write(502,'(4F20.10)') r1, r2, r3, r4
      if( kx == 0   .and. ky == L/2 ) write(503,'(4F20.10)') r1, r2, r3, r4
    end do ! kx 
    write(401,*)
  end do ! ky 
  close(401)
  close(501)
  close(502)
  close(503)

  open(401,file='sub_Chi_VBS_q.txt',status='unknown')
  open(501,file='sub_m2_VBS_pipi.txt',status='unknown')
  open(502,file='sub_m2_VBS_pi0.txt',status='unknown')
  open(503,file='sub_m2_VBS_0pi.txt',status='unknown')
  do ky = 0, L 
    do kx = 0, L 
      call calc_av_er(Nav,Chi_q_VBS(:,kx,ky,3),r1,r2)
      call calc_av_er(Nav,Chi_q_VBS(:,kx,ky,4),r3,r4)
      write(401,'(2I5,4F20.10)') kx, ky, r1, r2, r3, r4
      r1 = r1 / dble(N)
      r2 = r2 / dble(N)
      r3 = r3 / dble(N)
      r4 = r4 / dble(N)
      if( kx == L/2 .and. ky == L/2 ) write(501,'(4F20.10)') r1, r2, r3, r4
      if( kx == L/2 .and. ky == 0   ) write(502,'(4F20.10)') r1, r2, r3, r4
      if( kx == 0   .and. ky == L/2 ) write(503,'(4F20.10)') r1, r2, r3, r4
    end do ! kx 
    write(401,*)
  end do ! ky 
  close(401)
  close(501)
  close(502)
  close(503)
  ! 
  ! write correlation ratio
  ! 
  open(204,file='SS_correlation_ratio_1.txt',status='unknown')
  open(701,file='SS_correlation_length_1.txt',status='unknown')
  do i = 1, Nav
    r1 = Chi_q_SzSz(i,L/2  ,L/2  ) + 2d0*Chi_q_SpSm(i,L/2  ,L/2  )
    r2 = Chi_q_SzSz(i,L/2+1,L/2  ) + 2d0*Chi_q_SpSm(i,L/2+1,L/2  )
    r3 = Chi_q_SzSz(i,L/2  ,L/2+1) + 2d0*Chi_q_SpSm(i,L/2  ,L/2+1)
    r4 = Chi_q_SzSz(i,L/2-1,L/2  ) + 2d0*Chi_q_SpSm(i,L/2-1,L/2  )
    r5 = Chi_q_SzSz(i,L/2  ,L/2-1) + 2d0*Chi_q_SpSm(i,L/2  ,L/2-1)
    r6 = 0.25d0*(r2+r3+r4+r5)
    Chi_tmp2(i) = 1d0 - r6/r1
write(999,*) i, Chi_tmp2(i)
    Chi_tmp3(i) = dsqrt(r1/r6-1d0)/2d0/dacos(-1d0)
  end do ! i 
  call calc_av_er(Nav,Chi_tmp2,r1,r2)
  call calc_av_er(Nav,Chi_tmp3,r3,r4)
  write(204,'(2F20.10)') r1, r2
  write(701,'(I5,F10.5,2F20.10)') L, J2, r3, r4
  close(204)
  close(701)

  open(205,file='SS_correlation_ratio_2.txt',status='unknown')
  open(702,file='SS_correlation_length_2.txt',status='unknown')
  do i = 1, Nav
    r1 = Chi_q_SzSz(i,L/2  ,L/2  ) + 2d0*Chi_q_SpSm(i,L/2  ,L/2  )
    r2 = Chi_q_SzSz(i,L/2-1,L/2-1) + 2d0*Chi_q_SpSm(i,L/2-1,L/2-1)
    r3 = Chi_q_SzSz(i,L/2+1,L/2-1) + 2d0*Chi_q_SpSm(i,L/2+1,L/2-1)
    r4 = Chi_q_SzSz(i,L/2-1,L/2+1) + 2d0*Chi_q_SpSm(i,L/2-1,L/2+1)
    r5 = Chi_q_SzSz(i,L/2+1,L/2+1) + 2d0*Chi_q_SpSm(i,L/2+1,L/2+1)
    r6 = 0.25d0*(r2+r3+r4+r5)
    Chi_tmp2(i) = 1d0 - r6/r1
    Chi_tmp3(i) = dsqrt(r1/r6-1d0)/2d0/dacos(-1d0)
  end do ! i 
  call calc_av_er(Nav,Chi_tmp2,r1,r2)
  call calc_av_er(Nav,Chi_tmp3,r3,r4)
  write(205,'(2F20.10)') r1, r2
  write(702,'(I5,F10.5,2F20.10)') L, J2, r3, r4
  close(205)
  close(702)

  open(206,file='VBS_correlation_ratio_1.txt',status='unknown')
  open(703,file='VBS_correlation_length_1.txt',status='unknown')
  do i = 1, Nav
    r1 = Chi_q_VBS(i,L/2  ,0,1)
    r2 = Chi_q_VBS(i,L/2+1,0,1)
    r3 = Chi_q_VBS(i,L/2-1,0,1)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 1'
    r4 = 0.5d0*(r2+r3)
    Chi_tmp2(i) = 1d0 - r4/r1
    r1 = Chi_q_VBS(i,0,L/2  ,2)
    r2 = Chi_q_VBS(i,0,L/2+1,2)
    r3 = Chi_q_VBS(i,0,L/2-1,2)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 2'
    r4 = 0.5d0*(r2+r3)
    Chi_tmp3(i) = 1d0 - r4/r1
    r1 = Chi_q_VBS(i,L/2  ,0,1) + Chi_q_VBS(i,0,L/2  ,2)
    r2 = Chi_q_VBS(i,L/2+1,0,1) + Chi_q_VBS(i,0,L/2+1,2)
    r3 = Chi_q_VBS(i,L/2-1,0,1) + Chi_q_VBS(i,0,L/2-1,2)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 3'
    r4 = 0.5d0*(r2+r3)
    Chi_tmp4(i) = 1d0 - r4/r1
    Chi_tmp5(i) = dsqrt(r1/r4-1d0)/2d0/dacos(-1d0)
  end do ! i 
  call calc_av_er(Nav,Chi_tmp2,r1,r2)
  call calc_av_er(Nav,Chi_tmp3,r3,r4)
  call calc_av_er(Nav,Chi_tmp4,r5,r6)
  call calc_av_er(Nav,Chi_tmp5,r7,r8)
  write(206,'(6F20.10)') r5, r6, r1, r2, r3, r4
  write(703,'(I5,F10.5,2F20.10)') L, J2, r7, r8
  close(206)
  close(703)

  open(207,file='VBS_correlation_ratio_2.txt',status='unknown')
  open(704,file='VBS_correlation_length_2.txt',status='unknown')
  do i = 1, Nav
    r1 = Chi_q_VBS(i,L/2,  0,1)
    r2 = Chi_q_VBS(i,L/2,  1,1)
    r3 = Chi_q_VBS(i,L/2,L-1,1)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 4'
    r4 = 0.5d0*(r2+r3)
    Chi_tmp2(i) = 1d0 - r4/r1
    r1 = Chi_q_VBS(i,  0,L/2,2)
    r2 = Chi_q_VBS(i,  1,L/2,2)
    r3 = Chi_q_VBS(i,L-1,L/2,2)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 5'
    r4 = 0.5d0*(r2+r3)
    Chi_tmp3(i) = 1d0 - r4/r1
    r1 = Chi_q_VBS(i,L/2,  0,1) + Chi_q_VBS(i,  0,L/2,2)
    r2 = Chi_q_VBS(i,L/2,  1,1) + Chi_q_VBS(i,  1,L/2,2)
    r3 = Chi_q_VBS(i,L/2,L-1,1) + Chi_q_VBS(i,L-1,L/2,2)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 6'
    r4 = 0.5d0*(r2+r3)
    Chi_tmp4(i) = 1d0 - r4/r1
    Chi_tmp5(i) = dsqrt(r1/r4-1d0)/2d0/dacos(-1d0)
  end do ! i 
  call calc_av_er(Nav,Chi_tmp2,r1,r2)
  call calc_av_er(Nav,Chi_tmp3,r3,r4)
  call calc_av_er(Nav,Chi_tmp4,r5,r6)
  call calc_av_er(Nav,Chi_tmp5,r7,r8)
  write(207,'(6F20.10)') r5, r6, r1, r2, r3, r4
  write(704,'(I5,F10.5,2F20.10)') L, J2, r7, r8
  close(207)
  close(704)

  open(208,file='VBS_correlation_ratio_0.txt',status='unknown')
  open(705,file='VBS_correlation_length_0.txt',status='unknown')
  do i = 1, Nav
    r1 = Chi_q_VBS(i,L/2  ,0,1)
    r2 = Chi_q_VBS(i,L/2+1,0,1)
    r3 = Chi_q_VBS(i,L/2-1,0,1)
    r4 = Chi_q_VBS(i,L/2,  1,1)
    r5 = Chi_q_VBS(i,L/2,L-1,1)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 7'
    if( abs(r4-r5) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 8'
    r6 = 0.25d0*(r2+r3+r4+r5)
    Chi_tmp2(i) = 1d0 - r6/r1
    r1 = Chi_q_VBS(i,0,L/2  ,2)
    r2 = Chi_q_VBS(i,0,L/2+1,2)
    r3 = Chi_q_VBS(i,0,L/2-1,2)
    r4 = Chi_q_VBS(i,  1,L/2,2)
    r5 = Chi_q_VBS(i,L-1,L/2,2)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 9'
    if( abs(r4-r5) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 10'
    r6 = 0.25d0*(r2+r3+r4+r5)
    Chi_tmp3(i) = 1d0 - r6/r1
    r1 = Chi_q_VBS(i,L/2  ,0,1) + Chi_q_VBS(i,0,L/2  ,2)
    r2 = Chi_q_VBS(i,L/2+1,0,1) + Chi_q_VBS(i,0,L/2+1,2)
    r3 = Chi_q_VBS(i,L/2-1,0,1) + Chi_q_VBS(i,0,L/2-1,2)
    r4 = Chi_q_VBS(i,L/2,  1,1) + Chi_q_VBS(i,  1,L/2,2)
    r5 = Chi_q_VBS(i,L/2,L-1,1) + Chi_q_VBS(i,L-1,L/2,2)
    if( abs(r2-r3) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 11'
    if( abs(r4-r5) > 1d-4 ) stop 'Chi_q_VBS symmetry wrong 12'
    r6 = 0.25d0*(r2+r3+r4+r5)
    Chi_tmp4(i) = 1d0 - r6/r1
    Chi_tmp5(i) = dsqrt(r1/r6-1d0)/2d0/dacos(-1d0)
  end do ! i 
  call calc_av_er(Nav,Chi_tmp2,r1,r2)
  call calc_av_er(Nav,Chi_tmp3,r3,r4)
  call calc_av_er(Nav,Chi_tmp4,r5,r6)
  call calc_av_er(Nav,Chi_tmp5,r7,r8)
  write(208,'(6F20.10)') r5, r6, r1, r2, r3, r4
  write(705,'(I5,F10.5,2F20.10)') L, J2, r7, r8
  close(208)
  close(705)
  ! 
  ! write spin correlation function in real space for specific dirrections 
  ! 
  Chi_SzSz_10(:,:,:) = Chi_SzSz_10(:,:,:) * 0.25d0
  Chi_SzSz_01(:,:,:) = Chi_SzSz_01(:,:,:) * 0.25d0
  Chi_SzSz_11(:,:,:) = Chi_SzSz_11(:,:,:) * 0.25d0
  Chi_SpSm_10(:,:,:) = Chi_SpSm_10(:,:,:) * 0.25d0
  Chi_SpSm_01(:,:,:) = Chi_SpSm_01(:,:,:) * 0.25d0
  Chi_SpSm_11(:,:,:) = Chi_SpSm_11(:,:,:) * 0.25d0

  open(102,file='Chi_SS_10.txt',status='unknown')
  open(103,file='Chi_SS_01.txt',status='unknown')
  open(104,file='Chi_SS_11.txt',status='unknown')
  do i = 1, N
    do x = 0, L-1 
      call calc_av_er(Nav,Chi_SzSz_10(:,x,i),r1,r2)
      call calc_av_er(Nav,Chi_SpSm_10(:,x,i),r3,r4)
      Chi_tmp2(:) = Chi_SzSz_10(:,x,i)+2d0*Chi_SpSm_10(:,x,i)
      call calc_av_er(Nav,Chi_tmp2,r5,r6)
      if( x == 0 ) then 
        if( abs(r1) > 1d-5 ) stop 'Chi_SS_10 wrong 1'
        if( abs(r2) > 1d-5 ) stop 'Chi_SS_10 wrong 2'
        if( abs(r3) > 1d-5 ) stop 'Chi_SS_10 wrong 3'
        if( abs(r4) > 1d-5 ) stop 'Chi_SS_10 wrong 4'
        if( abs(r5) > 1d-5 ) stop 'Chi_SS_10 wrong 5'
        if( abs(r6) > 1d-5 ) stop 'Chi_SS_10 wrong 6'
        r1 = 0.25d0
        r3 = 0.25d0
        r5 = 0.75d0
      end if 
      write(102,'(2I5,6F20.10)') x, i, r5, r6, r1, r2, r3, r4

      call calc_av_er(Nav,Chi_SzSz_01(:,x,i),r1,r2)
      call calc_av_er(Nav,Chi_SpSm_01(:,x,i),r3,r4)
      Chi_tmp2(:) = Chi_SzSz_01(:,x,i)+2d0*Chi_SpSm_01(:,x,i)
      call calc_av_er(Nav,Chi_tmp2,r5,r6)
      if( x == 0 ) then 
        if( abs(r1) > 1d-5 ) stop 'Chi_SS_01 wrong 1'
        if( abs(r2) > 1d-5 ) stop 'Chi_SS_01 wrong 2'
        if( abs(r3) > 1d-5 ) stop 'Chi_SS_01 wrong 3'
        if( abs(r4) > 1d-5 ) stop 'Chi_SS_01 wrong 4'
        if( abs(r5) > 1d-5 ) stop 'Chi_SS_01 wrong 5'
        if( abs(r6) > 1d-5 ) stop 'Chi_SS_01 wrong 6'
        r1 = 0.25d0
        r3 = 0.25d0
        r5 = 0.75d0
      end if 
      write(103,'(2I5,6F20.10)') x, i, r5, r6, r1, r2, r3, r4
 
      call calc_av_er(Nav,Chi_SzSz_11(:,x,i),r1,r2)
      call calc_av_er(Nav,Chi_SpSm_11(:,x,i),r3,r4)
      Chi_tmp2(:) = Chi_SzSz_11(:,x,i)+2d0*Chi_SpSm_11(:,x,i)
      call calc_av_er(Nav,Chi_tmp2,r5,r6)
      if( x == 0 ) then 
        if( abs(r1) > 1d-5 ) stop 'Chi_SS_11 wrong 1'
        if( abs(r2) > 1d-5 ) stop 'Chi_SS_11 wrong 2'
        if( abs(r3) > 1d-5 ) stop 'Chi_SS_11 wrong 3'
        if( abs(r4) > 1d-5 ) stop 'Chi_SS_11 wrong 4'
        if( abs(r5) > 1d-5 ) stop 'Chi_SS_11 wrong 5'
        if( abs(r6) > 1d-5 ) stop 'Chi_SS_11 wrong 6'
        r1 = 0.25d0
        r3 = 0.25d0
        r5 = 0.75d0
      end if 
      write(104,'(2I5,6F20.10)') x, i, r5, r6, r1, r2, r3, r4

    end do ! x 
    write(102,*)
    write(103,*)
    write(104,*)
  end do ! i 
  close(102)
  close(103)
  close(104)
  ! 
  ! check energy 
  ! 
  Etot = 0d0
  do i0 = 1, 4 
    do i1 = 1, N 
      call bond_map(L,N,i0,i1,i2) 
      x1 = map1(1,i1)
      y1 = map1(2,i1)
      x2 = map1(1,i2)
      y2 = map1(2,i2)
      x = x2-x1 
      y = y2-y1 
      if( x < 0 ) x = x + L
      if( y < 0 ) y = y + L
      if( x < 0 .or. x > L-1 ) stop 'x is wrong 3' 
      if( y < 0 .or. y > L-1 ) stop 'y is wrong 3'
      if( x == 0 .and. y == 0 ) stop 'x == 0 .and. y == 0'
      if( i0 <= 2 ) then 
        Etot = Etot + J1 * (Chi_R_SzSz_av(x,y)+2d0*Chi_R_SpSm_av(x,y))
      else 
        Etot = Etot + J2 * (Chi_R_SzSz_av(x,y)+2d0*Chi_R_SpSm_av(x,y))
      end if 
    end do ! i1 
  end do ! i0 
  open(105,file='E_chk.txt',status='unknown')
  write(105,'(F20.12)') Etot*4d0 
  close(105)

  deallocate(map1,map2,chk)
  deallocate(Chi_tmp1,Chi_tmp2,Chi_tmp3,Chi_tmp4)
  deallocate(Chi_R_SzSz,Chi_R_SpSm)
  deallocate(Chi_q_SzSz,Chi_q_SpSm)
  deallocate(Chi_SzSz_10,Chi_SpSm_10)
  deallocate(Chi_SzSz_01,Chi_SpSm_01)
  deallocate(Chi_SzSz_11,Chi_SpSm_11)
  deallocate(Chi_R_SzSz_av,Chi_R_SpSm_av)
  deallocate(Chi_R_VBS,Chi_q_VBS)

end program 
