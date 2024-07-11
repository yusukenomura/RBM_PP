!
! initialize the weight for spin quantum projection
! algorithm taken from mVMC package
!
subroutine InitSPWeight(myrank) 
  use mod_main , only : spin_parity 
  use mod_PP   , only : N_SP_Gauss_Leg, SP_GL_Cos, SP_GL_Sin, SP_GL_CosCos, SP_GL_CosSin, SP_GL_SinSin, para_spin_proj 
  implicit none 
  integer , intent(in) :: myrank
  integer :: spidx
  real(8) , allocatable :: beta(:), weight(:)
  real(8) :: w 
  real(8) :: r1
  real(8) :: pi = dacos(-1d0)

  allocate( beta  (N_SP_Gauss_Leg) ); beta = 0d0 
  allocate( weight(N_SP_Gauss_Leg) ); weight = 0d0 

  allocate( SP_GL_Cos(N_SP_Gauss_Leg) ); SP_GL_Cos = 0d0
  allocate( SP_GL_Sin(N_SP_Gauss_Leg) ); SP_GL_Sin = 0d0
  allocate( SP_GL_CosCos(N_SP_Gauss_Leg) ); SP_GL_CosCos = 0d0
  allocate( SP_GL_CosSin(N_SP_Gauss_Leg) ); SP_GL_CosSin = 0d0
  allocate( SP_GL_SinSin(N_SP_Gauss_Leg) ); SP_GL_SinSin = 0d0

  call GaussLeg(N_SP_Gauss_Leg,0d0,pi,beta,weight)

  if( N_SP_Gauss_Leg == 1 ) then 
    beta(1) = 0d0
  else if( N_SP_Gauss_Leg == 2 ) then 
    beta(1) = -0.5d0*pi
    beta(2) =  0.5d0*pi
  end if

  if( myrank == 0 ) write(6,*) 
  r1 = 0d0 
  do spidx = 1, N_SP_Gauss_Leg 
    SP_GL_Cos(spidx) = sqrt(2d0)*cos(0.5d0*beta(spidx))    
    SP_GL_Sin(spidx) = sqrt(2d0)*sin(0.5d0*beta(spidx))    
    SP_GL_CosCos(spidx) = SP_GL_Cos(spidx)*SP_GL_Cos(spidx)
    SP_GL_CosSin(spidx) = SP_GL_Cos(spidx)*SP_GL_Sin(spidx)
    SP_GL_SinSin(spidx) = SP_GL_Sin(spidx)*SP_GL_Sin(spidx)

!! change S = 0 only !! 
    w = 0.5d0 * sin(beta(spidx)) * weight(spidx)
!! change S = 0 only !!
    para_spin_proj(spidx) = w

    r1 = r1 + w 
    if( myrank == 0 ) write(6,'(a,4F10.5)') ' chk for beta/pi, weight, w, r1 =', beta(spidx)/pi, weight(spidx), w, r1
 
  end do ! spidx 

  if( N_SP_Gauss_Leg == 1 ) then 
    para_spin_proj(1) = 1d0 
  else if( N_SP_Gauss_Leg == 2 ) then 
    if( spin_parity == 1 ) then 
      para_spin_proj(1) =  1d0 
      para_spin_proj(2) =  1d0 
    else 
      para_spin_proj(1) =  1d0 
      para_spin_proj(2) = -1d0 
    end if
  end if

  if( myrank == 0 ) then 
    write(6,'(a,F13.8)') ' chk for \int_0^pi sin(x) dx (should be 2) =', r1*2d0
    write(6,*) 
    write(6,'(a)') ' para_spin_proj, SP_GL_Cos, SP_GL_Sin'
    do spidx = 1, N_SP_Gauss_Leg
      write(6,'(I5,3F20.10)') spidx, para_spin_proj(spidx), SP_GL_Cos(spidx), SP_GL_Sin(spidx)
    end do ! spidx
    write(6,*) 
    write(6,*) 
  end if

  deallocate(beta,weight)
  return 

end subroutine 
!
! set parameters for gauss-legendre quadrature
! algorithm taken from mVMC package
!
subroutine GaussLeg(n,x1,x2,x,w)
  implicit none 
  integer, intent(in) :: n 
  real(8) , intent(in)  :: x1, x2
  real(8) , intent(out) :: x(n) 
  real(8) , intent(out) :: w(n) 
  integer ::  m, k, j, i
  real(8) ::  z1, z, xm, xl, pp, p3, p2, p1
  real(8) ::  pi = dacos(-1d0)

  m =(n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i = 1, m
    z = cos( pi * (dble(i-1)+0.75d0)/(dble(n)+0.5d0) )
    do k = 1, 100000
      p1 = 1.0d0
      p2 = 0.0d0 
      do j = 1, n
        p3 = p2
        p2 = p1
        p1 = ( (2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/dble(j)
      end do 
      pp = dble(n) * (z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z = z1 - p1/pp
      if( abs(z-z1) < 4.0d-14 ) exit 
    end do ! k  
    x(i) = xm - xl*z
    x(n+1-i)= xm + xl*z
    w(i) = 2.0d0*xl/((1.0d0-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do ! i  

  return

end subroutine 

