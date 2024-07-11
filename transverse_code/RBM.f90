!
!
!
function cCosh(x,y)
  implicit none
  real(8) , intent(in) :: x, y
  complex(8) :: cCosh

  cCosh = 0.5d0 * dcmplx(cosh(x)*cos(y),sinh(x)*sin(y))

end function
!
!
!
function cTanh(x,y)
  implicit none
  real(8) , intent(in) :: x, y
  complex(8) :: cTanh
  complex(8) :: c1, c2

  c1 = dcmplx(tanh(x),tan(y))
  c2 = dcmplx(1d0,tanh(x)*tan(y))
  cTanh = c1/c2

end function
!
! Subroutine for calculating theta 
!   see Eq.(S9) in Supplementary Materials, Carleo and Troyer, Science (2017)
!   theta(j) = \sum_i  W_ij * \sigma_i 
! 
subroutine calc_theta(sigma,thetac,thetar)
  use mod_main , only : N, N_MP_Trans
  use mod_RBM  , only : alphac, alphar, Mc, Mr, Ncopy, W_map, Wcirr, Wrirr, bcirr, brirr
  implicit none
  integer    , intent(in)  :: sigma(N)
  complex(8) , intent(out) :: thetac(Mc)
  real(8)    , intent(out) :: thetar(Mr)
  integer :: f, j, i, mpidx, spidx, icopy, iw
  integer :: Mtmpc, Mtmpr

  Mtmpc = 2 * Ncopy * alphac  
  Mtmpr = 2 * Ncopy * alphar  

!$omp parallel default(shared) private(mpidx,spidx,icopy,f,i,j,iw)
!$omp do
  do mpidx = 1, N_MP_Trans   ! for momentum projection
    j = (mpidx-1) * Mtmpc
    do spidx = 1, -1, -2     ! for spin-parity projection  
    do icopy = 1, Ncopy      ! loop for copies 
    do f     = 1, alphac     ! loop for independent neuron

      j = j + 1   
      thetac(j) = bcirr(f)

      do i = 1, N
        iw = W_map(i,icopy,mpidx)  
        thetac(j) = thetac(j) +  Wcirr(iw,f) * dble(sigma(i)) * dble(spidx)
      end do ! i 

    end do ! f 
    end do ! icopy 
    end do ! spidx 
    if( j /= mpidx*Mtmpc ) stop 'j /= mpidx*Mtmpc'
  end do ! mpidx 
!$omp end do

!$omp do
  do mpidx = 1, N_MP_Trans   ! for momentum projection
    j = (mpidx-1) * Mtmpr
    do spidx = 1, -1, -2     ! for spin-parity projection  
    do icopy = 1, Ncopy      ! loop for copies 
    do f     = 1, alphar     ! loop for independent neuron

      j = j + 1   
      thetar(j) = brirr(f)

      do i = 1, N
        iw = W_map(i,icopy,mpidx)  
        thetar(j) = thetar(j) +  Wrirr(iw,f) * dble(sigma(i)) * dble(spidx)
      end do ! i 

    end do ! f 
    end do ! icopy 
    end do ! spidx 
    if( j /= mpidx*Mtmpr ) stop 'j /= mpidx*Mtmpr'
  end do ! mpidx 
!$omp end do
!$omp end parallel

  return 
end subroutine 
!
! fast update of theta
!
subroutine update_theta(sigma,i1,i2,thetac,thetar)
  use mod_main , only : N, N_MP_Trans
  use mod_RBM  , only : alphac, alphar, Mc, Mr, Ncopy, W_map, Wcirr, Wrirr
  implicit none
  integer    , intent(in)    :: sigma(N)
  integer    , intent(in)    :: i1, i2
  complex(8) , intent(inout) :: thetac(Mc)
  real(8)    , intent(inout) :: thetar(Mr)
  integer :: f, j, mpidx, spidx, icopy, iw1, iw2
  integer :: Mtmpc, Mtmpr

  Mtmpc = 2 * Ncopy * alphac  
  Mtmpr = 2 * Ncopy * alphar  

!$omp parallel default(shared) private(mpidx,spidx,icopy,f,j,iw1,iw2)
!$omp do
  do mpidx = 1, N_MP_Trans   ! for momentum projection
    j = (mpidx-1) * Mtmpc
    do spidx = 1, -1, -2     ! for spin-parity projection  
    do icopy = 1, Ncopy      ! loop for copies 
    do f     = 1, alphac     ! loop for independent neuron
      j = j + 1   
      iw1 = W_map(i1,icopy,mpidx)  
      iw2 = W_map(i2,icopy,mpidx)  
      thetac(j) = thetac(j) - 2d0 * Wcirr(iw1,f) * dble(sigma(i1)) * dble(spidx) 
      thetac(j) = thetac(j) - 2d0 * Wcirr(iw2,f) * dble(sigma(i2)) * dble(spidx)
    end do ! f 
    end do ! icopy 
    end do ! spidx 
    if( j /= mpidx*Mtmpc ) stop 'j /= mpidx*Mtmpc'
  end do ! mpidx 
!$omp end do

!$omp do
  do mpidx = 1, N_MP_Trans   ! for momentum projection
    j = (mpidx-1) * Mtmpr
    do spidx = 1, -1, -2     ! for spin-parity projection  
    do icopy = 1, Ncopy      ! loop for copies 
    do f     = 1, alphar     ! loop for independent neuron
      j = j + 1   
      iw1 = W_map(i1,icopy,mpidx)  
      iw2 = W_map(i2,icopy,mpidx)  
      thetar(j) = thetar(j) - 2d0 * Wrirr(iw1,f) * dble(sigma(i1)) * dble(spidx)
      thetar(j) = thetar(j) - 2d0 * Wrirr(iw2,f) * dble(sigma(i2)) * dble(spidx)
    end do ! f 
    end do ! icopy 
    end do ! spidx 
    if( j /= mpidx*Mtmpr ) stop 'j /= mpidx*Mtmpr'
  end do ! mpidx 
!$omp end do
!$omp end parallel

  return 
end subroutine 
