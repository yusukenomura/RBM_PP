!
! Subroutine for calculating amplitude of RBM+PP wave function
! 
subroutine calc_amplitude_RBM_PP(thetac,thetar,pfaff,sgn_pfaff,psi_x_wo_pf,psi_x)
  use mod_main , only : N_QP_Full, N_MP_Trans, para_mp_trans, wf_av
  use mod_RBM  , only : alphac, alphar, Mc, Mr, Ncopy
  use mod_PP   , only : N_SP_Gauss_Leg, para_spin_proj
  implicit none
  complex(8) , intent(in)  :: thetac(Mc)                  ! thetac(j) = \sum_i  W_ij * \sigma_i
  real(8)    , intent(in)  :: thetar(Mr)                  ! thetar(j) = \sum_i  W_ij * \sigma_i
  real(8)    , intent(in)  :: pfaff(N_QP_Full)            ! pfaffian from pair-product part
  integer    , intent(in)  :: sgn_pfaff(N_MP_Trans)       ! sign for pair-product part (coming from order exchange of operators, gauge transformation)
  complex(8) , intent(out) :: psi_x_wo_pf(N_QP_Full)      ! <x|psi> without Pf  (|x> = |sigma_1, sigma_2, ..., sigma_N>) 
  complex(8) , intent(out) :: psi_x(0:N_MP_Trans*Ncopy*2) ! <x|psi>   (|x> = |sigma_1, sigma_2, ..., sigma_N>) 

  real(8) , external :: grnd 
  real(8) :: r1
  complex(8) :: c1
  complex(8) :: psi_x_wo_pf_tmp(N_QP_Full*Ncopy)
  integer :: j1, j2, f
  integer :: qpidx, mpidx, spidx
  integer :: icopy, itmp 

  !!!!!!!!!! specialized to spin-parity projection !!!!!!!!!!

  if( N_SP_Gauss_Leg /= 2 ) stop 'N_SP_Gauss_Leg /= 2'

!$omp parallel do default(shared) & 
!$omp private(mpidx,spidx,qpidx,icopy,itmp,j1,j2,f,r1,c1)   
  do itmp = 1, N_QP_Full*Ncopy   ! for quantum projection
    qpidx = (itmp-1)/Ncopy + 1  
    mpidx = (qpidx-1)/2 + 1 
    spidx = qpidx - (mpidx-1)*2 
    j1    = (itmp-1) * alphac
    j2    = (itmp-1) * alphar

    c1 = 0.5d0**alphac 
    do f = 1, alphac  ! loop for independent neuron
      j1 = j1 + 1
      c1 = c1 * dcmplx(cosh(dble(thetac(j1)))*cos(dimag(thetac(j1))),sinh(dble(thetac(j1)))*sin(dimag(thetac(j1))))
    end do ! f
    r1 = 1d0
    do f = 1, alphar  ! loop for independent neuron
      j2 = j2 + 1
      r1 = r1 * cosh(thetar(j2)) ! neglect irrelevant factor of 2
    end do ! f
    c1 = c1 * r1 * dble(sgn_pfaff(mpidx)) * para_mp_trans(mpidx) * para_spin_proj(spidx)
    psi_x_wo_pf_tmp(itmp) = c1
    psi_x(itmp) = c1 * pfaff(qpidx) 

    if( j1 /= itmp*alphac ) stop 'j1 /= itmp*alphac'
    if( j2 /= itmp*alphar ) stop 'j2 /= itmp*alphar'
  end do ! itmp
!$omp end parallel do

  psi_x(0) = 0d0
  do itmp = 1, N_MP_Trans*Ncopy*2
    psi_x(0) = psi_x(0) + psi_x(itmp) 
  end do ! itmp  

  itmp = 0
  do qpidx = 1, N_QP_Full
    psi_x_wo_pf(qpidx) = 0d0
    do icopy = 1, Ncopy
      itmp = itmp + 1
      psi_x_wo_pf(qpidx) =  psi_x_wo_pf(qpidx) + psi_x_wo_pf_tmp(itmp)
    end do ! icopy  
  end do ! qpidx  
  if( itmp /= N_QP_Full*Ncopy ) stop 'itmp wrong'


  if( grnd() < 1d-6 ) then 
    c1 = 0d0 
    do qpidx = 1, N_QP_Full
      c1 = c1 + pfaff(qpidx)*psi_x_wo_pf(qpidx)
    end do ! qpidx
    if( abs(c1-psi_x(0))/(abs(c1)+1d-200) > 1d-4 .and. abs(c1-psi_x(0)) > wf_av*1d-8 ) then
      write(6,'(a)') 'Warning psi calculation might be wrong'
      write(6,'(2E15.5)') c1
      write(6,'(2E15.5)') psi_x(0)
      write(6,'(2E15.5)') abs(c1-psi_x(0))
      write(6,'(2E15.5)') abs(c1-psi_x(0))/(abs(c1)+1d-200)
      write(6,*) 
    end if
  end if


  return 
end subroutine 
! 
! calculating average and error
! 
subroutine calc_av_er(N,data_,av,er) 
  implicit none 
  integer , intent(in)  :: N 
  real(8) , intent(in)  :: data_(N)
  real(8) , intent(out) :: av, er
  integer :: i

  av = 0d0
  er = 0d0
  do i = 1, N
    av = av + data_(i)
  end do ! i
  av = av / dble(N)

  do i = 1, N
    er = er + (data_(i)-av)**2
  end do ! i 
  er = er / dble(N)
  er = dsqrt(er)

  return 
end subroutine 
!
!
!
subroutine measure_time(t1,time) 
  implicit none 
  integer , intent(in)    :: t1
  real(8) , intent(inout) :: time  
  integer :: t2, t_rate, t_max

  call system_clock(t2,t_rate,t_max)
 
  if ( t2 < t1 ) then
    time = time + dble(t2-t1+t_max+1) / dble(t_rate)
  else
    time = time + dble(t2-t1) / dble(t_rate)
  endif

  return 
end subroutine 
