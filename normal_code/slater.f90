!
! set the entire slater matrix 
! algorithm taken from mVMC package
!
! Notice : We possess transposed matrix  
!
subroutine update_slater_elements() 
  use mod_main , only : N, N_QP_Full, N_MP_Trans 
  use mod_PP   , only : N_SP_Gauss_Leg, Slater_Elm, Slater, SP_GL_CosSin, SP_GL_CosCos, SP_GL_SinSin, & 
                        MP_Trans_Idx, MP_Trans_Sgn, Orbital_Idx, Orbital_Sgn
  implicit none 
  integer :: ri, ori, tri, sgni
  integer :: rj, orj, trj, sgnj
  integer :: qpidx, mpidx, spidx
  integer :: i1, i2
  real(8) :: cs, cc, ss
  real(8) :: slt_ij, slt_ji 
  real(8) :: rdum1
  
  qpidx = 0 
  do mpidx = 1, N_MP_Trans
  do spidx = 1, N_SP_Gauss_Leg
    qpidx = qpidx + 1 

    cs = SP_GL_CosSin(spidx)
    cc = SP_GL_CosCos(spidx)
    ss = SP_GL_SinSin(spidx)

    do ri = 1, N 
      ori  = ri 
      tri  = MP_Trans_Idx(ori,mpidx)
      sgni = MP_Trans_Sgn(ori,mpidx)
      do rj = 1, N 
        orj  = rj 
        trj  = MP_Trans_Idx(orj,mpidx)
        sgnj = MP_Trans_Sgn(orj,mpidx)

        i1 = Orbital_Idx(tri,trj)
        i2 = Orbital_Idx(trj,tri)
        slt_ij = Slater(i1) * dble( Orbital_Sgn(tri,trj)*sgni*sgnj )
        slt_ji = Slater(i2) * dble( Orbital_Sgn(trj,tri)*sgni*sgnj )
        Slater_Elm(rj,  ri,  qpidx) = -(slt_ij - slt_ji)*cs   ! up   - up
        Slater_Elm(rj+N,ri,  qpidx) =  slt_ij*cc + slt_ji*ss  ! up   - down
        Slater_Elm(rj  ,ri+N,qpidx) = -slt_ij*ss - slt_ji*cc  ! down - up
        Slater_Elm(rj+N,ri+N,qpidx) =  (slt_ij - slt_ji)*cs   ! down - down 
      end do ! rj 
    end do ! ri
    !  
    ! check 
    ! 
    do i2 = 1, 2*N 
    do i1 = 1, 2*N 
      rdum1 = Slater_Elm(i1,i2,qpidx) + Slater_Elm(i2,i1,qpidx)
      if( abs(rdum1) > 1.0d-10 ) stop 'Slater_Elm is not skew-symmetric' 
    end do ! i1 
    end do ! i2
 
  end do ! spidx  
  end do ! mpidx

  if( qpidx /= N_QP_Full ) stop 'qpidx /= N_QP_Full (slater)'

  return
end subroutine 
!
! calculate the derivative of Pfaffian with respect to fij parameters
!
subroutine calculate_derivative_Pfaffian(deriv,Ele_Idx,Pf_M,Inv_M)
  use mod_main , only : N, N_QP_Full, N_MP_Trans 
  use mod_PP   , only : N_Slater, N_SP_Gauss_Leg
  implicit none
  real(8) , intent(out) :: deriv(N_Slater,N_QP_Full) 
  integer , intent(in)  :: Ele_Idx(N)
  real(8) , intent(in)  :: Pf_M(N_QP_Full)
  real(8) , intent(in)  :: Inv_M(N,N,N_QP_Full)

  integer :: spidx, mpidx, qpidx

!$omp parallel do default(shared) private(qpidx,mpidx,spidx)
  do mpidx = 1, N_MP_Trans
    qpidx = (mpidx-1)*N_SP_Gauss_Leg 
    do spidx = 1, N_SP_Gauss_Leg
      qpidx = qpidx + 1 
      call calculate_derivative_Pfaffian_child(mpidx,spidx,deriv(:,qpidx),Ele_Idx,Pf_M(qpidx),Inv_M(:,:,qpidx))
    end do ! spidx 
  end do ! mpidx 
!$omp end parallel do 

  return 
end subroutine 
!
! calculate the derivative of Pfaffian with respect to fij parameters
!
subroutine calculate_derivative_Pfaffian_child(mpidx,spidx,deriv,Ele_Idx,Pf_M,Inv_M)
  use mod_main , only : N
  use mod_PP   , only : Npair, N_Slater, SP_GL_CosSin, SP_GL_CosCos, SP_GL_SinSin, & 
                        MP_Trans_Idx, MP_Trans_Sgn, Orbital_Idx, Orbital_Sgn
  implicit none
  integer , intent(in)  :: mpidx, spidx
  real(8) , intent(out) :: deriv(N_Slater)
  integer , intent(in)  :: Ele_Idx(N)
  real(8) , intent(in)  :: Pf_M
  real(8) , intent(in)  :: Inv_M(N,N)

  integer :: orbidx
  integer :: msi, msj, ri, rj, ori, orj, tri, trj, sgni, sgnj
  real(8) :: cs, cc, ss 
  integer :: Trans_Orb_Idx(N,N)
  integer :: Trans_Orb_Sgn(N,N)

  deriv(:) = 0d0 
  do msj = 1, N !! N = 2*Npair  
    rj = Ele_Idx(msj)
    orj = rj 
    trj  = MP_Trans_Idx(orj,mpidx)
    sgnj = MP_Trans_Sgn(orj,mpidx)
    do msi = 1, N 

      ri = Ele_Idx(msi)
      ori = ri 
      tri  = MP_Trans_Idx(ori,mpidx)
      sgni = MP_Trans_Sgn(ori,mpidx)

      Trans_Orb_Idx(msi,msj) = Orbital_Idx(tri,trj) 
      Trans_Orb_Sgn(msi,msj) = Orbital_Sgn(tri,trj)*sgni*sgnj 
      !if( Trans_Orb_Sgn(msi,msj,mpidx) /= 1 ) write(6,*) 'Warning: Trans_Orb_Sgn /= 1 '

    end do ! msi
  end do ! msj
  !
  ! calculating Tr(X^{-1}*dX/df_{msi,msj})=-2*alpha(sigma(msi),sigma(msj))(X^{-1})_{msi,msj}
  ! 1/2 factor Tahara paper Eq. (87) and 2 cancels out ?? 
  ! 
  ! Notice : The sign before InvM is changed because we possess transposed matrix
  ! 
  cs = Pf_M*SP_GL_CosSin(spidx)
  cc = Pf_M*SP_GL_CosCos(spidx)
  ss = Pf_M*SP_GL_SinSin(spidx)

  do msj = 1, Npair
  do msi = 1, Npair 
    !
    ! upup
    ! 
    orbidx = Trans_Orb_Idx(msi,msj)
    deriv(orbidx) = deriv(orbidx) &  
                  - Inv_M(msi,msj)*cs*Trans_Orb_Sgn(msi,msj)
    !
    ! updn
    ! 
    orbidx = Trans_Orb_Idx(msi,msj+Npair)
    deriv(orbidx) = deriv(orbidx) & 
                  + Inv_M(msi,msj+Npair)*cc*Trans_Orb_Sgn(msi,msj+Npair)
    !
    ! dnup
    ! 
    orbidx = Trans_Orb_Idx(msi+Npair,msj)
    deriv(orbidx) = deriv(orbidx) & 
                  - Inv_M(msi+Npair,msj)*ss*Trans_Orb_Sgn(msi+Npair,msj)
    !
    ! dndn
    ! 
    orbidx = Trans_Orb_Idx(msi+Npair,msj+Npair)
    deriv(orbidx) = deriv(orbidx) & 
                  + Inv_M(msi+Npair,msj+Npair)*cs*Trans_Orb_Sgn(msi+Npair,msj+Npair)
  end do ! msi  
  end do ! msj 

  return 
end subroutine 
