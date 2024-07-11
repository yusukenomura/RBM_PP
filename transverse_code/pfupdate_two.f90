!
! calculate new Pfaffian using fast update
! algorithm taken from mVMC package
! 
subroutine calculate_new_PfM_two(msa,msb,Tmp_Ele_Idx,InvM,Pf_old,Pf_new)
  use mod_main , only : N, N_QP_Full
  implicit none 
  integer , intent(in)  :: msa, msb
  integer , intent(in)  :: Tmp_Ele_Idx(N)
  real(8) , intent(in)  :: InvM(N,N,N_QP_Full) 
  real(8) , intent(in)  :: Pf_old(N_QP_Full)
  real(8) , intent(out) :: Pf_new(N_QP_Full)

  integer :: qpidx 

!$omp parallel do default(shared) private(qpidx) 
  do qpidx = 1, N_QP_Full
    call calculate_new_PfM_two_child(qpidx,msa,msb,Tmp_Ele_Idx,InvM(:,:,qpidx),Pf_old(qpidx),Pf_new(qpidx))
  end do ! qpidx
!$omp end parallel do 

  return 
end subroutine 
!
!
!
subroutine calculate_new_PfM_two_child(qpidx,msa,msb,Tmp_Ele_Idx,InvM,Pf_old,Pf_new)
  use mod_main , only : N
  use mod_PP   , only : Npair, Slater_Elm
  implicit none 
  integer , intent(in)  :: qpidx
  integer , intent(in)  :: msa, msb
  integer , intent(in)  :: Tmp_Ele_Idx(N)
  real(8) , intent(in)  :: InvM(N,N)
  real(8) , intent(in)  :: Pf_old
  real(8) , intent(out) :: Pf_new
  integer :: msi, msj
  integer :: rsa, rsb, rsi
  real(8) :: ratio, vec_ba, InvM_ab  
  real(8) :: vec_a(N), vec_b(N) 
  !real(8) :: vec_tmp(N)
  real(8) :: rtmp
  real(8) :: p_a, p_b, q_a, q_b, bMa

  rsa = Tmp_Ele_Idx(msa) + ((msa-1)/Npair)*N
  rsb = Tmp_Ele_Idx(msb) + ((msb-1)/Npair)*N

  do msi = 1, N
    rsi = Tmp_Ele_Idx(msi) + ((msi-1)/Npair)*N
    vec_a(msi) = Slater_Elm(rsi,rsa,qpidx)
    vec_b(msi) = Slater_Elm(rsi,rsb,qpidx)
  end do ! msi 

  vec_ba = vec_b(msa)
  InvM_ab = InvM(msb,msa) !!! InvM is transpose(true InvM)

  p_a = 0d0; p_b = 0d0
  q_a = 0d0; q_b = 0d0
  bMa = 0d0

  do msi = 1, N
    p_a = p_a + InvM(msi,msa) * vec_a(msi)
    p_b = p_b + InvM(msi,msb) * vec_a(msi)
    q_a = q_a + InvM(msi,msa) * vec_b(msi)
    q_b = q_b + InvM(msi,msb) * vec_b(msi)
  end do ! msi

  !vec_tmp(:) = 0d0    
  !call dgemm('T','N',N,1,N, 1.0d0,InvM,N,vec_a,N,0d0,vec_tmp,N)
  !call dgemm('N','N',N,1,N,-1.0d0,InvM,N,vec_a,N,0d0,vec_tmp,N)
  !
  !do msi = 1, N 
  !  bMa = bMa + vec_b(msi) * vec_tmp(msi)
  !end do ! msi

  do msi = 1, N 
    rtmp = 0d0
    do msj = 1, N 
      rtmp = rtmp + InvM(msj,msi) * vec_a(msj)
    end do ! msj 
    bMa = bMa + vec_b(msi) * rtmp
  end do ! msi

  ratio = InvM_ab*vec_ba + InvM_ab*bMa + p_a*q_b - p_b*q_a
  Pf_new = ratio*Pf_old

  return
end subroutine 
!
! calculate new inverse matrix using fast update
! algorithm taken from mVMC package
!
subroutine calculate_new_InvM_two(msa,msb,Tmp_Ele_Idx,InvM)
  use mod_main, only : N, N_QP_Full
  implicit none 
  integer , intent(in)    :: msa, msb
  integer , intent(in)    :: Tmp_Ele_Idx(N)
  real(8) , intent(inout) :: InvM(N,N,N_QP_Full) 

  integer :: qpidx 

!$omp parallel do default(shared) private(qpidx) 
  do qpidx = 1, N_QP_Full
    call calculate_new_InvM_two_child(qpidx,msa,msb,Tmp_Ele_Idx,InvM(:,:,qpidx))
  end do ! qpidx
!$omp end parallel do 
  
  return 
end subroutine 
!
!
!
subroutine calculate_new_InvM_two_child(qpidx,msa,msb,Tmp_Ele_Idx,InvM)
  use mod_main, only : N
  use mod_PP,   only : Npair, Slater_Elm
  implicit none 
  integer , intent(in)    :: qpidx
  integer , intent(in)    :: msa, msb
  integer , intent(in)    :: Tmp_Ele_Idx(N)
  real(8) , intent(inout) :: InvM(N,N) 

  integer :: msi, msj
  integer :: rsa, rsb, rsi
  real(8) :: vecP(N), vecQ(N) 
  real(8) :: vecS(N), vecT(N) 
  real(8) :: Det, InvDet, bMa
  real(8) :: a, b, c, d, e, f
  real(8) :: p_i, p_j, q_i, q_j, s_i, s_j, t_i, t_j

  rsa = Tmp_Ele_Idx(msa) + ((msa-1)/Npair)*N
  rsb = Tmp_Ele_Idx(msb) + ((msb-1)/Npair)*N
  !
  ! vecS, vecT are temporally used as
  ! vecS(i) = sltE(a,j), vecT(i) = sltE(b,j) 
  !
  do msi = 1, N
    rsi = Tmp_Ele_Idx(msi) + ((msi-1)/Npair)*N
    vecS(msi) = Slater_Elm(rsi,rsa,qpidx)
    vecT(msi) = Slater_Elm(rsi,rsb,qpidx)
  end do ! msi
  !
  ! Calculate vecP, vecQ 
  ! vecP(i)= sum_j InvM(i,j)*sltE(a,j) 
  ! vecQ(i)= sum_j InvM(i,j)*sltE(b,j)
  ! 
  !call dgemm('T','N',N,1,N,1.0d0,InvM,N,vecS,N,0d0,vecP,N)
  !call dgemm('T','N',N,1,N,1.0d0,InvM,N,vecT,N,0d0,vecQ,N)
  do msi = 1, N 
    vecP(msi) = 0d0
    vecQ(msi) = 0d0
    do msj = 1, N 
      vecP(msi) = vecP(msi) + InvM(msj,msi)*vecS(msj)
      vecQ(msi) = vecQ(msi) + InvM(msj,msi)*vecT(msj)
    end do ! msj
  end do ! msi
  ! 
  ! Set coefficients 
  ! 
  bMa = 0d0
  do msi = 1, N
    bMa = bMa + vecT(msi) * vecP(msi)
  end do ! msi 

  a = -vecP(msa)
  b =  vecP(msb)
  c =  vecQ(msa)
  d = -vecQ(msb)
  e = -bMa - vecT(msa)
  f = InvM(msb,msa)

  Det = a*d - b*c - e*f
  InvDet = 1.0/Det
  !
  ! Calculate vecS, vecT 
  ! vecS(i)= InvM(a,i)/D, vecT(i)= InvM(b,i)/D 
  !
  do msi = 1, N
    vecS(msi) = InvDet * InvM(msi,msa)
    vecT(msi) = InvDet * InvM(msi,msb)
  end do ! msi 
  !
  ! Update InvM 
  !
  do msi = 1, N
    p_i = vecP(msi)
    q_i = vecQ(msi)
    s_i = vecS(msi)
    t_i = vecT(msi)

    do msj = 1, N
      p_j = vecP(msj)
      q_j = vecQ(msj)
      s_j = vecS(msj)
      t_j = vecT(msj)

      InvM(msj,msi) = InvM(msj,msi)               & 
                    + a * (q_i * t_j - q_j * t_i)       &
                    + b * (q_i * s_j - q_j * s_i)       &
                    + c * (p_i * t_j - p_j * t_i)       &
                    + d * (p_i * s_j - p_j * s_i)       &
                    + e * Det * (s_i * t_j - s_j * t_i) &
                    + f * InvDet * (p_i * q_j - q_i * p_j)
    end do ! msj

    InvM(msa,msi) = InvM(msa,msi) - c*t_i - d*s_i - f*InvDet*q_i
    InvM(msb,msi) = InvM(msb,msi) - a*t_i - b*s_i + f*InvDet*p_i

  end do ! msi

  do msi = 1, N
    p_i = vecP(msi)
    q_i = vecQ(msi)
    s_i = vecS(msi)
    t_i = vecT(msi)

    InvM(msi,msa) = InvM(msi,msa) + c*t_i + d*s_i + f*InvDet*q_i
    InvM(msi,msb) = InvM(msi,msb) + a*t_i + b*s_i - f*InvDet*p_i
  end do ! msi

  InvM(msa,msb) = InvM(msa,msb) - f*InvDet
  InvM(msb,msa) = InvM(msb,msa) + f*InvDet

  return
end subroutine

