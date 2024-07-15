module mod_main
  real(8) :: J1, J2 
  real(8) :: wf_av
  real(8) :: shift
  real(8) :: shift_av
  integer :: Nupdate                           ! number of update between measurement
  integer :: spin_parity                       ! = +1 (S=even), -1 (S=odd)
  integer :: L                                 ! system size
  integer :: N                                 ! number of visible units
  integer :: Nchi                              ! number of elements for correlation function
  integer :: Nv                                ! number of variational parameters 
  integer :: Nv_RBM                            ! number of variational parameters for RBM part 
  integer :: N_MP_Trans                        ! number of quantum projection for translation and point group symmetry 
  integer :: N_QP_Full                         ! total number of quantum projection = N_SP_Gauss_Leg*N_MP_Trans
  logical , allocatable :: subB_or_not(:)      ! whether the site belong to sublattice B or not  
  complex(8) , allocatable :: para_mp_trans(:) ! factor for momentum projection 
end module 


module mod_RBM
  integer :: Mc                             ! number of independent theta angles 
  integer :: Mr                             ! number of independent theta angles 
  integer :: alphac                         ! number of independent hidden unit
  integer :: alphar                         ! number of complex independent hidden unit
  integer , parameter :: Ncopy = 1          ! number of copies of hidden unit
  complex(8) , allocatable :: Wcirr(:,:)    ! irreducible part of RBM interaction parameters
  complex(8) , allocatable :: bcirr(:)      ! irreducible part of RBM bias parameters
  complex(8) , allocatable :: Wcirr_av(:,:) ! average of Wirr for last several iterations
  complex(8) , allocatable :: bcirr_av(:)   ! average of birr for last several iterations 
  real(8) , allocatable :: Wrirr(:,:)       ! irreducible part of RBM interaction parameters
  real(8) , allocatable :: brirr(:)         ! irreducible part of RBM bias parameters
  real(8) , allocatable :: Wrirr_av(:,:)    ! average of Wirr for last several iterations
  real(8) , allocatable :: brirr_av(:)      ! average of birr for last several iterations 

  integer , allocatable :: W_map(:,:,:)     ! map for W interactions
end module

module mod_sample
  integer :: Nsample
  integer , allocatable :: x_sample(:,:)
  complex(8) , allocatable :: psi_x_sample(:,:) ! 
  real(8) , allocatable :: PfM_sample(:,:)
end module



module mod_PP 
  logical :: l_AP              ! if .true., we use anti-periodic boundary condition
  integer :: Npair             ! number of pair = number of electrons with up or down spin 
  integer :: N_SP_Gauss_Leg    ! number of points for the Gauss-Legendre quadrature 
  integer :: N_Slater          ! number of parameters in slater part
  integer :: Nfastup           ! number of fast update of Pfaffian 
  integer :: lwork_pfa         ! lwork in pfapack
  real(8) , allocatable :: Slater(:)                  ! variatonal parameters in slater part
  real(8) , allocatable :: Slater_av(:)               ! average of fij for last several iterations 
  real(8) , allocatable :: Slater_Elm(:,:,:)          ! entire slater matrix 
  real(8) , allocatable :: SP_GL_Cos(:), SP_GL_Sin(:) ! parameter in quantum projection [cos(beta/2) and sin(beta/2)]
  real(8) , allocatable :: SP_GL_CosSin(:)            ! parameter in quantum projection
  real(8) , allocatable :: SP_GL_CosCos(:)            ! parameter in quantum projection
  real(8) , allocatable :: SP_GL_SinSin(:)            ! parameter in qunatum projection 
  integer , allocatable :: Orbital_Idx(:,:)           ! index for orbital 
  integer , allocatable :: Orbital_Sgn(:,:)           ! sign for fij (used in anti-periodic case)
  integer , allocatable :: MP_Trans_Idx(:,:)          ! index for translational quantum projection
  integer , allocatable :: MP_Trans_Sgn(:,:)          ! sign for translational operation (used in anti-periodic case)
  real(8) , allocatable :: para_spin_proj(:)          ! factor for spin projection 

  real(8) :: slater_abs_av = 0.25d0
end module



module mod_pdposv
  integer , parameter :: N_rhs = 1   ! b has 1 column 
  integer , parameter :: Nb_rhs = 1  ! basic block size for column of b   
  integer , parameter :: ia = 1      ! row index in global(A) indicating the first row of sub(A)
  integer , parameter :: ja = 1      ! column index in global(A) indicating the first coulumn of sub(A)
  integer , parameter :: ib = 1      ! row index in global(B) indicating the first row of sub(B)
  integer , parameter :: jb = 1      ! column index in global(B) indicating the first coulumn of sub(B)
  integer , parameter :: rsrc = 0    ! process row over which the first row of matrix is distributed
  integer , parameter :: csrc = 0    ! process column over which the first column of matrix is distributed 
  integer , external  :: numroc      ! NUMber of Row Or Columns
  integer :: nprow, npcol            ! number of processes for row and column
  integer :: Nb                      ! basic block size

  integer :: mylda, myldb            ! local leading dimension of a and b
  integer :: myrow, mycol            ! index for row and column processes
  integer :: desca(9), descb(9)      ! descripter used in pdposv
  integer :: ictxt                   ! BLACS context
  integer :: mynumar, mynumac        ! row and column size of local piece of sub(A)  
  integer :: mynumbr, mynumbc        ! row and column size of local piece of sub(B) 
  integer , allocatable :: numar(:)  ! row sizes of pieces of sub(A)
  integer , allocatable :: numac(:)  ! column sizes of pieces of sub(A) 
  integer , allocatable :: numbr(:)  ! row sizes of pieces sub(B)
  integer , allocatable :: numbc(:)  ! column sizes of pieces of sub(B)
  integer , allocatable :: numtmp(:) ! for communication 
  integer :: maxnumar  
  integer :: maxnumbr 
  real(8) , allocatable :: amat(:,:), bmat(:,:)

  integer :: info

!!! Variables used in SR   
  real(8) , allocatable :: s_save(:,:)
  real(8) , allocatable :: gchk(:)
  integer :: nSmat
  integer :: kmax, kmin 
  integer :: offset
  real(8) :: Smax, Smin, Sdiag_av, Sdiag_av_RBM, Sdiag_av_PP
  real(8) :: Smat_scaling_factor
  real(8) :: lambda
  real(8) :: lambda_RBM_RBM, lambda_PP_PP, lambda_RBM_PP
  real(8) :: rcut

real(8) , allocatable :: s_save_gomi(:,:)
real(8) :: lam_chk
end module
