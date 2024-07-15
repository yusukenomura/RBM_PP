  ! 
  ! pre-define scalapack grid
  ! 
  if( mod(Nv,Nmpi) == 0 ) then 
    Nb = Nv / Nmpi 
  else 
    Nb = Nv / Nmpi + 1
  end if
  nprow = Nmpi
  npcol = 1 

  call sl_init(ictxt,nprow,npcol)
  call blacs_gridinfo(ictxt,nprow,npcol,myrow,mycol)
  if( myrow == -1 ) stop 'myrow == -1'
  if( mycol == -1 ) stop 'mycol == -1'

  mynumar = numroc(Nv,Nb,myrow,rsrc,nprow);
  mynumac = numroc(Nv,Nb,mycol,csrc,npcol);
  if( mynumac /= Nv ) stop 'mynumac /= Nv'

  mynumbr = mynumar
  mynumbc = numroc(N_rhs,Nb_rhs,mycol,csrc,npcol);
  if( mynumbc /= 1 ) stop 'mynumbc /= 1'

  numar(myrank) = mynumar
  numac(myrank) = mynumac
  numbr(myrank) = mynumbr
  numbc(myrank) = mynumbc

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  allocate( numtmp(0:Nmpi-1) ); numtmp = 0

  call MPI_ALLREDUCE(numar,numtmp,Nmpi,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  numar(:) = numtmp(:) 

  call MPI_ALLREDUCE(numac,numtmp,Nmpi,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  numac(:) = numtmp(:) 

  call MPI_ALLREDUCE(numbr,numtmp,Nmpi,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  numbr(:) = numtmp(:) 

  call MPI_ALLREDUCE(numbc,numtmp,Nmpi,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  numbc(:) = numtmp(:) 

  i1 = 0; i2 = 0 
  do i = 0, Nmpi-1
    if( numar(i) /= numbr(i) ) stop 'numar(i) /= numbr(i)'
    if( numac(i) /= Nv ) stop 'numac(i) /= Nv'
    if( numbc(i) /= 1   ) stop 'numbc(i) /= 1'
    if( numar(i) < 0 ) stop 'numar(i) < 0'
    if( numbr(i) < 0 ) stop 'numbr(i) < 0'
    i1 = i1 + numar(i)
    i2 = i2 + numbr(i)
  end do ! i  
  if( i1 /= Nv ) stop 'numar is wrong'
  if( i2 /= Nv ) stop 'numbr is wrong'

  maxnumar = maxval(numar)
  maxnumbr = maxval(numbr)
  if( maxnumar /= maxnumbr ) stop 'maxnumar /= maxnumbr'
  if( maxnumar /= Nb ) stop 'maxnumar /= Nb'

  deallocate( numtmp )
  call blacs_gridexit( ictxt )

