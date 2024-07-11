# Copyright (c) 2018 Yusuke Nomura

#======= please edit this part ========
FC = mpiifort -O3 -qopenmp 
#FC = mpif90 -qopenmp -debug full -check all -warn all -traceback -g 
FFLAGS =
LIBBLAS = -mkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
LIBPFA  = libpfapack.a
#======================================

.SUFFIXES: .f .f90
.f.o:
	${FC} ${FFLAGS} -c $<

.f90.o:
	${FC} ${FFLAGS} -c $<

TARGET = RBM_solver_SR.x
OBJECTS = mt.o main.o calc_phys_quantity.o make_sample.o sub.o RBM.o SR.o interface_for_PP.o matrix.o pfupdate_two.o spin_projection.o slater.o parameter.o module.o 
MOD_FILES = mod_main.mod mod_rbm.mod mod_pp.mod mod_pdposv.mod


${TARGET} : ${OBJECTS} 
	${FC} ${FFLAGS} -o $@ ${OBJECTS} ${LIBBLAS} ${LIBPFA} 

main.o: main.f90 module.o
sub.o: sub.f90 module.o 
RBM.o: RBM.f90 module.o 
make_sample.o: make_sample.f90 module.o
calc_phys_quantity.o: calc_phys_quantity.f90 module.o 
interface_for_PP.o: interface_for_PP.f90 module.o
matrix.o: matrix.f90 module.o
pfupdate_two.o: pfupdate_two.f90 module.o 
spin_projection.o: spin_projection.f90 module.o 
parameter.o: parameter.f90 module.o 
slater.o: slater.f90 module.o 


clean: 
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES} *genmod* 
