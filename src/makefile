all: all-redirect

BUILD_DIR = /opt/teems-solver
PETSC_DIR = $(BUILD_DIR)/lib/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CC          = mpicc
MPICH_DIR    = $(BUILD_DIR)/lib/mpi
FC          = $(MPICH_DIR)/bin/mpif90 -fopenmp
F77          = $(MPICH_DIR)/bin/mpif77 -fopenmp
F77FLAGS   = -I. -I/usr/local/include
CFLAGS    = -I. -I${PETSC_DIR}/native/include -I ${PETSC_DIR}/include -fopenmp
FCFLAGS     =
LDFLAGS     = -fopenmp

LIBS        =  ${PETSC_KSP_LIB} $(BUILD_DIR)/lib/ma48/src/.libs/libma48.so $(BUILD_DIR)/lib/ma51/src/.libs/libma51.so

.PHONY : clean
	-rm -f *.o 
	-rm -f hsl

all-redirect: hsl
OBJECTS = main.o ha_cgeiof.o ha_cgefparse.o ha_cgetab.o ha_newfparse.o ha_newmfparse.o ha_mp48.o hsl_mp48d.o ddeps.o hsl_mc66d.o ddeps90.o 
hsl: $(OBJECTS)
	${CLINKER} -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -Ofast -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -Ofast -c $<
%.o: %.f
	$(F77) $(F77FLAG) -Ofast -c $<

# dependencies
hsl_mp48d.o: ddeps.o
hsl_mc66d.o: ddeps90.o
ha_mp48.o: hsl_mc66d.o hsl_mp48d.o

