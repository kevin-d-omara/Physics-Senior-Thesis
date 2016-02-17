#CFLAGS   = -O2 -w
#CFLAGS	 = -O2  # -g -fbounds-check -fbacktrace
#CFLAGS = -g -fbounds-check -Wall -fbacktrace -finit-real=nan
CFLAGS = -O2 #   -openmp
LFLAGS   = #-openmp
DPLIBS   =
COMPILER77 = gfortran
COMPILER90 = gfortran
#COMPILER77 = ifort
#COMPILER90 = ifort
#LIBS = -L/Users/jtstaker/LIBRARIES/lapack-3.4.2 -llapack -L/Users/jtstaker/LIBRARIES/lapack-3.4.2/BLAS -lblas
#LIBS = -L/home/jtstaker/Libraries/lapack-3.4.1 -llapack -L/home/jtstaker/Libraries/BLAS -lblas
#LIBS = -L/opt/intel/mkl/lib/ia32 -lmkl_lapack95 -lmkl_blas95 -lpthread -L/opt/intel/lib/ia32 -liomp5 -L/home/jtstaker/Dvi ropbox/JohnsonResearch/lapack-3.4.0 -llapack -lblas
#LIBS = -L/opt/intel/mkl/9.1.021/lib/em64t  -lmkl_lapack  -lmkl_em64t -lmkl -lvml -lguide  -lpthread
#LIBS = -L/Users/joshuast/Dropbox/JohnsonResearch/lapack-3.4.0/SRC -lmkl_lapack -lmkl_ia32  -lmkl -lvml -lguide -lpthread
#LIBS = -L/Users/joshuast/Dropbox/JohnsonResearch/lapack-3.4.0/SRC -llapack -lblas
#LIBS = -L/home/jtstaker/Downloads/lapack-3.4.0 -llapack -lblas
#LIBS = -L/usr/local/lib -llapack -lrefblas

HF_OBJECTS =  \
phf_modules_full.o \
phfmaster.o \
phf_input.o  \
phf_denslib.o \
phf_tbme.o \
phf_applyh.o \
rotlib.o \
zcholeskylib.o \
phf_diag.o \
eiglib.o  \
libra.o  \
hfampFAST.o \
hfampPAIRFAST.o \
hfamp_menu.o \
gaujactest.o \
gauleg.o  \
phflib.o \
easyProjection.o \
zgesv.f 


%.o : %.f90
	$(COMPILER90) $(CFLAGS) -c $(<F)  

%.o : %.f
	$(COMPILER77) $(CFLAGS) -c $(<F)

#.f.o:;  (cd $(<D); $(COMPILER) -c  $(CFLAGS) $(<F))

ephf.x : $(HF_OBJECTS)
	$(COMPILER90) -o ephf.x $(HF_OBJECTS) $(LFLAGS) $(DPLIBS) $(LIBS)

ephf.nol.x : $(HF_OBJECTS)
	    $(COMPILER90) -o ephf.x $(HF_OBJECTS) $(LFLAGS) $(DPLIBS)
	
clean :
	rm -v *.o *.mod
