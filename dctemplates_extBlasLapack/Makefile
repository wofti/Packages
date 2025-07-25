# vars used also for F2CLIBS/libF77.a  and  F2CLIBS/libI77.a
SHELL = /bin/sh
CC_F2CLIBS = gcc
CFLAGS_F2CLIBS = -O3 -DNON_UNIX_STDIO
export

# vars used only to build libiteratortemplates.a and xtest
CC = gcc
CFLAGS = -O3

# Fortran compat libs
f2cLIBS = F2CLIBS/libF77.a F2CLIBS/libI77.a

templates_o = CG.o BiCG.o CGS.o BiCGSTAB.o GMRES.o QMR.o Cheby.o  \
            Jacobi.o SOR.o WT_utils.o
#            dlapack.o dblas.o 

xtest_o = Tester.o MatGen.o dspdchk.o dnsychk.o MatVec.o PSolve.o \
        Utils.o templates_dlaran.o

dl_lapack_blas_o = dlapack.o dblas.o


# compile
.f.o:;  $(CC) $(CFLAGS) -c $<

# make libiteratortemplates.a
libiteratortemplates.a: $(templates_o)
	ar rucs libiteratortemplates.a $(templates_o)
	@echo "to make libdl_lapack_blas.a type: make libdl_lapack_blas.a"
	@echo "target libdl_lapack_blas.a also makes libdlblas.a and libdllapack.a"
	@echo "for testing with libdl_lapack_blas.a type:  make xtest ; ./xtest"
	@echo "for testing with external BLAS + LAPACK type:  make xtest_ext ; ./xtest_ext"

# make libF77.a and libI77.a 
F2CLIBS/libF77.a:
	make -C F2CLIBS/libF77

F2CLIBS/libI77.a:
	make -C F2CLIBS/libI77


# make lib containing double and long int version of lapack and blas
libdl_lapack_blas.a: $(dl_lapack_blas_o)
	ar rucs libdl_lapack_blas.a $(dl_lapack_blas_o)
	ar rucs libdlblas.a   dblas.o
	ar rucs libdllapack.a dlapack.o

# make the test prog
xtest:   $(xtest_o) libdl_lapack_blas.a libiteratortemplates.a $(f2cLIBS)
	$(CC) $(CFLAGS) -o $@ $(xtest_o) libdl_lapack_blas.a libiteratortemplates.a $(f2cLIBS) -lm

xtest_ext: $(xtest_o) libiteratortemplates.a $(f2cLIBS)
	$(CC) $(CFLAGS) -o $@ $(xtest_o) libiteratortemplates.a $(f2cLIBS) -lm -lblas -llapack

# clean up
clean:
	rm -f libiteratortemplates.a *.o core a.out xtest* test.results *.trace $(f2cLIBS)
	rm -f F2CLIBS/libF77/*.o F2CLIBS/libI77/*.o 
	rm -f libdlblas.a
	rm -f libdllapack.a
	rm -f libdl_lapack_blas.a
