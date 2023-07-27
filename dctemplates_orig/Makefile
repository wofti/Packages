# vars used also for F2CLIBS/libF77.a  and  F2CLIBS/libI77.a
SHELL = /bin/sh
CC_F2CLIBS = icc
CFLAGS_F2CLIBS = -O2 -DNON_UNIX_STDIO
export

# vars used only to build libiteratortemplates.a and xtest
CC = icc
CFLAGS = -O2

f2cLIBS = F2CLIBS/libF77.a F2CLIBS/libI77.a

templates_o = Tester.o MatGen.o dspdchk.o dnsychk.o \
            CG.o BiCG.o CGS.o BiCGSTAB.o GMRES.o QMR.o Cheby.o  \
            Jacobi.o SOR.o Utils.o \
            MatVec.o PSolve.o \
            dlapack.o dblas.o 
xtest_o = Tester.o MatVec.o


# compile
.f.o:;  $(CC) $(CFLAGS) -c $<

# make libiteratortemplates.a
libiteratortemplates.a: $(templates_o) $(f2cLIBS)
	ar rucs libiteratortemplates.a $(templates_o)
	echo "for testing type:  make xtest ; ./xtest"

# make libF77.a and libI77.a 
F2CLIBS/libF77.a:
	make -C F2CLIBS/libF77

F2CLIBS/libI77.a:
	make -C F2CLIBS/libI77
	

# make the test prog
xtest:   $(xtest_o) libiteratortemplates.a $(f2cLIBS)
	$(CC) $(CFLAGS) -o $@ $(xtest_o) libiteratortemplates.a $(f2cLIBS) -lm

# clean up
clean:
	rm -f libiteratortemplates.a *.o core a.out xtest* test.results *.trace $(f2cLIBS)
	rm -f F2CLIBS/libF77/*.o F2CLIBS/libI77/*.o 
