include ../../Makefile.inc

LIBO = ap.o linalg.o alglibinternal.o alglibmisc.o blas.o ialglib.o

libmatrix.a : $(LIBO)
	$(AR) libmatrix.a $^
	
ap.o : ap.cpp ap.h
	$(CMP) -o $@ $<

linalg.o : linalg.cpp linalg.h
	$(CMP) -o $@ $<

alglibinternal.o : alglibinternal.cpp alglibinternal.h
	$(CMP) -o $@ $<

alglibmisc.o : alglibmisc.cpp alglibmisc.h
	$(CMP) -o $@ $<
	
blas.o : blas.cpp blas.h
	$(CMP) -o $@ $<

ialglib.o : ialglib.cpp ialglib.h
	$(CMP) -o $@ $<
