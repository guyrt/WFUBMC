include ../Makefile.inc

READER = reader.o linkage_reader.o arff_reader.o binary_bed_reader.o

libreader.a : $(READER)
	$(AR) libreader.a $^

reader.o : reader.cpp reader.h
	$(CMP) -o $@ $<

linkage_reader.o : linkage_reader.cpp linkage_reader.h reader.h ../param/param_reader.h
	$(CMP) -o $@ $<
	
binary_bed_reader.o : binary_bed_reader.cpp binary_bed_reader.h reader.h ../param/param_reader.h
	$(CMP) -o $@ $<

arff_reader.o : arff_reader.cpp arff_reader.h reader.h reader.cpp ../param/param_reader.h
	$(CMP) -o $@ $<
