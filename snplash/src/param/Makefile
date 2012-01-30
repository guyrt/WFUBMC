include ../Makefile.inc

PARAMS = param_reader.o engine_param_reader.o

libparams.a : $(PARAMS)
	$(AR) libparams.a $^
	
param_reader.o : param_reader.cpp param_reader.h
	$(CMP) -o $@ $<

engine_param_reader.o : engine_param_reader.cpp engine_param_reader.h param_reader.o
	$(CMP) -o $@ $<
