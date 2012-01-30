include ../Makefile.inc

LR_O = ../engine/utils/lr.o ../engine/utils/float_ops.o ../engine/utils/statistics.o ../engine/utils/vecops.o

ENGINE_PARAM_O = ../param/engine_param_reader.o ../param/param_reader.o
ENGINE_PARAM_H = $(subst .o,.h,${ENGINE_PARAM_O})

test : make_libs test_lr.o 
	$(LINK) -lboost_unit_test_framework-mt  -o test_snpadt test_lr.o $(LR_O) ../engine/libs/libmatrix.a 
	
make_libs : 
#	cd ..; $(MAKE)

test_lr.o : engine/utils/test_lr_engine.cpp ../engine/utils/lr.hh ../engine/utils/lr.cpp
	$(CMP) -o $@ $< -lboost_unit_test_framework-mt 
