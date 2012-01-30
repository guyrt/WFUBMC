include ../../Makefile.inc

BAGGING = bagging.o

#Root of project.
ROOT = ../..

LOCAL_ENGINE_INCLUDE = $(foreach inc,${ENGINE_INCLUDE},${ROOT}/${inc})

libbagging.a : $(BAGGING)
	$(AR) libbagging.a $^
	
bagging.o : bagging.cpp bagging.h $(ROOT)/engine/combinable.h $(LOCAL_ENGINE_INCLUDE)
	$(CMP) -o $@ $<
