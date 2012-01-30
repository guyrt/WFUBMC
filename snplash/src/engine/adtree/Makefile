include ../../Makefile.inc

ADTREE = adtree.o condition.o ad_rule.o ad_tree_data.o precondition.o

#Root of project.
ROOT = ../..

LOCAL_ENGINE_INCLUDE = $(foreach inc,${ENGINE_INCLUDE},${ROOT}/${inc})

libadtree.a : $(ADTREE)
	$(AR) libadtree.a $^
	
adtree.o : adtree.cpp adtree.h $(ROOT)/engine/combinable.h $(LOCAL_ENGINE_INCLUDE)
	$(CMP) -o $@ $<

condition.o : condition.cpp condition.h 
	$(CMP) -o $@ $<
	
ad_rule.o : ad_rule.cpp ad_rule.h condition.h condition.cpp
	$(CMP) -o $@ $<

precondition.o : precondition.cpp precondition.h condition.h condition.cpp
	$(CMP) -o $@ $<

ad_tree_data.o : ad_tree_data.cpp ad_tree_data.h
	$(CMP) -o $@ $<

