#
# What the #$@^ did you do to this Makefile, Tommy!?!?
#
# Answer: don't hate.  It's undergoing some much needed changes.  Here's
# a cheat sheet:
#
#  The variable $@ is the thing being build (the thing before a :)
#  The variable $< is the first variable in the list of dependencies.
#  The variable $^ is a list of all dependencies.
#
#

include Makefile.inc

snplash : $(ALLO) $(SNPGWA) $(QSNPGWA) $(LD) $(DANDEL) $(INTERTWOLOG) $(ADTREE) $(MACHLEARN) $(UTILS) $(OUTPUT_HELPERS) param/libparams.a reader/libreader.a engine/libs/libmatrix.a engine/adtree/libadtree.a engine/bagging/libbagging.a
	$(LINK) -o $@ $^
	
param/libparams.a : force_look
	cd param; $(MAKE) 

reader/libreader.a : force_look
	cd reader; $(MAKE)
	
engine/libs/libmatrix.a : force_look
	cd engine/libs; $(MAKE)
	
engine/adtree/libadtree.a : force_look
	cd engine/adtree; $(MAKE)
	
engine/bagging/libbagging.a : force_look
	cd engine/bagging; $(MAKE)

engine/snp_data.o : engine/snp_data.cpp engine/snp_data.hh param/param_reader.h reader/linkage_reader.h
	$(CMP) -o $@ $<

engine/data_plugin.o : engine/data_plugin.cpp engine/snp_data.cpp engine/snp_data.hh param/param_reader.h reader/linkage_reader.h engine/data_plugin.h engine/data_plugin.cpp
	$(CMP) -o $@ $<

engine/randwh.o : engine/randwh.cpp engine/randwh.h
	$(CMP) -o $@ $<

engine/crossval/cross_val.o : engine/crossval/cross_val.cpp engine/crossval/cross_val.h engine/combinable.h $(ENGINE_INCLUDE)
	$(CMP) -o $@ $<

engine/ld/ld.o : engine/ld/ld.cpp engine/ld/ld.h $(ENGINE_INCLUDE) engine/em/em.h
	$(CMP) -o $@ $<

engine/dandelion/dandelion.o : engine/dandelion/dandelion.cpp engine/dandelion/dandelion.hh engine/em/em.h engine/output/dandelion_out.hh
	$(CMP) -o $@ $<

engine/em/em.o : engine/em/em.cpp engine/em/em.h engine/em/haplotype.h
	$(CMP) -o $@ $<

engine/em/haplotype.o :  engine/em/haplotype.cpp  engine/em/haplotype.h
	$(CMP) -o  $@ $<

engine/intertwolog/intertwolog.o : engine/intertwolog/intertwolog.cpp engine/intertwolog/intertwolog.hh engine/output/intertwolog_out.hh $(ENGINE_INCLUDE)
	$(CMP) -o $@ $<

engine/utils/lr.o : engine/utils/lr.cpp engine/utils/lr.hh engine/utils/vecops.hh
	$(CMP) -o $@ $<

engine/output/dandelion_out.o : engine/output/dandelion_out.cpp engine/output/dandelion_out.hh engine/utils/zaykin.hh $(ENGINE_INCLUDE)
	$(CMP) -o $@ $<

engine/output/dprime_out.o : engine/output/dprime_out.cpp engine/output/dprime_out.h $(ENGINE_INCLUDE)
	$(CMP) -o $@ $<

engine/output/intertwolog_out.o : engine/output/intertwolog_out.cpp engine/output/intertwolog_out.hh
	$(CMP) -o $@ $<

engine/output/qsnpgwa_out.o : engine/output/qsnpgwa_out.cpp engine/output/qsnpgwa_out.hh engine/output/output.h param/engine_param_reader.h engine/output/snpinfo_out.o
	$(CMP) -o $@ $<

engine/output/snpgwa_out.o : engine/output/snpgwa_out.cpp engine/output/snpgwa_out.hh engine/output/output.h param/engine_param_reader.h engine/output/snpinfo_out.o engine/snpgwa/genostats.o engine/snpgwa/haplostats.o engine/snpgwa/popstats.o
	$(CMP) -o $@ $<
	
engine/ouput/snpinfo_out.o : engine/ouput/snpinfo_out.cpp engine/ouput/snpinfo_out.hh
	$(CMP) -o $@ $<

engine/output/output.o : engine/output/output.cpp engine/output/output.h
	$(CMP) -o $@ $<

engine/qsnpgwa/qsnpgwa.o : engine/qsnpgwa/qsnpgwa.cpp engine/qsnpgwa/qsnpgwa.hh $(ENGINE_INCLUDE) engine/output/qsnpgwa_out.hh
	$(CMP) -o $@ $<

engine/qsnpgwa/cont_genostats.o : engine/qsnpgwa/cont_genostats.cpp engine/qsnpgwa/cont_genostats.hh engine/output/qsnpgwa_out.hh engine/utils/lr.hh
	$(CMP) -o $@ $<

engine/qsnpgwa/cont_popstats.o : engine/qsnpgwa/cont_popstats.cpp engine/qsnpgwa/cont_popstats.hh engine/output/qsnpgwa_out.hh
	$(CMP) -o $@ $<

engine/snpgwa/genostats.o : engine/snpgwa/genostats.cpp engine/snpgwa/genostats.hh engine/utils/lr.hh engine/output/snpgwa_out.hh logger/log.hh
	$(CMP) -o $@ $<

engine/snpgwa/haplostats.o : engine/snpgwa/haplostats.cpp engine/snpgwa/haplostats.hh  engine/output/snpgwa_out.hh engine/utils/zaykin.hh logger/log.hh engine/em/em.o
	$(CMP) -o $@ $<

engine/snpgwa/popstats.o: engine/snpgwa/popstats.cpp engine/snpgwa/popstats.hh engine/output/snpgwa_out.hh logger/log.hh
	$(CMP) -o $@ $<

engine/snpgwa/snpgwa.o : engine/snpgwa/snpgwa.cpp engine/snpgwa/snpgwa.h engine/output/snpgwa_out.hh $(ENGINE_INCLUDE) logger/log.hh engine/snpgwa/haplostats.o engine/snpgwa/genostats.o engine/snpgwa/popstats.o
	$(CMP) -o $@ $<

engine/utils/linear_regression.o : engine/utils/linear_regression.cpp engine/utils/linear_regression.hh engine/utils/exceptions.h
	$(CMP) -o $@ $<

engine/utils/statistics.o : engine/utils/statistics.cpp engine/utils/statistics.h engine/utils/exceptions.h
	$(CMP) -o $@ $<

engine/utils/stringutils.o : engine/utils/stringutils.cpp engine/utils/stringutils.h
	$(CMP) -o $@ $<
	
engine/utils/vecops.o : engine/utils/vecops.cpp engine/utils/vecops.hh
	$(CMP) -o $@ $<

engine/utils/zaykin.o : engine/utils/zaykin.cpp engine/utils/zaykin.hh engine/em/em.h engine/utils/vecops.hh
	$(CMP) -o $@ $<
	
engine/utils/float_ops.o : engine/utils/float_ops.cpp engine/utils/float_ops.hh
	$(CMP) -o $@ $<

logger/log.o : logger/log.cpp logger/log.hh
	$(CMP) -o $@ $<

snpadt.o : snpadt.cpp force_look
	$(CMP) $<



.PHONY: clean docs
clean :
	find . -name "*.o" | xargs /bin/rm -f
	find . -name "*.a" | xargs /bin/rm -f

docs :
	doxygen docs/doxy_config.txt
	
# Used to force entry into subdirs.
force_look : 
	true
