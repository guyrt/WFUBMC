# Version information

SVNDEV_T = $(shell svnversion -n .)
SVNDEV = -D'SVN_REV="$(SVNDEV_T)"'
SNPLASH_V = -D'SNPLASH_VERSION="1.1.$(SVNDEV_T)"'

SNPGWA_V = -D'SNPGWA_VERSION="4.3.$(SVNDEV_T)"'
QSNPGWA_V = -D'QSNPGWA_VERSION="3.1.$(SVNDEV_T)"'
I2L_V = -D'INTERTWOLOG_VERSION="1.1.$(SVNDEV_T)"'
DPRIME_V = -D'DPRIME_VERSION="4.0.$(SVNDEV_T)"'


VERSION_INFO = $(SVNDEV) $(SNPLASH_V) 

#CMP = icpc -O3 -ipo -no-prec-div -static -openmp -c -g -Wall -wd981 -wd383 $(VERSION_INFO)
#LINK = icpc -O3 -ipo -no-prec-div -static -openmp -g -Wall -wd981 -wd383
#AR = xiar rvs

#CMP = g++ -O3 -fopenmp -g -c -static -Wall -Wunused $(VERSION_INFO)
#LINK = g++ -O3 -fopenmp -g -Wall -static -Wunused

CMP = g++ -O3 -fopenmp -g -c -Wall -Wunused $(VERSION_INFO)
LINK = g++ -O3 -fopenmp -g -Wall -Wunused
AR = ar rvs 

MAKE = make

ALLO = snpadt.o engine/data_plugin.o engine/randwh.o engine/snp_data.o engine/utils/statistics.o engine/em/em.o engine/em/haplotype.o engine/output/output.o  
UTILS = engine/utils/stringutils.o engine/utils/lr.o engine/utils/linear_regression.o logger/log.o engine/utils/vecops.o engine/utils/float_ops.o
SNPGWA = engine/snpgwa/genostats.o engine/snpgwa/snpgwa.o engine/snpgwa/popstats.o engine/output/snpgwa_out.o engine/snpgwa/haplostats.o engine/utils/zaykin.o
QSNPGWA = engine/qsnpgwa/qsnpgwa.o engine/qsnpgwa/cont_popstats.o engine/output/qsnpgwa_out.o engine/qsnpgwa/cont_genostats.o

OUTPUT_HELPERS = engine/output/snpinfo_out.o

DANDEL = engine/dandelion/dandelion.o engine/output/dandelion_out.o
LD = engine/ld/ld.o engine/output/dprime_out.o
INTERTWOLOG = engine/intertwolog/intertwolog.o engine/output/intertwolog_out.o
MACHLEARN = engine/crossval/cross_val.o

ENGINE_INCLUDE = engine/randwh.h engine/snp_data.hh engine/data_plugin.h param/engine_param_reader.h param/param_reader.h engine/engine.h

