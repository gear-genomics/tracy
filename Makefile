DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
SDSL_ROOT ?= ${PWD}/src/sdslLite/
EBROOTHTSLIB ?= ${PWD}/src/htslib/
JLIB ?= ${PWD}/src/jlib/

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX=g++
CXXFLAGS += -std=c++14 -isystem ${JLIB} -isystem ${EBROOTHTSLIB} -isystem ${SDSL_ROOT}/include -pedantic -W -Wall
LDFLAGS += -L${EBROOTHTSLIB} -L${EBROOTHTSLIB}/lib -L${SDSL_ROOT}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time

ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif
ifeq (${SDSL_ROOT}, ${PWD}/src/sdslLite/)
	SUBMODULES += .sdsl
endif


# External sources
SDSLSOURCES = $(wildcard src/xxsds/lib/*.cpp)
SOURCES = $(wildcard src/*.cpp) $(wildcard src/*.h)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
PBASE=$(shell pwd)

# Targets
BUILT_PROGRAMS = src/tracy
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.sdsl: $(SDSLSOURCES)
	if [ -r src/xxsds/install.sh ]; then cd src/xxsds/ && ./install.sh ${PBASE}/src/sdslLite && cd ../../ && touch .sdsl; fi

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoheader && autoconf && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

src/tracy: ${SUBMODULES} ${SOURCES}
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	if [ -r src/xxsds/install.sh ]; then rm -rf src/sdslLite/; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
