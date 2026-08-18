# To build against HTSlib system libraries
#
#   make HTSLIBINCDIR=/usr/include HTSLIBLIBDIR=/usr/lib all
#
PWD = $(shell pwd)
HTSLIBINCDIR ?= ${PWD}/src/htslib/
HTSLIBLIBDIR ?= ${PWD}/src/htslib/
BOOSTINCDIR ?=
BOOSTLIBDIR ?=

# Debug or static build
DEBUG ?= 0
STATIC ?= 0

# Submodules
SDSL_ROOT ?= ${PWD}/src/sdslLite/
JLIB ?= ${PWD}/src/jlib/

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX ?= g++
CXXFLAGS += -std=c++17 -isystem ${JLIB} -isystem ${HTSLIBINCDIR} -isystem ${SDSL_ROOT}/include -pedantic -W -Wall
LDFLAGS += -L${HTSLIBLIBDIR} -L${SDSL_ROOT}/lib -lboost_iostreams -lboost_filesystem -lboost_program_options -lboost_date_time -ldl -lpthread

# Boost location
ifneq (${BOOSTINCDIR},)
	CXXFLAGS += -isystem ${BOOSTINCDIR}
endif
ifneq (${BOOSTLIBDIR},)
	LDFLAGS += -L${BOOSTLIBDIR} -Wl,-rpath,${BOOSTLIBDIR}
endif

# Flags for static compile
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -lhts -lz -llzma -lbz2 -ldeflate
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${HTSLIBLIBDIR}
endif

# Flags for debugging, profiling and releases
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${HTSLIBINCDIR}, ${PWD}/src/htslib/)
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
	if [ -r src/xxsds/install.sh ]; then cd src/xxsds/build && ./clean.sh && cmake -DCMAKE_INSTALL_PREFIX=${PBASE}/src/sdslLite -DSDSL_BUILD_TESTS=OFF -DSDSL_BUILD_EXAMPLES=OFF -DSDSL_BUILD_TUTORIAL=OFF .. && make -j 8 sdsl && make install && cd ../../../ && touch .sdsl; fi

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoreconf -i && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

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
