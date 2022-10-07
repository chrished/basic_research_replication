SYS_ARCH = $(shell sh -c 'uname -m 2>/dev/null')
SYS_OS = $(shell sh -c 'uname -s 2>/dev/null')

ifeq ($(SYS_ARCH),i686)
MEXSUFFIX = mexglx
MEXARCH = glnx86
endif
ifeq ($(SYS_ARCH),x86_64)
MEXSUFFIX = mexa64
MEXARCH = glnxa64
endif

MEX = mex
CC = gcc
CXX = g++
CXXOPTIMFLAGS = -O3
LDOPTIMFLAGS = -O3
MEXFLAGS = -cxx -v CC="$(CC)" CXX="$(CXX)" CXXOPTIMFLAGS="$(CXXOPTIMFLAGS)" LDOPTIMFLAGS="$(LDOPTIMFLAGS)"

all: eigsim.${MEXSUFFIX} firmsim.${MEXSUFFIX}

eigsim.${MEXSUFFIX}: eigsim.cpp
	$(MEX) $(MEXFLAGS) $(LIBS) eigsim.cpp -output eigsim.${MEXSUFFIX}

firmsim.${MEXSUFFIX}: firmsim.cpp
	$(MEX) $(MEXFLAGS) $(LIBS) firmsim.cpp -output firmsim.${MEXSUFFIX}

.SUFFIXES: .${MEXSUFFIX}

clean:
	rm -f eigsim.${MEXSUFFIX}
	rm -f firmsim.${MEXSUFFIX}

