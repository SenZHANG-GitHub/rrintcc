
## -O3: optimization (-On, n=0,1,2,3...)
## -O3 might be problematic for scientific calculation!! Use -O2 instead for stable version
## -I. : specify current folder as the .h files folder
## -Wall print warning messages; -w : turn off all warnings
## -static: On systems that support dynamic linking, this prevents linking with the shared libraries.
## -lx:	link to library libx.so 

SYS = UNIX

#CXX_WIN = x86_64-w64-mingw32-c++.exe
#CXX_UNIX = g++

# Never enter this, this RInside-based program is only for Linux
ifeq ($(SYS),WIN)
 CXX = $(CXX_WIN)
endif

ifeq ($(SYS),UNIX)
 CXX = $(shell $(R_HOME)/bin/R CMD config CXX)
endif

##########################################################################
## Settings for RInside

## comment this out if you need a different version of R, 
## and set set R_HOME accordingly as an environment variable
R_HOME := 		$(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 		$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := 		$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

##########################################################################

OUTPUT = rrintcc_BOOST

CXX_FLAGS += -w -I. -m64 --std=c++11 $(shell $(R_HOME)/bin/R CMD config CPPFLAGS) $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) $(shell $(R_HOME)/bin/R CMD config CXXFLAGS)
LIB += $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

SRC = rrintcc_BOOST.cpp utility.cpp stats.cpp 
HDR = utility.h stats.h
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT)

$(OUTPUT) :
	$(CXX) $(CXX_FLAGS) -o $(OUTPUT) $(OBJ) $(LIB)

$(OBJ) : $(HDR)

.cpp.o :
	$(CXX) $(CXX_FLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	rm -f *.o *~
