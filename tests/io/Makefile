CXX = g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl
#CRYPTOLIB = -lcrypto


HTSSRC = ../../htslib
# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
CPPFLAGS += -I"$(realpath $(HTSSRC))"
LIBS := "$(realpath $(HTSSRC))/libhts.a" $(LIBS)
else
$(info HTSSRC not defined, assuming systemwide installation)
LIBS += -lhts
endif


NGSLIB = ../../include
CPPFLAGS += -I"$(realpath $(NGSLIB))"


#modied from htslib makefile
FLAGS = -O3 -fPIC
FLAGS2 = $(CPPFLAGS) $(FLAGS) $(LDFLAGS)
CXXFLAGS := $(FLAGS2) $(CXXFLAGS)


#for compiling with ZSTD which is used for .bgen file format
ifeq ($(WITH_ZSTD),1)
LIBS += -lzstd
CXXFLAGS += -D__ZSTD__
endif

NGSLIB_SRC = ../../src
CXXNGS = $(wildcard $(NGSLIB_SRC)/*.cpp) $(wildcard $(NGSLIB_SRC)/io/*.cpp)
CXXSRC = $(wildcard *.cpp)

OBJ = $(CXXSRC:.cpp=.o) $(CXXNGS:.cpp=.o)


PROGRAMS = testfasta

all: $(PROGRAMS)

#-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp


$(PROGRAMS): $(OBJ)
	$(CXX) $(FLAGS) -o testfasta *.o $(LIBS)

clean:
	rm -f *.o


