#Set inclusion paths here (if boost, bamtools, or args are installed outside your path)
INCLUDE_DIRS=-I/usr/local/include/bamtools
#Set library paths here (if boost or bamtools are installed outside your path)
LIBRARY_PATHS=
#Set to 0 if you encounter linker errors regarding strings from the bamtools library
ABI=1
#Provide full paths here to .a archives for libraries which should be statically linked
STATIC_LIBS=
#List of remaining libraries that will be dynamically linked
LIBS=-lbamtools -lboost_filesystem -lboost_regex -lboost_system -lz

CC=g++
STDLIB=-std=c++14
CFLAGS=-Wall -I. $(STDLIB) -D_GLIBCXX_USE_CXX11_ABI=$(ABI) -O3

rnaseqc: BED.cpp Expression.cpp GTF.cpp RNASeQC.cpp Metrics.cpp
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $(LIBRARY_PATHS) -o $@ $^ $(STATIC_LIBS) $(LIBS)

.PHONY: clean

clean:
	rm rnaseqc
