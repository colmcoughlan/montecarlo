# Set compiler to g++.
link=g++-5
# Set options for the compiler
INC=-I/Users/admin/git/quickfits
LINKS=-L/Users/admin/git/quickfits -lcfitsio -lquickfits -fopenmp -lfftw3 -lfftw3_threads
OPTS=-O3

all:
	${link} ${INC} -o conv src/conv.cpp ${LINKS} ${OPTS} 

clean:
	rm conv
