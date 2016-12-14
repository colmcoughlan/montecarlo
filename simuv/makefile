# Set compiler to g++.
link=g++-5
# Set options for the compiler
INC=-I/home/coughlan/astro/quickfits
LINKS=-L/home/coughlan/astro/quickfits -lcfitsio -lquickfits -lcfitsio -lgsl -fopenmp -lm -lblas
OPTS=-O3

all:
	${link} ${INC} -o simuv src/simuv.cpp ${LINKS} ${OPTS} 

clean:
	rm simuv
