# Set compiler to g++.
link=g++-5
# Set options for the compiler
INC=-I/Users/admin/git/quickfits
LINKS=-L/Users/admin/git/quickfits -lcfitsio -lquickfits -lgsl -fopenmp
OPTS=-O3

all:
	${link} ${INC} -o simuv src/simuv.cpp ${LINKS} ${OPTS} 

clean:
	rm simuv
