# Set compiler to g++.
link=g++-5
# Set options for the compiler
INC=-I/Users/admin/git/quickfits
LINKS=-L/Users/admin/git/quickfits -L/opt/local/lib/ -lcfitsio -lquickfits
OPTS=-O3

all:
	${link} ${INC} -o mcaverage src/mcaverage.cpp ${LINKS} ${OPTS} 

clean:
	rm mcaverage
