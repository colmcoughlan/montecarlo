# Usage

To use just call the executable: ./conv.

This software convolves a radio FITS map with a specified elliptical beam

# Installation instructions

Dependencies: a C compiler, cfitsio and quickfits (https://github.com/colmcoughlan/quickfits).

On Ubuntu cfitsio can be found in the cfitsio-dev package. On a mac the best way to install it might be through homebrew (brew install cfitsio).
More mac instructions can be found at https://heasarc.gsfc.nasa.gov/fitsio/fitsio_macosx.html. You can read about homebrew at http://brew.sh/.

If you're using linux you probably already have a c compiler (gcc). On a mac you can install this through homebrew with "brew install gcc5", to install gcc version 5.

Once you have the dependencies installed, make sure they are on your library path. On linux this involves adding the directories to $LD_LIBRARY_PATH.

To build, edit the makefiles such that:

  1. link=\<your C compiler>
  2. INC=-I/<path to quickfits.h directory>
  3. LINKS=-L/<path to quickfits.a directory> -lcfitsio -lquickfits -fopenmp

Note this assumes that cfitsio etc. have been installed onto your path already. If this is not the case, you should specify additional -I and 
-L paths to their location.

Finally, run make -f makefile to generate the executable.
