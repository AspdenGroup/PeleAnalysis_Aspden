#!/bin/bash

set -x

# build the surface decimation tool
# reduces a surface millions of points down to a specified number in a sensible way
# this bit might be a bit tricky; you may find some library dependencies...

export TOP=${PWD}

# we're going to follow the instructions in here
#cat ${TOP}/qslim-2.1/README.txt

cd ${TOP}/qslim-2.1/libgfx

./configure

cd src; make -j 8

cd ../../mixkit

./configure

cd src; make -j 8

# doesn't work, but don't need it to
# cd ../../tools/qslim
# make all

# this is what we really need
cd ${TOP}/MyQSlim

/bin/cp ../GNUmakefileMyQSlim GNUmakefile

make -j 8

