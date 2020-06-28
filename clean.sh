#!/bin/sh

## Delete python temporary files
rm -rf python/CDSK.egg*
rm -rf python/build
rm -rf python/dist
rm -rf python/tmp
rm -rf python/var


## Delete R temporary files
#rm -f R/CDSK/NAMESPACE
rm -f R/CDSK/man/*.Rd
rm -f R/*.tar.gz
rm -f R/CDSK/src/RcppExports.cpp
rm -f R/CDSK/R/RcppExports.R
rm -f R/CDSK/src/*.o
rm -f R/CDSK/src/*.so

