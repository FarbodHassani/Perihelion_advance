#!/bin/bash
rm it
gfortran -O6 -ffpe-summary=none -o it Schw2D_GR_discrete.f90 subs.f90 zeroin16.f 
./it
#now=$(date +"%Y.%m.%d.%H.%M.%S") &&
#mkdir ./Results/$now &&
#cp schw_perihelion.dat schw_txy.dat sample.h ./Results/$now 
