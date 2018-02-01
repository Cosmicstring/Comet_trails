#!/bin/bash

#cp /home/andrija/Documents/Documents/MercuryPlay/My\ C++\ files/small.in /home/andrija/Documents/Documents/MercuryPlay/mercury6/
#cp /home/andrija/Documents/Documents/MercuryPlay/My\ C++\ files/big.in /home/andrija/Documents/Documents/MercuryPlay/mercury6/
gfortran -c mercury6_2.for -w
gfortran -o mercury6 mercury6_2.o -w
gfortran -c element6.for -w
gfortran -o element6 element6.o -w
gfortran -c close6.for -w
gfortran -o close6 close6.o -w
rm -f *.dmp 
rm -f *.aei 
rm -f *.tmp 
rm -f *.out 

