#!/bin/bash

g++ -o Generate_SOHO_Image.o -c Generate_SOHO_Image_utils.cpp
g++ -o Generate_SOHO_Image Generate_SOHO_Image.cpp Generate_SOHO_Image.o

