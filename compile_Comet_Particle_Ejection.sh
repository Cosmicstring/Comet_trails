#!/bin/bash

g++ -o Comet_Particle_Ejection.o -c Comet_Particle_Ejection_utils.cpp
g++ -o Comet_Particle_Ejection Comet_Particle_Ejection.cpp Comet_Particle_Ejection.o
