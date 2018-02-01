#!/bin/bash

g++ -o Particle_Orbit.o -c Particle_Orbit_utils.cpp
g++ -o Particle_Orbit Particle_Orbit.cpp Particle_Orbit.o

