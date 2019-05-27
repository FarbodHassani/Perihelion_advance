#!/bin/bash
rm RK
g++ Newtonian_leap_frog.cpp libnew.cpp -o RK

./RK