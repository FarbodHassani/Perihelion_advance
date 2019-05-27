#!/bin/bash
rm RK
g++ Newtonian_RK2_midpoint.cpp libnew.cpp -o rk2

./rk2