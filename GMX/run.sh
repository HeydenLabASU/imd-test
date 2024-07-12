#!/bin/bash

rm run-NPT.cpt run-NPT.gro run-NPT.trr run-NPT.edr run-NPT.log run-NPT.tpr run-NPT.xvg imd.gro mdout.mdp

gmx grompp -f run-NPT.mdp -c equi-NPT+posres.gro -p topol.top -imd imd.gro -o run-NPT.tpr

gmx mdrun -nt 1 -deffnm run-NPT -imdport 8888 -imdwait
