#!/bin/bash

rm run.cpt run.gro run.trr run.edr run.log run.tpr run.xvg imd.gro mdout.mdp

gmx grompp -f run.mdp -c struct.gro -p topol.top -imd imd.gro -o run.tpr

gmx mdrun -nt 1 -deffnm run -imdport 8888 -imdwait
