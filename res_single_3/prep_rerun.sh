#!/bin/bash

rm -f \#* no_sc* sc_no_sigma* sc_sigma* 

grompp -f params_no_sc.mdp -c struct_single.gro -p top.top -maxwarn 1 -o no_sc.tpr
mdrun -deffnm no_sc -rerun traj.xtc -reprod -nt 1 -v

#grompp -f params_sc_no_sigma.mdp -c struct_single.gro -p top.top -maxwarn 1 -o sc_no_sigma.tpr
#mdrun -deffnm sc_no_sigma -rerun traj.xtc -reprod -nt 1 -v

grompp -f params_sc_sigma.mdp -c struct_single.gro -p top.top -maxwarn 1 -o sc_sigma.tpr
mdrun -deffnm sc_sigma -rerun traj.xtc -reprod -nt 1 -v

rm \#*

