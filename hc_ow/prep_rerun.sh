#!/bin/bash

rm -f \#* no_sc* sc* ener*
rm -f *.xvg

grompp -f params_no_sc.mdp -c struct_single.gro -p top.top -maxwarn 2 -o no_sc.tpr
mdrun -deffnm no_sc -rerun traj.trr -reprod -nt 1 -v 

grompp -f params_sc.mdp -c struct_single.gro -p top.top -maxwarn 2 -o sc.tpr
mdrun -deffnm sc -rerun traj.trr -reprod -nt 1 -v 

grompp -f params_no_sc.mdp -c struct_single.gro -p top_charge.top -maxwarn 2 -o no_sc_charge.tpr
mdrun -deffnm no_sc_charge -rerun traj.trr -reprod -nt 1 -v 

grompp -f params_sc.mdp -c struct_single.gro -p top_charge.top -maxwarn 2 -o sc_charge.tpr
mdrun -deffnm sc_charge -rerun traj.trr -reprod -nt 1 -v

rm \#*

