#!/bin/bash

rm -f *.tpr \#* no_sc* sc_no_sigma* sc_sigma* traj* confout* ener*

#grompp -f params_em.mdp -c struct_single.gro -p top_run.top -o em.tpr -maxwarn 1
#mdrun -s em.tpr -c struct_single.gro

rm *.pdb

#grompp -f params.mdp -c struct_single.gro -p top.top -o run.tpr -maxwarn 2
#mdrun -s run.tpr -v
cd ../hc_ow
./prep.sh
cd ../hc_ow_charge
cp ../hc_ow/traj.trr .

rm \#*
