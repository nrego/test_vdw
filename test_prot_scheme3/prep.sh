#!/bin/bash

rm -f *.tpr \#* no_sc* sc_no_sigma* sc_sigma* traj* confout* ener*

#grompp -f params_em.mdp -c struct_single.gro -p top_run.top -o em.tpr -maxwarn 1
#mdrun -s em.tpr -c struct_single.gro

rm *.pdb
cp ../test_prot_scheme2/traj.xtc .


rm -f \#*
