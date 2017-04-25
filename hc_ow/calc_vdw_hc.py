## Do vdw calculation from md trajectory##

from __future__ import division, print_function

import MDAnalysis
import numpy as np

import argparse
try:
    u_for_prev = u_for
except:
    pass
univ = MDAnalysis.Universe('no_sc.tpr', 'traj.xtc')
from mdtools import dr

alc_indices = np.arange(1)
atm_indices = np.arange(univ.atoms.n_atoms)

excls = {
         0:(0)
}

# dictionary of alchemical atom type transformations.
#    keyed by alchemical index of an atom that is to be transformed
#    valued by tuple (Astate_idx, Bstate_idx)
#  NOTE: This is topology specific!!!
alc_types = {
    0: ('DUM_HC', 'HC')
}

type_lookup = {
    'HC': 0,
    'DUM_HC': 1,
    'OW_spc': 2,
    'HW_spc': 1
}
#format: 24 atom types, atomtype i is a tuple with (sigma, eps)
atmtypes = [
(2.64953e-01,  6.56888e-02), # HC
(0.00000e+00,  0.00000e+00), # DUM_HC, HW_spc
(3.16557e-01,  6.50629e-01)  # OW_spc
]


# set combined sigma to this if it's smaller
sc_sigma = 0.3
sc_sigma6 = sc_sigma**6
sc_alpha = 0.5

n_atmtype = len(atmtypes)
# Generate VdW lookup table
c6_lut = np.zeros(n_atmtype**2)
c12_lut = np.zeros(n_atmtype**2)
sig_lut = np.zeros(n_atmtype**2)
sig6_lut = np.zeros(n_atmtype**2)
for i, payload_i in enumerate(atmtypes):
    sig_i, eps_i = payload_i
    for j, payload_j in enumerate(atmtypes):
        sig_j, eps_j = payload_j
        eps = np.sqrt(eps_i*eps_j)

        idx = i*n_atmtype + j

        sig = (sig_i + sig_j) / 2.0

        c6 = 4*eps*sig**6
        c12 = 4*eps*sig**12
        c6_lut[idx] = c6
        c12_lut[idx] = c12


        sig_lut[idx] = sig
        sig_6 = c12/c6
        if c6 == 0 or c12 == 0 or sig < sc_sigma:
            sig_6 = sc_sigma**6
        sig6_lut[idx] = sig_6

lmbda = 0.0
for_lmbdas = [0.5, 1.0]

n_frames = univ.trajectory.n_frames
my_diffs = np.zeros((len(for_lmbdas), n_frames, 2))

for window_idx, lmbda_for in enumerate(for_lmbdas):
    

    for i_frame in range(n_frames):
        univ.trajectory[i_frame]
        my_diffs[window_idx, i_frame, 0] = univ.trajectory.time
        univ.atoms.positions = univ.atoms.positions / 10.0
        # Calculate VdW energy differences between lambdas
        u_lmbda = 0.0
        u_for = 0.0
        for i in alc_indices:
            # all atoms separated by more than nrexcl bonds (i.e. not excluded)
            # Note: If i and j are both in alc_indices, skip if j !> i
            incl_indices = np.setdiff1d(atm_indices, excls[i])

            atm_i = univ.atoms[i]
            # from tpr file, should be A state topology 
            type_i = type_lookup[atm_i.type]
            name_i_a, name_i_b = alc_types[i]
            type_i_a = type_lookup[name_i_a]
            type_i_b = type_lookup[name_i_b]

            assert type_i == type_i_a

            for j in incl_indices:

                atm_j = univ.atoms[j]
                if j in alc_indices:
                    if j < i:
                        continue
                    name_j_a, name_j_b = alc_types[j]
                    type_j_a = type_lookup[name_j_a]
                    type_j_b = type_lookup[name_j_a]
                else:
                    type_j_a = type_j_b = type_lookup[atm_j.type]        

                lut_idx_a = type_i_a * n_atmtype + type_j_a
                lut_idx_b = type_i_b * n_atmtype + type_j_b

                r_ij_sq = np.sum((atm_i.position - atm_j.position)**2)
                if r_ij_sq >= 1:
                    continue

                # state A params for i
                c6_a = c6_lut[lut_idx_a]
                c12_a = c12_lut[lut_idx_a]
                sig_a = sig_lut[lut_idx_a]
                sig6_a = sig6_lut[lut_idx_a]
                #if sig_a < sc_sigma:
                #    sig_a = sc_sigma
                #    sig6_a = sc_sigma6

                c6_b = c6_lut[lut_idx_b]
                c12_b = c12_lut[lut_idx_b]
                sig_b = sig_lut[lut_idx_b]
                sig6_b = sig6_lut[lut_idx_b]  

                #if sig_b < sc_sigma:
                #    sig_b = sc_sigma
                #    sig6_b = sc_sigma6 

                denom_lmbda_a = (sc_alpha*sig6_a*lmbda + r_ij_sq**3)
                denom_for_a = (sc_alpha*sig6_a*lmbda_for + r_ij_sq**3)

                denom_lmbda_b = (sc_alpha*sig6_b*(1-lmbda) + r_ij_sq**3)
                denom_for_b = (sc_alpha*sig6_b*(1-lmbda_for) + r_ij_sq**3)

                u_lmbda += (1-lmbda) * ((c12_a/denom_lmbda_a**2) - (c6_a/denom_lmbda_a)) + (lmbda) * ( (c12_b/denom_lmbda_b**2) - (c6_b/denom_lmbda_b))
                u_for += (1-lmbda_for) * ((c12_a/denom_for_a**2) - (c6_a/denom_for_a)) + (lmbda_for) * ( (c12_b/denom_for_b**2) - (c6_b/denom_for_b))


        my_diffs[window_idx, i_frame, 1] = u_for
        print("frame {}".format(i_frame))
        print("u_for {}".format(u_for))