## Do vdw calculation from md trajectory##

from __future__ import division, print_function

import MDAnalysis
import numpy as np

import argparse

parser = argparse.ArgumentParser('find VdW for scheme3')
parser.add_argument('-s', '--tprfile', metavar='TPR', required=True, type=str,
                    help='TPR file')
parser.add_argument('-f', '--xtcfile', metavar='XTC', required=True, type=str,
                    help='XTC file')
parser.add_argument('--lmbda', metavar='LMBDA', default=0.0, type=float,
                    help='Initial lambda value')
parser.add_argument('--sc-alpha', metavar='SC_ALPHA', default=0.5, type=float,
                    help='Alpha value for SC LJ (default: 0.5)')
parser.add_argument('--sc-sigma', metavar='SC_SIGMA', default=0.3, type=float,
                    help='SC minimum sigma (default: 0.3 nm)')
args = parser.parse_args()

# set combined sigma to this if it's smaller
sc_sigma = args.sc_sigma
sc_sigma6 = sc_sigma**6
sc_alpha = args.sc_alpha

lmbda = args.lmbda

if lmbda == 1.0:
    for_lmbdas = [0.9]
elif lmbda == 0.0:
    for_lmbdas = [0.1]
else:
    for_lmbdas = [lmbda-0.1, lmbda+0.1]

tprfile = args.tprfile
xtcfile = args.xtcfile

univ = MDAnalysis.Universe(tprfile, xtcfile)

alc_indices = np.arange(877, 888)
atm_indices = np.arange(univ.atoms.n_atoms)


excls = {
         877:(870, 872, 873, 874, 875, 876, 877, 878, 
            879, 880, 884, 888),
         878:(870, 872, 873, 874, 875, 876, 877, 878, 
            879, 880, 881, 882, 883, 884, 885, 886, 887, 888),
         879:(872, 874, 875, 876, 877, 878, 879, 880, 
            881, 882, 883, 884, 885, 886, 887),
         880:(872, 874, 875, 876, 877, 878, 879, 880, 
            881, 882, 883, 884, 885, 886, 887),
         881:(874, 878, 879, 880, 881, 882, 883, 884),
         882:(874, 878, 879, 880, 881, 882, 883, 884),
         883:(874, 878, 879, 880, 881, 882, 883, 884),
         884:(872, 874, 875, 876, 877, 878, 879, 880, 
            881, 882, 883, 884, 885, 886, 887),
         885:(874, 878, 879, 880, 884, 885, 886, 887),
         886:(874, 878, 879, 880, 884, 885, 886, 887),
         887:(874, 878, 879, 880, 884, 885, 886, 887)
}


# 14 pair list for each excluded atom
pairs = {
    877: (870, 873, 879, 880, 884, 888), # HB3
    878: (870, 873, 888),
    879: (872, 875, 876, 877, 881, 882, 883, 885, 886, 887),
    880: (872, 875, 876, 877, 885, 886, 887),
    881: (874, 879, 884),
    882: (874, 879, 884),
    883: (874, 879, 884),
    884: (872, 875, 876, 877, 881, 882, 883),
    885: (874, 879, 880),
    886: (874, 879, 880),
    887: (874, 879, 880)
}
# dictionary of alchemical atom type transofmations.
#    keyed by alchemical index of an atom that is to be transformed
#    valued by tuple (Astate_idx, Bstate_idx)
#  NOTE: This is topology specific!!!
alc_types = {
    877: ('DUM_HC', 'HC'),
    878: ('CT', 'DUM_CT'),
    879: ('HC', 'DUM_HC'),
    880: ('CT', 'DUM_CT'),
    881: ('HC', 'DUM_HC'),
    882: ('HC', 'DUM_HC'),
    883: ('HC', 'DUM_HC'),
    884: ('CT', 'DUM_CT'),
    885: ('HC', 'DUM_HC'),
    886: ('HC', 'DUM_HC'),
    887: ('HC', 'DUM_HC')
}

type_lookup = {
    'N3': 0,
    'H': 1,
    'CT': 2,
    'HP': 3,
    'HC': 4,
    'C': 5,
    'O': 6,
    'N': 7,
    'H1': 8,
    'SH': 9,
    'HS': 10,
    'OH': 11,
    'HO': 12,
    'CA': 13,
    'HA': 14,
    'O2': 15,
    'CC': 5,
    'NB': 16,
    'CR': 17,
    'H5': 18,
    'NA': 7,
    'CW': 17,
    'H4': 19,
    'DUM_HC': 20,
    'DUM_CT': 21,
    'DUM': 22,
    'OW_spc': 23,
    'HW_spc': 20
}
#format: 24 atom types, atomtype i is a tuple with (sigma, eps)
atmtypes = [
(3.25000e-01,  7.11280e-01), # N3
(1.06908e-01,  6.56888e-02), # H
(3.39967e-01,  4.57730e-01), # CT
(1.95998e-01,  6.56888e-02), # HP
(2.64953e-01,  6.56888e-02), # HC
(3.39967e-01,  3.59824e-01), # C
(2.95992e-01,  8.78640e-01), # O
(3.25000e-01,  7.11280e-01), # N
(2.47135e-01,  6.56888e-02), # H1
(3.56359e-01,  1.04600e+00), # SH
(1.06908e-01,  6.56888e-02), # HS
(3.06647e-01,  8.80314e-01), # OH
(0.00000e+00,  0.00000e+00), # HO
(3.39967e-01,  3.59824e-01), # CA
(2.59964e-01,  6.27600e-02), # HA
(2.95992e-01,  8.78640e-01), # O2
(3.25000e-01,  7.11280e-01), # NB
(3.39967e-01,  3.59824e-01), # CR
(2.42146e-01,  6.27600e-02), # H5
(2.51055e-01,  6.27600e-02), # H4
(0.00000e+00,  0.00000e+00), # DUM_HC
(0.00000e+00,  0.00000e+00), # DUM_CT
(0.00000e+00,  0.00000e+00), # DUM
(3.16557e-01,  6.50629e-01)  # OW_spc
]


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

fudge_vdw = 0.5

n_frames = univ.trajectory.n_frames
#n_frames = 1
my_diffs = np.zeros((n_frames, len(for_lmbdas)+1))

for window_idx, lmbda_for in enumerate(for_lmbdas):
    print("foreign lambda: {}".format(lmbda_for))
    for i_frame in range(n_frames):
        univ.trajectory[i_frame]
        if my_diffs[i_frame, 0] == 0.0:
            my_diffs[i_frame, 0] = univ.trajectory.time
        elif my_diffs[i_frame, 0] != univ.trajectory.time:
            print("WARNING: DIFFERENT TIMES: {}  {}".format(my_diffs[i_frame,0], univ.trajectory.time))
        univ.atoms.positions = univ.atoms.positions / 10.0
        # Calculate VdW energy differences between lambdas
        u_lmbda = 0.0
        u_for = 0.0
        for i in alc_indices:
            # all atoms separated by more than nrexcl bonds (i.e. not excluded)
            # Note: If i and j are both in alc_indices, skip if j !> i
            incl_indices = np.setdiff1d(atm_indices, excls[i])

            # 1-4 pairs
            pair_indices = pairs[i]

            assert np.intersect1d(incl_indices, pair_indices).size == 0, "Double counting some pairs as 14 pairs!!"

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
                    type_j_b = type_lookup[name_j_b]
                else:
                    type_j_a = type_j_b = type_lookup[atm_j.type]  

                #print("j: {}".format(j))
                #print("  type a: {}, type b: {}".format(type_j_a, type_j_b))      

                lut_idx_a = type_i_a * n_atmtype + type_j_a
                lut_idx_b = type_i_b * n_atmtype + type_j_b

                r_ij_sq = np.sum((atm_i.position - atm_j.position)**2)
                #print("  r_ij_sq: {}".format(r_ij_sq))
                if r_ij_sq >= 1:
                    continue

                # state A params for i
                c6_a = c6_lut[lut_idx_a]
                c12_a = c12_lut[lut_idx_a]
                sig_a = sig_lut[lut_idx_a]
                sig6_a = sig6_lut[lut_idx_a]

                c6_b = c6_lut[lut_idx_b]
                c12_b = c12_lut[lut_idx_b]
                sig_b = sig_lut[lut_idx_b]
                sig6_b = sig6_lut[lut_idx_b]  

                denom_lmbda_a = (sc_alpha*sig6_a*lmbda + r_ij_sq**3)
                denom_for_a = (sc_alpha*sig6_a*lmbda_for + r_ij_sq**3)

                denom_lmbda_b = (sc_alpha*sig6_b*(1-lmbda) + r_ij_sq**3)
                denom_for_b = (sc_alpha*sig6_b*(1-lmbda_for) + r_ij_sq**3)

                this_u_lmbda = (1-lmbda) * ((c12_a/denom_lmbda_a**2) - (c6_a/denom_lmbda_a)) + (lmbda) * ( (c12_b/denom_lmbda_b**2) - (c6_b/denom_lmbda_b))
                this_u_for = (1-lmbda_for) * ((c12_a/denom_for_a**2) - (c6_a/denom_for_a)) + (lmbda_for) * ( (c12_b/denom_for_b**2) - (c6_b/denom_for_b))
                #print("  u_lmbda contrib: {}".format(this_u_lmbda))
                #print("  u_for contrib: {}".format(this_u_for))
                u_lmbda += this_u_lmbda
                u_for += this_u_for

            for j in pair_indices:

                atm_j = univ.atoms[j]
                if j in alc_indices:
                    if j < i:
                        continue
                    name_j_a, name_j_b = alc_types[j]
                    type_j_a = type_lookup[name_j_a]
                    type_j_b = type_lookup[name_j_b]
                else:
                    type_j_a = type_j_b = type_lookup[atm_j.type]  

                #print("j (14 pair): {}".format(j))
                #print("  type a: {}, type b: {}".format(type_j_a, type_j_b))      

                lut_idx_a = type_i_a * n_atmtype + type_j_a
                lut_idx_b = type_i_b * n_atmtype + type_j_b

                r_ij_sq = np.sum((atm_i.position - atm_j.position)**2)
                #print("  r_ij_sq: {}".format(r_ij_sq))
                if r_ij_sq >= 1:
                    continue

                # state A params for i
                c6_a = c6_lut[lut_idx_a]
                c12_a = c12_lut[lut_idx_a]
                sig_a = sig_lut[lut_idx_a]
                sig6_a = sig6_lut[lut_idx_a]

                c6_b = c6_lut[lut_idx_b]
                c12_b = c12_lut[lut_idx_b]
                sig_b = sig_lut[lut_idx_b]
                sig6_b = sig6_lut[lut_idx_b]  

                denom_lmbda_a = (sc_alpha*sig6_a*lmbda + r_ij_sq**3)
                denom_for_a = (sc_alpha*sig6_a*lmbda_for + r_ij_sq**3)

                denom_lmbda_b = (sc_alpha*sig6_b*(1-lmbda) + r_ij_sq**3)
                denom_for_b = (sc_alpha*sig6_b*(1-lmbda_for) + r_ij_sq**3)

                this_u_lmbda = fudge_vdw * ((1-lmbda) * ((c12_a/denom_lmbda_a**2) - (c6_a/denom_lmbda_a)) + (lmbda) * ( (c12_b/denom_lmbda_b**2) - (c6_b/denom_lmbda_b)))
                this_u_for = fudge_vdw * ((1-lmbda_for) * ((c12_a/denom_for_a**2) - (c6_a/denom_for_a)) + (lmbda_for) * ( (c12_b/denom_for_b**2) - (c6_b/denom_for_b)))
                #print("  u_lmbda contrib: {}".format(this_u_lmbda))
                #print("  u_for contrib: {}".format(this_u_for))
                u_lmbda += this_u_lmbda
                u_for += this_u_for


        my_diffs[i_frame, window_idx+1] = u_for - u_lmbda
        print("frame {}".format(i_frame))
        print("delta u {}".format(u_for - u_lmbda))


print("saving data")

headerstr = 'time (ps)' + ''.join(['     for_lmbda={}'.format(for_lmbda) for for_lmbda in for_lmbdas])
np.savetxt('diffs_{}.dat'.format(lmbda), my_diffs, header=headerstr)



