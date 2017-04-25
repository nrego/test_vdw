
from mdtools import dr


fe_ds = dr.loadXVG('sc_no_sigma.xvg')
fe_dat = fe_ds.data[lmbda_for]

plt.plot(fe_dat, '-o', label=r'$\lambda={}$ (Gromacs FE)'.format(lmbda))
plt.plot(my_diffs[:,0], my_diffs[:,1], '-o', label=r'$\lambda={}$ (calculated directly)'.format(lmbda))