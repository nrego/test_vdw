fe_ds = dr.loadXVG("sc_sigma.xvg")

for i, lmbda in enumerate(for_lmbdas):
    fe_dat = fe_ds.data[lmbda]
    my_dat = my_diffs[i]
    plt.plot(fe_dat, '-o', markersize=10, label=r'Gromacs FE output')
    plt.plot(my_dat[:,0], my_dat[:,1], '-o', linewidth=4, label='my code')
    
    plt.legend()
    plt.xlabel('time (ps)')
    plt.ylabel(r'$\Delta U$ (kJ/mol)')
    plt.title(r'SC LJ $(\lambda_1={})$'.format(lmbda))
    plt.show()

for i, lmbda in enumerate(for_lmbdas):
    fe_dat = fe_ds.data[lmbda]
    my_dat = my_diffs[i]
    plt.plot(my_dat[:,0], np.abs(fe_dat-my_dat[:,1]), '-o', linewidth=4, markersize=10, label='abs difference')
    
    plt.legend()
    plt.xlabel('time (ps)')
    plt.ylabel(r'$\Delta\Delta U$ (kJ/mol)')
    plt.title(r'SC LJ $(\lambda_1={})$'.format(lmbda))
    plt.show()