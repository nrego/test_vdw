fe_ds = dr.loadXVG("sc_sigma.xvg")

fname = 'img_lam_{}_sc{}'
for i, lmbda in enumerate(for_lmbdas):
    plt.clf()
    fe_dat = fe_ds.data[lmbda]
    my_dat = my_diffs[i]
    plt.plot(fe_dat, '-o', markersize=10, label=r'Gromacs FE output')
    plt.plot(my_dat[:,0], my_dat[:,1], '-o', linewidth=4, label='my code')
    
    plt.legend(loc=0)
    plt.xlabel('time (ps)')
    plt.ylabel(r'$\Delta U$ (kJ/mol)')
    plt.title(r'SC LJ $(\lambda_1={})$'.format(lmbda))
    #plt.ylim(-10,2)
    #plt.show()
    plt.savefig(fname.format(int(lmbda*10), '.png'), bbox_inches='tight')

for i, lmbda in enumerate(for_lmbdas):
    plt.clf()
    fe_dat = fe_ds.data[lmbda]
    my_dat = my_diffs[i]
    diffs = np.abs( (fe_dat-my_dat[:,1]) / my_dat[:,1])
    diffs = np.abs((fe_dat-my_dat[:,1]))
    plt.plot(diffs, '-o', linewidth=4, markersize=10, label='abs difference')
    
    plt.legend(loc=0)
    plt.xlabel('time (ps)')
    #plt.ylabel(r'$\frac{\Delta\Delta U}{\Delta U}$')
    plt.ylabel(r'$\Delta U$ (kJ/mol)')
    plt.title(r'SC LJ $(\lambda_1={})$'.format(lmbda))
    plt.savefig(fname.format(int(lmbda*10), '_abs_diff.png'), bbox_inches='tight')
    #plt.show()