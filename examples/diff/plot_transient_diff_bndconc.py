import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# input data
L = 1.00  # length of domain
D = 0.10  # diffusion coefficient
R = 2.00  # reaction rate
cL = 1.00  # left boundary concentration
cR = 0.50  # right boundary flux
ci = 1.00  # initial concentration

# time stepping
dt = 0.01
ts_total = 101
ts_every = 20

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# time stepping calculation
t_arr = np.arange(0, ts_total, ts_every) * dt
num_step = t_arr.size

# color maps
col_num, col_ana = [], []
cmap_num, cmap_ana = plt.get_cmap('Blues'), plt.get_cmap('Grays')
for shade in np.linspace(0.3, 0.7, num_step):
    col_num.append(mcolors.to_hex(cmap_num(shade)))
    col_ana.append(mcolors.to_hex(cmap_ana(shade)))

# numerical solution
for i in range(num_step):

    # get data
    num_df = pd.read_csv(f'fvchem_fvm1d/examples/output_transient_diff_bndconc/c_face_{ts_every * i}.csv')
    x_num, c_num = num_df['x'].values, num_df['u'].values

    # plot results
    if i == num_step - 1:
        plt.plot(x_num, c_num, color=col_num[i], alpha=0.7, marker='o', label='FVChem')
    else:
        plt.plot(x_num, c_num, color=col_num[i], alpha=0.7, marker='o')

# analytical solution
for i in range(num_step):

    # steady component
    x_ana = np.linspace(0, L, 100)
    c_ss = cL + (cR - cL)*(x_ana/L) + (R/2/D)*x_ana*(L - x_ana)

    # transient component
    c_tr = np.zeros_like(x_ana)
    for n in range(1, 1000):
        bn = (2/n/np.pi) * ((ci - cL) - (-1)**n * (ci - cR)) - 2*R*L**2/(D*n**3 * np.pi**3) * (1 - (-1)**n)
        c_tr += bn * np.sin(n * np.pi * x_ana / L) * np.exp(-D * (n * np.pi / L)**2 * t_arr[i])

    # total concentration
    c_ana = c_ss + c_tr

    # plot results
    if i == 0:
        plt.plot(x_ana, np.full_like(x_ana, ci), color=col_ana[i], linestyle='--', linewidth=1)
    elif i == num_step - 1:
        plt.plot(x_ana, c_ana, color=col_ana[i], linestyle='--', linewidth=1, label='Analytical')
    else:
        plt.plot(x_ana, c_ana, color=col_ana[i], linestyle='--', linewidth=1)

# plot results
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.show()
