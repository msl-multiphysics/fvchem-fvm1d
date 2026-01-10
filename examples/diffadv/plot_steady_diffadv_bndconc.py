import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# input data
L = 1.00  # length of domain
D = 0.10  # diffusion coefficient
u = 2.00  # advection velocity
cL = 1.00  # left boundary concentration
cR = 0.50  # right boundary concentration

# numerical solution
num_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_diffadv_bndconc/c_face.csv')
x_num, c_num = num_df['x'].values, num_df['u'].values
idx = np.argsort(x_num)
x_num, c_num = x_num[idx], c_num[idx]

# analytical solution
x_ana = np.linspace(0, L, 100)
c_ana = cL + (cR - cL)*(np.exp(u*x_ana/D) - 1)/(np.exp(u*L/D) - 1)

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(x_ana, c_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(x_num, c_num, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.show()
