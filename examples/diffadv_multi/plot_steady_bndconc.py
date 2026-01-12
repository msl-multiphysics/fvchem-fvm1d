import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### input data

# length of domain
L = 1.00

# component 1
D1 = 0.10  # diffusion coefficient
u1 = 2.00  # advection velocity
c1L = 1.00  # left boundary concentration
c1R = 0.50  # right boundary concentration

# component 2
D2 = 0.20  # diffusion coefficient
u2 = 10.00  # advection velocity
c2L = 1.00  # left boundary concentration
c2R = 0.50  # right boundary concentration

### numerical solution
num1_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_bndconc/c1_face.csv')
num2_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_bndconc/c2_face.csv')
x_num1, c_num1 = num1_df['x'].values, num1_df['u'].values
x_num2, c_num2 = num2_df['x'].values, num2_df['u'].values

### analytical solution
x_ana = np.linspace(0, L, 100)
c1_ana = c1L + (c1R - c1L)*(np.exp(u1*x_ana/D1) - 1)/(np.exp(u1*L/D1) - 1)
c2_ana = c2L + (c2R - c2L)*(np.exp(u2*x_ana/D2) - 1)/(np.exp(u2*L/D2) - 1)

### plot component 1

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(x_ana, c1_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(x_num1, c_num1, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.show()

### plot component 2

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(x_ana, c2_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(x_num2, c_num2, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.show()
