import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# input data
L = 1.00  # length of domain
D = 0.10  # diffusion coefficient
R = 2.00  # reaction rate
cL = 1.00  # left boundary concentration
NR = 0.50  # right boundary flux
b = -1.00  # non-constant source term 1
c0 = 2.00  # non-constant source term 2

### concentration profile

# numerical solution
num_df = pd.read_csv('fvchem_fvm1d/examples/output_scl1d_function/c_face.csv')
x_num, c_num = num_df['x'].values, num_df['u'].values

# analytical solution
x_ana = np.linspace(0, L, 100)
l = np.sqrt(b/D)
c1 = cL - c0
c2 = -(NR/(D*l) + c1*np.sinh(l*L))/np.cosh(l*L)
c_ana = c0 + c1*np.cosh(l*x_ana) + c2*np.sinh(l*x_ana)

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

### source term profile

# numerical solution
num_df = pd.read_csv('fvchem_fvm1d/examples/output_scl1d_function/r_face.csv')
x_num, r_num = num_df['x'].values, num_df['u'].values
idx = np.argsort(x_num)
x_num, r_num = x_num[idx], r_num[idx]

# analytical solution
x_ana = np.linspace(0, L, 100)
r_ana = b * (c_ana - c0)

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(x_ana, r_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(x_num, r_num, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.xlabel('Distance [m]')
plt.ylabel('Source term [mol m$^{-3}$ s$^{-1}$]')
plt.legend()
plt.show()
