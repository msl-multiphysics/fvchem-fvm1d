import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# input data
L = 1.00  # length of domain
D = 0.10  # diffusion coefficient
R = 2.00  # reaction rate
k = 0.20  # mass transfer coefficient
cextL = 1.00  # left external concentration
cR = 1.00  # right boundary concentration

# numerical solution
num_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_bndmtrn/c_face.csv')
x_num, c_num = num_df['x'].values, num_df['u'].values

# analytical solution
a_mat = np.array([[1, k], [-L/D, 1]])
b_vec = np.array([k*cextL, cR + (R*L**2)/(2*D)])
c1, c2 = np.linalg.solve(a_mat, b_vec)
x_ana = np.linspace(0, L, 100)
c_ana = (-R/2/D)*x_ana**2 - (c1/D)*x_ana + c2

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
