import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# input data
LA = 0.50  # length of domain
DA = 0.20  # diffusion coefficient
RA = 2.00  # reaction rate
LB = 1.00  # length of domain
DB = 0.10  # diffusion coefficient
RB = 1.00  # reaction rate
cL = 1.00  # left boundary concentration
NR = 0.50  # right boundary flux

# numerical solution
numA_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_diff_itrcont/ca_face.csv')
numB_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_diff_itrcont/cb_face.csv')
xA_num, cA_num = numA_df['x'].values, numA_df['u'].values
xB_num, cB_num = numB_df['x'].values, numB_df['u'].values
idxA = np.argsort(xA_num)
idxB = np.argsort(xB_num)
xA_num, cA_num = xA_num[idxA], cA_num[idxA]
xB_num, cB_num = xB_num[idxB], cB_num[idxB]

# analytical solution
c3 = cL
c2 = NR - RB*(LA + LB)
c1 = c2 + LA*(RB - RA)
c4 = -RA/(2*DA)*LA**2 - c1/DA*LA + c3 + RB/(2*DB)*LA**2 + c2/DB*LA
xA_ana = np.linspace(0, LA, 100)
cA_ana = (-RA/2/DA)*xA_ana**2 - (c1/DA)*xA_ana + c3
xB_ana = np.linspace(LA, LA+LB, 100)
cB_ana = (-RB/2/DB)*xB_ana**2 - (c2/DB)*xB_ana + c4

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(xA_ana, cA_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(xB_ana, cB_ana, color='#515151', linestyle='--', linewidth=1)
plt.plot(xA_num, cA_num, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.plot(xB_num, cB_num, color='#2070B4', alpha=0.7, marker='o')
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.show()
