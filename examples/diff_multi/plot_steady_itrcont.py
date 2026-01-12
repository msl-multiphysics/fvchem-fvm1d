import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

### input data

# length of domain
L = 1.00  

# diffusion matrix
D11 = 1.00
D12 = 0.50
D21 = 0.00
D22 = 2.00

# reaction matrix
R11 = 2.00
R12 = 0.50
R21 = 0.50
R22 = 1.00

# boundary conditions
c1L = 1.00
c2L = 0.50
c1R = 2.00
c2R = 1.00

### numerical solution
num1a_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_itrcont/c1_a_face.csv')
num1b_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_itrcont/c1_b_face.csv')
num2a_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_itrcont/c2_a_face.csv')
num2b_df = pd.read_csv('fvchem_fvm1d/examples/output_steady_itrcont/c2_b_face.csv')
x1a_num, c1a_num = num1a_df['x'].values, num1a_df['u'].values
x1b_num, c1b_num = num1b_df['x'].values, num1b_df['u'].values
x2a_num, c2a_num = num2a_df['x'].values, num2a_df['u'].values
x2b_num, c2b_num = num2b_df['x'].values, num2b_df['u'].values
x1_num = np.concatenate([x1a_num, x1b_num])
c1_num = np.concatenate([c1a_num, c1b_num])
x2_num = np.concatenate([x2a_num, x2b_num])
c2_num = np.concatenate([c2a_num, c2b_num])

### analytical solution

# matrix functions via eigendecomposition
def matfunc_via_eig(M: np.ndarray, f):
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    fw = f(w)  # elementwise on eigenvalues (possibly complex)
    return V @ np.diag(fw) @ Vinv

def matsqrt(M: np.ndarray):
    return matfunc_via_eig(M, np.sqrt)

def matsin(M: np.ndarray):
    return matfunc_via_eig(M, np.sin)

def matcos(M: np.ndarray):
    return matfunc_via_eig(M, np.cos)

# build matrices
D = np.array([[D11, D12],
              [D21, D22]], dtype=float)

A = np.array([[R11, R12],
              [R21, R22]], dtype=float)

cL = np.array([c1L, c2L], dtype=float)
cR = np.array([c1R, c2R], dtype=float)

# M = D^{-1} A
M = np.linalg.solve(D, A)

# B = sqrt(M)
B = matsqrt(M)

# Precompute matrix sin/cos terms
sinBL = matsin(B * L)
cosBL = matcos(B * L)

sinBL_inv = np.linalg.inv(sinBL)

# evaluate analytical solution
x_ana = np.linspace(0, L, 200)
c_ana = np.zeros((len(x_ana), 2), dtype=complex)

for i, x in enumerate(x_ana):
    sinBx = matsin(B * x)
    cosBx = matcos(B * x)
    c = cosBx @ cL + sinBx @ (sinBL_inv @ (cR - cosBL @ cL))
    c_ana[i, :] = c

# for real inputs, the result should typically be real
c_ana = np.real_if_close(c_ana, tol=1e6).astype(float)
c1_ana = c_ana[:, 0]
c2_ana = c_ana[:, 1]

### plot component 1

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(x_ana, c1_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(x1_num, c1_num, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.tight_layout()
plt.show()

### plot component 2

# plot formatting
plt.figure(figsize=(4/3*3.5, 3.5))
plt.rcParams['font.family'] = 'IBM Plex Sans'
plt.rcParams['font.size'] = 12
plt.rc('legend', fontsize=10)

# plot results
plt.plot(x_ana, c2_ana, color='#515151', linestyle='--', linewidth=1, label='Analytical')
plt.plot(x2_num, c2_num, color='#2070B4', alpha=0.7, marker='o', label='FVChem')
plt.xlabel('Distance [m]')
plt.ylabel('Concentration [mol m$^{-3}$]')
plt.legend()
plt.tight_layout()
plt.show()
