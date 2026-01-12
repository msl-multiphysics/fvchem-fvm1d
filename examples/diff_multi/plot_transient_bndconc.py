import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

### input data

# length of domain
L = 1.00  

# diffusion matrix (full)
D11 = 1.00
D12 = 0.50
D21 = 0.00
D22 = 2.00

# reaction matrix (full)
R11 = 2.00
R12 = 0.50
R21 = 0.50
R22 = 1.00

# boundary conditions (Dirichlet)
c1L = 1.00
c2L = 0.50
c1R = 2.00
c2R = 1.00

# initial conditions (constant)
c1i = 1.00
c2i = 2.00

# time stepping
dt = 0.001
ts_total = 101
ts_every = 20

# analytical series controls
Nx = 100
N_modes = 1000

### matrix functions via eigendecomposition
def matfunc_via_eig(M: np.ndarray, f):
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    return V @ np.diag(f(w)) @ Vinv

def matsqrt(M: np.ndarray):
    return matfunc_via_eig(M, np.sqrt)

def matsin(M: np.ndarray):
    return matfunc_via_eig(M, np.sin)

def matcos(M: np.ndarray):
    return matfunc_via_eig(M, np.cos)

def matexp(M: np.ndarray):
    return matfunc_via_eig(M, np.exp)

### build matrices / vectors
D = np.array([[D11, D12],
              [D21, D22]], dtype=float)

A = np.array([[R11, R12],
              [R21, R22]], dtype=float)

cL = np.array([c1L, c2L], dtype=float)
cR = np.array([c1R, c2R], dtype=float)
ci = np.array([c1i, c2i], dtype=float)

### time array + colors (same style)
t_arr = np.arange(0, ts_total, ts_every) * dt
num_step = t_arr.size

col_num, col_ana = [], []
cmap_num, cmap_ana = plt.get_cmap('Blues'), plt.get_cmap('Grays')
for shade in np.linspace(0.3, 0.7, num_step):
    col_num.append(mcolors.to_hex(cmap_num(shade)))
    col_ana.append(mcolors.to_hex(cmap_ana(shade)))

### steady Dirichlet solution c_s(x) using matrix sin/cos
# D c'' + A c = 0
M = np.linalg.solve(D, A)  # D^{-1}A
B = matsqrt(M)

sinBL = matsin(B * L)
cosBL = matcos(B * L)

sinBL_inv = np.linalg.inv(sinBL)

x_ana = np.linspace(0, L, Nx)
c_s = np.zeros((Nx, 2), dtype=complex)

for i, x in enumerate(x_ana):
    sinBx = matsin(B * x)
    cosBx = matcos(B * x)
    c_s[i, :] = cosBx @ cL + sinBx @ (sinBL_inv @ (cR - cosBL @ cL))

c_s = np.real_if_close(c_s, tol=1e6).astype(float)  # (Nx,2)

### transient coefficients for w(x,t) = c - c_s
# w(0,t)=w(L,t)=0, w(x,0)=ci - c_s(x)
w0 = ci[None, :] - c_s  # (Nx,2)

k = np.arange(1, N_modes + 1)
S = np.sin(np.pi * x_ana[:, None] * k[None, :] / L)  # (Nx, N_modes)

# modal coefficients a_n = (2/L) ∫ w0(x) sin(nπx/L) dx
a0 = np.zeros((N_modes, 2), dtype=complex)
for j in range(N_modes):
    sj = S[:, j]
    a0[j, 0] = (2.0 / L) * np.trapz(w0[:, 0] * sj, x_ana)
    a0[j, 1] = (2.0 / L) * np.trapz(w0[:, 1] * sj, x_ana)

def c_analytical_at_time(t: float) -> np.ndarray:
    """
    Returns c(x,t) on x_ana as (Nx,2) float array.
    """
    # evolve each mode by a matrix exponential
    a_t = np.zeros_like(a0)
    for j in range(N_modes):
        lam = (np.pi * (j + 1) / L) ** 2
        G = A - lam * D
        a_t[j, :] = matexp(G * t) @ a0[j, :]

    w_t = S @ a_t  # (Nx,2)
    w_t = np.real_if_close(w_t, tol=1e6).astype(float)
    return c_s + w_t

### plotting function
def plot_component(comp_idx: int, comp_name: str, num_file_pattern: str):
    # plot formatting (same as your sample)
    plt.figure(figsize=(4/3*3.5, 3.5))
    plt.rcParams['font.family'] = 'IBM Plex Sans'
    plt.rcParams['font.size'] = 12
    plt.rc('legend', fontsize=10)

    # numerical solution snapshots
    for i in range(num_step):
        # read numerical snapshot
        num_df = pd.read_csv(num_file_pattern.format(ts=ts_every * i))
        x_num, c_num = num_df['x'].values, num_df['u'].values

        # plot numerical
        if i == num_step - 1:
            plt.plot(x_num, c_num, color=col_num[i], alpha=0.7, marker='o', label='FVChem')
        else:
            plt.plot(x_num, c_num, color=col_num[i], alpha=0.7, marker='o')

    # analytical solution snapshots
    for i in range(num_step):
        t = t_arr[i]

        if i == 0:
            # plot initial condition as dashed constant line
            plt.plot(x_ana, np.full_like(x_ana, ci[comp_idx]),
                     color=col_ana[i], linestyle='--', linewidth=1)
        else:
            c_t = c_analytical_at_time(t)
            y = c_t[:, comp_idx]

            if i == num_step - 1:
                plt.plot(x_ana, y, color=col_ana[i], linestyle='--', linewidth=1,
                         label='Analytical')
            else:
                plt.plot(x_ana, y, color=col_ana[i], linestyle='--', linewidth=1)

    # finalize plot
    plt.xlabel('Distance [m]')
    plt.ylabel('Concentration [mol m$^{-3}$]')
    plt.legend()
    plt.tight_layout()
    plt.show()

### plot results

# component 1:
plot_component(
    comp_idx=0,
    comp_name="c1",
    num_file_pattern="fvchem_fvm1d/examples/output_transient_bndconc/c1_face_{ts}.csv"
)

# component 2:
plot_component(
    comp_idx=1,
    comp_name="c2",
    num_file_pattern="fvchem_fvm1d/examples/output_transient_bndconc/c2_face_{ts}.csv"
)
