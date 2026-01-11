import os

import numpy as np
import pandas as pd

# settings
x_min = 0.0
x_max = 1.0
num_cell = 20
r_ref = 2.0  # reference reaction rate

# create input file directory
os.makedirs('fvchem_fvm1d/examples/input_scl1d_file', exist_ok=True)

# cell data
dx = (x_max - x_min) / num_cell
x_cell = 0.5*dx + np.arange(num_cell)*dx
r_cell = r_ref * x_cell  # reaction rate
id_cell = np.arange(num_cell)
cell_df = pd.DataFrame({'cid': id_cell, 'x': x_cell, 'u': r_cell})
cell_df.to_csv('fvchem_fvm1d/examples/input_scl1d_file/r_cell.csv', index=False)

# face data
x_face = np.arange(num_cell+1)*dx
r_face = r_ref * x_face  # reaction rate
id_face = -np.arange(1, num_cell+2)
face_df = pd.DataFrame({'fid': id_face, 'x': x_face, 'u': r_face})
face_df.to_csv('fvchem_fvm1d/examples/input_scl1d_file/r_face.csv', index=False)
