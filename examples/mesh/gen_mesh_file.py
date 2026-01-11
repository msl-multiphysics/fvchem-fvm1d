import os

import numpy as np
import pandas as pd

# settings
x_min = 0.0
dx_list = [0.2, 0.5, 0.3]  # lengths of segments
nc_list = [10, 20, 15]  # number of cells

# create input file directory
os.makedirs('fvchem_fvm1d/examples/input_mesh_file', exist_ok=True)

# face data
face_x_list = []
for i in range(len(nc_list)):
    start_x = x_min + sum(dx_list[:i])
    end_x = x_min + sum(dx_list[:i+1])
    face_x_seg = np.linspace(start_x, end_x, nc_list[i] + 1)
    if i > 0:
        face_x_seg = face_x_seg[1:]  # avoid duplicate face at segment boundary
    face_x_list.append(face_x_seg)
face_x = np.concatenate(face_x_list)
num_face = face_x.size
face_id = -np.arange(1, num_face + 1)
pd.DataFrame({'fid': face_id, 'x': face_x}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_face.csv', index=False)

# cell data
cell_x = 0.5 * (face_x[:-1] + face_x[1:])
cell_dx = np.diff(face_x)
num_cell = cell_x.size
cell_id = np.arange(num_cell)
pd.DataFrame({'cid': cell_id, 'x': cell_x, 'dx': cell_dx}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_cell.csv', index=False)

# cell-cell data
cc_nid0 = np.concatenate([[-1], np.arange(num_cell - 1)])
cc_nid1 = np.concatenate([np.arange(1, num_cell), [-num_face]])
cc_dist0 = np.concatenate([[cell_x[0] - face_x[0]], cell_x[1:] - cell_x[:-1]])
cc_dist1 = np.concatenate([cell_x[1:] - cell_x[:-1], [face_x[-1] - cell_x[-1]]])
pd.DataFrame({'cid': cell_id, 'nid_0': cc_nid0, 'nid_1': cc_nid1, 'dist_0': cc_dist0, 'dist_1': cc_dist1}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_cell_cell.csv', index=False)

# cell-face data
cf_fid0 = -np.arange(1, num_cell + 1)
cf_fid1 = -np.arange(2, num_cell + 2)
cf_dist0 = cell_x - face_x[:-1]
cf_dist1 = face_x[1:] - cell_x
cf_norm0 = -np.ones(num_cell)
cf_norm1 = np.ones(num_cell)
pd.DataFrame({'cid': cell_id, 'fid_0': cf_fid0, 'fid_1': cf_fid1, 'dist_0': cf_dist0, 'dist_1': cf_dist1, 'norm_0': cf_norm0, 'norm_1': cf_norm1}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_cell_face.csv', index=False)

# face-cell data
fc_cid0 = np.concatenate([[-1], np.arange(num_cell)])
fc_cid1 = np.concatenate([np.arange(num_cell), [-num_cell]])
fc_dist0 = np.concatenate([[0], cell_x - face_x[:-1]])
fc_dist1 = np.concatenate([face_x[1:] - cell_x, [0]])
pd.DataFrame({'fid': face_id, 'cid_0': fc_cid0, 'cid_1': fc_cid1, 'dist_0': fc_dist0, 'dist_1': fc_dist1}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_face_cell.csv', index=False)

# region data
reg_id = np.zeros_like(cell_id)
start_idx = 0
for i in range(len(nc_list)):
    end_idx = start_idx + nc_list[i]
    reg_id[start_idx:end_idx] = i
    start_idx = end_idx
pd.DataFrame({'rid': reg_id, 'cid': cell_id}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_reg.csv', index=False)

# boundary data
bnd_id = np.arange(2*len(nc_list))
bnd_cid = np.zeros_like(bnd_id, dtype=int)
bnd_loc = np.zeros_like(bnd_id, dtype=int)
for i in range(len(nc_list)):
    # left boundary of segment
    bnd_cid[2*i] = sum(nc_list[:i])
    bnd_loc[2*i] = 0  # left
    # right boundary of segment
    bnd_cid[2*i + 1] = sum(nc_list[:i + 1]) - 1
    bnd_loc[2*i + 1] = 1  # right
pd.DataFrame({'bid': bnd_id, 'cid': bnd_cid, 'loc': bnd_loc}).to_csv('fvchem_fvm1d/examples/input_mesh_file/mesh_bnd.csv', index=False)
