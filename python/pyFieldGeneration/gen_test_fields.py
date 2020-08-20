# FIXME: Update generation of test fields...

import numpy as np
import IDSimPy.preprocessing.field_generation as fg


# define simple linear scalar field:
grid_points = [[0, 2, 5, 15], [0, 2, 10], [0, 2, 5, 7, 10]]
X, Y,Z = np.meshgrid(grid_points[0], grid_points[1], grid_points[2], indexing='ij')
S = X+Y+Z
fields = [{'name': 'test_field', 'data': S}]
dat = {"grid_points": grid_points, "meshgrid": [X,Y,Z], "fields": fields}
fg.write_3d_scalar_fields_to_hdf5(dat,'../testfields/test_linear_scalar_field_01.h5')


# define simple linear vector field:
grid_points = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20], [-10, 0, 10], [0, 10]]
X,Y,Z = np.meshgrid(grid_points[0], grid_points[1], grid_points[2], indexing='ij')
S_x = X+Y+Z
S_y1 = np.zeros(np.shape(X)) + 5.0
S_z1 = np.zeros(np.shape(X)) + 1.0
S_y2 = np.zeros(np.shape(X)) + 15.0
S_z2 = np.zeros(np.shape(X)) + 11.0


fields = [
    {'name': 'test_vectorfield_1', 'data': [S_x, S_y1, S_z1]},
    {'name': 'test_vectorfield_2', 'data': [S_x, S_y2, S_z2]}
]

dat = {"grid_points": grid_points, "meshgrid": [X, Y, Z], "fields": fields}
fg.write_3d_vector_fields_to_hdf5(dat, '../testfields/test_linear_vector_field_01.h5')