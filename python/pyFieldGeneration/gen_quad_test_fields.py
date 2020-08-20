import numpy as np
import IDSimPy.preprocessing.field_generation as fg

# define 2d axial symmetric pressure field and translate it to 3d cartesian:
grid_r = [np.linspace(0, 0.01, 30)]
grid_z = [np.linspace(0, 0.1, 200)]
#grid_r = np.arange(0,20)
#grid_z = 10

R, Z = np.meshgrid(grid_r, grid_z, indexing='ij')
P = R*Z
P[:, 1:10] = 0

X, Y, Z, P3 = fg.transform_2d_axial_to_3d(R, Z, P)
grid_points = [X[0, :, 0], Y[:, 0, 0], Z[0, 0, :]]
grid_points_len = [len(X[0, :, 0]), len(Y[:, 0, 0]), len(Z[0, 0, :])]

print(P3.shape)
fields = [{'name': 'radial_pressure', 'data': P3}]
dat = {"grid_points": grid_points, "fields": fields}
fg.write_3d_scalar_fields_to_hdf5(dat, '../testfields/quad_dev_pressure_radial.h5')


# define simple 3d pressure field for testing / development:
dx = 0.4
dyz = 0.1
grid_points = [np.linspace(-dx, dx, 300), np.linspace(-dyz, dyz, 40), np.linspace(-dyz, dyz, 40)]

X, Y, Z = np.meshgrid(grid_points[0], grid_points[1], grid_points[2], indexing='ij')
print(np.shape(X))
S = np.zeros(np.shape(X))+10.0#(dx-X)*10
fields = [{'name': 'pressure', 'data': S}]
dat = {"grid_points": grid_points, "meshgrid": [X, Y, Z], "fields": fields}
fg.write_3d_scalar_fields_to_hdf5(dat, '../testfields/quad_dev_pressure_3d.h5')

# define simple 3d vector flow field for testing / development:
S_x = ((dyz-np.abs(Y)) * (dyz-np.abs(Z)))*150/(dyz**2)
S_y = np.zeros(np.shape(X))
S_z = np.zeros(np.shape(X))

fields = [{'name': 'velocity', 'data': [S_x, S_y, S_z]}]
dat = {"grid_points": grid_points, "meshgrid": [X, Y, Z], "fields": fields}
fg.write_3d_vector_fields_to_hdf5(dat, '../testfields/quad_dev_flow_3d.h5')


# define simple 3d vector electrical field for testing / development:
S_x = 0.0+np.zeros(np.shape(X))
S_y = Y * -10.0
S_z = Z * -10.0

fields = [{'name': 'electric field', 'data': [S_x, S_y, S_z]}]
dat = {"grid_points": grid_points, "meshgrid": [X, Y, Z], "fields": fields}
fg.write_3d_vector_fields_as_vtk_point_data(dat, '../testfields/quad_dev_field.h5')

