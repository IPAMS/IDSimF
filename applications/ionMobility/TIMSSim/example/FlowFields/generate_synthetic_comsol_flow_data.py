import numpy as np 
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def write_vector(vector, fh):
    for val in vector:
        fh.write(str(val)+'\t')
    fh.write('\n')

def write_comsol_txt_format(grid_x, grid_y, grid_z, data, filename):
    
    with open(filename, 'w') as fh: 
        fh.write('% Grid\n')
        write_vector(grid_x, fh)
        write_vector(grid_y, fh)
        if grid_z is not None: 
            write_vector(grid_z, fh)
        
        fh.write('% Data\n')
        if grid_z is not None:
            print(np.shape(data))
            print(f'lx: {len(grid_x)} ly: {len(grid_y)} lz: {len(grid_z)} ')
            for i_z in range(len(grid_z)):
                for i_y in range(len(grid_y)):
                    write_vector(data[:,i_y,i_z], fh)
        else:
            print(np.shape(data))
            print(f'lx: {len(grid_x)} ly: {len(grid_y)} ')
            for i_y in range(len(grid_y)):
                write_vector(data[:,i_y], fh)
        


def quadratic_yz_ramp_x(nx_points, ny_points, nz_points, x_length = 0.01, y_spread=0.004, z_spread=0.005):
    grid_x = np.linspace(0, x_length, nx_points)
    grid_y = np.linspace(-y_spread, y_spread, ny_points)
    grid_z = np.linspace(-z_spread, z_spread, nz_points)
    xx, yy, zz = np.meshgrid(grid_x, grid_y, grid_z, indexing='ij')

    data = (xx/x_length) * (1- ((yy/y_spread)**2.0 + (zz/z_spread)**2.0)/2.0)

    return(grid_x, grid_y, grid_z, data)

def quadratic_r_ramp_x(nx_points, nr_points,  x_length = 0.01, r_length=0.03):
    grid_x = np.linspace(0, x_length, nx_points)
    grid_r = np.linspace(0, r_length, nr_points)
    xx, rr = np.meshgrid(grid_x, grid_r, indexing='ij')

    data = (xx/x_length) * (1- (rr/r_length)**2.0)

    return(grid_x, grid_r, data)

def isosurf_plot(grid_x, grid_y, grid_z, values, isomin, isomax):
    xx, yy, zz = np.meshgrid(grid_x, grid_y, grid_z, indexing='ij')

    fig= go.Figure(data=go.Isosurface(
        x=xx.flatten(),
        y=yy.flatten(),
        z=zz.flatten(),
        value=values.flatten(),
        isomin=isomin,
        isomax=isomax
    ))

    fig.show()

def contour_plot(grid_x, grid_y, values, levels=10):
    xx, yy = np.meshgrid(grid_x, grid_y, indexing='ij')

    plt.contourf(xx, yy, values, levels=levels)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('r')
    plt.show()

def prepare_tims_flow_field():
    # prepare synthetic flow files for TIMS exampe
    grid_x, grid_r, data = quadratic_r_ramp_x(80, 80, 80/1000.0, 80/1000.0)

    p_2d = data * 2000 # Pa
    T_2d = (data * 100) + 298 # K
    vx_2d = data * -70 # m/s

    write_comsol_txt_format(grid_x, grid_r, None, p_2d, 'TIMS_2dAxi_p.txt')
    write_comsol_txt_format(grid_x, grid_r, None, T_2d, 'TIMS_2dAxi_T.txt')
    write_comsol_txt_format(grid_x, grid_r, None, vx_2d, 'TIMS_2dAxi_vx.txt')

    grid_x, grid_y, grid_z, data = quadratic_yz_ramp_x(50, 40, 40,  80/1000.0, 80/1000.0, 80/1000.0)

    p_3d = data * 2000 # Pa
    T_3d = (data * 100) + 298 # K
    vx_3d = data * -70 # m/s
    vy_3d = data * 10 # m/s    
    vz_3d = data * 0 # m/s        

    write_comsol_txt_format(grid_x, grid_y, grid_z, p_3d, 'TIMS_3d_p.txt')
    write_comsol_txt_format(grid_x, grid_y, grid_z, T_3d, 'TIMS_3d_T.txt')
    write_comsol_txt_format(grid_x, grid_y, grid_z, vx_3d,'TIMS_3d_vx.txt')
    write_comsol_txt_format(grid_x, grid_y, grid_z, vy_3d,'TIMS_3d_vy.txt')
    write_comsol_txt_format(grid_x, grid_y, grid_z, vz_3d,'TIMS_3d_vz.txt')



def main():
    prepare_tims_flow_field()
    # flow fields have to be converted to PA with "comsol_to_pa.lua" from SIMION



if __name__ == "__main__":
    main()