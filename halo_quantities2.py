import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
import os
import sys
import ytree
import numpy as np
"""
Run this second:

python halo_quantities2.py Natural DD0274

(or python halo_quantities2.py Tempest8 DD0274)
"""

def get_rockstar_data(rstar_fn, halo_id):
    """
    Use ytree to get all of the halo centroids, virial radii, and redshift info; store in a dict
    """
    # load up dataset and get appropriate TreeNode
    a = ytree.load(rstar_fn)
    t = a[a["Orig_halo_ID"] == halo_id][0]

    redshift_arr = t['prog', 'redshift']
    x_arr = t['prog', 'x'].in_units('unitary')
    y_arr = t['prog', 'y'].in_units('unitary')
    z_arr = t['prog', 'z'].in_units('unitary')
    rvir_arr = t['prog', 'Rvir'].convert_to_units('kpc')

    return {'redshift_arr':redshift_arr, 'x_arr':x_arr, 'y_arr':y_arr, 'z_arr':z_arr, 'rvir_arr':rvir_arr}

def read_rockstar_center(rockstar_data, ds):
    """
    Interpolate halo center from rockstar merger tree
    """
    redshift = ds.current_redshift
    redshift_arr = rockstar_data['redshift_arr']
    x = np.interp(redshift, redshift_arr, rockstar_data['x_arr'].in_units('unitary'))
    y = np.interp(redshift, redshift_arr, rockstar_data['y_arr'].in_units('unitary'))
    z = np.interp(redshift, redshift_arr, rockstar_data['z_arr'].in_units('unitary'))

    # Construct YTArray with correct units of original dataset (i.e., unitary)
    #arr = np.array([x,y,z]) * rockstar_data['x_arr'][0].uq
    arr = ds.arr([x,y,z], 'unitary')
    return arr

if __name__ == '__main__':
    run = sys.argv[1]
    fn = sys.argv[2]
    cwd = os.getcwd()
    path = '/Users/chummels/scratch/Tempest/new'
    full_fn = os.path.join(path, run, fn, fn)
    data_ds = yt.load(full_fn)
    #hc = HaloCatalog(data_ds=data_ds, finder_method='hop', finder_kwargs={"dm_only":False})
    #hc.create()
    # Find centroid of Natural halo from rockstar files
    rockstar_data = get_rockstar_data(os.path.join(path, 'tree_27.dat'), 27)
    c = read_rockstar_center(rockstar_data, data_ds)
    offset = data_ds.quan(150, 'kpc').to('unitary')
    x_min = c[0]-offset
    x_max = c[0]+offset
    y_min = c[1]-offset
    y_max = c[1]+offset
    z_min = c[2]-offset
    z_max = c[2]+offset
    print(x_min)
    print(x_max)
    print(y_min)
    print(y_max)
    print(z_min)
    print(z_max)
    # Instantiate a catalog using those two paramter files
    halos_ds = yt.load(os.path.join(run, fn, 'catalog', 'catalog.0.h5'))
    hc = HaloCatalog(data_ds=data_ds, halos_ds=halos_ds, 
                    output_dir='.')
    hc.add_filter('quantity_value', 'particle_mass', '>', 1E10, 'Msun')
    # 0.49243498  0.48236668  0.50483814
    hc.add_filter('quantity_value', 'particle_position_x', '>', x_min.v, 'unitary')
    hc.add_filter('quantity_value', 'particle_position_x', '<', x_max.v, 'unitary')
    hc.add_filter('quantity_value', 'particle_position_y', '>', y_min.v, 'unitary')
    hc.add_filter('quantity_value', 'particle_position_y', '<', y_max.v, 'unitary')
    hc.add_filter('quantity_value', 'particle_position_z', '>', z_min.v, 'unitary')
    hc.add_filter('quantity_value', 'particle_position_z', '<', z_max.v, 'unitary')
    hc.add_callback("iterative_center_of_mass", inner_ratio=0.05)
    #hc.add_callback("sphere")
    #hc.add_callback("sphere_field_max_recenter", 'density')
    hc.add_recipe("calculate_virial_quantities", ["radius", "matter_mass"])
    hc.create()
    os.rename(os.path.join(cwd, 'catalog.0.h5'), os.path.join(cwd, run, fn, 'catalog.0.h5'))
