import yt
from yt_astro_analysis.halo_analysis.api import HaloCatalog
import os
import ytree
import numpy as np
import trident
import sys
import h5py as h5
"""
Run this third:

python halo_quantities3.py Natural
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

def read_rockstar_rvir(rockstar_data, ds):
    """
    Interpolate halo virial radius from rockstar merger tree
    """
    redshift = ds.current_redshift
    redshift_arr = rockstar_data['redshift_arr']
    rvir_arr = rockstar_data['rvir_arr']
    rvir = np.interp(redshift, redshift_arr, rvir_arr)
    rvir = ds.quan(rvir, 'kpccm').in_units('kpc')
    return rvir

def distance(x1,y1,z1, x2,y2,z2):
    return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5

@yt.particle_filter(requires=["particle_type"], filtered_type='all')
def stars(pfilter, data):
    filter = data[(pfilter.filtered_type, "particle_type")] == 2
    return filter

@yt.particle_filter(requires=["particle_type"], filtered_type='all')
def DM(pfilter, data):
    filter = data[(pfilter.filtered_type, "particle_type")] != 2
    return filter

if __name__ == '__main__':
    run = sys.argv[1]
    fn = sys.argv[2]
    #path = '/mnt/c/scratch/sciteam/chummels'
    path = '/Users/chummels/scratch/Tempest/new'
    cwd = os.getcwd()
    data_ds = yt.load(os.path.join(path, run, fn, fn))
    halo_ds = yt.load(os.path.join(cwd, run, fn, 'catalog.0.h5'))
    
    rockstar_data = get_rockstar_data(os.path.join(path,'tree_27.dat'), 27)
    c = read_rockstar_center(rockstar_data, data_ds)
    rvir = read_rockstar_rvir(rockstar_data, data_ds)
    f = open(os.path.join(cwd, run, fn, 'README'), 'w')
    ad = halo_ds.all_data()
    dists = distance(c[0], c[1], c[2], ad['particle_position_x'].convert_to_units('unitary'), ad['particle_position_y'].convert_to_units('unitary'), ad['particle_position_z'].convert_to_units('unitary'))
    winner = np.argmin(dists)
    center = [ad['particle_position_x'][winner].in_units('unitary'), ad['particle_position_y'][winner].in_units('unitary'), ad['particle_position_z'][winner].in_units('unitary')]

    r200 = ad['radius_200'][winner].to('kpc')
    m200 = ad['matter_mass_200'][winner].to('Msun')

    print("Original Halo Position and rvir")
    print(c)
    print(rvir)
    print('HaloID = %s' % ad['particle_identifier'][winner])
    print('distance = %g kpc' % dists[winner].in_units('kpc'))
    print('R200 = %g kpc' % r200)
    print('M200 = %g Msun' % m200)
    #print('Mvir = %g Msun' % ad['particle_mass'][winner].in_units('Msun'))
    #print('Rvir = %g kpc' % ad['virial_radius'][winner].in_units('kpc'))
    print('center = %s' % center)

    data_ds.add_particle_filter('stars')
    data_ds.add_particle_filter('DM')
    trident.add_ion_fields(data_ds, ['O VI'])
    sp = data_ds.sphere(center, r200)
    sp_core = data_ds.sphere(center, (10, 'kpc'))
    cgm = sp - sp_core
    m_tot = sp.quantities.total_mass()
    print(m_tot.to('Msun'))
    m_gas = np.sum(sp['cell_mass']).to('Msun')
    m_star = np.sum(sp[('stars', 'particle_mass')]).to('Msun')
    m_DM = np.sum(sp[('DM', 'particle_mass')]).to('Msun')
    m_CGM = np.sum(cgm[('gas', 'cell_mass')]).to('Msun')
    m_HI = np.sum(cgm[('gas', 'H_p0_mass')]).to('Msun')
    m_OVI = np.sum(cgm[('gas', 'O_p5_mass')]).to('Msun')
    cold = cgm.cut_region(["(obj['temperature'] < 1e4)"])
    cool = cgm.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e5)"])
    warm = cgm.cut_region(["(obj['temperature'] > 1e5) & (obj['temperature'] < 1e6)"])
    hot = cgm.cut_region(["(obj['temperature'] > 1e6)"])
    m_cold = np.sum(cold[('gas', 'cell_mass')]).to('Msun')
    m_cool = np.sum(cool[('gas', 'cell_mass')]).to('Msun')
    m_warm = np.sum(warm[('gas', 'cell_mass')]).to('Msun')
    m_hot = np.sum(hot[('gas', 'cell_mass')]).to('Msun')
    m_tot = m_gas+m_star+m_DM
    print("m_gas = %g Msun" % m_gas)
    print("m_star = %g Msun" % m_star)
    print("m_DM = %g Msun" % m_DM)
    print("m_tot = %g Msun" % m_tot)
    print("m_CGM = %g Msun" % m_CGM)
    print("m_cold = %g Msun" % m_cold)
    print("m_cool = %g Msun" % m_cool)
    print("m_warm = %g Msun" % m_warm)
    print("m_hot = %g Msun" % m_hot)
    print("m_HI = %g Msun" % m_HI)
    print("m_OVI = %g Msun" % m_OVI)

    f.write("Original Halo Position and rvir\n")
    f.write("%s\n" % c)
    f.write("%s\n" % rvir)
    f.write('HaloID = %s\n' % ad['particle_identifier'][winner])
    f.write('distance = %g kpc\n' % dists[winner].to('kpc'))
    f.write('R200 = %g kpc\n' % r200)
    f.write('M200 = %g Msun\n' % m200)
    #f.write('Mvir = %g Msun\n' % ad['particle_mass'][winner].to('Msun'))
    #f.write('Rvir = %g kpc\n' % ad['virial_radius'][winner].to('kpc'))
    f.write('center = %s\n' % center)

    f.write("m_gas = %g Msun\n" % m_gas)
    f.write("m_star = %g Msun\n" % m_star)
    f.write("m_DM = %g Msun\n" % m_DM)
    f.write("m_tot = %g Msun\n" % m_tot)
    f.write("m_CGM = %g Msun\n" % m_CGM)
    f.write("m_cold = %g Msun\n" % m_cold)
    f.write("m_cool = %g Msun\n" % m_cool)
    f.write("m_warm = %g Msun\n" % m_warm)
    f.write("m_hot = %g Msun\n" % m_hot)
    f.write("m_HI = %g Msun\n" % m_HI)
    f.write("m_OVI = %g Msun\n" % m_OVI)
    f.close()

    d = [r200, m200, m_gas, m_star, m_DM, m_tot, m_CGM, m_cold, m_cool, m_warm, m_hot, m_HI, m_OVI]
    f = h5.File(os.path.join(cwd, run, fn, 'halo.h5'), 'w')
    f.create_dataset('d', data=d)
    f.close()

    #for ax in 'xyz':
    #    p = yt.ProjectionPlot(data_ds, ax, 'H_number_density', center=center, width=(200, 'kpc'), data_source=sp)
    #    p.set_zlim('H_number_density', 1e11,1e23)
    #    p.annotate_sphere(center, (10, 'kpc'), circle_args={'color':'white', 'alpha':0.5, 'linestyle':'dashed', 'linewidth':2})
    #    p.annotate_marker(center, coord_system='data')
    #    p.save("%s/%s/" % (run, fn))
    #prof = yt.ProfilePlot(sp, "radius", ["cell_mass", 'H_p0_mass', 'O_p5_mass'],
    #                      weight_field=None,
    #                      accumulation=False)
    #prof.set_unit("radius", "kpc")
    #prof.set_log('radius', False)
    #prof.set_unit("cell_mass", "Msun")
    #prof.set_unit("H_p0_mass", "Msun")
    #prof.set_unit("O_p5_mass", "Msun")
    #prof.set_ylim('H_p0_mass', 1e1, 1e9)
    #prof.set_ylim('O_p5_mass', 1e1, 1e5)
    #prof.save("%s/" % run)
    #prof = yt.ProfilePlot(cgm, "temperature", ["cell_mass"])
    #prof.set_xlim(5e3, 5e7)
    #prof.set_unit("cell_mass", "Msun")
    #prof.set_ylim('cell_mass', 1e2, 1e7)
    #prof.save("%s/" % run)
    #phas = yt.PhasePlot(cgm, "density", "temperature", "cell_mass", weight_field=None)
    #phas.set_xlim(1e-30, 1e-23)
    #phas.set_ylim(1e4, 1e8)
    #phas.set_unit('cell_mass', 'Msun')
    #phas.set_zlim('cell_mass', 1e2, 3e7)
    #phas.save("%s/" % run)
