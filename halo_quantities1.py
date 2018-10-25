import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
import os
import sys
"""
Run this first:

python halo_quantities1.py Natural DD0274

(or python halo_quantities1.py Tempest8 DD0274)

"""

run = sys.argv[1]
fn = sys.argv[2]
cwd = os.getcwd()
path = '/Users/chummels/scratch/Tempest/new'
full_fn = os.path.join(path, run, fn, fn)
data_ds = yt.load(full_fn)
hc = HaloCatalog(data_ds=data_ds, finder_method='hop', finder_kwargs={"dm_only":False})
hc.create()
dest_dir = os.path.join(cwd, run, fn, 'catalog')
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)
os.rename(os.path.join(cwd, 'halo_catalogs/catalog/catalog.0.h5'), os.path.join(cwd, run, fn, 'catalog', 'catalog.0.h5'))
