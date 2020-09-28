import os
import glob
import numpy  as np
import pandas as pd
import tables as tb
import h5py

from IC.invisible_cities.reco.psf_functions    import create_psf
from IC.invisible_cities.reco.psf_functions    import hdst_psf_processing
from IC.invisible_cities.reco.psf_functions    import add_empty_sensors_and_normalize_q
from IC.invisible_cities.reco.psf_functions    import add_variable_weighted_mean

import IC.invisible_cities.core.core_functions as     coref
import IC.invisible_cities.io.dst_io         as     dstio

from IC.invisible_cities.database              import load_db
from IC.invisible_cities.io.kdst_io      import psf_writer

from IC.invisible_cities.evm.nh5     import PSFfactors
from IC.invisible_cities.io.table_io import make_table

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"]          = 10, 10
plt.rcParams["font.size"]               = 16
plt.rcParams["figure.max_open_warning"] = 100

db = load_db.DataSiPM('next100', -1)

# Input and output path
#psf_path   = '/Users/halmamol/NEXT/files/NEXT100/psf/'
psf_path   = '/dipc/halmazan/files/psf/'
psf_filename = 'next100.kr83m.psf.h5'
out_psf  = psf_path + psf_filename

psf_files = glob.glob(psf_path+'psf_*.h5')

df_psf = pd.read_hdf(psf_files[0], 'PSF/PSFs')

evts    = []
factors = []

for i in range(1, 100):
    print(f'Running file {i} ')
    filein = psf_path + f'psf_{i}.h5'
    try:
        df = pd.read_hdf(filein, 'PSF/PSFs')
    except:
        #print(f'{filein} not found')
        continue
        
    if len(df.nevt.values) != 2440000:
        print(i, len(df.nevt.values))
    else:
    #if len(df.nevt.values) == 2440000:
        evts   .append(df.nevt  .values)
        factors.append(df.factor.values)

print('Taking array of events ')
evts_all    = np.array(evts)
factors_all = np.array(factors)

print('Estimating new events and factors ')
new_evts = evts_all.sum(axis=0).astype('int32')
tmp_factors = (evts_all * factors_all).sum(axis=0)
#new_factors = tmp_factors / new_evts
new_factors = tmp_factors / new_evts

print('Introducing them')
df_psf.nevt   = new_evts
df_psf.factor = new_factors

print(f'Creating file {out_psf}')
with tb.open_file(out_psf, 'w') as outfile:
    print(f'Taking array of events ')
    psf_table = make_table(outfile,
                           group       = "PSF",
                           name        = "PSFs",
                           fformat     = PSFfactors,
                           description = "XYZ dependent point spread functions",
                           compression = 'ZLIB4')

    row = psf_table.row
    for xr, yr, zr, x, y, z, f, ne in zip(df_psf.xr.values, df_psf.yr.values, df_psf.zr.values, 
                                          df_psf.x.values, df_psf.y.values, df_psf.z.values, 
                                          df_psf.factor.values, df_psf.nevt.values):
        row["xr"    ] = xr
        row["yr"    ] = yr
        row["zr"    ] = zr
        row["x"     ] = x
        row["y"     ] = y
        row["z"     ] = z
        row["factor"] = f
        row["nevt"  ] = ne
        row.append()
        
psf = pd.read_hdf(out_psf, 'PSF/PSFs')
print(psf.head())
