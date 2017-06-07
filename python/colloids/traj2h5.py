import sys, os.path, argparse, shutil
import numpy as np
from pandas import DataFrame
import trackpy as tp
from colloids import experiment as xp
from colloids.progressbar import ProgressBar

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert coordinates and trajectory data to a HDF5 format compatible with trackpy.PandasHDFStoreSingleNode.')
    parser.add_argument('trname', help='Trajectory file name')
    parser.add_argument('--out', help='The file.h5 where to save the data. By default the same name as trajectory file, one level below in the directory tree.')
    args = parser.parse_args()
    trname = args.trname
    print(trname)
    x = xp.Experiment(trname)
    
    if args.out is None:
        args.out = os.path.join(
            os.path.split(x.path)[0], 
            os.path.splitext(x.trajfile)[0]+'.h5'
            )
    print('to %s'%args.out)
    assert x.size < 2**15, "Too many time steps, use larger data types"
    assert x.nb_trajs < 2**31, "Too many trajectories, use larger data types"
    
    pro = ProgressBar(x.size)
    with tp.PandasHDFStoreSingleNode(args.out) as s:
         for t, name in x.enum():
            pro.animate(t)
            d = {x:c for x,c in zip(
                'xyz', 
                np.loadtxt(name, skiprows=2, unpack=True, dtype=np.float16)
                )}
            d['frame'] = np.full(d['x'].shape, t, dtype=np.int16)
            d['particle'] = x.p2tr(t).astype(np.int32)
            s.put(DataFrame(d))
    pro.animate(x.size)
