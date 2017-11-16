import sys, os, os.path, argparse, shutil
import numpy as np
from pandas import DataFrame
import trackpy as tp
from colloids import experiment as xp
from colloids.progressbar import ProgressBar
from colloids.particles import Linker

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a HDF5 format compatible with trackpy.PandasHDFStoreSingleNode to coordinates and trajectory data.')
    parser.add_argument('h5name', help='HDF5 file name')
    parser.add_argument('--out', help='The folder where to save the data. By default the same name as HDF5 file, without extension.')
    args = parser.parse_args()
    h5name = args.h5name
    print(h5name)
    
    if args.out is None:
        args.out = os.path.splitext(h5name)[0]
    print('to %s'%args.out)
    os.mkdir(args.out)
    basename = os.path.join(args.out, os.path.split(os.path.basename(h5name))[0])
    
    with tp.PandasHDFStoreSingleNode(h5name) as s:
        nbzeros = len('%d'%(len(s)-1))
        datfile = +'_t%0{:d}d.dat'.format(nbzeros)
        p2trfile = +'_t%0{:d}d.p2tr'.format(nbzeros)
        pro = ProgressBar(len(s))
        for t, frame in s:
            pos = np.vstack((np.zeros(2,3), frame.as_matrix(['x','y','z'])))
            pos[0] = 1, len(pos)-2, 1
            np.savetxt(datfile%t, pos, fmt='%g')
            np.savetxt(p2trfile%t, frame.asmatrix(['particle'], fmt='%d')
            #if t==0:
            #    linker = Linker(len(frame))
            #else:
            #    linker.loadFrame(frame.asmatrix(['particle'])
            pro.animate(t)
        #TODO save the linker to a traj file

