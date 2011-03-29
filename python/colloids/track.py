#    Copyright 2009 Mathieu Leocmach
#
#    This file is part of Colloids.
#
#    Colloids is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Colloids is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Colloids.  If not, see <http://www.gnu.org/licenses/>.
#
from __future__ import with_statement #for python 2.5, useless in 2.6
import numpy as np
import os.path, subprocess, shlex
from colloids import lif, vtk

def bestParams(inputPath, outputPath, radMins=np.arange(2.5, 5, 0.5), radMaxs=np.arange(6, 32, 4), t=0, serie=None):
    if serie==None:
        serie = lif.Reader(inputPath).chooseSerieIndex()
    s = lif.Reader(inputPath).getSeries()[serie]
    out_head = outputPath+'_thr0_radMin%(m)g_radMax%(M)g_'
    if s.getNbFrames>1:
        dt = s.getTimeLapse()
        dt_min, dt_sec = divmod(dt, 60)
        dt_sec, dt_msec = divmod(dt_sec, 1)
        if dt_min > 0:
            out_head += '%dmin_'%dt_min
            dt_msec=0
        if dt_sec > 1:
            out_head += '%ds'%dt_sec
            if dt_msec > 0.01:
                out_head += '%d'%np.floor(dt_msec*100)
            out_head +='_'
    out_noext = out_head+'t%(t)0'+('%d'%len('%d'%s.getNbFrames()))+'d'
    for radMin in radMins:
	for radMax in radMaxs:
            if not os.path.exists(
                (out_noext+'.intensity')%{'m':radMin, 'M':radMax, 't':t}
                ):
                p = subprocess.Popen(
                    shlex.split(
                        (
                            ('/home/mathieu/test/bin/tracker -L -i %(lif)s -o %(out)s --serie %(serie)d'%{
                                'lif':inputPath,
                                'out':out_head,
                                'serie':serie,
                                }) +
                            ' --radiusMin %(m)g --radiusMax %(M)g --threshold 0 --onlyTimeStep %(t)d --exportIntensities --fftPlanning 32'
                        )%{
                            'm':radMin,
                            'M':radMax, 't':t
                            }
                        ),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                    )
                out, err = p.communicate()
                if not p.returncode==0:
                    print 'At radMin=%(m)g and radMax=%(M)g'%{'m':radMin, 'M':radMax}
                    print out
                    print err
                    return
    hists = np.asarray([[
	np.histogram(
		np.loadtxt(
                    (out_noext+'.intensity')%{'m':radMin, 'M':radMax, 't':t}
                    ),
		bins=range(300)
		)[0]
	for radMax in radMaxs] for radMin in radMins])
    return hists

def getBandPassMask(serie, radMin, radMax):
    maskshape = serie.getFrameShape()
    maskshape[-1] = maskshape[-1]/2 + 1
    mask = np.zeros(maskshape, dtype=bool)
    freqs = np.repeat(np.expand_dims(serie.getFrameShape(),1),2, axis=1)
    freqs[:,0] /= float(2*radMin)
    freqs[:,1] /= float(2*radMax)
    freqs[-1] *= serie.getZXratio()
    xs = (np.repeat(np.arange(maskshape[0]/2, dtype=float)[:, np.newaxis],2, axis=1)/freqs[0])**2
    ys = (np.repeat(np.arange(maskshape[1]/2, dtype=float)[:, np.newaxis],2, axis=1)/freqs[1])**2
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            z = freqs[-1] * np.where(x+y<1, np.sqrt(1 - x - y), 0)
            mask[i, j, z[1]:z[0]] = True
            mask[-i, j, z[1]:z[0]] = True
            mask[i, -j, z[1]:z[0]] = True
            mask[-i, -j, z[1]:z[0]] = True
    return mask

def saveIntensityVTK(fname):
	v=vtk.Polydata()
	v.points = np.loadtxt(fname, skiprows=2)
	v.scalars = [('Intensity', np.loadtxt(os.path.splitext(fname)[0]+'.intensity'))]
	n = fname.split('_t')
	n[-2]+='_test'
	n[-1] = n[-1][:-3]+'vtk'
	v.save('_t'.join(n))
