#
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
import numpy as np
import  matplotlib.pyplot as plt
from os.path import splitext
from xml.dom.minidom import parse

name2column = {'q4':0, 'q6':1, 'q8':2, 'q10':3, 'w4':4, 'w6':5, 'w8':6, 'w10':7}

def draw_cloud(filename, plane=('q6', 'w6'), cloudType=('raw', 'raw'), timeAveraged=False):
    """
Draw a bond orientational order 2D histogram for various subsets of particles :
    all, crystalline, liquid, MRC0

cloudType can be 'raw', 'cg' or 'surf'
    """
    base, ext = splitext(filename)
    base, digits = base.split('_t')
    #load data
    #coords = np.loadtxt(filename, skiprows=2)
    if timeAveraged:
        cloudExt = '.boo'
    else:
        cloudExt = '.cloud'
    spaceCloud = np.loadtxt('_space_t'.join([base, digits])+cloudExt)
    #keep only the particles having data (inside)
    inside2 = np.where(spaceCloud[:,name2column['q6']]>0)
    cloud = np.loadtxt('_t'.join([base, digits])+cloudExt)[inside2]
    spaceCloud = spaceCloud[inside2]
    surfCloud = np.loadtxt('_surf_t'.join([base, digits])+cloudExt)[inside2]
    #coords = coords[inside2]
    subsets = {
    'all':np.arange(cloud.shape[0]),
    'liquid':np.where(spaceCloud[:,name2column['q6']]<0.25)[0],
    'crystal':np.where(spaceCloud[:,name2column['q6']]>0.35)[0],
    'MRCO':np.where(np.bitwise_and(
            spaceCloud[:,name2column['q6']]<0.35,
            spaceCloud[:,name2column['q6']]>0.25
            ))[0],
        'FCC':np.where(np.bitwise_and(
            spaceCloud[:,name2column['q6']]>0.35,
            spaceCloud[:,name2column['w4']]<0
            ))[0],
        'HCP':np.where(np.bitwise_and(
            spaceCloud[:,name2column['q6']]>0.35,
            spaceCloud[:,name2column['w4']]>0
            ))[0]
    }
    cloudTypes = {
        'raw': cloud,
        'cg': spaceCloud,
        'surf': surfCloud
        }
    x = cloudTypes[cloudType[0]][:, name2column[plane[0]]]
    y = cloudTypes[cloudType[1]][:, name2column[plane[1]]]
    ranges = [[min([0, x.min()]), x.max()],[min([0, y.min()]), y.max()]]
    
    i=230
    for name, sel in subsets.iteritems():
        if(len(sel)>0):
            H, xedges, yedges = np.histogram2d(x[sel], y[sel], bins=[100,100], range=ranges)
        else:
            H=np.zeros((100,100))
        i+=1
        plt.subplot(i).set_aspect('equal')
        plt.imshow(
            np.log(np.rot90(H)/H.max()+1)*255,
            extent=ranges[0]+ranges[1],
            aspect='auto'
            )
        #plt.axis(ranges[0]+ranges[1])
        plt.xlabel(plane[0])
        plt.ylabel(plane[1])
        plt.title(name)
    plt.show()
