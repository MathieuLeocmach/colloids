import numpy as np
import numexpr

def lab2rgb(im):
    l, a, b = np.dsplit(im, 3)
    y = numexpr.evaluate('(l + 16) / 116.0')
    x = numexpr.evaluate('a/500.0 + y')
    z = numexpr.evaluate('y - b/200.0')
    xyz = numexpr.evaluate(
        'where(y**3>0.008856, y**3, (y - 16.0 / 116.0) / 7.787)*ref/100',
        {
            'y':np.dstack((x,y,z)), 
            'ref': np.array([[[95.047, 100, 108.883]]]), #Observer= 2 deg Illuminant= D65
            }
        )
    rgb = np.inner(xyz, np.array([
        [3.2406, -1.5372, -0.4986],
        [-0.9689, 1.8758, 0.0415],
        [0.0557, -0.2040, 1.0570]
        ]))
    rgb = numexpr.evaluate('where(rgb > 0.0031308, 1.055 * rgb**(1/2.4) - 0.055, 12.92 * rgb)')
    #clip colors
    return np.minimum(np.maximum(rgb, 0), 1)

def msh2rgb(im):
    m, s, h = np.dsplit(im, 3)
    l = numexpr.evaluate('m * cos(s)')
    a = numexpr.evaluate('m * sin(s) * cos(h)')
    b = numexpr.evaluate('m * sin(s) * sin(h)')
    return lab2rgb(np.dstack((l,a,b)))
    
def colorscale(im, centre=88):
    """Take real values from 0 to 1 and return rgb colors on the diverging color scale blue-white-red. Taken from Paraview."""
    x = np.mod(im, 1)
    m = numexpr.evaluate('2*(80-centre)*abs(x-0.5)+centre')
    #m = numexpr.evaluate('88-16*abs(x-0.5)')
    s = numexpr.evaluate('1.08 * 2*abs(0.5-x)')
    h = numexpr.evaluate('where(x<0.5, -1.1 - 2*x*0.561, 1.061 - 2*(x-0.5)*0.561)')
    return msh2rgb(np.dstack((m,s,h)))
    
