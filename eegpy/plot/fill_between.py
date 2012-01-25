#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for plotting topography-maps of values distributed over the head surface.

"""
#################
# Module-Import #
#################

#eegpy-modules
import eegpy
from eegpy.misc import FATALERROR, debug

#Third-party
try:
    import numpy as N
    from scipy.optimize import leastsq, fmin#_bfgs#_powell
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

try:
    import pylab as P
    from matplotlib import cm
    from matplotlib.patches import Polygon
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

from matplotlib.mlab import griddata
#from griddata import griddata



def fill_between(x, y1, y2, ax=P.gca(), alpha=0.5, **kwargs):
    """Plot distribution to a head surface, derived from some sensor locations.

    The sensor locations are first projected onto the best fitting sphere and
    finally projected onto a circle (by simply ignoring the z-axis).

    :Parameters:
      x:
      ,y1,y2
      ax: mpl axes
        axes to plot to. Standard is pylab.
      view: one of 'top' and 'rear'
        Defines from where the head is viewed.
      **kwargs:
        All additional arguments will be passed to `pylab.imshow()`.


    :Returns:
      (map, head, sensors)
        The corresponding matplotlib objects are returned if plotted, ie.
        if plothead is set to `False`, `head` will be `None`.

          map
            The colormap that makes the actual plot, a
            matplotlib.image.AxesImage instance.
          head
            What is returned by `plot_head_outline()`.
          sensors
            The dots marking the electrodes, a matplotlib.lines.Line2d
            instance.
    """
    # add x,y2 in reverse order for proper polygon filling
    verts = zip(x,y1) + [(x[i], y2[i]) for i in range(len(x)-1,-1,-1)]
    poly = Polygon(verts, **kwargs)
    poly.set_alpha(alpha)
    ax.add_patch(poly)
    ax.autoscale_view()
    return poly




if __name__=="__main__":
    x = N.arange(0, 2, 0.01)
    y1 = N.sin(2*N.pi*x)
    y2 = N.sin(4*N.pi*x) + 2
    ax = P.gca()
    
    poly = fill_between(x, y1, y2, facecolor="g")
    poly.set_alpha(0.5)
    P.show()
    
