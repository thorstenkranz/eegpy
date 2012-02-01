#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for plotting topography-maps and connectivity-diagrams 
on a schematic head surface.


"""
__docformat__ = "restructuredtext"

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
except ImportError:
    raise FATALERROR('SciPy or NumPy not found!\nPlease visit www.scipy.org or numeric.scipy.org for more information.')

from matplotlib.mlab import griddata
#from griddata import griddata



def plot_head_topography(topography, sensorlocations, plotsensors=False,
                       resolution=51, masked=True, plothead=True,
                       plothead_kwargs=None,  ax=P, view='top', **kwargs):
    """Plot distribution to a head surface, derived from some sensor locations.

The sensor locations are first projected onto the best fitting sphere and
finally projected onto a circle (by simply ignoring the z-axis).

Parameters:
    topography: array
        A vector of some values corresponding to each sensor.
    sensorlocations: (nsensors x 3) array
        3D coordinates of each sensor. The order of the sensors has to match
        with the `topography` vector.
    plotsensors: bool
        If True, sensor will be plotted on their projected coordinates.
        No sensor are shown otherwise.
    plothead: bool
        If True, a head outline is plotted.
    plothead_kwargs: dict
        Additional keyword arguments passed to `plot_head_outline()`.
    resolution: int
        Number of surface samples along both x and y-axis.
    masked: bool
        If True, all surface sample extending to head outline will be
        masked.
    ax: mpl axes
        axes to plot to. Standard is pylab.
    view: one of 'top' and 'rear'
        Defines from where the head is viewed.
    kwargs:
        All additional arguments will be passed to `P.imshow()`.


Returns:
  (map, head, sensors)
    The corresponding matplotlib objects are returned if plotted, ie.
    if plothead is set to `False`, head will be `None`.

      map
        The colormap that makes the actual plot, a
        matplotlib.image.AxesImage instance.
      head
        What is returned by `plot_head_outline()`.
      sensors
        The dots marking the electrodes, a matplotlib.lines.Line2d
        instance.
    """
    # give sane defaults
    if plothead_kwargs is None:
        plothead_kwargs = {}
    
    #assert sensorlocations is an numpy-arrays and swap x and y coordinates
    #swapping is necessary as the eeglab .ced files have another interpretation of this
    sensorlocations = N.array(sensorlocations)
    tmp = sensorlocations[:,1].copy()
    sensorlocations[:,1] = sensorlocations[:,0]
    sensorlocations[:,0]=tmp[:]
    sensorlocations[:,0]*=-1
    del tmp


    # error function to fit the sensor locations to a sphere
    def err(params):
        r, cx, cy, cz = params
        #print r,cx,cy,cz
        rv = (sensorlocations[:, 0] - cx) ** 2 + (sensorlocations[:, 1] - cy) ** 2 + (sensorlocations[:, 2] - cz) ** 2 - r ** 2
        rv = abs(rv.sum())
        #print "rv: ",rv
        return rv
    

    # initial guess of sphere parameters (radius and center)
    params = N.array([1.0, 0.0, 0.0, 0.0])

    # do fit
    r, cx, cy, cz = fmin(err,params,disp=0)#leastsq(err, params)#
    #print "Results of fit:", r, cx,cy,cz

    # project the sensor locations onto the sphere
    sphere_center = N.array((cx, cy, cz))
    sproj = sensorlocations - sphere_center
    sproj = r * sproj / N.c_[N.sqrt(N.sum(sproj ** 2, axis=1))]
    sproj += sphere_center
    #print "sproj.shape:",sproj.shape

    if view == 'top':
        # Generate a grid and interpolate using the griddata module
        x = N.arange(cx - r, cx + r, 2. * r / resolution)
        y = N.arange(cy - r, cy + r, 2. * r / resolution)
        x, y = P.meshgrid(x, y)
        # fit topology onto xy projection of sphere
        topo = griddata(sproj[:, 0], sproj[:, 1],\
                N.ravel(N.array(topography)), x, y)
        # mask values outside the head
        if masked:
            notinhead = N.greater_equal((x - cx) ** 2 + (y - cy) ** 2,
                                        (1.0 * r) ** 2)
            topo = N.ma.masked_where(notinhead, topo)
    elif view=='rear':
        # Generate a grid and interpolate using the griddata module
        x = N.arange(cx - r, cx + r, 2. * r / resolution)
        z = N.arange(cz - r, cz + r, 2. * r / resolution)
        x, z = P.meshgrid(x, z)
        # fit topology onto xz projection of sphere and discard all frontal electrodes
        topography = N.ravel(N.array(topography))
        topography = topography[sproj[:,1]<cy]
        sproj = sproj[sproj[:,1]<cy]
        #print "sproj.shape:",sproj.shape
        topo = griddata(sproj[:, 0], sproj[:, 2],topography, x, z)
        # mask values outside the head
        if masked:
            notinhead = N.greater_equal((x - cx) ** 2 + (z - cz) ** 2,
                                        (1.0 * r) ** 2)
            topo = N.ma.masked_where(notinhead, topo)
    else:
        raise ValueError("view must be one of 'top' and 'rear'")

    # show surface
    map = ax.imshow(topo, origin="lower", extent=(-r, r, -r, r), **kwargs)
    ax.axis('off')

    if plothead:
        # plot scaled head outline
        head = plot_head_outline(scale=r, shift=(cx, cy), view=view, **plothead_kwargs)
    else:
        head = None

    if plotsensors:
        sensors = plot_sensors(sproj,ax,"wo",view=view)
    else:
        sensors = None
    
    if view == 'top':
        ax.xlim((cx-(r*1.2),cx+(r*1.2)))
        ax.ylim((cy-(r*1.2),cy+(r*1.2)))
    elif view=='rear':
        ax.xlim((cx-(r*1.2),cx+(r*1.2)))
        ax.ylim((cz-(r*0.4),cz+(r*1.2)))
    
    return map, head, sensors

def plot_head_connectivity(connectivity, sensorlocations, plotsensors=False,
                       vmin=None, vmax=None, cm=cm.jet, plothead=True,
                       plothead_kwargs=None,  ax=P, view='top', **kwargs):
    """Plot connectivity on a head surface, derived from some sensor locations.

        The sensor locations are first projected onto the best fitting sphere and
        finally projected onto a circle (by simply ignoring the z-axis).

        :param connectivity: Connectivity matrix
        :type connectivity: matrix
        :param sensorlocations: array (nsensors x 3), 3D coordinates of each sensor. The order of the sensors has to match with the `connectivity` matrix.
        :param plotsensors: bool; if True, sensor will be plotted on their projected coordinates. No sensors are shown otherwise.
        :param plothead: bool; If True, a head outline is plotted.
        :param plothead_kwargs: Additional keyword arguments passed to `plot_head_outline()`.
        :param vmin,vmax: Minimum and maximum value to be used for graphics.
        :param ax: matplotlib axes to plot to. Standard is pylab.
        :param view: one of 'top' and 'rear'; Defines from where the head is viewed.
        :param kwargs: All additional arguments will be passed to `P.imshow()`.
        :returns: (map, head, sensors)

        The corresponding matplotlib objects are returned if plotted, i.e.,
        if plothead is set to `False`, head will be `None`.

          map
            The colormap that makes the actual plot, a
            matplotlib.image.AxesImage instance.
          head
            What is returned by :py:meth:`plot_head_outline`.
          sensors
            The dots marking the electrodes, a matplotlib.lines.Line2d
            instance.

        ..seealso: :py:meth:`plot_head_topography`
    """
    #Some assertions:
    assert len(connectivity.shape)==2 and connectivity.shape[0]==connectivity.shape[1], "connectivity must be a quadratic matrix"
    assert connectivity.shape[0] == sensorlocations.shape[0], "connectivity and sensorlocations must have same length"
    if vmin!=None and vmax != None:
        assert  vmin<vmax, "vmin(=%f) must be smaller than vmax(=%f)" % (vmin,vmax)
    # give sane defaults
    if plothead_kwargs is None:
        plothead_kwargs = {}
    
    #assert sensorlocations is an numpy-arrays and swap x and y coordinates
    #swapping is necessary as the eeglab .ced files have another interpretation of this
    sensorlocations = N.array(sensorlocations)
    tmp = sensorlocations[:,1].copy()
    sensorlocations[:,1] = sensorlocations[:,0]
    sensorlocations[:,0]=tmp[:]
    sensorlocations[:,0]*=-1
    del tmp


    # error function to fit the sensor locations to a sphere
    def err(params):
        r, cx, cy, cz = params
        #print r,cx,cy,cz
        rv = (sensorlocations[:, 0] - cx) ** 2 + (sensorlocations[:, 1] - cy) ** 2 + (sensorlocations[:, 2] - cz) ** 2 - r ** 2
        rv = abs(rv.sum())
        #print "rv: ",rv
        return rv
    

    # initial guess of sphere parameters (radius and center)
    params = N.array([1.0, 0.0, 0.0, 0.0])

    # do fit
    r, cx, cy, cz = fmin(err,params,disp=0)#leastsq(err, params)#
    #print "Results of fit:", r, cx,cy,cz

    # project the sensor locations onto the sphere
    sphere_center = N.array((cx, cy, cz))
    sproj = sensorlocations - sphere_center
    sproj = r * sproj / N.c_[N.sqrt(N.sum(sproj ** 2, axis=1))]
    sproj += sphere_center
    #print "sproj.shape:",sproj.shape
    
    #vmin, vmax: give sane defaults
    #first, make copy of connectivity  and set diagonal elements to zero
    conn = connectivity.copy()
    conn *= N.ones((conn.shape[0]),"d")-N.diag(N.ones((conn.shape[0]),"d"))
    if vmin==None:
        vmin=conn.min()
    if vmax==None:
        vmax=conn.max()
    #Now transform values of conn to be between 0 and 1
    conn = (conn-vmin)/(vmax-vmin)
    conn[conn>1] = N.ones(conn[conn>1].shape,"d") 
    conn[conn<0] = N.zeros(conn[conn<0].shape,"d")
    
    if view == 'top':
        #
        fig = ax.gca()
        for i in range(connectivity.shape[0]):
            for j in range(i):
                dx=sproj[j,0]-sproj[i,0]
                dy=sproj[j,1]-sproj[i,1]
                if conn[i,j] != conn[j,i]:
                    if conn[i,j]>0.0:
                        #ax.arrow(sproj[i,0],sproj[i,1],dx,dy,lw=conn[i,j]*5+1,ec=cm(conn[i,j]),head_width=N.sqrt(dx**2+dy**2)*conn[i,j]*0.03,zorder=100-conn[i,j])
                        arr1 = P.Arrow(sproj[i,0], sproj[i,1], dx, dy, width=(conn[i,j]*5+1)/30,ec=cm(conn[i,j]),fc=cm(conn[i,j]),zorder=100+conn[i,j])
                        fig.add_patch(arr1)
                    if conn[j,i]>0.0:
                        #ax.arrow(sproj[i,0],sproj[i,1],dx,dy,lw=conn[j,i]*5+1,ec=cm(conn[j,i]),head_width=N.sqrt(dx**2+dy**2)*conn[j,i]*0.03,zorder=100-conn[j,i])
                        arr1 = P.Arrow(sproj[i,0], sproj[i,1], dx, dy, width=(conn[j,i]*5+1)/30,ec=cm(conn[j,i]),fc=cm(conn[j,i]),zorder=100+conn[j,i])
                        fig.add_patch(arr1)
                else:
                    if conn[i,j]>0.0:
                        ax.arrow(sproj[i,0],sproj[i,1],sproj[j,0]-sproj[i,0],sproj[j,1]-sproj[i,1],lw=conn[i,j]*5+1,ec=cm(conn[i,j]),zorder=100-conn[i,j])
    #elif view=='rear':
    #    pass
    else:
        raise ValueError("view must be one of 'top' and 'rear'")

    # show surface
    #map = ax.imshow(topo, origin="lower", extent=(-r, r, -r, r), **kwargs)
    ax.axis('off')
    ax.axis('equal')

    if plothead:
        # plot scaled head outline
        head = plot_head_outline(scale=r, shift=(cx, cy), view=view, **plothead_kwargs)
    else:
        head = None

    if plotsensors:
        sensors = plot_sensors(sproj,ax,"wo",view=view)
    else:
        sensors = None
    
    if view == 'top':
        ax.xlim((cx-(r*1.2),cx+(r*1.2)))
        ax.ylim((cy-(r*1.2),cy+(r*1.2)))
    elif view=='rear':
        ax.xlim((cx-(r*1.2),cx+(r*1.2)))
        ax.ylim((cz-(r*0.4),cz+(r*1.2)))
    
    return map, head, sensors

def plot_connectivity_legend(num_lines=5, vmin=0, vmax=1,ax=P,cm=cm.jet):
    """Plot a legend that explains arrows in connectivity plot"""
    #TODO: Make it handle symmetric measures as well
    vals = N.linspace(vmin,vmax,num_lines)
    pvals = N.linspace(0,1,num_lines)
    fig = ax.gca()
    fig.axis('off')
    for i in range(num_lines):
        #ax.xlim([0,1])
        ax.ylim([0,1])
        ax.xlim([0,1])
        arr1 = ax.Arrow(0.1, i*0.1 ,0.3,0, width=(pvals[i]*5+1)/50,ec=cm(pvals[i]),fc=cm(pvals[i]))
        fig.add_patch(arr1)
        ax.text(0.5,i*0.1,"%.3f"%vals[i], va="center")
        
        
def plot_head_outline(scale=1, shift=(0, 0), color='k', linewidth='3', ax=P, view='top', **kwargs):
    """Plots a simple outline of a head viewed from the top.

    The plot contains schematic representations of the nose and ears. The
    size of the head is basically a unit circle for nose and ears attached
    to it.

    :Parameters:
      scale: float
        Factor to scale the size of the head.
      shift: 2-tuple of floats
        Shift the center of the head circle by these values.
      color: matplotlib color spec
        The color the outline should be plotted in.
      linewidth: int
        Linewidth of the head outline.
      ax: mpl axes
        axes to plot to. Standard is pylab.
      view: one of 'top' and 'rear'
        Defines from where the head is viewed.
      kwargs:
        All additional arguments are passed to `P.plot()`.

    :Returns:
      Matplotlib lines2D object
        can be used to tweak the look of the head outline.
    """

    rmax = 0.5
    # factor used all the time
    fac = 2 * N.pi * 0.01
    
    if view=='top':
        # Koordinates for the ears
        EarX1 =  -1 * N.array(
                [.497, .510, .518, .5299,
                .5419, .54, .547, .532, .510,
                rmax * N.cos(fac * (54 + 42))])
        EarY1 = N.array(
                [.0655, .0775, .0783, .0746, .0555,
                -.0055, -.0932, -.1313, -.1384,
                rmax * N.sin(fac * (54 + 42))])
        EarX2 = N.array(
                [rmax * N.cos(fac * (54 + 42)),
                .510, .532, .547, .54, .5419,
                .5299, .518, .510, .497] )
        EarY2 = N.array(
                [rmax * N.sin(fac * (54 + 42)),
                -.1384, -.1313, -.0932, -.0055,
                .0555, .0746, .0783, .0775, .0655] )
    
        # Coordinates for the Head
        HeadX1 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 2)), (21,))
        HeadY1 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 2)), (21,))
        HeadX2 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 28)), (21,))
        HeadY2 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 28)), (21,))
        HeadX3 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 54)), (43,))
        HeadY3 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 54)), (43,))
    
        # Coordinates for the Nose
        NoseX = N.array([.18 * rmax, 0, -.18 * rmax])
        NoseY = N.array([rmax - 0.004, rmax * 1.15, rmax - 0.004])
    
        # Combine to one
        X = N.concatenate((EarX2, HeadX1, NoseX, HeadX2, EarX1, HeadX3))
        Y = N.concatenate((EarY2, HeadY1, NoseY, HeadY2, EarY1, HeadY3))
    
        X *= 2 * scale
        Y *= 2 * scale
        X += shift[0]
        Y += shift[1]
    elif view=='rear':
        # Koordinates for the ears
        EarX1 =  -1 * N.array(
                [.497, .510, .518, .5299,
                .5419, .54, .537, .525, .510,
                rmax * N.cos(fac * (54 + 42))])
        EarY1 = N.array(
                [.0655, .0775, .0783, .0746, .0555,
                -.0055, -.0932, -.1713, -.1784,
                rmax * N.sin(fac * (54 + 42))])
        EarX2 = N.array(
                [rmax * N.cos(fac * (54 + 42)),
                .510, .525, .537, .54, .5419,
                .5299, .518, .510, .497] )
        EarY2 = N.array(
                [rmax * N.sin(fac * (54 + 42)),
                -.1784, -.1713, -.0932, -.0055,
                .0555, .0746, .0783, .0775, .0655] )
    
        # Coordinates for the Head
        HeadX1 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 2)), (23,))
        HeadY1 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 2)), (23,))
        HeadX2 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 26)), (23,))
        HeadY2 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 26)), (23,))
        HeadX3 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 54)), (17,))
        HeadY3 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 54)), (17,))
        HeadX4 = N.fromfunction(
                lambda x: rmax * N.cos(fac * (x + 80)), (17,))
        HeadY4 = N.fromfunction(
                lambda y: rmax * N.sin(fac * (y + 80)), (17,))
    
        # Coordinates for the Neck
        NeckX = N.array([rmax * N.cos(fac * 70),rmax * N.cos(fac * 80)])
        NeckY = N.array([rmax * -1.1, rmax * -1.1])
    
        # Combine to one
        X = N.concatenate((EarX2, HeadX1, HeadX2, EarX1, HeadX3, NeckX, HeadX4))
        Y = N.concatenate((EarY2, HeadY1, HeadY2, EarY1, HeadY3, NeckY, HeadY4))
    
        X *= 2 * scale
        Y *= 2 * scale
        X += shift[0]
        Y += shift[1]
    else:
        raise ValueError("view must be one of 'top' and 'rear'")

    return ax.plot(X, Y, color=color, linewidth=linewidth)

def plot_sensors(sproj,ax,format = 'wo',view='top'):
    # plot projected sensor locations
    # reorder sensors so the ones below plotted first
    # TODO: please fix with more elegant solution
    zenum = [x[::-1] for x in enumerate(sproj[:, 2].tolist())]
    zenum.sort()
    indx = [ x[1] for x in zenum ]
    if view=='top':
        return ax.plot(sproj[indx, 0], sproj[indx, 1], format,zorder=10e50)
    elif view=='rear':
        return ax.plot(sproj[indx, 0], sproj[indx, 2], format,zorder=10e50)
    else:
        raise ValueError("view must be one of 'top' and 'rear'")
    

if __name__=="__main__":
    vals = N.array([9,4,7])
    #slocs = N.array([[1.5,0.5,0.5],[-0.5,0.5,0.5],[0.5,-0.5,0.5]])
    #plot_head_topography(vals, slocs, plotsensors=False, masked=True, plothead=True,view='rear')
    from eegpy.plot.sensorlocations import SensorLocations
    sloc = SensorLocations("sample_locs/Standard-10-20-Cap19.ced")
    sensors = sloc.get_all_cartesian3d()
    connectivity = N.random.standard_normal((sensors.shape[0],sensors.shape[0]))
    print connectivity.min(), connectivity.max()
    plot_head_connectivity(connectivity, sensors, plotsensors=True, plothead=True, vmin=0)
    #plot_head_outline(view='rear')
    P.show()
    
