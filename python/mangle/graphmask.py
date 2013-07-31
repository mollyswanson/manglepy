"""
module for plotting mangle polygon masks using matplotlib's basemap
geographic mapping module. Contains the class:`Skymap` class which
extends the class:`Basemap` class to make it cleaner
for use with astronomical data and the following functions:

:func:`main`: wrapper for plot_mangle_map that allows this script to
be called from the unix command line

:func:`plot_mangle_map`: creates Skymap and draw polygons in it

:func:`polysplit`: takes a list of nan-separated polygons (as in
mangle's .list format) and splits into separate polygons

:func:`read_poly_list`: read in a graphics polygon file in
mangle's .list format

:func:`draw_weighted_polygons`: function that draws the polygons
using matplotlib's PolyCollection (works w/o Basemap)

:func:`ispolycw`: determines whether the polygon points are listed in a clockwise
order (when plotted with azimuth increasing to the left as when looking at the sky)

:func:`ispolyccw`: determines whether the polygon points are listed in a counterclockwise
order (when plotted with azimuth increasing to the left as when looking at the sky)

:func:`expecting`: generic function to mimic matlab's nargout - counts the
number of output arguments expected from within a function

:func:`get_graphics_polygons`: takes a set of mangle polygons and returns
sets of points along the edge of each polygon for plotting

:func:`trim_mask`: trims the input polygons to only lie within
a given region
"""

import matplotlib
#if calling from command line, use a backend without a gui so the image is written straight to a file
if __name__=='__main__':
    matplotlib.use('Agg')    
from pylab import *
import subprocess
import os
import tempfile
import inspect
from matplotlib.collections import PolyCollection
import numpy.ma as ma
import matplotlib.cm as cm

__author__ = "Molly Swanson"
__version__ = "1.0"
__email__  = "mwswanson@cfa.harvard.edu"

#Names of cyclindrical and pseudocyclindrical projections in Basemap
_cylproj = ['cyl','merc','mill','gall']
_pseudocyl = ['moll','robin','sinu','mbtfpq','vandg','hammer']

#Try to load Basemap - if it isn't installed, run without Basemap and just plot RA/dec as x and y values
try:
    from mpl_toolkits.basemap import Basemap as Basemap
except ImportError as e:
    print "WARNING: could not import Basemap:",e
    print "Polygons will be plotted using azimuth and elevation values for x and y."
    print "Visit http://matplotlib.github.com/basemap/ to download and install Basemap."
    # TO DO: check basemap version and only print out below if applicable
    #print "basemap versions earlier than 1.0.3 (not yet released as of 30 Nov 2011) require patch - do"
    #print "  import mpl_toolkits.basemap"
    #print "  mpl_toolkits.basemap.__file__"
    #print "to find the location of the __init__.py file to patch, and replace /path-to/ with the path in the following."
    #print "then do"
    #print "sudo patch /path-to/__init__.py -i basemapinit.patch -o __init__.py.new"
    #print "sudo mv __init__.py.new /path-to/__init__.py"
    useBasemap = False
else:
    useBasemap = True

try:
    import mangle
except ImportError as e:
    useManglepy = False
    useBasemap = False
else:
    useManglepy = True

def main(argv=None):
    """
    Running graphmask.py as a script from the unix command line will take
    a polygon file and plot the polygons to an output graphics file, e.g.,
    $$ python graphmask.py mypolys.pol mypolys.eps
    will generate an .eps file illustrating the polygons in mypolys.pol
    in an all-sky plot.
    Flags:
    -r --range lets you define the azimuth and elevation boundaries of the plot:
               -r azmin azmax elmin elmax
               e.g. python graphmask.py mypolys.pol mypolys.eps -r 330 350 -20 -10
    -a --autoscale will determine the boundaries automatically from the range
               covered by the polygons in the file
               e.g. python graphmask.py -a mypolys.pol mypolys.eps
    You can also specify any keyword arguments used by the plot_mangle_map function
               to customize the plot as keyword=value
               e.g. python graphmask.py mypolys.pol mypolys.eps facecolor='red' plottitle='My Pretty Polygons' projection='hammer'
    """
    import sys
    import argparse
    if argv is None:
        argv=sys.argv
    #define command line parser to interpret command line arguments
    parser = argparse.ArgumentParser(description='call plot_mangle_map from command line.')
    parser.add_argument('infile', help='name of input polygon file in polygon or list format')
    parser.add_argument('outfile', help='name of output image to generate - can be  png, pdf, ps, eps or svg')
    parser.add_argument('-r','--range', type=float, nargs=4,help='min azimuth,max azimuth, min elevation, max elevation of the plot')
    parser.add_argument('-a','--autoscale', action='store_true',default=False,help='use this to automatically determine the plot range from the polygon file')
    #parse the command line for the arguments above.  Others are saved as kwargstrings to be processed later.
    args,kwargstrings = parser.parse_known_args()
    #convert range from arguments into minaz,maxaz,minel,maxel
    if args.range is not None:
        minaz=args.range[0]
        maxaz=args.range[1]
        minel=args.range[2]
        maxel=args.range[3]
    else:
        minaz=None
        maxaz=None
        minel=None
        maxel=None
    #build keyword argument dictionary out of any additional command line arguments that contain equals signs
    kwargs={}
    for s in kwargstrings:
        try:
            key,value=s.split('=',1)
            #ignore an unrecognized argument without an equals sign
        except ValueError:
            print "Ignoring keyword argument "+s+" because it does not contain an equals sign."
            print "To pass keyword arguments, use key=value on command line."
        #check to see if value given is an integer, float, or True/False/None.  Otherwise pass it along as a string
        else:
            try:
                value=int(value)
            except ValueError:
                try:
                    value=float(value)
                except ValueError:
                    if value=='True':
                        value=True
                    elif value=='False':
                        value=False
                    elif value=='None':
                        value=None
            kwargs[key]=value
    #call plot_mangle_map to draw the polygons
    print "plotting polygons from " + args.infile + " ..."
    plot_mangle_map(args.infile,outfilename=args.outfile,autoscale=args.autoscale,minaz=minaz,maxaz=maxaz,minel=minel,maxel=maxel,**kwargs)
    print "plotted polygons in " + args.outfile + "."


def plot_mangle_map(polys, outfilename=None,graphicsfilename=None,cmap='gray_r',plottitle=None,autoscale=False,enlarge_border=.1,bgcolor='auto',drawgrid=True,gridlinewidth=.1,pointsper2pi=30,minaz=None,maxaz=None,minel=None,maxel=None,cenaz=None,cenel=None,**kwargs):
    """
    Function to draw a set of mangle polygons
    Input: polys
           The polygons to draw. Can be the name of a mangle file in
           polygon (.pol or .ply) format, .fits format, or .list
           format (NaN-separated vectors), or an instance of the
           mangle.Mangle class.  If polys is not in .list format, a
           temporary .list format file is created.

           Note that using mangle to convert the polygons into a .list
           format is the slowest part of the process so starting with
           a .list format file is the fastest way to plot.  Use
           graphicsfilename keyword argument to save a copy of the .list
           file created to speed up future plotting.
           
    Return value(s): - matplotlib.collections.PolyCollection object containing
                       the polygons that were plotted
                     - (optional) Skymap instance defining the map projection
                     e.g. : p=plot_mangle_map('mypolys.pol')
                            # p is the plotted PolyCollection
                            p,m=plot_mangle_map('mypolys.pol')
                            # p is the plotted PolyCollection, m is the Skymap instance

    Optional keyword arguments:

    outfilename: name of image file to save.  if None, just draw on screen.
                 Accepts any format valid for matplotlib's savefig
                 (usually png, pdf, ps, eps, and svg)

    graphicsfilename: name of graphics (.list format) file to save to speed
                      up subsequent plotting.  Should have a .list suffix.
                      e.g. : p,m=plot_mangle_map('mypolys.pol',graphicsfilename='mypolys.list')
                             #this takes a long time, and after looking at the plot you decide
                             #it would be much better to plot the polygons in red.
                             p,m=plot_mangle_map('mypolys.list',facecolor='r')
                             #this plots much faster.

                       Note that the .list file saved is projection-specific
                       (and plot-range-specific) since the polygons are trimmed
                       to the map projection region plotted.  If you change the
                       projection or the plot range, you might need to use the
                       original polygon file again.                      

    pointsper2pi: number of points along each circle edge used in making the
                  .list file with mangle's poly2poly.  default is 30.

    cmap: colormap used to plot the weight values of the polygons.  Can be any
          valid matplotlib colormap, or a string defining one, e.g. cmap=cm.jet
          or cmap='jet'.  Default is 'gray_r' which plots black/grayscale polygons
          on a white background.

    plottitle: a title for the plot.

    autoscale: if autoscale is True, the plot range will be determined by the
               min and max azimuth and elevation of the polygons. Default is False.

    enlarge_border: fraction to enlarge the range determined by autoscale to give
                    a reasonable border around the plot.  Default is .1 (10%)

    bgcolor: background color of the map.  If 'auto', use the color at the zero
             value of the colormap.  Default is 'auto'.

    drawgrid: if True, draw grid lines for azimuth and elevation (latitude and longitude)
              and labels for them.

    gridlinewidth: line width for the grid lines.  Use gridlinewidth=0 to draw
                   axis labels but no grid lines.

    minaz,maxaz,minel,maxel: if these are defined, use these values to define the
                             corner regions of the plot.

    cenaz: central azimuth value to use for map projection (equiv to Basemap's lon_0)

    cenel: central elevation to use for map projeciton (equiv to Basemap's lat_0)

    Additional keyword arguments are passed to the Skymap class, which accepts
    any keyword arguments recognized by the Basemap class.  Commonly used ones are:

    projection: the map projection used by Basemap.  print the basemap module variable
                ``supported_projections`` to see a list of the options.

    width, height: the width and height of the region to be plotted,
                   in units of degrees on the sky.  cenaz,cenel,width,height
                   can be used instead of minaz,maxaz,minel,maxel to define
                   the plotted region.

    Further additional keyword arguments not accepted by Skymap/Basemap are passed
    to the draw_weighted_polygons function, which accepts

    holes: can be 'cw'==clockwise,'ccw'==counterclockwise, or None.  If 'cw' ('ccw'),
           polygons with the edge points ordered clockwise (counterclockwise) are
           treated as holes and plotted in the background color.  Mangle polygons
           that have holes in them will be written as separate polygons in the .list
           format, with the exterior border in counterclockwise order and the holes
           in clockwise order, so 'cw' is the default.
           If you know your polygons have no holes, using None can speed up
           the plotting.
           Note that clockwise and counterclockwise are defined as when looking
           up at the sky, so it is reversed from what geographers would expect
           from looking down at the surface of the Earth.

    emptyweight: value of polygon weight indicating a polygon that should not be plotted.
                 Default is 0, but emptyweight can be set to something else to allow for
                 weights whose allowed value range includes 0.

    draw_colorbar: if True, draw a colorbar showing the color scale of the weight
                   values plotted.  Default is True unless all weight values are
                   equal. If you want to control how the colorbar is drawn,
                   use draw_colorbar=False and then use the returned
                   PolyCollection object as your ScalarMappable object.

    vmin,vmax: minimum and maximum values for the colormap scale.

    Still further keywords are passed on to the PolyCollection plotter, which can
    used to customize the appearance of the plotted polygons, e.g.,
    facecolor='red',edgecolor='black',linewidth=0.5.

    Examples:
     import graphmask
     #plot overlapping SDSS plates in a small example patch of sky,
     #using a transparent face color to illustrate the overlap.
     graphmask.plot_mangle_map('dr9plates.pol',minaz=117.5,maxaz=125.5,minel=30, maxel=36.5,facecolor=[0, 0, 1, .2],linewidth=0)
     #plot a boss geometry file starting from a .fits file, and save a graphics .list file for faster plotting later
     p,m=graphmask.plot_mangle_map('boss_geometry_2011_05_20.fits',graphicsfilename='boss_geometry_2011_05_20.list',cenaz=270,gridlinewidth=.1,projection='moll')
    
    """
    #convert colormap string to colormap
    if type(cmap) is str:
        try:
            cmap=vars(cm)[cmap]
        except:
            raise KeyError(cmap + ' is not a valid matplotlib colormap.')
    if bgcolor=='auto':
        bgcolor=cmap(0)
    #if autoscaling, find az, el range from min, max in polygon file
    azrangemin=360
    shiftmin=0
    if autoscale:
        azel,weight=read_poly_list(polys,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi)
        for shift in arange(0,360,10):
            az_shifted=(azel[:,0]+shift) % 360 -shift
            minaz=nanmin(az_shifted)
            maxaz=nanmax(az_shifted)            
            azrange=maxaz-minaz
            if azrange<azrangemin:
                azrangemin=azrange
                shiftmin=shift
                
        az_shifted=(azel[:,0]+shiftmin) % 360 -shiftmin
        minaz=nanmin(az_shifted)
        maxaz=nanmax(az_shifted)  
        minel=nanmin(azel[:,1])
        maxel=nanmax(azel[:,1])
        azrange=maxaz-minaz
        elrange=maxel-minel
        minaz=minaz-enlarge_border*azrange
        maxaz=maxaz+enlarge_border*azrange
        minel=minel-enlarge_border*elrange
        maxel=maxel+enlarge_border*elrange
        if maxel-minel>=180:
            minel=-89.99
            maxel=89.99
            cenel=0
        if maxaz-minaz>=360:
            if cenaz is None:
                cenaz=mean([maxaz,minaz])
            minaz=cenaz-179.99
            maxaz=cenaz+179.99
    #if Basemap is installed, make a Skymap instance in which to draw the plot
    if useBasemap:
        #list of keywords to pass to the map instance
        mapkws=inspect.getargspec(Basemap.__init__).args+inspect.getargspec(Skymap.__init__).args
        mapkwargs=kwargs.copy()
        for key in kwargs.keys():
            if key in mapkws:
                #if key is a keyword for Skymap or Basemap, delete it from kwargs
                val=kwargs.pop(key)
            else:
                #if key is not a keyword for Skymap or Basemap, delete it from mapkwargs
                val=mapkwargs.pop(key)        
        m=Skymap(minaz=minaz,maxaz=maxaz,minel=minel,maxel=maxel,cenaz=cenaz,cenel=cenel,**mapkwargs)
        azel,weight=m.get_graphics_polygons(polys,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi)
        if bgcolor is not None:
            m.drawmapboundary(fill_color=bgcolor)
        else:
            m.drawmapboundary()
        if drawgrid:
            m.draw_grid(linewidth=gridlinewidth)
        if azel is not None:
            p=draw_weighted_polygons(m.world2pix(azel),weight=weight,cmap=cmap,bgcolor=bgcolor,**kwargs)
        else:
            p=None
    #if Basemap is not installed, draw the polygons using az,el as x,y
    else:
        corners_defined=(None not in [minaz,maxaz,minel,maxel])
        if not corners_defined:
            if cenaz is None:
                cenaz=180.
            minaz=cenaz-179.99
            maxaz=cenaz+179.99
            minel=-89.99
            maxel=89.99
        else:
            if cenaz is None:
                cenaz=mean([minaz,maxaz])            
        axis([minaz,maxaz,minel,maxel])
        if bgcolor is not None:
            gca().set_axis_bgcolor(bgcolor)
        gca().invert_xaxis()
        azel,weight=get_graphics_polygons(polys,trimmer_poly=[minaz,maxaz,minel,maxel],cenaz=cenaz,graphicsfilename=graphicsfilename)
        if azel is not None:
            p=draw_weighted_polygons(azel,weight=weight,holes='ccw',cmap=cmap,bgcolor=bgcolor,**kwargs)
        else:
            p=None
        m=None
    if plottitle is not None:
        title(plottitle)
    if outfilename is not None:
        savefig(outfilename)
    else:
        show()
    howmany=expecting()
    if howmany<2:
        return p
    else:
        return p,m

def polysplit(xy):
    """
    Function that splits a NaN-separated array of polygon points into a list of polygons suitable for
    plotting with matplotlib's PolyCollection.
    input argument: Nx2 numpy array of az,el points along polygon outlines, with the breaks between the
    different polygons denoted by a row with two NaN's. This is mangle's .list format.
    output: a list of az,el arrays with each list entry containing the outline of one polygon.  This is
    the format required for plotting the polygons with matplotlib's PolyCollection.
    """
    polyend=nonzero(all(isnan(xy),axis=1))[0]
    polystart=append(0,polyend[:-1]+1)
    polys=[xy[slice(*polyspan)] for polyspan in zip(polystart,polyend)]
    return polys

def read_poly_list(polys,pointsper2pi=30,graphicsfilename=None):
    """
    Function to read in graphics polygons from a .list format mangle file.  If input is a .pol/.ply
    polygon file, a .fits file, or a mangle.Mangle instance, a (temporary) graphics .list file will
    be generated.

    Argument:
    polys: .list/.pol/.ply/.fits file or mangle.Mangle instance containing polygons
    
    Keyword arguments:
    pointsper2pi: used when generating a .list file - the number of points per 2pi along the edges
                  of the polygon used for defining the polygon outlines.

    graphicsfilename: file in which to save the .list file to use for future plotting.
                      if graphicsfilename=None, the .list file will be deleted after it is read in.

    Returns:
    azel: NaN-separated numpy array of az,el points defining the boundaries of the polygons.

    weight (optional): array of weight values for each polygon

    midpoints (optional): array of midpoints for each polygon 
    """
    jpol = tempfile.NamedTemporaryFile('rw+b')
    #if mangle.Mangle instance, write polygons to temporary file
    if useManglepy:
        if isinstance(polys,mangle.Mangle):
            polys.writeply(jpol.name)
            filename=jpol.name
        else:
            filename=polys
    else:
        filename=polys
    #try loading file - will be successful if file is .list format
    try:
        azel=loadtxt(filename)
    #if file couldn't be read as .list file, run poly2poly to make a .list file and read it in
    except:
        jlist = tempfile.NamedTemporaryFile('rw+b')
        if graphicsfilename is None:
            graphicsfilename=jlist.name                    
        call = ' '.join(('poly2poly -q -ol'+str(pointsper2pi),filename,graphicsfilename))
        subprocess.check_call(call,shell=True)
        azel=loadtxt(graphicsfilename)
        weights=loadtxt(graphicsfilename+'.weight')
        try:
            os.remove(jlist.name+'.weight')
        except:
            pass
        jlist.close()
    #read in polygon weights from separate .list.weight file
    else:
        weights=loadtxt(filename+'.weight')
    if len(shape(weights))==1:
        weight=array([weights[1]])
        #also try to read in midpoints from .list.weight file (not present if .list file created with mangle2.2 or earlier)
        try:
            midpoints=array([weights[2:4]])
        except IndexError:
            midpoints=None
    else:
        weight=weights[:,1]
        try:
            midpoints=weights[:,2:4]
        except IndexError:
            midpoints=None
    howmany=expecting()
    if howmany<2:
        return azel
    if howmany==2:
        return azel,weight
    else: 
        return azel,weight,midpoints

def draw_weighted_polygons(inpolys,weight=None,holes='cw',emptyweight=0,cmap='gray_r',facecolor=None,draw_colorbar=True, edgecolor='face',linewidth=.8, bgcolor='auto', vmin=None,vmax=None, **kwargs):
    """
    function that draws polygons using matplotlib's PolyCollection with a
    ScalarMappable weight value that can be used to draw a colorbar.

    Can be used for polygons defined in map x,y coordinates in a Skymap/Basemap instance,
    or with az,el points to plot the polygons using az,el as x,y.

    Argument:
    inpolys: a NaN-separated array of x,y points to plot along the edge of a polygon.
    
    Keyword arguments:
    weight: an array of weights with one weight per polygon, to be colored according
            to the colormap.
    
    holes: can be 'cw'==clockwise,'ccw'==counterclockwise, or None.  If 'cw' ('ccw'),
           polygons with the edge points ordered clockwise (counterclockwise) are
           treated as holes and plotted in the background color.  Mangle polygons
           that have holes in them will be written as separate polygons in the .list
           format, with the exterior border in counterclockwise order and the holes
           in clockwise order, so 'cw' is the default.
           If you know your polygons have no holes, using None can speed up
           the plotting.
           Note that clockwise and counterclockwise are defined as when looking
           up at the sky, so it is reversed from what geographers would expect
           from looking down at the surface of the Earth.

           If you are using draw_weighted_polygons to plot az,el points directly, you
           should use holes='ccw' since the orientation is reversed if a map projection
           is not applied.

    emptyweight: value of polygon weight indicating a polygon that should not be plotted.
                 Default is 0, but emptyweight can be set to something else to allow for
                 weights whose allowed value range includes 0.

    draw_colorbar: if True, draw a colorbar showing the color scale of the weight
                   values plotted.  Default is True unless all weight values are
                   equal. If you want to control how the colorbar is drawn,
                   use draw_colorbar=False and then use the returned
                   PolyCollection object as your ScalarMappable object.

    vmin,vmax: minimum and maximum values for the colormap scale.

    Further keywords are passed on to the PolyCollection plotter, which can
    used to customize the appearance of the plotted polygons, e.g.,
    facecolor='red',edgecolor='black',linewidth=0.5.
    """
    #convert colormap string into colormap
    if type(cmap) is str:
        try:
            cmap=vars(cm)[cmap]
        except:
            raise KeyError(cmap + ' is not a valid matplotlib colormap.')
    #split NaN-separated polygon array into list of separate polygons
    polys=polysplit(inpolys)
    if weight is not None:
        weights=weight
    else:
        weights=ones(len(polys))
    #make masked array, masking out polygons that are holes or whose values are equal to the emptyweight value
    if holes is not None:
        if holes=='cw':
            polyholes=array(map(ispolycw,polys))
        elif holes=='ccw':
            polyholes=array(map(ispolyccw,polys))
        else:
            raise ValueError('holes value should be "cw" to denote clockwise holes, "ccw" to denote counterclockwise holes, or None')
        if emptyweight is not None:
            masked_weight=ma.masked_array(weights,mask=polyholes | (weights==emptyweight),fill_value=emptyweight)
        else:
            masked_weight=ma.masked_array(weights,mask=polyholes)
    else:
        if emptyweight is not None:
            masked_weight=ma.masked_array(weights,(weights==emptyweight),fill_value=emptyweight)
        else:
            masked_weight=weights
    #set the background color and set the colormap to plot "bad" polygons as the background color
    if bgcolor is None:
        bgcolor=gca().get_axis_bgcolor()
    elif bgcolor=='auto':
        bgcolor=cmap(0)
    cmap.set_bad(bgcolor)
    #plot the polygons as a PolyCollection
    p=PolyCollection(polys,cmap=cmap,edgecolor=edgecolor,facecolor=facecolor,linewidth=linewidth, **kwargs)
    #if facecolor=None, use the weight values for the color
    if facecolor is None:
        #if no boundary values were specified for the colormap, determine automatically from range of weights
        if None in [vmin,vmax]:
            if vmin is None:
                vmin=masked_weight.min()
            if vmax is None: 
                vmax=masked_weight.max()
            #if all weight values are the same, don't draw a colorbar
            if vmin==vmax:
                vmin=0
                draw_colorbar=False
            else:
                if vmin>0:
                    vmin=max(0,vmin-.1*(vmax-vmin))
                else:
                    vmin=vmin-.1*(vmax-vmin)
        #apply the colorbar limits and the ScalarMappable weights array to the PolyCollection object
        p.set_clim(vmin=vmin,vmax=vmax)
        p.set_array(masked_weight)
    else:
        draw_colorbar=False
    #add the PolyCollection to the current axis
    gca().add_collection(p)
    if draw_colorbar:
        colorbar(p)
    return p

def ispolyccw(poly):
    """
    Tests if a series of points around the edges of a polygon is wound counterclockwise.
    Returns True for a counterclockwise polygon, False for a clockwise polygon.
    (when looking up at the sky).  This applies to x,y map projection coordinates - if used on
    a polygon in az,el coordinates, the definitions of clockwise and counterclockwise are reversed.
    """
    tot=sum([cross(*vecs) for vecs in zip(poly,poly[1:])])
    return tot<0

def ispolycw(poly):
    """
    Tests if a series of points around the edges of a polygon is wound counterclockwise.
    Returns True for a clockwise polygon, False for a counterclockwise polygon.
    (when looking up at the sky).  This applies to x,y map projection coordinates - if used on
    a polygon in az,el coordinates, the definitions of clockwise and counterclockwise are reversed.
    """
    tot=sum([cross(*vecs) for vecs in zip(poly,poly[1:])])
    return tot>0
    ## Below is an alternative test that applies to az,el points
    #tot=array([0.0,0.0,0.0])
    #p=poly*pi/180.
    #cosp=cos(p)
    #sinp=sin(p)
    #x=cosp[:,1]*cosp[:,0]
    #y=cosp[:,1]*sinp[:,0]
    #z=sinp[:,1]
    #p=vstack([x,y,z]).transpose()
    #mid=mean(p,0)
    #tot=sum([cross(*vecs) for vecs in zip(p-mid,p[1:]-mid)],0)
    #return dot(tot,mid)<0

## {{{ http://code.activestate.com/recipes/284742/ (r4)
import inspect,dis
def expecting():
    """
    Call within a function to find out how many output arguments the caller of the function
    is expecting.  This is to emulate matlab's nargout.
    """
    f = inspect.currentframe()
    f = f.f_back.f_back
    c = f.f_code
    i = f.f_lasti
    bytecode = c.co_code
    instruction = ord(bytecode[i+3])
    if instruction == dis.opmap['UNPACK_SEQUENCE']:
        howmany = ord(bytecode[i+4])
        return howmany
    elif instruction == dis.opmap['POP_TOP']:
        return 0
    return 1
## end of http://code.activestate.com/recipes/284742/ }}}

def get_graphics_polygons(polys,trimmer_poly=None,cenaz=180,projection='cyl',pointsper2pi=30,graphicsfilename=None):
    """
    Takes polygons from a polygon file (.list,.ply/.pol, or .fits format) or mangle.Mangle
    instance and returns a NaN-separated array of az,el points along the edges of the polygons.
    Polygons can be trimmed to only return the polygons that lie within a given region defined
    by a set of trimmer polygons.

    Argument:
    polys: .list/.pol/.ply/.fits file or mangle.Mangle instance containing polygons
    
    Keyword arguments:
    trimmer_poly: .list/.pol/.ply/.fits file or mangle.Mangle instance containing the trimmer
                  polygons defining the boundary of the region to be plotted.  This would
                  typically be a rectangle spanning the min and max az,el to be plotted. For
                  cylindrical/pseudocylindrical projections, it is recommended to cut out the
                  north and south poles with small holes, and define the trimmer polygons such
                  that there is a cut along the meridian 180 degrees opposite the central
                  azimuth of the map so no polygons cross this meridian.
    
    pointsper2pi: used when generating a .list file - the number of points per 2pi along the edges
                  of the polygon used for defining the polygon outlines.

    graphicsfilename: file in which to save the .list file to use for future plotting.
                      if graphicsfilename=None, the .list file will be deleted after it is read in.

    cenaz: the central meridian of the plot.  If 'projection' is in the list of cylindrical
           or pseudocylindrical projections (true for the default projection 'cyl'), azimuth
           points will be adjusted so they all lie between cenaz-180 and cenaz+180.

    projection: the Basemap map projection.  (ignore if you are not using Basemap)
    """
    trim=True
    if trimmer_poly is None:
        trim=False
    if useManglepy:
        if not isinstance(polys,mangle.Mangle):
            if polys[-5:]=='.list':
                trim=False
            if polys[-5:]=='.fits':
                polys=mangle.Mangle(polys)
    else:
        trim=False
    if trim:
        trimmed_polys=trim_mask(polys,trimmer_poly)
        if trimmed_polys.npoly>0:
            azel,weight,midpoints=read_poly_list(trimmed_polys,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi)
        else:
            print "WARNING: no polygons in range to be plotted"
            azel=None
            weight=None
            midpoints=None
    else:
        azel,weight,midpoints=read_poly_list(polys,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi) 
    if (projection in _cylproj) or (projection in _pseudocyl):
        #azimuth has to be between cenaz-180 and cenaz+180
        azel[:,0]=(azel[:,0]-(cenaz-180))%360+(cenaz-180)
    howmany=expecting()
    if howmany<2:
        return azel
    if howmany==2:
        return azel,weight
    else: 
        return azel,weight,midpoints

if useBasemap:
    class Skymap(Basemap):
        """
        Class that extends the Basemap class to provide appropriate default behavior and
        additional functions for plotting maps of the sky rather than maps of the Earth.

        Key differences from Basemap:
        - keyword argument Celestial is set to True, so plots are drawn with aziumuth increasing
          to the left, and rsphere is set to 180/pi so that xy map coordinates are (approximately)
          in units of degrees on the sky rather than meters on the surface of the Earth. 
            
        - Instead of specifying corners with llcrnrlon, llcrnrlat, urcrnrlon, and urcrnrlat,
          use minaz,minel,maxaz,maxel for the minimum and maximum azimuth and elevation.

        - Likewise, use cenaz and cenel for the central azimuth and elevation rather than
          lat_0 and lon_0.

        - cenaz, cenal, elmin, elmax, azmin, and azmax are defined as class variables.

        - default projection is McBride-Thomas Flat Polar Quartic ('mbtfpq') for all-sky plots,
          gnomonic ('gnom') for small regions (spanning < 90 degrees), and cylindrical ('cyl')
          for larger regions. 

        - if defining a range to plot around a central point using length and width, the length
          and width should be given in units of degrees.

        - meridian labeling automatically uses celestial conventions, with aziumuth values
          ranging from 0 to 360 degrees, increasing to the left.

        - drawparallels and drawmeridians methods write labels on the left and bottom of the
          plot respectively by default. 

        - converting to and from az,el to x,y map coordinates is done by the world2pix and
          pix2world methods rather than calling the Skymap instance itself

        - A set of mangle polygons to use for trimming input polygons to plot in this particular
          map region (i.e, to deal with meridian crossing, poles, etc.) is generated and saved
          as a mangle.Mangle instance in the class variable trimmer_poly.

        Methods included:

        __init__: wrapper for Basemap's __init__ to set up a map with celestial defaults

        drawmeridians: wraps Basemap's drawmeridians using celestial conventions for labels

        drawparallels: wraps Basemap's drawparallels using celestial conventions for lablels

        mlabelzero22pi: formatting function for celestial meridians, ranging from 0 to 360

        draw_grid: draw a grid of meridian and parallel lines

        pix2world: convert from x,y map coordinates to az, el

        world2pix: convert from az, el to x, y map coordinates

        make_trimmer_poly: 

        get_graphics_polygons:
        """
        def __init__(self,
                     minaz=None,
                     maxaz=None,
                     minel=None,
                     maxel=None,
                     cenaz=None,
                     cenel=None,
                     width=None,
                     height=None,
                     projection=None,
                     **kwargs):
            """
    Sets up a skymap (extension of basemap) with the specified map projection to use
    for plotting on the celestial sphere.
    
    Keyword arguments:
    
    minaz,maxaz,minel,maxel: if these are defined, use these values to define the
    corner regions of the plot.
    
    cenaz: central azimuth value to use for map projection (equiv to Basemap's lon_0)
    
    cenel: central elevation to use for map projeciton (equiv to Basemap's lat_0)

    Additional keyword arguments passed along to the Basemap __init__ function, so
    any keyword arguments recognized by the Basemap class can be used. Common ones are:

    projection: the map projection used by Basemap.  print the basemap module variable
                ``supported_projections`` to see a list of the options.

    width, height: the width and height of the region to be plotted,
                   in units of degrees on the sky.  cenaz,cenel,width,height
                   can be used instead of minaz,maxaz,minel,maxel to define
                   the plotted region.            
    """
            # check if corners or range is defined
            corners_defined=(None not in [minaz,maxaz,minel,maxel])
            range_defined=(None not in [width,height,cenaz,cenel])

            #use corners to define map projection region
            if corners_defined:
                #default projection is gnomonic for regions spanning less than 90 degrees,
                #cylindrical for larger regions
                if projection is None:
                    if maxaz-minaz<90 and maxel-minel<90:
                        projection='gnom'
                    else:
                        projection='cyl'
                #if cenel and cenaz not defined, use center of az,el range
                if cenel is None:
                    cenel = mean([minel,maxel])
                if cenaz is None:
                    cenaz = mean([minaz,maxaz])
                #convert azimuth values into longitude values
                if (projection in _pseudocyl) or (projection in _cylproj):
                    lon_0=cenaz                
                    llcrnrlon=minaz
                    urcrnrlon=maxaz
                else:
                    lon_0=-cenaz                
                    llcrnrlon=-maxaz
                    urcrnrlon=-minaz
                #call Basemap's init with celestial defaults
                Basemap.__init__(self,llcrnrlon=llcrnrlon,llcrnrlat=minel,urcrnrlon=urcrnrlon,urcrnrlat=maxel,lon_0=lon_0,lat_0=cenel,projection=projection, rsphere=180./pi, celestial=True,**kwargs)
            elif range_defined:
                #default projection is gnomonic for regions spanning less than 90 degrees,
                #cylindrical for larger regions
                if projection is None:
                    if width<90 and height<90:
                        projection='gnom'
                    else:
                        projection='cyl'
                #convert azimuth values into longitude values
                if (projection in _pseudocyl) or (projection in _cylproj):
                    lon_0=-cenaz
                else:
                    lon_0=-cenaz
                #call Basemap's init with celestial defaults
                Basemap.__init__(self,lon_0=lon_0,lat_0=cenel,width=width,height=height, projection=projection,rsphere=180./pi,celestial=True,**kwargs)
            else:
                #if no corners or range defined, plot as all-sky
                #default projection is mptfpq
                if projection is None:
                    projection='mbtfpq'
                if cenaz is None:
                    cenaz=180.
                if cenel is None:
                    cenel=0.
                #convert azimuth values into longitude values
                if (projection in _pseudocyl) or (projection in _cylproj):
                    lon_0=cenaz
                else:
                    lon_0=-cenaz
                #call Basemap's init with celestial defaults
                Basemap.__init__(self,lon_0=lon_0,lat_0=cenel,projection=projection, rsphere=180./pi,celestial=True,**kwargs)

            #save cenaz, cenel as class variables
            self.cenel=self.projparams['lat_0']
            if (projection in _pseudocyl) or (projection in _cylproj):
                self.cenaz=self.projparams['lon_0']
            else:
                self.cenaz=-self.projparams['lon_0']

            #save elmin, elmax, azmin, azmax as class variables
            self.elmin=self.latmin
            self.elmax=self.latmax
            if (projection in _pseudocyl) or (projection in _cylproj):
                self.azmin=self.lonmin
                self.azmax=self.lonmax
            else:
                self.azmin=(-self.lonmax+360) % 360
                self.azmax=(-self.lonmin+360) % 360
            if self.azmin>=self.azmax:
                self.azmin-=360

            #make a trimmer set of mangle polygons for to use to clip polygons to this map projection
            self.make_trimmer_poly()
        
        def mlabelzero22pi(lon):
            """
            Formats a longitude value using celestial conventions, with values ranging
            from 0 to 360 degrees.  Use fmt=mlabelzero22pi as keyword argument to drawmeridians
            """
            fmt='%g'
            if rcParams['text.usetex']:
                lonlabstr = r'${%s\/^{\circ}}$'%fmt
            else:
                lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
                lonlab = lonlabstr%lon
            return lonlab

        def drawmeridians(self,meridians,fmt=mlabelzero22pi,labels=[0,0,0,1], linewidth=1.,**kwargs):
            """
            Wrapper for Basemap's drawmeridians with celestial defaults.
            Call the same way as Basemap version.
            """
            if self.projection in _cylproj:
                meridians1=(meridians-(self.cenaz-180))%360+(self.cenaz-180)
                Basemap.drawmeridians(self,meridians1,fmt=fmt,labels=[0,0,0,0],linewidth=linewidth,**kwargs)
                Basemap.drawmeridians(self,meridians,fmt=fmt,labels=labels,linewidth=0,**kwargs)
            else:
                Basemap.drawmeridians(self,meridians,fmt=fmt,labels=labels,linewidth=linewidth,**kwargs)

        def drawparallels(self,parallels,labels=[1,0,0,0], linewidth=1.,**kwargs):
            """
            Wrapper for Basemap's drawparallels with celestial defaults.
            Call the same way as Basemap version.
            """
            Basemap.drawparallels(self,parallels,labels=labels,linewidth=linewidth,**kwargs)
        
        def draw_grid(self,linewidth=1.,**kwargs):
            """
            Draw grid on skymap by calling drawmeridians and drawparallels functions.
            Automatically determines grid spacing.
            Keyword arguments are passed along to drawmeridians and drawparallels.
            """
            if self.projection in _pseudocyl:
                if self.projection == 'vandg':
                    self.drawparallels(np.arange(-90.,91.,30.),linewidth=linewidth,labels=[0,0,0,0],**kwargs)
                else:
                    self.drawparallels(np.arange(-90.,91.,30.),linewidth=linewidth,**kwargs)
                if self.projection == 'robin':
                    self.drawmeridians(np.arange(0.,360.,60.),linewidth=linewidth,**kwargs)
                elif self.projection == 'mbtfpq':
                    self.drawmeridians(np.arange(0.,360.,90.),linewidth=linewidth,**kwargs)
                else:
                    self.drawmeridians(np.arange(0.,360.,60.),linewidth=linewidth,labels=[0,0,0,0],**kwargs)
            elif self.azmax-self.azmin>359:
                self.drawparallels(np.arange(-90.,91.,30.),labels=[1, 0, 0, 0],linewidth=linewidth,**kwargs)
                self.drawmeridians(np.arange(0.,360.,60.),labels=[0, 0, 0, 1],linewidth=linewidth,**kwargs)
            else:
                autolocate_az=matplotlib.ticker.AutoLocator()
                autolocate_el=matplotlib.ticker.AutoLocator()
                self.drawparallels(autolocate_el.bin_boundaries(self.elmin,self.elmax),linewidth=linewidth,**kwargs)
                self.drawmeridians(autolocate_az.bin_boundaries(self.azmin,self.azmax),linewidth=linewidth,**kwargs)


        def pix2world(self,xy):
            """
            function to convert x/y map projection coordinates into az/el
            coordinates.

            input argument should be Nx2 numpy array of x/y coordinates
            (1st column x, 2nd column y)

            returns Nx2 numpy array of az/el coordinates
            (1st column az, 2nd column el)

            NaN values in input are retained as NaN values in output
            """
            xyout=vstack([Basemap.__call__(self,xy[:,0],xy[:,1],inverse=True)]).transpose()
            xyout[isnan(xy)]=nan
            return xyout

        def world2pix(self,xy):
            """
            function to convert az/el coordinates into  x/y map projection
            coordinates.

            input argument should be Nx2 numpy array of az/el coordinates
            (1st column az, 2nd column el)

            returns Nx2 numpy array of x/y coordinates
            (1st column x, 2nd column y)

            NaN values in input are retained as NaN values in output
            """
            xyout=vstack([Basemap.__call__(self,xy[:,0],xy[:,1],inverse=False)]).transpose()
            xyout[isnan(xy)]=nan
            return xyout

        def make_trimmer_poly(self):
            """
            function to create a set of mangle polygons to use for trimming polygons to
            plot in this map projection.

            For projections 'aeqd','npaeqd','spaeqd','stere','npstere','spstere','ortho',
            the trimmer polygons cover the whole sky minus a hole around the antipode
            of the projection.  Hole has a radius of 90.01 degrees for ortho projection
            (i.e., the unseen hemisphere), or .01 degrees for the others.

            For projections 'gnom','cass','poly','lcc','eqdc','aea','laea','nplaea','splaea',
            trimmer polygons cover a rectangular region with great circle edges around
            the map region.

            For projections 'omerc','tmerc', no trimmer polygons are defined.

            The resulting trimmer polygons are saved as a mangle.Mangle instance in
            the class variable trimmer_poly.
            """
            
            mapproj=self.projection
            cenaz=self.cenaz
            cenel=self.cenel
            if (mapproj in _cylproj) or (mapproj in _pseudocyl):
                trim=array([[cenaz, cenaz+179.99,0,89.99],
                            [cenaz-179.99,cenaz,0,89.99],
                            [cenaz, cenaz+179.99,-89.99,0],
                            [cenaz-179.99,cenaz,-89.99,0]])
                trimtype='r'            
            elif mapproj in ['aeqd','npaeqd','spaeqd','stere','npstere','spstere','ortho']:
                #make a mangle file with a single polygon defining a hemisphere 
                trim=array([[cenaz,cenel, 90]])
                jtrim = tempfile.NamedTemporaryFile('rw+b')
                jlist = tempfile.NamedTemporaryFile('rw+b')
                savetxt(jtrim.name,trim,fmt='%.12f')
                call = ' '.join(('poly2poly -q -ic -ol4',jtrim.name,jlist.name))
                subprocess.check_call(call,shell=True)
                jtrim.close()
                #load points on equator of projection
                eq=loadtxt(jlist.name)
                jlist.close()
                os.remove(jlist.name+'.weight')
                #make polygons of 4 quadrants covering the whole sky minus a hole around the antipode of the projection
                #hole has radius of .01 deg for aeqd and stereo projections, 90.01 for ortho projection (i.e., the unseen hemisphere)
                if mapproj=='ortho':
                    radius=-90.01
                else:
                    radius=-0.01
                trim=array([[eq[0,0], eq[0,1], 90.,  eq[1,0], eq[1,1], 90., cenaz+180., -cenel, radius],
                            [eq[0,0], eq[0,1], -90.,  eq[1,0], eq[1,1], 90., cenaz+180., -cenel, radius],
                            [eq[0,0], eq[0,1], 90.,  eq[1,0], eq[1,1], -90., cenaz+180., -cenel, radius],
                            [eq[0,0], eq[0,1], -90.,  eq[1,0], eq[1,1], -90., cenaz+180., -cenel, radius]])
                trimtype='c'
            elif mapproj in ['gnom','cass','poly','lcc','eqdc','aea','laea','nplaea','splaea']:
                #define a square with great circle edges around the map region
                xmin=self.xmin
                ymin=self.ymin
                xmax=self.xmax
                ymax=self.ymax
                xmid=mean([xmin,xmax])
                ymid=mean([ymin,ymax])                
                v=array([[xmin,ymin],[xmin,ymid],[xmid,ymid],[xmid,ymin],
                         [xmid,ymin],[xmid,ymid],[xmax,ymid],[xmax,ymin],
                         [xmin,ymid],[xmin,ymax],[xmid,ymax],[xmid,ymid],
                         [xmid,ymid],[xmid,ymax],[xmax,ymax],[xmax,ymid]])
                trim=reshape(self.pix2world(v),(-1,8))
                trimtype='v'
            elif mapproj in ['omerc','tmerc']:
                trim=None
                trimtype=None
            else:
                raise ValueError('projection '+mapproj+' not supported by Skymap')
            jtrim = tempfile.NamedTemporaryFile('rw+b')
            jpol = tempfile.NamedTemporaryFile('rw+b',suffix='.pol')
            savetxt(jtrim.name,trim,fmt='%.12f')
            call = ' '.join(('poly2poly -q -i'+trimtype,jtrim.name,jpol.name))
            subprocess.check_call(call,shell=True)
            jtrim.close()
            if useManglepy:
                trimmer_poly=mangle.Mangle(jpol.name,keep_ids=True)
                self.trimmer_poly=trimmer_poly
            else:
                self.trimmer_poly=None

        def get_graphics_polygons(self,polys,pointsper2pi=30,graphicsfilename=None):
            """
    Takes polygons from a polygon file (.list,.ply/.pol, or .fits format) or mangle.Mangle
    instance and returns a NaN-separated array of az,el points along the edges of the polygons.
    Polygons are trimmed to only return the polygons that lie within a given region defined
    by the trimmer polygons of this Skymap instance.

    Argument:
    polys: .list/.pol/.ply/.fits file or mangle.Mangle instance containing polygons
    
    Keyword arguments:    
    pointsper2pi: used when generating a .list file - the number of points per 2pi along the edges
                  of the polygon used for defining the polygon outlines.

    graphicsfilename: file in which to save the .list file to use for future plotting.
                      if graphicsfilename=None, the .list file will be deleted after it is read in.
    """
            
            azel,weight,midpoints=get_graphics_polygons(polys,trimmer_poly=self.trimmer_poly,projection=self.projection,
                                                        cenaz=self.cenaz,pointsper2pi=pointsper2pi,graphicsfilename=graphicsfilename)
            howmany=expecting()
            if howmany<2:
                return azel
            if howmany==2:
                return azel,weight
            else: 
                return azel,weight,midpoints

def trim_mask(polys,trimmer_poly):
    """
    wrapper for mangle unix shell script trim_mask.sh, which trims input polygons to lie
    within a given region defined by a set of trimmer polygons.

    Arguments:
    polys: .list/.pol/.ply/.fits file or mangle.Mangle instance containing polygons
    
    trimmer_poly: .list/.pol/.ply/.fits file or mangle.Mangle instance containing the trimmer
                  polygons defining the boundary of the region to be plotted. 
 
    returns mangle.Mangle instance containing the trimmed polygons.
    """
    import shutil
    jpol = tempfile.NamedTemporaryFile('rw+b')
    jpols = tempfile.NamedTemporaryFile('rw+b')
    if isinstance(trimmer_poly,mangle.Mangle):
        trimmer_poly.writeply(jpol.name)
    elif isinstance(trimmer_poly,list):
        cenaz=mean(trimmer_poly[0:2])
        cenel=mean(trimmer_poly[2:4])
        azrange=trimmer_poly[1]-trimmer_poly[0]
        elrange=trimmer_poly[3]-trimmer_poly[2]
        trim=array([[cenaz, cenaz+.5*azrange,cenel,cenel+.5*elrange],
                    [cenaz-.5*azrange, cenaz,cenel,cenel+.5*elrange],
                    [cenaz, cenaz+.5*azrange,cenel-.5*elrange,cenel],
                    [cenaz-.5*azrange, cenaz,cenel-.5*elrange,cenel]])
        trimtype='r'          
        jtrim = tempfile.NamedTemporaryFile('rw+b')
        savetxt(jtrim.name,trim,fmt='%.12f')
        call = ' '.join(('poly2poly -q -i'+trimtype,jtrim.name,jpol.name))
        subprocess.check_call(call,shell=True)
        jtrim.close()        
    else:
        shutil.copyfile(trimmer_poly,jpols.name)        
    if isinstance(polys,mangle.Mangle):
        polys.writeply(jpols.name)
    else:
        shutil.copyfile(polys,jpols.name)
    jrast = tempfile.NamedTemporaryFile('rw+b',suffix='.pol')
    pwd=os.getcwd()
    jdir = tempfile.mkdtemp()
    os.chdir(jdir)
    call = ' '.join(('trim_mask.sh',jpols.name,jpol.name,jrast.name))
    try:
        x=subprocess.check_output(call,shell=True,stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as x:
        os.chdir(pwd)
        print x.output
        raise
    except:
        os.chdir(pwd)
        raise
    else:
        os.chdir(pwd)
    os.rmdir(jdir)
    trimmed_polys=mangle.Mangle(jrast.name,keep_ids=True)
    jrast.close()    
    jpol.close()
    jpols.close()
    return trimmed_polys

if __name__=='__main__':
    sys.exit(main())
