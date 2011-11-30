import matplotlib
if __name__=='__main__':
    matplotlib.use('Agg')    
from pylab import *
import subprocess
import os
import tempfile
import inspect
import mangle
from matplotlib.collections import PolyCollection
import numpy.ma as ma
import matplotlib.cm as cm
_cylproj = ['cyl','merc','mill','gall']
_pseudocyl = ['moll','robin','sinu','mbtfpq','vandg','hammer']

try:
    from mpl_toolkits.basemap import Basemap as Basemap
except ImportError as e:
    print "WARNING: could not import Basemap:",e
    print "Polygons will be plotted using azimuth and elevation values for x and y."
    print "Visit http://matplotlib.github.com/basemap/ to download and install Basemap."
    print "basemap versions earlier than 1.0.3 (not yet released as of 30 Nov 2011) require patch - do"
    print "  import mpl_toolkits.basemap"
    print "  mpl_toolkits.basemap.__file__"
    print "to find the location of the __init__.py file to patch, and replace /path-to/ with the path in the following."
    print "then do"
    print "sudo patch /path-to/__init__.py -i basemapinit.patch -o __init__.py.new"
    print "sudo mv __init__.py.new /path-to/__init__.py"
    useBasemap = False
else:
    useBasemap = True

def main(argv=None):
    import sys
    import argparse
    if argv is None:
        argv=sys.argv
    parser = argparse.ArgumentParser(description='call plot_mangle_map from command line.')
    parser.add_argument('infile', help='name of input polygon file in polygon or list format')
    parser.add_argument('outfile', help='name of output image to generate - can be  png, pdf, ps, eps or svg')
    parser.add_argument('-r','--range', type=float, nargs=4,help='min azimuth,max azimuth, min elevation, max elevation of the plot')
    parser.add_argument('-a','--autoscale', action='store_true',default=False,help='use this to automatically determine the plot range from the polygon file')
    args,kwargstrings = parser.parse_known_args()
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
        
    kwargs={}
    for s in kwargstrings:
        try:
            key,value=s.split('=',1)
        except ValueError:
            print "Ignoring keyword argument "+s+" because it does not contain an equals sign."
            print "To pass keyword arguments, use key=value on command line."
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
    plot_mangle_map(args.infile,outfilename=args.outfile,autoscale=args.autoscale,minaz=minaz,maxaz=maxaz,minel=minel,maxel=maxel,**kwargs)


def plot_mangle_map(polys, outfilename=None,graphicsfilename=None,cmap='gray_r',plottitle=None,autoscale=False,enlarge_border=.1,bgcolor='auto',drawgrid=True,gridlinewidth=.1,pointsper2pi=30,minaz=None,maxaz=None,minel=None,maxel=None,cenaz=None,cenel=None,**kwargs):
    if type(cmap) is str:
        try:
            cmap=vars(cm)[cmap]
        except:
            raise KeyError(cmap + ' is not a valid matplotlib colormap.')
    if bgcolor=='auto':
        bgcolor=cmap(0)
    if autoscale:
        azel,weight=read_poly_list(polys)
        minaz=nanmin(azel[:,0])
        maxaz=nanmax(azel[:,0])
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
        p=draw_weighted_polygons(m.world2pix(azel),weight=weight,cmap=cmap,bgcolor=bgcolor,**kwargs)
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
        p=draw_weighted_polygons(azel,weight=weight,holes='ccw',cmap=cmap,bgcolor=bgcolor,**kwargs)
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
    polyend=nonzero(all(isnan(xy),axis=1))[0]
    polystart=append(0,polyend[:-1]+1)
    polys=[xy[slice(*polyspan)] for polyspan in zip(polystart,polyend)]
    return polys

def read_poly_list(polys,pointsper2pi=30,graphicsfilename=None):
    jpol = tempfile.NamedTemporaryFile('rw+b')
    if isinstance(polys,mangle.Mangle):
        polys.writeply(jpol.name)
        filename=jpol.name
    else:
        filename=polys
    try:
        azel=loadtxt(filename)
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
    else:
        weights=loadtxt(filename+'.weight')
    if len(shape(weights))==1:
        weight=array([weights[1]])
        midpoints=array([weights[2:4]])
    else:
        weight=weights[:,1]
        midpoints=weights[:,2:4]
    howmany=expecting()
    if howmany<2:
        return azel
    if howmany==2:
        return azel,weight
    else: 
        return azel,weight,midpoints

def draw_weighted_polygons(inpolys,weight=None,holes='cw',emptyweight=0,cmap='gray_r',facecolor=None,draw_colorbar=True, edgecolor='face',linewidth=.8, bgcolor='auto', vmin=None,vmax=None, **kwargs):
    if type(cmap) is str:
        try:
            cmap=vars(cm)[cmap]
        except:
            raise KeyError(cmap + ' is not a valid matplotlib colormap.')
    polys=polysplit(inpolys)
    if weight is not None:
        weights=weight
    else:
        weights=ones(len(polys))            
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
    if bgcolor is None:
        bgcolor=gca().get_axis_bgcolor()
    elif bgcolor=='auto':
        bgcolor=cmap(0)
    cmap.set_bad(bgcolor)
    p=PolyCollection(polys,cmap=cmap,edgecolor=edgecolor,facecolor=facecolor,linewidth=linewidth, **kwargs)
    if facecolor is None:
        if None in [vmin,vmax]:
            if vmin is None:
                vmin=masked_weight.min()
            if vmax is None: 
                vmax=masked_weight.max()
            if vmin==vmax:
                vmin=0
                draw_colorbar=False
            else:
                if vmin>0:
                    vmin=max(0,vmin-.1*(vmax-vmin))
                else:
                    vmin=vmin-.1*(vmax-vmin)
        p.set_clim(vmin=vmin,vmax=vmax)
        p.set_array(masked_weight)
    else:
        draw_colorbar=False        
    gca().add_collection(p)
    if draw_colorbar:
        colorbar(p)
    return p

def ispolyccw(poly):
    tot=sum([cross(*vecs) for vecs in zip(poly,poly[1:])])
    return tot<0

def ispolycw(poly):
    tot=sum([cross(*vecs) for vecs in zip(poly,poly[1:])])
    return tot>0
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
    """Return how many values the caller is expecting.  Emulates matlab's nargout."""
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
    trim=True
    if trimmer_poly is None:
        trim=False        
    if not isinstance(polys,mangle.Mangle):
        if polys[-5:]=='.list':
            trim=False            
    if trim:
        trimmed_polys=trim_mask(polys,trimmer_poly)
        azel,weight,midpoints=read_poly_list(trimmed_polys,graphicsfilename=graphicsfilename)
    else:
        azel,weight,midpoints=read_poly_list(polys,graphicsfilename=graphicsfilename) 
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
            corners_defined=(None not in [minaz,maxaz,minel,maxel])
            range_defined=(None not in [width,height,cenaz,cenel])
            
            if corners_defined:
                if projection is None:
                    if maxaz-minaz<90 and maxel-minel<90:
                        projection='gnom'
                    else:
                        projection='cyl'
                if cenel is None:
                    cenel = mean([minel,maxel])
                if cenaz is None:
                    cenaz = mean([minaz,maxaz])            
                if (projection in _pseudocyl) or (projection in _cylproj):
                    lon_0=cenaz                
                    llcrnrlon=minaz
                    urcrnrlon=maxaz
                else:
                    lon_0=-cenaz                
                    llcrnrlon=-maxaz
                    urcrnrlon=-minaz
                Basemap.__init__(self,llcrnrlon=llcrnrlon,llcrnrlat=minel,urcrnrlon=urcrnrlon,urcrnrlat=maxel,lon_0=lon_0,lat_0=cenel,projection=projection, rsphere=180./pi, celestial=True,**kwargs)
            elif range_defined:
                if projection is None:
                    if width<90 and height<90:
                        projection='gnom'
                    else:
                        projection='cyl'
                if (projection in _pseudocyl) or (projection in _cylproj):
                    lon_0=cenaz
                else:
                    lon_0=-cenaz
                Basemap.__init__(self,lon_0=lon_0,lat_0=cenel,width=width,height=height, projection=projection,rsphere=180./pi,celestial=True,**kwargs)
            else:
                if projection is None:
                    projection='mbtfpq'
                if cenaz is None:
                    cenaz=180.
                if cenel is None:
                    cenel=0.
                if (projection in _pseudocyl) or (projection in _cylproj):
                    lon_0=cenaz
                else:
                    lon_0=-cenaz          
                Basemap.__init__(self,lon_0=lon_0,lat_0=cenel,projection=projection, rsphere=180./pi,celestial=True,**kwargs)

            self.cenel=self.projparams['lat_0']
            if (projection in _pseudocyl) or (projection in _cylproj):
                self.cenaz=self.projparams['lon_0']
            else:
                self.cenaz=-self.projparams['lon_0']
            
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

            self.make_trimmer_poly()
        
        def mlabelzero22pi(lon):
            fmt='%g'
            if rcParams['text.usetex']:
                lonlabstr = r'${%s\/^{\circ}}$'%fmt
            else:
                lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
                lonlab = lonlabstr%lon
            return lonlab

        def drawmeridians(self,meridians,fmt=mlabelzero22pi,labels=[0,0,0,1], linewidth=1.,**kwargs):
            if self.projection in _cylproj:
                meridians1=(meridians-(self.cenaz-180))%360+(self.cenaz-180)
                Basemap.drawmeridians(self,meridians1,fmt=fmt,labels=[0,0,0,0],linewidth=linewidth,**kwargs)
                Basemap.drawmeridians(self,meridians,fmt=fmt,labels=labels,linewidth=0,**kwargs)
            else:
                Basemap.drawmeridians(self,meridians,fmt=fmt,labels=labels,linewidth=linewidth,**kwargs)

        def drawparallels(self,parallels,labels=[1,0,0,0], linewidth=1.,**kwargs):
            Basemap.drawparallels(self,parallels,labels=labels,linewidth=linewidth,**kwargs)
        
        def draw_grid(self,linewidth=1.,**kwargs):
            if self.projection in _pseudocyl:
                if self.projection != 'vandg':
                    self.drawparallels(np.arange(-90.,91.,30.),linewidth=linewidth,**kwargs)
                if self.projection == 'robin':
                    self.drawmeridians(np.arange(0.,360.,60.),linewidth=linewidth,**kwargs)
                elif self.projection == 'mbtfpq':
                    self.drawmeridians(np.arange(0.,360.,90.),linewidth=linewidth,**kwargs)
            elif self.azmax-self.azmin>359:
                self.drawparallels(np.arange(-90.,91.,30.),labels=[1, 0, 0, 0],linewidth=linewidth,**kwargs)
                self.drawmeridians(np.arange(0.,360.,60.),labels=[0, 0, 0, 1],linewidth=linewidth,**kwargs)
            else:
                autolocate_az=matplotlib.ticker.AutoLocator()
                autolocate_el=matplotlib.ticker.AutoLocator()
                self.drawparallels(autolocate_el.bin_boundaries(self.elmin,self.elmax),linewidth=linewidth,**kwargs)
                self.drawmeridians(autolocate_az.bin_boundaries(self.azmin,self.azmax),linewidth=linewidth,**kwargs)

        def world2pix(self,xy):
            xyout=vstack([Basemap.__call__(self,xy[:,0],xy[:,1],inverse=False)]).transpose()
            xyout[isnan(xy)]=nan
            return xyout

        def pix2world(self,xy):
            xyout=vstack([Basemap.__call__(self,xy[:,0],xy[:,1],inverse=True)]).transpose()
            xyout[isnan(xy)]=nan
            return xyout

        def make_trimmer_poly(self):
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
            trimmer_poly=mangle.Mangle(jpol.name,keep_ids=True)
            self.trimmer_poly=trimmer_poly

        def get_graphics_polygons(self,polys,pointsper2pi=30,graphicsfilename=None):
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
