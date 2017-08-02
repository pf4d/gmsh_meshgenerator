from pylab import *
import subprocess
import os

from colored           import fg, attr
from shapely.geometry  import Polygon
from shapely.geometry  import Point as shapelyPoint
from shapely.ops       import cascaded_union


def get_text(text, color='white', atrb=0, cls=None):
  """
  Returns text ``text`` from calling class ``cls`` for printing at a later time.

  :param text: the text to print
  :param color: the color of the text to print
  :param atrb: attributes to send use by ``colored`` package
  :param cls: the calling class
  :type text: string
  :type color: string
  :type atrb: int
  :type cls: object
  """
  if cls is not None:
    color = cls.color()
  if atrb != 0:
    text = ('%s%s' + text + '%s') % (fg(color), attr(atrb), attr(0))
  else:
    text = ('%s' + text + '%s') % (fg(color), attr(0))
  return text


def print_text(text, color='white', atrb=0, cls=None):
  """
  Print text ``text`` from calling class ``cls`` to the screen.

  :param text: the text to print
  :param color: the color of the text to print
  :param atrb: attributes to send use by ``colored`` package
  :param cls: the calling class
  :type text: string
  :type color: string
  :type atrb: int
  :type cls: object
  """
  if cls is not None:
    color = cls.color()
  if atrb != 0:
    text = ('%s%s' + text + '%s') % (fg(color), attr(atrb), attr(0))
  else:
    text = ('%s' + text + '%s') % (fg(color), attr(0))
  print text



class MeshGenerator(object):
  """
  generate a mesh.
  """
  def __init__(self, x, y, fn, direc):
    """
    Generate a mesh with DataInput object <dd>, output filename <fn>, and
    output directory <direc>.
    """
    #self.color = 'grey_46'
    self.color = '43'
    s    = "::: INITIALIZING MESHGENERATOR :::"
    print_text(s, self.color)
    self.fn          = fn
    self.direc       = direc
    self.x, self.y   = x, y
    if not os.path.exists(direc):
      os.makedirs(direc)
    self.f           = open(direc + fn + '.geo', 'w')
    self.fieldList   = []  # list of field indexes created.
    self.contourlist = []  # list of contours
    self.loop_a      = []  # list of loop-list strings
    self.dim         = 2   # default to 2D mesh

  def create_contour(self, var, zero_cntr, skip_pts, distance_check):
    """
    Create a contour of the data field with index <var> of <dd> provided at
    initialization.  <zero_cntr> is the value of <var> to contour, <skip_pts>
    is the number of points to skip in the contour, needed to prevent overlap.
    """
    s    = "::: creating contour skipping %i point(s) :::"
    print_text(s % skip_pts, self.color)

    skip_pts = skip_pts + 1

    # create contour :
    fig = figure()
    ax  = fig.add_subplot(111)
    ax.set_aspect('equal')
    cnt = ax.contour(self.x, self.y, var, [zero_cntr])
    close(fig)

    # remove skip points and last point to avoid overlap :
    for c in cnt.allsegs[0]:
      s    = "    - contour created, length %s nodes -"
      print_text(s % shape(c)[0], self.color)
      cont = self.remove_skip_points(c, skip_pts)
      cont = self.eliminate_intersections(cont, distance_check)
      self.contourlist.append(cont)

  def num_contours(self):
    """
    Return the integer number of contours for this mesh.
    """
    return len(self.contourlist)

  def remove_skip_points(self, contour, skip_pts):
    """
    remove every other <skip_pts> node from <contour>.
    """
    # remove skip points and last point to avoid overlap :
    cont = contour[::skip_pts,:][:-1,:]
    s    = "    - contour points skipped, new length %s nodes -"
    print_text(s % shape(cont)[0], self.color)
    return cont

  def add_contour(self, cont_array):
    """
    This is an alternative to the create_contour method that allows you to
    manually specify contour points.
    Inputs:
    cont_array : A numpy array of contour points (i.e. array([[1,2],[3,4],...]))
    """
    s = "::: manually setting contour with %s nodes:::"
    print_text(s % shape(cont_array)[0], self.color)
    self.contourlist.append(cont_array)

  def plot_contour(self):
    """
    Plot the contour created with the "create_contour" method.
    """
    s = "::: plotting contour :::"
    print_text(s, self.color)
    
    cmap = get_cmap('jet')
    fig  = figure()
    ax   = fig.add_subplot(111)
    
    colors = [cmap(x) for x in np.linspace(0,1,len(self.contourlist))]
    for i,(cnt, col) in enumerate(zip(self.contourlist, colors)):
      ax.plot(cnt[:,0], cnt[:,1], '-', color=col, lw = 3.0, label=i)
    
    ax.set_aspect('equal')
    ax.set_title("contour")
    leg = ax.legend(bbox_to_anchor=(1.2,0.5), loc='center right', ncol=1)
    leg.get_frame().set_alpha(0.0)
    leg.get_frame().set_color('w')
    
    plt.tight_layout(rect=[0.005,0.01,0.88,0.995])
    show()

  def eliminate_intersections(self, contour, dist=10):
    """
    Eliminate intersecting boundary elements. <dist> is an integer specifiying
    how far forward to look to eliminate intersections.  If any intersections
    are found, this method is called recursively until none are found.
    """
    s    = "::: eliminating intersections :::"
    print_text(s, self.color)

    class Point:
      def __init__(self,x,y):
        self.x = x
        self.y = y

    def ccw(A,B,C):
      return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

    def intersect(A,B,C,D):
      return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

    flag = ones(len(contour))
    intr = False
    for ii in range(len(contour)-1):

      A = Point(*contour[ii])
      B = Point(*contour[ii+1])

      for jj in range(ii, min(ii + dist, len(contour)-1)):

        C = Point(*contour[jj])
        D = Point(*contour[jj+1])

        if intersect(A,B,C,D) and ii!=jj+1 and ii+1!=jj:
          s    = "    - intersection found between node %i and %i -"
          print_text(s % (ii+1, jj), 'red')
          flag[ii+1] = 0
          flag[jj]   = 0
          intr       = True

    counter  = 0
    new_cont = zeros((int(sum(flag)),2))
    for ii,fl in enumerate(flag):
      if fl:
        new_cont[counter,:] = contour[ii,:]
        counter += 1

    s    = "    - eliminated %i nodes -"
    print_text(s % sum(flag == 0), self.color)
    # call again if required :
    if intr:
      self.eliminate_intersections(new_cont, dist)
    return new_cont

  def restart(self):
    """
    clear all contents from the .geo file.
    """
    self.f.close
    self.f = open(self.direc + self.fn + '.geo', 'w')
    s = 'Reopened \"' + self.direc + self.fn + '.geo\".'
    print_text(s, self.color)

  def write_gmsh_contour(self, lc=1.0, boundary_extend=True):
    """
    write the contour created with create_contour to the .geo file with mesh
    spacing <lc>.  If <boundary_extend> is true, the spacing in the interior
    of the domain will be the same as the distance between nodes on the contour.
    """
    #FIXME: sporadic results when used with ipython, does not stops writing the
    #       file after a certain point.  calling restart() then write again
    #       results in correct .geo file written.  However, running the script
    #       outside of ipython works.
    s    = "::: writing gmsh contour to \"%s%s.geo\" :::"
    print_text(s % (self.direc, self.fn), self.color)
      
    f    = self.f
      
    # write the file to .geo file :
    f.write("// Mesh spacing\n")
    f.write("lc = %f;\n\n" % lc)

    ctr  = 0
    lp_a = []
    for j,c in enumerate(self.contourlist):
      
      pts = size(c[:,0])
      
      f.write("// Points\n")
      for i in range(pts):
        f.write("Point(%i) = {%f,%f,0,lc};\n" % (ctr+i, c[i,0], c[i,1]))
      
      f.write("\n// Lines\n")
      loop = "{"
      for i in range(pts-1):
        f.write("Line(%i) = {%i,%i};\n" % (ctr+i, ctr+i, ctr+i+1))
        loop += "%i," % (ctr+i)
      loop += "%i}" % (ctr+pts-1)
      f.write("Line(%i) = {%i,%i};\n" % (ctr+pts-1, ctr+pts-1, ctr))
      
      f.write("// Line loop\n")
      f.write("Line Loop(%i) = %s;\n\n" % (j, loop))
      self.loop_a.append(loop)

      lp_a.append(j)
      ctr += pts
      
    f.write("// Surface\n")
    srf_loop = "{"
    for i in lp_a[:-1]:
      srf_loop += "%i," % i
    srf_loop += "%i}" % lp_a[-1]
    f.write("Plane Surface(0) = %s;\n\n" % srf_loop)
      
    if not boundary_extend:
      f.write("Mesh.CharacteristicLengthExtendFromBoundary = 0;\n\n")

  def extrude(self, h, n_layers):
    """
    Extrude the mesh <h> units with <n_layers> number of layers.
    """
    s    = "::: extruding gmsh contour %i layers :::" % n_layers
    print_text(s, self.color)
    
    self.dim = 3   # now we have a 3D mesh
    self.f.write("Extrude {0,0,%i}{Surface{0};Layers{%i};}\n\n" % (h, n_layers))

  def add_box(self, field, vin, xmin, xmax, ymin, ymax, zmin, zmax):
    """
    add a box to the mesh.  e.g. for Byrd Glacier data:

      add_box(10000, 260000, 620000, -1080000, -710100, 0, 0)

    """
    f  = self.f
    fd = str(field)

    f.write("Field[" + fd + "]      =  Box;\n")
    f.write("Field[" + fd + "].VIn  =  " + float(vin)  + ";\n")
    f.write("Field[" + fd + "].VOut =  lc;\n")
    f.write("Field[" + fd + "].XMax =  " + float(xmax) + ";\n")
    f.write("Field[" + fd + "].XMin =  " + float(xmin) + ";\n")
    f.write("Field[" + fd + "].YMax =  " + float(ymax) + ";\n")
    f.write("Field[" + fd + "].YMin =  " + float(ymin) + ";\n")
    f.write("Field[" + fd + "].ZMax =  " + float(zmax) + ";\n")
    f.write("Field[" + fd + "].ZMin =  " + float(zmin) + ";\n\n")

    self.fieldList.append(field)

  def add_edge_attractor(self, field, contour_index, NNodesByEdge=10):
    """
    """
    fd = str(field)
    f  = self.f

    loop = self.loop_a[contour_index]

    f.write("Field[" + fd + "]              = Attractor;\n")
    f.write("Field[" + fd + "].EdgesList    = " + loop + ";\n")
    f.write("Field[" + fd + "].NNodesByEdge = %i;\n\n" % NNodesByEdge)

  def add_threshold(self, field, ifield, lcMin, lcMax, distMin, distMax):
    """
    """
    fd = str(field)
    f  = self.f

    f.write("Field[" + fd + "]         = Threshold;\n")
    f.write("Field[" + fd + "].IField  = " + str(ifield)  + ";\n")
    f.write("Field[" + fd + "].LcMin   = " + str(lcMin)   + ";\n")
    f.write("Field[" + fd + "].LcMax   = " + str(lcMax)   + ";\n")
    f.write("Field[" + fd + "].DistMin = " + str(distMin) + ";\n")
    f.write("Field[" + fd + "].DistMax = " + str(distMax) + ";\n\n")

    self.fieldList.append(field)

  def finish(self):
    """
    figure out background field and close the .geo file.
    """
    f     = self.f
    flist = self.fieldList

    # get a string of the fields list :
    l = ""
    for i,j in enumerate(flist):
      l += str(j)
      if i != len(flist) - 1:
        l += ', '

    # make the background mesh size the minimum of the fields :
    if len(flist) > 0:
      idx = str(max(flist)+1)
      f.write("Field[" + idx + "]            = Min;\n")
      f.write("Field[" + idx + "].FieldsList = {" + l + "};\n")
      f.write("Background Field    = " + idx + ";\n\n")
    else:
      f.write("Background Field = 0;\n\n")

    s = 'finished, closing \"' + self.direc + self.fn + '.geo\".'
    print_text(s, self.color)
    f.close()

  def close_file(self):
    """
    close the .geo file down for further editing.
    """
    s    = '::: finished, closing \"' + self.direc + self.fn + '.geo\" :::'
    print_text(s, self.color)
    self.f.close()

  def create_mesh(self):
    """
    create the mesh to file.
    """
    cmd = 'gmsh -%i %s%s.geo' % (self.dim, self.direc, self.fn)
    s = "\nExecuting :\n\n\t%s\n\n" % cmd
    print_text(s, self.color)
    subprocess.call(cmd.split())

  def check_dist(self, dist=1.0):
    """
    remove points in contour that are not a linear distance of at least
    <dist> from previous point.
    """
    s    = "::: ensuring that no nodes are less than %f from each other :::"
    print_text(s % dist, self.color)
    lin_dist = lambda p1, p2: sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    new_contourlist = []
    for coords in self.contourlist:
      mask = ones(len(coords), dtype=bool)

      i = 0
      while(i < len(coords)-1):
        p1 = xycoords[i]
        j = i + 1
        while(j < len(coords) and \
              lin_dist(p1, coords[j]) < dist):
          mask[j] = 0
          j += 1
        i = j

      # fix end of array
      i = -1
      while(len(coords) + i >= 0 and (not mask[i] or \
            lin_dist(coords[0], coords[i]) < dist)):
        mask[i] = 0
        i -= 1

      new_contourlist.append(coords[mask])

      # print results
      s    = "    - removed %i points closer than %f to one another -"
      print_text(s % (len(mask) - sum(mask), dist), self.color)


    self.longest_cont = new_contourlist

  def intersection(self, new_contour):
    """
    Take the geometric intersection of current coordinates with <new_contour>.
    """
    p2 = Polygon(zip(new_contour[:,0], new_contour[:,1]))
      
    exterior = zip(self.contourlist[0][:,0], self.contourlist[0][:,1])

    if len(self.contourlist) > 1:
      interior = []
      for c in self.contourlist[1:]:
        interior.append(zip(c[:,0], c[:,1]))
      p1 = Polygon(exterior, interior)
    else:
      p1 = Polygon(exterior)
  
    # calculate the intersection :
    p3 = p1.intersection(p2)
 
    # convert the Shapely exterior into the format for GMSH : 
    exterior  = array(zip(p3.exterior.xy[:][0], p3.exterior.xy[:][1]))[1:]
    
    s    = "::: intersection contour created, length %s nodes :::"
    print_text(s % shape(exterior)[0], self.color)

    # create a new list of contours from the intersection and its interior :
    self.contourlist = [exterior]
    interiors        = []
    for interior in p3.interiors:
      int_contour = array(zip(interior.xy[:][0], interior.xy[:][1]))[1:]
      self.contourlist.append(int_contour)
    
      s    = "    - interior contour, length %s nodes -"
      print_text(s % shape(interior)[0], self.color)
        


  def extend_edge(self, r):
    """
    Extends a 2d contour out from points labeled in self.edge by a distance
    <r> (radius) in all directions.
    """
    xycoords = self.longest_cont

    polygons = []
    for i, v in enumerate(xycoords):
      polygons.append(shapelyPoint(v[0],v[1]).buffer(r))

    # union of our original polygon and convex hull
    p1 = cascaded_union(polygons)
    p3 = cascaded_union(p1)

    xycoords_buf = array(zip(p3.exterior.xy[:][0], p3.exterior.xy[:][1]))
    self.longest_cont = xycoords_buf

  def convert_msh_to_xml(self):
    """
    convert <mshfile> .msh file to .xml file <xmlfile> via dolfin-convert.
    """
    msh = self.direc + self.fn + '.msh'
    xml = self.direc + self.fn + '.xml'

    cmd = 'dolfin-convert ' + msh + ' ' + xml
    s   = "\nExecuting :\n\n\t %s\n\n" % cmd
    print_text(s, self.color)
    subprocess.call(cmd.split())
    
    cmd = 'gzip -f ' + xml
    s   = "\nExecuting :\n\n\t %s\n\n" % cmd
    print_text(s, self.color)
    subprocess.call(cmd.split())

  def convert_msh_to_xdmf(self):
    """
    convert <mshfile> .msh file to .xdmf file <xdmffile>.
    """
    self.convert_msh_to_xml(mshfile, xdmffile)
    
    xdmf      = self.direc + self.fn + '.xdmf'
    xml       = self.direc + self.fn + '.xml.gz'
    mesh      = Mesh(xml)
    mesh_file = XDMFFile(mesh.mpi_comm(), xdmf)
    mesh_file.write(mesh)


class linear_attractor(object):
  r"""
  Create an attractor object which refines with min and max cell radius 
  :math:`l_{min}`, :math:`l_{max}` over data field :math:`f`.  The 
  :math:`f_{max}` parameter specifies a max value for which to apply the 
  minimum cell size such that if :math:`f_i` is less than :math:`f_{max}`,
  the cell size in this region will be :math:`l_{max}`.  If *inv* = ``True``
  the object refines on the inverse of the data field :math:`f`.

  .. math::

     h_i = \begin{cases}
             l_{min},     & f_i > f_{max} \\
             l_{max},     & f_i < f_{max} \\
             f_i,         & otherwise \\
           \end{cases}

  Args:

    :spline: the iterpolator which is evaluated at x- and y-coordinates
    :f:      the :class:`~numpy.array` being interpolated
    :f_max:  maximum value of *f* for *l_max* and *l_min*
    :l_min:  minimum cell size
    :l_max:  maximum cell size
    :inv:    boolean, invert *f* for attraction

  """
  def __init__(self, spline, f, f_max, l_min, l_max, inv=True):
    """
    Refine the mesh off of data field <f> using spline <spline> with the
    cell radius defined as :

               {l_min,     f_i > f_max
    cell_h_i = {l_max,     f_i < f_max
               {f_i,       otherwise

    If <inv> is True, refine off of the inverse of <f> instead.

    """
    self.spline   = spline
    self.field    = f
    self.l_min    = l_min
    self.l_max    = l_max
    self.f_max    = f_max
    self.inv      = inv

  def op(self, x, y, z, entity):
    """
    Method which evaluates this linear attractor.

    Args:

      :x:      the x-coordinate
      :y:      the y-coordinate
      :z:      the z-coordinate (not used)
      :entity: not used

    Returns:
    
      :lc:     characteristic radius for a given cell at *(x,y)*
    """
    l_min = self.l_min
    l_max = self.l_max
    f     = self.field
    v     = self.spline(x,y)[0][0]
    if self.inv:
      if v < self.f_max:
        lc = l_max - (l_max - l_min) / f.max() * v
      else:
        lc = l_min
    else:
      if v < self.f_max:
        lc = l_min + (l_max - l_min) / f.max() * v
      else:
        lc = l_max
    return lc

class static_attractor(object):
  """
  """
  def __init__(self, spline, c, inv=False):
    """
    Refine the mesh off of data field <spline> with the cell radius
    defined as :

    cell_h_i = c * spline(x,y)

    """
    self.spline = spline
    self.c      = c
    self.inv    = inv

  def op(self, x, y, z, entity):
    """
    """
    if not self.inv:
      lc = self.c * self.spline(x,y)[0][0]
    else:
      lc = self.c * 1/self.spline(x,y)[0][0]
    return lc


class min_field(object):
  """
  Return the minimum of a list of attactor operator fields <f_list>.
  """
  def __init__(self, f_list):
    self.f_list = f_list

  def op(self, x, y, z, entity):
    l = []
    for f in self.f_list:
      l.append(f(x,y,z,entity))
    return min(l)


class max_field(object):
  """
  Return the minimum of a list of attactor operator fields <f_list>.
  """
  def __init__(self, f_list):
    self.f_list = f_list

  def op(self, x, y, z, entity):
    l = []
    for f in self.f_list:
      l.append(f(x,y,z,entity))
    return max(l)


class MeshRefiner(object):

  def __init__(self, di, fn, gmsh_file_name):
    """
    Creates a 2D or 3D mesh based on contour .geo file <gmsh_file_name>.
    Refinements are done on DataInput object <di> with data field index <fn>.
    """
    self.color = '43'
    s    = "::: initializing MeshRefiner on \"%s.geo\" :::" % gmsh_file_name
    print_text(s, self.color)

    self.field  = di.data[fn].T
    print_min_max(self.field, 'refinement field [m]')

    self.spline = RectBivariateSpline(di.x, di.y, self.field, kx=1, ky=1)

    #load the mesh into a GModel
    self.m = GModel.current()
    self.m.load(gmsh_file_name + '.geo')

    # set some parameters :
    GmshSetOption("Mesh", "CharacteristicLengthFromPoints", 0.0)
    GmshSetOption("Mesh", "CharacteristicLengthExtendFromBoundary", 0.0)
    GmshSetOption("Mesh", "Smoothing", 100.0)

  def add_linear_attractor(self, f_max, l_min, l_max, inv):
    """
    Refine the mesh with the cell radius defined as :

               {l_min,     field_i > f_max
    cell_h_i = {l_max,     field_i < f_max
               {field_i,   otherwise

    If <inv> is True, refine off of the inverse of <field> instead.

    """
    # field, f_max, l_min, l_max, hard_cut=false, inv=true
    a   = linear_attractor(self.spline, self.field, f_max, l_min, l_max,
                           inv=inv)
    aid = self.m.getFields().addPythonField(a.op)
    return a,aid

  def add_static_attractor(self, c=1, inv=False):
    """
    Refine the mesh with the cell radius defined as :

    cell_h_i = c * field_i

    returns a tuple, static_attractor object and id number.

    """
    # field, f_max, l_min, l_max, hard_cut=false, inv=true
    a   = static_attractor(self.spline, c, inv)
    aid = self.m.getFields().addPythonField(a.op)
    return a,aid

  def add_min_field(self, op_list):
    """
    Create a miniumum field of attactor operator lists <op_list>.
    """
    mf  = min_field(op_list)
    mid = self.m.getFields().addPythonField(mf.op)
    return mid

  def set_background_field(self, idn):
    """
    Set the background field to that of field index <idn>.
    """
    self.m.getFields().setBackgroundFieldId(idn)

  def finish(self, gui=True, dim=3, out_file_name='mesh'):
    """
    Finish and create the .msh file.  If <gui> is True, run the gui program,
    Otherwise, create the .msh file with dimension <dim> and filename
    <out_file_name>.msh.
    """
    self.out_file_name = out_file_name

    #launch the GUI
    if gui:
      print_text("::: opening GUI :::", self.color)
      FlGui.instance().run()

    # instead of starting the GUI, we could generate the mesh and save it
    else:
      s    = "::: writing %s.msh :::" % out_file_name
      print_text(s, self.color)
      self.m.mesh(dim)
      self.m.save(out_file_name + ".msh")
  
  def convert_msh_to_xml(self):
    """
    convert ``self.out_file_name``.msh file to .xml file 
    ``self.out_file_name``.xml via dolfin-convert.
    """
    msh = self.out_file_name + '.msh'
    xml = self.out_file_name + '.xml'

    cmd = 'dolfin-convert ' + msh + ' ' + xml
    s   = "\nExecuting :\n\n\t %s\n\n" % cmd
    print_text(s, self.color)
    subprocess.call(cmd.split())
    
    cmd = 'gzip -f ' + xml
    s   = "\nExecuting :\n\n\t %s\n\n" % cmd
    print_text(s, self.color)
    subprocess.call(cmd.split())

  def convert_msh_to_xdmf(self):
    """
    convert ``self.out_file_name``.msh file to .xdmf file 
    ``self.out_file_name``.xdmf via dolfin-convert.
    """
    self.convert_msh_to_xml()

    xdmf      = self.out_file_name + '.xdmf'
    xml       = self.out_file_name + '.xml.gz'
    mesh      = Mesh(xml)
    mesh_file = XDMFFile(mesh.mpi_comm(), xdmf)
    mesh_file.write(mesh)



