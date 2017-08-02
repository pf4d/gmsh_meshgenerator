from gmsh_meshgenerator import MeshGenerator
from pylab              import *

msh_name = 'cellular_mesh'                # name of all outputs
out_dir  = 'meshes/'                      # directory for all outputs
L        = 1.0                            # characteristic length of domain :
x        = linspace(0, L, 100)            # x-coordinates
y        = linspace(0, L, 100)            # y-coordinates
X,Y      = meshgrid(x,y)                  # matrix of coordinates
S        = sin(4*pi*X/L) * sin(4*pi*Y/L)  # mathematical function to contour

# the MeshGenerator instance initialized with the matrix coordinates,
# output file name, and output directory :
m = MeshGenerator(x, y, msh_name, out_dir)

# add the exterior contour, a box :
m.add_contour(array([[0,0],[0,L],[L,L],[L,0]]))

# create and add contours from S.  Parameter "skip_pts" are the number of 
# contour nodes to remove from the resulting contour, if a lower-resolution
# contour is desired.  Parameter "distance_check" is the distance from each node
# to check for edge intersections which may arise from removing points; 
# increase the distance if overlaps occur :
m.create_contour(S, zero_cntr= 0.2, skip_pts=0, distance_check = 10)
m.create_contour(S, zero_cntr=-0.2, skip_pts=0, distance_check = 10)

# plot the resulting contours :
m.plot_contour()

# a new contour to take the intersection from :
lmin = 0.10
lmax = 0.90
outer_contour = array([[lmin,lmin],[lmin,lmax],[lmax,lmax],[lmax,lmin]])

# this function replaces the current set of contours with the intersection 
# of the current contours with "outer_contour" :
m.intersection(outer_contour)

# plot the modified contour, to be sure.  The contours are labeled by
# "contour_index" used by m.add_edge_attractor() below :
m.plot_contour()

# write the gmsh contour to the "msh_name".geo file with characteristic 
# cell diameter "lc".  If "boundary_extend"=True, the edge size of the contour 
# is extrapolated into the interior :
m.write_gmsh_contour(lc=0.1, boundary_extend=False)

# get the number of contours for the mesh :
num_ctrs = m.num_contours()

# add identical edge attractors to all contours with sequential identifier.
# the parameter "NNodesByEdge" adds extra nodes to the edge for the purpose
# of refinement, increase this if your edges are much longer than your 
# cell refinement "lcMin" set by m.add_threshold() below :
for i in range(num_ctrs):
  m.add_edge_attractor(field=i, contour_index=i, NNodesByEdge=10)

# for each edge attractor, add a threshold for mesh refinement in the vicinity 
# of the edge :
for i in range(num_ctrs):
  # parameters are respectively: field, ifield, lcMin, lcMax, distMin, distMax
  m.add_threshold(num_ctrs+i, i, 0.001, 0.1, 0, 0.25)

# finialize the creation of the "msh_name".geo file :
m.finish()

# instruct gmsh to create the mesh, saving to a "msh_name".msh 
# file in "out_dir" :
m.create_mesh()

# convert the "msh_name".msh file to a "msh_name".xml.gz file used by FEniCS :
m.convert_msh_to_xml()



