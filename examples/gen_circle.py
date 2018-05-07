from gmsh_meshgenerator import MeshGenerator
from pylab              import *

msh_name = 'circle_mesh'
out_dir  = 'meshes/'

x = linspace(-1.0, 1.0, 100)
y = linspace(-1.0, 1.0, 100)

X,Y = meshgrid(x,y)

S = 1 - sqrt(X**2 + Y**2)

m = MeshGenerator(x, y, msh_name, out_dir)

m.create_contour(S, zero_cntr=1e-16, skip_pts=0)
m.create_contour(S, zero_cntr=0.5,   skip_pts=0)
#m.plot_contour()

m.write_gmsh_contour(lc=0.1, boundary_extend=False)

m.add_edge_attractor(field=0, contour_index=0, NNodesByEdge=10)
m.add_edge_attractor(field=1, contour_index=1, NNodesByEdge=10)

#field, ifield, lcMin, lcMax, distMin, distMax
m.add_threshold(2, 0, 0.001, 0.1, 0, 0.25)
m.add_threshold(3, 1, 0.001, 0.1, 0, 0.25)

m.finish()

m.create_mesh()
m.convert_msh_to_xml()


