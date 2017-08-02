from gmsh_meshgenerator import MeshGenerator
from pylab              import *

msh_name = 'square_mesh'
out_dir  = 'meshes/'

m = MeshGenerator(None, None, msh_name, out_dir)

m.add_contour(array([[0,0],[0,1],[1,1],[1,0]]))
m.write_gmsh_contour(lc=1.0, boundary_extend=False)

m.add_edge_attractor(field=0, contour_index=0, NNodesByEdge=100)

#field, ifield, lcMin, lcMax, distMin, distMax
m.add_threshold(1, 0, 0.001, 0.1, 0, 0.25)
m.finish()

m.create_mesh()
m.convert_msh_to_xml()






