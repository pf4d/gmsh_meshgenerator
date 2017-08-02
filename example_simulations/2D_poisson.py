from fenics import *

mesh = Mesh("../meshes/cellular_mesh.xml.gz")

# variational problem :
Q   = FunctionSpace(mesh, 'CG', 1)
u   = TrialFunction(Q)
v   = TestFunction(Q)

bc  = DirichletBC(Q, 0, lambda x, on_boundary : on_boundary)
f   = Constant(1)

a   = inner(grad(u), grad(v)) * dx
L   = f*v*dx

# solve variational problem
uh = Function(Q)

A = assemble(a)
b = assemble(L)

solve(a == L, uh, bc)

from pylab                   import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib              import colors, ticker

# generate figure :
v   = uh.compute_vertex_values(mesh)
x   = mesh.coordinates()[:,0]
y   = mesh.coordinates()[:,1]
t   = mesh.cells()

fig = figure(figsize=(8,7))
ax  = fig.add_subplot(111)
cm  = get_cmap('viridis')

c  = ax.tricontourf(x, y, t, v, 10, cmap=cm)
tp = ax.triplot(x, y, t, '-', color='k', lw=0.1, alpha=0.5)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.axis('equal')
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
  
divider = make_axes_locatable(gca())
cax  = divider.append_axes('right', "5%", pad="3%")
cbar = colorbar(c, cax=cax, format='%.1e') 
tight_layout()
savefig('2D_poisson.pdf')



