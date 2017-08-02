from fenics import *

mesh = Mesh("../meshes/cellular_mesh.xml.gz")
    
xmin = MPI.min(mpi_comm_world(), mesh.coordinates()[:,0].min())
xmax = MPI.max(mpi_comm_world(), mesh.coordinates()[:,0].max())
ymin = MPI.min(mpi_comm_world(), mesh.coordinates()[:,1].min())
ymax = MPI.max(mpi_comm_world(), mesh.coordinates()[:,1].max())

# Sub domain for inflow (right)
class Right(SubDomain):
  def inside(self, x, on_boundary):
    return abs(x[0] - xmax) < DOLFIN_EPS and on_boundary
right = Right()

# Sub domain for inflow (right)
class Left(SubDomain):
  def inside(self, x, on_boundary):
    return abs(x[0] - xmin) < DOLFIN_EPS and on_boundary
left = Left()

# Sub domain for inflow (right)
class Top(SubDomain):
  def inside(self, x, on_boundary):
    return abs(x[1] - ymax) < DOLFIN_EPS and on_boundary
top = Top()

# Sub domain for inflow (right)
class Bottom(SubDomain):
  def inside(self, x, on_boundary):
    return abs(x[1] - ymin) < DOLFIN_EPS and on_boundary
bottom = Bottom()

sub_domains = FacetFunction('size_t', mesh, 0)
right.mark(sub_domains,  1)
left.mark(sub_domains,   2)
top.mark(sub_domains,    3)
bottom.mark(sub_domains, 4)

# Taylor-Hood element
Q1e   = FiniteElement("CG", mesh.ufl_cell(), 1)
Q2e   = FiniteElement("CG", mesh.ufl_cell(), 2)
QTHe  = MixedElement([Q2e, Q2e, Q1e])
Q     = FunctionSpace(mesh, 'CG', 1)
W     = FunctionSpace(mesh, QTHe)

# variational problem :
U   = TrialFunction(W)
Phi = TestFunction(W)

u_x,   u_y,   p = U
phi_x, phi_y, q = Phi

u = as_vector([u_x,   u_y])
v = as_vector([phi_x, phi_y])

# no penetration boundary condition for velocity :
u_n   = Constant(0.0)

# inflow and outflow boundary conditions :
u_r   = Constant(( -1, 0.0))
u_l   = Constant((  1, 0.0))
u_t   = Constant((0.0,   1))
u_b   = Constant((0.0,  -1))

# relavent measures :
ds     = Measure("ds", subdomain_data = sub_domains)
dG_0   = ds(0)
dG_r   = ds(1)
dG_l   = ds(2)
dG_t   = ds(3)
dG_b   = ds(4)
dG_in  = dG_r + dG_l
dG_out = dG_t + dG_b
dG_ext = dG_in + dG_out

# constants :
gamma = Constant(1e2)
h     = CellSize(mesh)
n     = FacetNormal(mesh)
I     = Identity(2)
eta   = Constant(1.0)
f     = Constant((0.0,0.0))
beta  = Constant(10.0)

def epsilon(u): return 0.5*(grad(u) + grad(u).T)
def sigma(u,p): return 2*eta*epsilon(u) - p*I

t   = dot(sigma(u,p), n)
s   = dot(sigma(v,q), n)

B_o = + inner(sigma(u,p),grad(v))*dx - div(u)*q*dx

B_g = - dot(n,t) * dot(v,n) * dG_0 \
      - dot(u,n) * dot(s,n) * dG_0 \
      + gamma/h * dot(u,n) * dot(v,n) * dG_0 \
      + beta * dot(u, v) * dG_0 \
      - inner(dot(sigma(u,p), n), v) * dG_ext \
      - inner(dot(sigma(v,q), n), u) * dG_ext \
      + gamma/h * inner(v,u) * dG_ext

F   = + dot(f,v) * dx \
      + gamma/h * u_n * dot(v,n) * dG_0 \
      - inner(dot(sigma(v,q), n), u_r) * dG_r \
      + gamma/h * inner(v,u_r) * dG_r \
      - inner(dot(sigma(v,q), n), u_l) * dG_l \
      + gamma/h * inner(v,u_t) * dG_t \
      - inner(dot(sigma(v,q), n), u_t) * dG_t \
      + gamma/h * inner(v,u_l) * dG_l \
      - inner(dot(sigma(v,q), n), u_b) * dG_b \
      + gamma/h * inner(v,u_b) * dG_b \

# solve variational problem
wh = Function(W)

A = assemble(B_o + B_g)
b = assemble(F)

solver = LUSolver('mumps')
solver.set_operator(A)

solver.solve(wh.vector(), b)
u0, u1, ph = wh.split(True)

from pylab                   import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib              import colors, ticker

# calculate array componets :
v0  = u0.compute_vertex_values(mesh)
v1  = u1.compute_vertex_values(mesh)
v   = sqrt(v0**2 + v1**2 + 1e-16)
v0  = v0 / v
v1  = v1 / v
x   = mesh.coordinates()[:,0]
y   = mesh.coordinates()[:,1]
t   = mesh.cells()

# generate velocity figure :
fig = figure(figsize=(8,7))
ax  = fig.add_subplot(111)

v[v > 2.0] = 2.0
cm = get_cmap('viridis')
c  = ax.tricontourf(x, y, t, v, cmap=cm)
q  = ax.quiver(x, y, v0, v1, pivot='middle',
                             color='k',
                             scale=60,
                             width=0.0015,
                             headwidth=4.0, 
                             headlength=4.0, 
                             headaxislength=4.0)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.axis('equal')
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_xticklabels([])
ax.set_yticklabels([])
  
divider = make_axes_locatable(gca())
cax  = divider.append_axes('right', "5%", pad="3%")
cbar = fig.colorbar(c, cax=cax, format='%.1f') 
tight_layout()
savefig('2Dstokes_nitsche_u.pdf')

# generate pressure figure :
v  = ph.compute_vertex_values(mesh)

fig = figure(figsize=(8,7))
ax  = fig.add_subplot(111)

c  = ax.tricontourf(x, y, t, v, 10, cmap=cm)
tp = ax.triplot(x, y, t, '-', color='k', lw=0.1, alpha=0.5)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.axis('equal')
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_xticklabels([])
ax.set_yticklabels([])
  
divider = make_axes_locatable(gca())
cax  = divider.append_axes('right', "5%", pad="3%")
cbar = colorbar(c, cax=cax) 
tight_layout()
savefig('2Dstokes_nitsche_p.pdf')



