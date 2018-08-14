from dolfin import *
from numpy import sqrt

parameters["allow_extrapolation"] = True

def kjell(self,u):
    return pow(inner(u, u), 0.5)

def epsilon(u):
    'Return strain-rate tensor'
    return 0.5*(grad(u) + transpose(grad(u)))

class STRESS:
    'Computation of stress for given u, p'

    def __init__(self, u, p, nu, mesh):
        bm = BoundaryMesh(mesh, 'exterior')
        self.bmV = VectorFunctionSpace(bm, 'CG', 1)

        # Compute stress tensor
        sigma = (2*nu*epsilon(u)) - (p*Identity(len(u)))

        # Compute stress on surface
        n = FacetNormal(mesh)
        F = -(sigma*n)

        # Compute normal and tangential components
        Fn = inner(F,n)       # scalar-valued
        Ft = F - (Fn*n) # vector-valued

        # Integrate against piecewise constants on the boundary
        scalar  = FunctionSpace(mesh, 'DG', 0)
        vector  = VectorFunctionSpace(mesh, 'CG', 1)
        scaling = FacetArea(mesh) # Normalise the computed stress relative to the size of the element

        v  = TestFunction(scalar)
        w  = TestFunction(vector)

        # Create functions
        self.Fn  = Function(scalar)
        self.Ftv = Function(vector)
        self.Ft  = Function(scalar)

        self.Ln  = 1/scaling*v*Fn*ds
        self.Ltv = 1/(2*scaling)*inner(w, Ft)*ds
        self.Lt  = 1/scaling*inner(v, kjell(self,self.Ftv))*ds

    def __call__(self, u):
        'Return  shear stress (Ft, Ftv)'

        # Assemble vectors
        assemble(self.Ltv, tensor=self.Ftv.vector())
        self.Ftv_bm = interpolate(self.Ftv, self.bmV)

        return self.Ftv_bm
