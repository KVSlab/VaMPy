# This file is modified and simplified from Womersely


from dolfin import (UserExpression, Mesh, MeshFunction, SubsetIterator, MPI, ds,
assemble, Constant, sqrt, FacetNormal, as_vector, SpatialCoordinate)

import numpy as np

from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps, romberg
from scipy.special import jn
import math

def x_to_r2(x, c, n):
    """Compute r**2 from a coordinate x, center point c, and normal vector n.

        r is defined as the distance from c to x', where x' is
        the projection of x onto the plane defined by c and n.
        """
    # Steps:
    # rv = x - c
    # rvn = rv . n
    # rp = rv - (rv . n) n
    # r2 = ||rp||**2

    rv = x-c
    rvn = rv.dot(n)
    rp = rv - rvn*n
    r2 = rp.dot(rp)

    return r2

def compute_radius(mesh, facet_domains, ind, center):
    d = len(center)
    it = SubsetIterator(facet_domains, ind)
    geom = mesh.geometry()
    #maxr2 = -1.0
    maxr2 = 0
    for i, facet in enumerate(it):
        ent = facet.entities(0)
        for v in ent:
            p = geom.point(v)
            r2 = sum((p[j] - center[j])**2 for j in range(d))
            maxr2 = max(maxr2, r2)
    r = MPI.max(MPI.comm_world, sqrt(maxr2))
    return r


def compute_boundary_geometry_acrn(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)

    d = mesh.geometry().dim()
    x = SpatialCoordinate(mesh)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0, name="one")*dsi)
    #assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"
    if A == 0:
        return None

    # Compute barycenter by integrating x components over all facets
    c = [assemble(x[i]*dsi) / A for i in range(d)]

    # Compute average normal (assuming boundary is actually flat)
    n = FacetNormal(mesh)
    ni = np.array([assemble(n[i]*dsi) for i in range(d)])
    n_len = np.sqrt(sum([ni[i]**2 for i in range(d)])) # Should always be 1!?
    normal = ni/n_len

    # Compute radius by taking max radius of boundary points
    # (assuming boundary points are on exact geometry)
    # r = compute_radius(mesh, facet_domains, ind, c)
    # This old estimate is a few % lower because of boundary discretization errors
    r = np.sqrt(A / math.pi)

    return A, c, r, normal


def compute_area(mesh, ind, facet_domains):
    # Some convenient variables
    assert facet_domains is not None
    dsi = ds(ind, domain=mesh, subdomain_data=facet_domains)

    # Compute area of boundary tesselation by integrating 1.0 over all facets
    A = assemble(Constant(1.0, name="one")*dsi)
    assert A > 0.0, "Expecting positive area, probably mismatch between mesh and markers!"
    return A


def fourier_coefficients(x, y, T, N):
    '''From x-array and y-spline and period T, calculate N complex Fourier coefficients.'''
    omega = 2*np.pi/T
    ck = []
    ck.append(1/T*simps(y(x), x))
    for n in range(1,N):
        c = 1/T*simps(y(x)*np.exp(-1j*n*omega*x), x)

        # Clamp almost zero real and imag components to zero
        if 1:
            cr = c.real
            ci = c.imag
            if abs(cr) < 1e-14: cr = 0.0
            if abs(ci) < 1e-14: ci = 0.0
            c = cr + ci*1j

        ck.append(2*c)
    return ck


class ParabolicComponent(UserExpression):
    # Subclassing the expression class restricts the number of arguments, args
    # is therefore a dict of arguments.
    def __init__(self, radius, center, normal, normal_component, period, nu, element, Q=None,
                 V=None):
        # Spatial args
        self.radius = radius
        self.center = center
        self.normal = normal
        self.normal_component = normal_component

        # Temporal args
        self.period = period
        if Q is not None:
            assert V is None, "Cannot provide both Q and V!"
            self.Qn = Q
            self.N = len(self.Qn)
        elif V is not None:
            self.Vn = V
            self.N = len(self.Vn)
        else:
            raise ValueError("Invalid transient data type, missing argument 'Q' or 'V'.")
        
        self.ns = np.arange(1, self.N)
        self.omega = 2 * np.pi / self.period
        
        # Physical args
        self.nu = nu

        # Internal state
        self.t = None
        self.scale_value = 1.0

        # Precomputation
        #self._precompute_bessel_functions()
        self._all_r_dependent_coeffs = {}

        super().__init__(element=element)

    def _precompute_r_dependent_coeffs(self, y): #y ? 
        pir2 = np.pi * self.radius**2
        # Compute intermediate terms for womersley function
        #r_dependent_coeffs = np.zeros(self.N, dtype=np.complex)
        r_dependent_coeffs = np.zeros(self.N)
        if hasattr(self, 'Vn'):
            #r_dependent_coeffs[0] = (self.Vn[0]/2.0) * (1 - y**2)
            r_dependent_coeffs[0] = self.Vn[0] * (1 - y**2)
            for n in self.ns:
                r_dependent_coeffs[n] = self.Vn[n] 
        elif hasattr(self, 'Qn'):
            r_dependent_coeffs[0] = (2*self.Qn[0]/pir2) * (1 - y**2)
            for n in self.ns:
                r_dependent_coeffs[n] = self.Qn[n] 
        else:
            raise ValueError("Missing Vn or Qn!")
        return r_dependent_coeffs

    def _get_r_dependent_coeffs(self, y):
        "Look for cached coeffs."
        key = y
        r_dependent_coeffs = self._all_r_dependent_coeffs.get(key)
        if r_dependent_coeffs is None:
            # Cache miss! Compute coeffs for this coordinate the first time.
            r_dependent_coeffs = self._precompute_r_dependent_coeffs(y)
            self._all_r_dependent_coeffs[key] = r_dependent_coeffs
        return r_dependent_coeffs

    def set_t(self, t):
        self.t = float(t) % self.period
        self._expnt = np.exp((self.omega * self.t * 1j) * self.ns)

    def eval(self, value, x):
        # Compute or get cached complex coefficients that only depend on r
        y = np.sqrt(x_to_r2(x, self.center, self.normal)) / self.radius
        coeffs = self._get_r_dependent_coeffs(y)

        # Multiply complex coefficients for x with complex exponential functions in time
        par = (coeffs[0] + np.dot(coeffs[1:], self._expnt)).real

        # Scale by negative normal direction and scale_value
        value[0] = -self.normal_component * self.scale_value * par

def make_parabolic_bcs(t, Q, mesh, nu, area, center, radius, normal,
                       element, scale_to=None, coeffstype="Q",
                       N=1001, num_fourier_coefficients=20, **NS_namespace):
    """Generate a list of expressions for the components of a profile."""
    # Compute transient profile as interpolation of given coefficients
    period = max(t)
    transient_profile = UnivariateSpline(t, Q, s=0, k=1)

    # Compute fourier coefficients of transient profile
    timedisc = np.linspace(0, period, N)

    Cn = fourier_coefficients(timedisc, transient_profile, period, num_fourier_coefficients)

    # Create Expressions for each direction
    expressions = []
    for ncomp in normal:
        if coeffstype == "Q":
            Q = Cn
            V = None
        elif coeffstype == "V":
            V = Cn
            Q = None
        expressions.append(ParabolicComponent(radius, center, normal, ncomp, period, nu,
                                              element, Q=Q, V=V))

    return expressions

