from __future__ import print_function
# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.

from dolfin import Expression, Mesh, MeshFunction, error, SubsetIterator, MPI, ds, assemble, Constant, sqrt,FacetNormal, as_vector, mpi_comm_world, SpatialCoordinate

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
            r2 = sum((p[j] - center[j])**2 for j in xrange(d))
            maxr2 = max(maxr2, r2)
    r = MPI.max(mpi_comm_world(), sqrt(maxr2))
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
    c = [assemble(x[i]*dsi) / A for i in xrange(d)]

    # Compute average normal (assuming boundary is actually flat)
    n = FacetNormal(mesh)
    ni = np.array([assemble(n[i]*dsi) for i in xrange(d)])
    n_len = np.sqrt(sum([ni[i]**2 for i in xrange(d)])) # Should always be 1!?
    normal = ni/n_len

    # Compute radius by taking max radius of boundary points
    # (assuming boundary points are on exact geometry)
    #r = compute_radius(mesh, facet_domains, ind, c)
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


class WomersleyComponent2(Expression):
    # Subclassing the expression class restricts the number of arguments, args is therefore a dict of arguments.
    def __init__(self, args): # TODO: Document args properly
        Expression.__init__(self)

        # Spatial args
        self.radius = args["radius"]
        self.center = args["center"]
        self.normal = args["normal"]
        self.normal_component = args["normal_component"]

        # Temporal args
        self.period = args["period"]
        self.Qn = args["Qn"]

        # Physical args
        self.nu = args["nu"]

        # Internal state
        self._t = None
        self._scale_value = 1.0

        # Precomputation
        self._values = np.zeros(10000)
        self._ys = np.linspace(0.0, 1.0, len(self._values))
        self._precompute_bessel_functions()
        self._precompute_y_coeffs()

    def _precompute_bessel_functions(self):
        '''Calculate the Bessel functions of the Womersley profile'''
        self.omega = 2 * np.pi / self.period
        self.N = len(self.Qn)
        self.ns = np.arange(1, self.N)

        # Allocate for 0...N-1
        alpha = np.zeros(self.N, dtype=np.complex)
        self.beta = np.zeros(self.N, dtype=np.complex)
        self.jn0_betas = np.zeros(self.N, dtype=np.complex)
        self.jn1_betas = np.zeros(self.N, dtype=np.complex)

        # Compute vectorized for 1...N-1 (keeping element 0 in arrays to make indexing work out later)
        alpha[1:] = self.radius * np.sqrt(self.ns * (self.omega / self.nu))
        self.beta[1:] = alpha[1:] * np.sqrt(1j**3)
        self.jn0_betas[1:] = jn(0, self.beta[1:])
        self.jn1_betas[1:] = jn(1, self.beta[1:])

    def _precompute_y_coeffs(self):
        "Compute intermediate terms for womersley function."
        n = len(self._values)
        self._y_coeffs = np.zeros((n,self.N), dtype=np.complex)
        pir2 = np.pi * self.radius**2
        for i, y in enumerate(self._ys):
            self._y_coeffs[i,0] = (2*self.Qn[0]/pir2) * (1 - y**2)
            for n in self.ns:
                tmp1 = self.Qn[n] / pir2
                tmp2 = 1.0 - jn(0, self.beta[n]*y) / self.jn0_betas[n]
                tmp3 = 1.0 - 2.0*self.jn1_betas[n] / (self.beta[n]*self.jn0_betas[n])
                self._y_coeffs[i,n] = tmp1 * (tmp2 / tmp3)

    def set_t(self, t):
        # Compute time dependent coeffs once
        self._t = float(t) % self.period
        self._expnt = np.exp((self.omega * self._t * 1j) * self.ns)

        # Compute values for this time for each y
        n = len(self._values)
        for i in xrange(n):
            # Multiply complex coefficients for x with complex exponential functions in time
            wom = (self._y_coeffs[i,0] + np.dot(self._y_coeffs[i,1:], self._expnt)).real

            # Scale by negative normal direction and scale_value
            self._values[i] = -self.normal_component * self._scale_value * wom

    def eval(self, value, x):
        y = np.sqrt(x_to_r2(x, self.center, self.normal)) / self.radius
        nm = len(self._values) - 1
        yi = max(0, min(nm, int(round(y*nm))))
        value[0] = self._values[yi]


class WomersleyComponent1(Expression):
    # Subclassing the expression class restricts the number of arguments, args
    # is therefore a dict of arguments.
    def init(self, args):
        # Spatial args
        self.radius = args["radius"]
        self.center = args["center"]
        self.normal = args["normal"]
        self.normal_component = args["normal_component"]

        # Temporal args
        self.period = args["period"]
        if "Q" in args:
            assert "V" not in args, "Cannot provide both Q and V!"
            self.Qn = args["Q"]
            self.N = len(self.Qn)
        elif "V" in args:
            self.Vn = args["V"]
            self.N = len(self.Vn)
        else:
            error("Invalid transient data type, missing argument 'Q' or 'V'.")

        # Physical args
        self.nu = args["nu"]

        # Internal state
        self.t = None
        self.scale_value = 1.0

        # Precomputation
        self._precompute_bessel_functions()
        self._all_r_dependent_coeffs = {}

    def _precompute_bessel_functions(self):
        '''Calculate the Bessel functions of the Womersley profile'''
        self.omega = 2 * np.pi / self.period
        self.ns = np.arange(1, self.N)

        # Allocate for 0...N-1
        alpha = np.zeros(self.N, dtype=np.complex)
        self.beta = np.zeros(self.N, dtype=np.complex)
        self.jn0_betas = np.zeros(self.N, dtype=np.complex)
        self.jn1_betas = np.zeros(self.N, dtype=np.complex)

        # Compute vectorized for 1...N-1 (keeping element 0 in arrays to make indexing work out later)
        alpha[1:] = self.radius * np.sqrt(self.ns * (self.omega / self.nu))
        self.beta[1:] = alpha[1:] * np.sqrt(1j**3)
        self.jn0_betas[1:] = jn(0, self.beta[1:])
        self.jn1_betas[1:] = jn(1, self.beta[1:])

    def _precompute_r_dependent_coeffs(self, y):
        pir2 = np.pi * self.radius**2
        # Compute intermediate terms for womersley function
        r_dependent_coeffs = np.zeros(self.N, dtype=np.complex)
        if hasattr(self, 'Vn'):
            #r_dependent_coeffs[0] = (self.Vn[0]/2.0) * (1 - y**2)
            r_dependent_coeffs[0] = self.Vn[0] * (1 - y**2)
            for n in self.ns:
                r_dependent_coeffs[n] = self.Vn[n] * (self.jn0_betas[n] - jn(0,
                                        self.beta[n]*y)) / (self.jn0_betas[n] - 1.0)
        elif hasattr(self, 'Qn'):
            r_dependent_coeffs[0] = (2*self.Qn[0]/pir2) * (1 - y**2)
            for n in self.ns:
                bn = self.beta[n]
                j0bn = self.jn0_betas[n]
                j1bn = self.jn1_betas[n]
                r_dependent_coeffs[n] = (self.Qn[n] / pir2) * (j0bn - jn(0,
                                            bn*y)) / (j0bn - (2.0/bn)*j1bn)
        else:
            error("Missing Vn or Qn!")
        return r_dependent_coeffs

    def _get_r_dependent_coeffs(self, y):
        "Look for cached womersley coeffs."
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
        wom = (coeffs[0] + np.dot(coeffs[1:], self._expnt)).real

        # Scale by negative normal direction and scale_value
        value[0] = -self.normal_component * self.scale_value * wom


def make_womersley_bcs(t, Q, mesh, nu, area, center, radius, normal,
                       v_degree, scale_to=None, coeffstype="Q",
                       N=1001, num_fourier_coefficients=20, **NS_namespace):
    """Generate a list of expressions for the components of a Womersley profile."""
    # Compute transient profile as interpolation of given coefficients
    period = max(t)
    transient_profile = UnivariateSpline(t, Q, s=0, k=1)

    # Compute fourier coefficients of transient profile
    timedisc = np.linspace(0, period, N)

    #from matplotlib.pyplot import plot, savefig
    #plot(timedisc, transient_profile(timedisc))
    #savefig("series.png")
    #sys.exit(0)
    Cn = fourier_coefficients(timedisc, transient_profile, period, num_fourier_coefficients)

    # Create Expressions for each direction
    expressions = []
    for ncomp in normal:
        args = {
            "radius": radius,
            "center": center,
            "normal": normal,
            "normal_component": ncomp,
            "period": period,
            "nu": nu,
            }
        args[coeffstype] = Cn
        expressions.append(WomersleyComponent1(degree=v_degree))
        expressions[-1].init(args)

    return expressions
