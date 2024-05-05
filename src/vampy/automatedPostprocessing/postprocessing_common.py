import argparse
from time import time

import numpy as np
from dolfin import sym, parameters, MPI, assemble, interpolate, Measure, FacetNormal, VectorFunctionSpace, \
    BoundaryMesh, Function, TestFunction, FunctionSpace, grad, inner, sqrt, TrialFunction, LUSolver, \
    FunctionAssigner, Mesh, ds

try:
    parameters["allow_extrapolation"] = True
except NameError:
    pass


def read_command_line():
    """Read arguments from commandline"""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Automated post-processing for vascular modeling.")

    parser.add_argument('-c', '--case',
                        type=str,
                        default="simulation/results_folder/model/data/1/Solutions",
                        help="Path to simulation results.",
                        metavar="PATH")

    parser.add_argument('-nu', '--nu',
                        type=float,
                        default=3.3018e-3,
                        help="Kinematic viscosity used in simulation. Measured in [mm^2/ms].")

    parser.add_argument('-r', '--rho',
                        type=float,
                        default=1060,
                        help="Fluid density used in simulation. Measured in [kg/m^3].")

    parser.add_argument('-T', '--T',
                        type=float,
                        default=951,
                        help="Duration of one cardiac cycle. Measured in [ms].")

    parser.add_argument('-dt', '--dt',
                        type=float,
                        default=0.0951,
                        help="Time step of simulation. Measured in [ms].")

    parser.add_argument('-vd', '--velocity-degree',
                        type=int,
                        default=1,
                        help="Degree of velocity element.")

    parser.add_argument('-pd', '--pressure-degree',
                        type=int,
                        default=1,
                        help="Degree of pressure element.")

    parser.add_argument('-sf', '--save-frequency',
                        type=int,
                        default=5,
                        help="Frequency of saving velocity to file.")

    parser.add_argument('-pf', '--probe-frequency',
                        type=int,
                        default=100,
                        help="Frequency of saving probes to file.")

    parser.add_argument('-ta', '--times-to-average',
                        type=float,
                        default=[],
                        nargs="+",
                        help="Time(s) during cardiac cycle to average, in the interval [0,T). Measured in [ms].")

    parser.add_argument('-sc', '--start-cycle',
                        type=int,
                        default=2,
                        help="Start post-processing from this cardiac cycle.")

    parser.add_argument('-ss', '--sample-step',
                        type=int,
                        default=1,
                        help="Step size that determines how many times data is sampled.")

    parser.add_argument('-ac', '--average-over-cycles',
                        default=False,
                        action='store_true',
                        help="Computes average over all cycles if True.")

    args = parser.parse_args()

    return args.case, args.nu, args.rho, args.dt, args.velocity_degree, args.pressure_degree, args.probe_frequency, \
           args.T, args.save_frequency, args.times_to_average, args.start_cycle, args.sample_step, \
           args.average_over_cycles


def rate_of_strain(strain, u, v, mesh, h):
    """
    Computes rate of strain

    Args:
        strain (Function): Function to save rate of strain to
        u (Function): Function for velocity field
        v (Function): Test function for velocity
        mesh: Mesh to compute strain rate on
        h (float): Cell diameter of mesh
    """
    dx = Measure("dx", domain=mesh)
    eps = sym(grad(u))
    f = sqrt(inner(eps, eps))
    x = assemble(inner(f, v) / h * dx)
    strain.vector().set_local(x.get_local())
    strain.vector().apply("insert")


def rate_of_dissipation(dissipation, u, v, mesh, h, nu):
    """
    Computes rate of dissipation

    Args:
        dissipation (Function): Function to save rate of dissipation to
        u (Function): Function for velocity field
        v (Function): Test function for velocity
        mesh: Mesh to compute dissipation rate on
        h (float): Cell diameter of mesh
        nu (float): Viscosity
    """
    dx = Measure("dx", domain=mesh)
    eps = sym(grad(u))
    f = 2 * nu * inner(eps, eps)
    x = assemble(inner(f, v) / h * dx)
    dissipation.vector().set_local(x.get_local())
    dissipation.vector().apply("insert")


class InterpolateDG:
    """
    interpolate DG function from the domain to the boundary. FEniCS built-in function interpolate does not work
    with DG function spaces. This class is a workaround for this issue. Basically, for each facet, we find the
    mapping between the dofs on the boundary and the dofs on the domain. Then, we copy the values of the dofs on the
    domain to the dofs on the boundary. This is done for each subspaces of the DG vector function space.
    """

    def __init__(self, V: VectorFunctionSpace, V_sub: VectorFunctionSpace, mesh: Mesh, boundary_mesh: Mesh) -> None:
        """
        Initialize the interpolator

        Args:
            V (VectorFunctionSpace): function space on the domain
            V_sub (VectorFunctionSpace): function space on the boundary
            mesh (Mesh): whole mesh
            boundary_mesh (Mesh): boundary mesh of the whole mesh
        """
        assert V.ufl_element().family() == "Discontinuous Lagrange", "V must be a DG space"
        assert V_sub.ufl_element().family() == "Discontinuous Lagrange", "V_sub must be a DG space"
        self.V = V
        self.v_sub = Function(V_sub)
        self.Ws = [V_sub.sub(i).collapse() for i in range(V_sub.num_sub_spaces())]
        self.ws = [Function(Wi) for Wi in self.Ws]
        self.w_sub_copy = [w_sub.vector().get_local() for w_sub in self.ws]
        self.sub_dofmaps = [W_sub.dofmap() for W_sub in self.Ws]
        self.sub_coords = [Wi.tabulate_dof_coordinates() for Wi in self.Ws]
        self.mesh = mesh
        self.sub_map = boundary_mesh.entity_map(self.mesh.topology().dim() - 1).array()
        self.mesh.init(self.mesh.topology().dim() - 1, self.mesh.topology().dim())
        self.f_to_c = self.mesh.topology()(self.mesh.topology().dim() - 1, self.mesh.topology().dim())
        self.dof_coords = V.tabulate_dof_coordinates()
        self.fa = FunctionAssigner(V_sub, self.Ws)

    def __call__(self, u_vec: np.ndarray) -> Function:
        """interpolate DG function from the domain to the boundary"""

        for k, (coords_k, vec, sub_dofmap) in enumerate(zip(self.sub_coords, self.w_sub_copy, self.sub_dofmaps)):
            for i, facet in enumerate(self.sub_map):
                cells = self.f_to_c(facet)
                # Get closure dofs on parent facet
                sub_dofs = sub_dofmap.cell_dofs(i)
                closure_dofs = self.V.sub(k).dofmap().entity_closure_dofs(
                    self.mesh, self.mesh.topology().dim(), [cells[0]])
                copy_dofs = np.empty(len(sub_dofs), dtype=np.int32)

                for dof in closure_dofs:
                    for j, sub_coord in enumerate(coords_k[sub_dofs]):
                        if np.allclose(self.dof_coords[dof], sub_coord):
                            copy_dofs[j] = dof
                            break
                sub_dofs = sub_dofmap.cell_dofs(i)
                vec[sub_dofs] = u_vec[copy_dofs]

            self.ws[k].vector().set_local(vec)

        self.fa.assign(self.v_sub, self.ws)

        return self.v_sub


class SurfaceProjector:
    """
    Project a function contains surface integral onto a function space V
    """

    def __init__(self, V: FunctionSpace) -> None:
        """
        Initialize the surface projector

        Args:
            V (FunctionSpace): function space to project onto
        """
        u = TrialFunction(V)
        v = TestFunction(V)
        a_proj = inner(u, v) * ds
        # keep_diagonal=True & ident_zeros() are necessary for the matrix to be invertible
        self.A = assemble(a_proj, keep_diagonal=True)
        self.A.ident_zeros()
        self.u_ = Function(V)
        self.solver = LUSolver(self.A)

    def __call__(self, f: Function) -> Function:
        v = TestFunction(self.u_.function_space())
        self.b_proj = inner(f, v) * ds
        self.b = assemble(self.b_proj)
        self.solver.solve(self.u_.vector(), self.b)
        return self.u_


class STRESS:
    """
    Computes the stress on a given mesh based on provided velocity and pressure fields.
    Note that the pressure term is unused since it disappears in the process of extracting the tangential component.

    FIXME: Currently works for P1P1, but not for higher order finite elements (e.g. P2P1)
    """

    def __init__(self, u: Function, V_dg: VectorFunctionSpace, V_sub: VectorFunctionSpace, nu: float, mesh: Mesh,
                 boundary_mesh: Mesh) -> None:
        """Initializes the StressComputer.

        Args:
            u (Function): The velocity field.
            nu (float): The kinematic viscosity.
            mesh (Mesh): The mesh on which to compute stress.
        """

        self.projector = SurfaceProjector(V_dg)
        self.interpolator = InterpolateDG(V_dg, V_sub, mesh, boundary_mesh)

        # Compute stress tensor
        sigma = 2 * nu * sym(grad(u))

        # Compute stress on surface
        n = FacetNormal(mesh)
        F = -(sigma * n)

        # Compute normal and tangential components
        Fn = inner(F, n)  # scalar-valued
        self.Ft = F - (Fn * n)  # vector-valued

        # Integrate against piecewise constants on the boundary
        # scalar = FunctionSpace(mesh, 'DG', 0)
        # vector = VectorFunctionSpace(mesh, 'CG', 1)
        # scaling = FacetArea(mesh)  # Normalise the computed stress relative to the size of the element
        #
        # v = TestFunction(scalar)
        # w = TestFunction(vector)
        #
        # # Create functions
        # self.Fn = Function(scalar)
        # self.Ftv = Function(vector)
        # self.Ft = Function(scalar)
        #
        # self.Ln = 1 / scaling * v * Fn * boundary_ds
        # self.Ltv = 1 / (2 * scaling) * inner(w, Ft) * boundary_ds
        # self.Lt = 1 / scaling * inner(v, self.norm_l2(self.Ftv)) * boundary_ds

    def __call__(self):
        """
        Compute stress for given velocity field u and pressure field p

        Returns:
            Ftv_mb (Function): Shear stress
        """

        self.Ftv = self.projector(self.Ft)
        self.Ftv_bd = self.interpolator(self.Ftv.vector().get_local())

        return self.Ftv_bd

        # Assemble vectors
        assemble(self.Ltv, tensor=self.Ftv.vector())
        self.Ftv_bm = interpolate(self.Ftv, self.bmV)

        return self.Ftv_bm

    def norm_l2(self, u):
        """
        Compute norm of vector u in expression form
        Args:
            u (Function): Function to compute norm of

        Returns:
            norm (Power): Norm as expression
        """
        return pow(inner(u, u), 0.5)


def get_dataset_names(data_file, num_files=100000, step=1, start=0, print_info=True,
                      vector_filename="/velocity/vector_%d"):
    """
    Read velocity fields datasets and extract names of files

    Args:
        data_file (HDF5File): File object of velocity
        num_files (int): Number of files
        step (int): Step between each data dump
        start (int): Step to start on
        print_info (bool): Prints info about data if true
        vector_filename (str): Name of velocity files

    Returns:
        names (list): List of data file names
    """
    check = True
    # Find start file
    t0 = time()
    while check:
        if data_file.has_dataset(vector_filename % start):
            check = False
            start -= step

        start += step

    # Get names
    names = []
    for i in range(num_files):

        index = start + i * step
        if data_file.has_dataset(vector_filename % index):
            names.append(vector_filename % index)

    t1 = time()

    # Print info
    if MPI.rank(MPI.comm_world) == 0 and print_info:
        print()
        print("=" * 6 + " Timesteps to average over " + "=" * 6)
        print("Length on data set names:", len(names))
        print("Start index:", start)
        print("Wanted num files:", num_files)
        print("Step between files:", step)
        print("Time used:", t1 - t0, "s")
        print()

    return names
