import argparse
from time import time

from dolfin import parameters, MPI, assemble, interpolate, Measure, FacetNormal, Identity, VectorFunctionSpace, \
    BoundaryMesh, Function, FacetArea, TestFunction, FunctionSpace, grad, inner, sqrt

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


def epsilon(u):
    """
    Computes the strain-rate tensor
    Args:
        u (Function): Velocity field

    Returns:
        epsilon (Function): Strain rate tensor of u
    """

    return 0.5 * (grad(u) + grad(u).T)


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
    eps = epsilon(u)
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
    eps = epsilon(u)
    f = 2 * nu * inner(eps, eps)
    x = assemble(inner(f, v) / h * dx)
    dissipation.vector().set_local(x.get_local())
    dissipation.vector().apply("insert")


class STRESS:
    """Computes the stress on a given mesh based on provided velocity and pressure fields."""

    def __init__(self, u, p, nu, mesh):
        """Initializes the StressComputer.

        Args:
            u (Function): The velocity field.
            p (Function): The pressure field.
            nu (float): The kinematic viscosity.
            mesh (Mesh): The mesh on which to compute stress.
        """
        boundary_ds = Measure("ds", domain=mesh)
        boundary_mesh = BoundaryMesh(mesh, 'exterior')
        self.bmV = VectorFunctionSpace(boundary_mesh, 'CG', 1)

        # Compute stress tensor
        sigma = (2 * nu * epsilon(u)) - (p * Identity(len(u)))

        # Compute stress on surface
        n = FacetNormal(mesh)
        F = -(sigma * n)

        # Compute normal and tangential components
        Fn = inner(F, n)  # scalar-valued
        Ft = F - (Fn * n)  # vector-valued

        # Integrate against piecewise constants on the boundary
        scalar = FunctionSpace(mesh, 'DG', 0)
        vector = VectorFunctionSpace(mesh, 'CG', 1)
        scaling = FacetArea(mesh)  # Normalise the computed stress relative to the size of the element

        v = TestFunction(scalar)
        w = TestFunction(vector)

        # Create functions
        self.Fn = Function(scalar)
        self.Ftv = Function(vector)
        self.Ft = Function(scalar)

        self.Ln = 1 / scaling * v * Fn * boundary_ds
        self.Ltv = 1 / (2 * scaling) * inner(w, Ft) * boundary_ds
        self.Lt = 1 / scaling * inner(v, self.norm_l2(self.Ftv)) * boundary_ds

    def __call__(self):
        """
        Compute stress for given velocity field u and pressure field p

        Returns:
            Ftv_mb (Function): Shear stress
        """

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


def get_dataset_names(data_file, num_files=100000, step=1, start=1, print_info=True,
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
