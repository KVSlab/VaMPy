from argparse import ArgumentParser
from time import time

from dolfin import parameters, MPI, assemble, interpolate, Measure, FacetNormal, Identity, VectorFunctionSpace, BoundaryMesh, Function, FacetArea, TestFunction, FunctionSpace, grad, inner, sqrt

try:
    from dolfin import *

    try:
        parameters["allow_extrapolation"] = True
    except NameError:
        pass
except ImportError:
    pass


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--case',
                        type=str,
                        default="simulation/results_folder/model/data/1/Solutions",
                        help="Path to simulation results (PostProc folder for executing compute_hemodynamic_indices.py and h5_files folder for executing compute_differences_h5_files.py): mesh.h5, u.h5, p.h5, nu.h5",
                        metavar="PATH")    
    parser.add_argument('--nu',
                        type=float,
                        default=3.3018e-3,
                        help="Kinematic viscosity used in simulation, measured in [mm^2/ms]")
    parser.add_argument('--rheology_model',
                        type=str,
                        default="Newtonian",
                        help="Set the rheology model used in the simulation and computation of WSS")
    parser.add_argument('--rho',
                        type=float,
                        default=1060,
                        help="Fluid density used in simulation, measured in [kg/m^3]")
    parser.add_argument('--dt',
                        type=float,
                        default=0.1,
                        help="Time step of simulation, measured in [ms]")
    parser.add_argument('--velocity-degree',
                        type=int,
                        default=1,
                        help="Degree of velocity element")
    parser.add_argument('--probe-frequency',
                        type=int,
                        default=100,
                        help="Frequency of saving probes to file")
    parser.add_argument('--T',
                        type=float,
                        default=1000,
                        help="Duration of one cardiac cycle. Measured in [ms].")
    parser.add_argument('--save-frequency',
                        type=int,
                        default=10,
                        help="Frequency of saving velocity to file.")
    parser.add_argument('--start-cycle',
                        type=int,
                        default=1,
                        help="Start post-processing from this cardiac cycle.")
    parser.add_argument('--sample-step',
                        type=int,
                        default=1,
                        help="Step size that determines how many times data is sampled.")
    parser.add_argument('--average-over-cycles',
                        default=True,
                        action='store_true',
                        help="Computes average over all cycles if True.")

    args = parser.parse_args()

    return args.case, args.nu, args.rheology_model, args.rho, args.dt, args.velocity_degree, args.probe_frequency, args.T, args.save_frequency, args.start_cycle, args.sample_step, args.average_over_cycles



def epsilon(u):
    """
    Computes the strain-rate tensor
    Args:
        u (Function): Velocity field

    Returns:
        epsilon (Function): Strain rate tensor of u
    """

    return 0.5 * (grad(u) + grad(u).T)


class STRESS:
    def __init__(self, u, p, nu, rheology_model, mesh):
        boundary_ds = Measure("ds", domain=mesh)
        boundary_mesh = BoundaryMesh(mesh, 'exterior')
        self.bmV = VectorFunctionSpace(boundary_mesh, 'CG', 1)

        print("STRESS computation based on the rheology_model={}!".format(rheology_model))

        # Compute stress tensor
        sigma = (2 * nu * epsilon(u)) - (p * Identity(len(u)))
        print("sigma=", sigma)
        print("pass")
        # Compute stress on surface
        n = FacetNormal(mesh)
        F = -(sigma * n)

        print("pass 2")

        # Compute normal and tangential components
        Fn = inner(F, n)  # scalar-valued
        Ft = F - (Fn * n)  # vector-valued

        print("pass 3")

        # Integrate against piecewise constants on the boundary
        scalar = FunctionSpace(mesh, 'DG', 0)
        vector = VectorFunctionSpace(mesh, 'CG', 1)
        scaling = FacetArea(mesh)  # Normalise the computed stress relative to the size of the element

        v = TestFunction(scalar)
        w = TestFunction(vector)

        print("pass 4")

        # Create functions
        self.Fn = Function(scalar)
        self.Ftv = Function(vector)
        self.Ft = Function(scalar)


        print("pass 5")

        self.Ln = 1 / scaling * v * Fn * boundary_ds
        self.Ltv = 1 / (2 * scaling) * inner(w, Ft) * boundary_ds
        self.Lt = 1 / scaling * inner(v, self.norm_l2(self.Ftv)) * boundary_ds

        print("pass 6")

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

