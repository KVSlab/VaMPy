"""Top-level package for VaMPy."""

from importlib.metadata import metadata

# Imports from post-processing
from .automatedPostprocessing import (
    compute_flow_and_simulation_metrics,
    compute_hemodynamic_indices,
    compute_velocity_and_pressure,
    postprocessing_common,
    visualize_probes,
)

# Imports from pre-processing
from .automatedPreprocessing import (
    DisplayData,
    ImportData,
    NetworkBoundaryConditions,
    ToolRepairSTL,
    automated_preprocessing,
    preprocessing_common,
    simulate,
    visualize,
    vmtk_pointselector,
)

# Imports from simulation scripts
from .simulation import Artery, Atrium, Probe, Womersley, simulation_common

meta = metadata("vampy")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]

__all__ = [
    "compute_flow_and_simulation_metrics",
    "compute_hemodynamic_indices",
    "compute_velocity_and_pressure",
    "postprocessing_common",
    "visualize_probes",
    "automated_preprocessing",
    "preprocessing_common",
    "DisplayData",
    "ImportData",
    "NetworkBoundaryConditions",
    "simulate",
    "ToolRepairSTL",
    "visualize",
    "vmtk_pointselector",
    "Artery",
    "Atrium",
    "simulation_common",
    "Probe",
    "Womersley",
]
