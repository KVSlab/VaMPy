"""Top-level package for VaMPy."""
from importlib.metadata import metadata

# Imports from post-processing
from .automatedPostProcessing import compute_flow_and_simulation_metrics
from .automatedPostProcessing import compute_hemodynamic_indices
from .automatedPostProcessing import compute_velocity_and_pressure
from .automatedPostProcessing import postprocessing_common
from .automatedPostProcessing import visualize_probes
# Imports from pre-processing
from .automatedPreProcessing import DisplayData
from .automatedPreProcessing import ImportData
from .automatedPreProcessing import NetworkBoundaryConditions
from .automatedPreProcessing import ToolRepairSTL
from .automatedPreProcessing import automatedPreProcessing
from .automatedPreProcessing import preprocessing_common
from .automatedPreProcessing import simulate
from .automatedPreProcessing import visualize
from .automatedPreProcessing import vmtkpointselector
# Imports from simulation scripts
from .simulation import Artery
from .simulation import Atrium
from .simulation import Probe
from .simulation import Womersley
from .simulation import simulation_common

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
    "automatedPreProcessing",
    "preprocessing_common",
    "DisplayData",
    "ImportData",
    "NetworkBoundaryConditions",
    "simulate",
    "ToolRepairSTL",
    "visualize",
    "vmtkpointselector",
    "Artery",
    "Atrium",
    "simulation_common",
    "Probe",
    "Womersley",
]
