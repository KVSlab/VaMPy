"""Top-level package for VaMPy."""
from importlib.metadata import metadata

from .automatedPostProcessing import compute_flow_and_simulation_metrics
from .automatedPostProcessing import compute_flow_and_simulation_metrics
from .automatedPostProcessing import compute_velocity_and_pressure
from .automatedPostProcessing import postprocessing_common
from .automatedPostProcessing import visualize_probes
from .automatedPreProcessing import run_pre_processing

meta = metadata("vampy")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]

__all__ = [
    "visualize_probes",
    "compute_velocity_and_pressure",
    "compute_flow_and_simulation_metrics",
    "compute_flow_and_simulation_metrics",
    "postprocessing_common"
]
