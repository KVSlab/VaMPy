# Copyright (c) 2025 Simula Research Laboratory
# SPDX-License-Identifier: GPL-3.0-or-later
"""Top-level package for VaMPy."""
from importlib.metadata import metadata

# Imports from post-processing
from .automatedPostprocessing import compute_flow_and_simulation_metrics
from .automatedPostprocessing import compute_hemodynamic_indices
from .automatedPostprocessing import compute_velocity_and_pressure
from .automatedPostprocessing import postprocessing_common
from .automatedPostprocessing import visualize_probes
# Imports from pre-processing
try:
    from .automatedPreprocessing import DisplayData
    from .automatedPreprocessing import ImportData
    from .automatedPreprocessing import NetworkBoundaryConditions
    from .automatedPreprocessing import repair_tools
    from .automatedPreprocessing import automated_preprocessing
    from .automatedPreprocessing import preprocessing_common
    from .automatedPreprocessing import simulate
    from .automatedPreprocessing import visualize
    from .automatedPreprocessing import vmtk_pointselector
except ModuleNotFoundError:
    print("WARNING: morphMan is not installed, automated preprocessing is not available")
# Imports from simulation scripts
try:
    from .simulation import Artery
    from .simulation import Atrium
    from .simulation import Probe
    from .simulation import Womersley
    from .simulation import simulation_common
except ModuleNotFoundError:
    print("WARNING: Oasis is not installed, running CFD is not available")

try:
    from .simulation import Probe
    from .simulation import Womersley
    from .simulation import simulation_common
    from .simulation import MovingAtrium
except ModuleNotFoundError:
    print("WARNING: OasisMove is not installed, running moving domain simulations (MovingAtrium) is not available")


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
    "repair_tools",
    "visualize",
    "vmtk_pointselector",
    "Artery",
    "Atrium",
    "simulation_common",
    "Probe",
    "Womersley",
    'MovingAtrium'
]
