try:
    from . import Artery, Atrium, Probe, Womersley, simulation_common
except ModuleNotFoundError:
    print("WARNING: Oasis is not installed, running CFD is not available")

try:
    from . import MovingAtrium, Probe, Womersley, simulation_common
except ModuleNotFoundError:
    print(
        "WARNING: OasisMove is not installed, running moving domain CFD is not available"
    )
