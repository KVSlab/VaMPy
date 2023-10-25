try:
    from . import Artery
    from . import Atrium
    from . import simulation_common
    from . import Womersley
    from . import Probe
except ModuleNotFoundError:
    print("WARNING: Oasis is not installed, running CFD is not available")

try:
    from . import simulation_common
    from . import Womersley
    from . import Probe
    from . import MovingAtrium
except ModuleNotFoundError:
    print("WARNING: OasisMove is not installed, running moving domain CFD is not available")
