try :
    from . import Artery
    from . import Atrium
    from . import simulation_common
    from . import Womersley
    from . import Probe
except ModuleNotFoundError :
    print("oasis is not installed, running CFD is not available")