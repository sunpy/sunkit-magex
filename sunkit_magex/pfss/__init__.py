# Import sunkit_magex.pfss sub-modules to have them available through sunkit_magex.pfss.{name}
try:
    import sunkit_magex.pfss.analytic
except ModuleNotFoundError:
    # If sympy isn't installed
    pass
import sunkit_magex.pfss.coords
import sunkit_magex.pfss.fieldline
# Import this to register map sources
import sunkit_magex.pfss.map
import sunkit_magex.pfss.sample_data
import sunkit_magex.pfss.tracing
import sunkit_magex.pfss.utils

from .input import Input
from .output import Output
from .pfss import pfss

__all__ = ['Input', 'Output', 'pfss']
