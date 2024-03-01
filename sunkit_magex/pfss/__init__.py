# Import this to register map sources
import sunkit_magex.pfss.map as _
from sunkit_magex.pfss import coords
from sunkit_magex.pfss import fieldline
from sunkit_magex.pfss import sample_data
from sunkit_magex.pfss import tracing
from sunkit_magex.pfss import utils
from sunkit_magex.pfss.input import Input
from sunkit_magex.pfss.output import Output
from sunkit_magex.pfss.pfss import pfss

__all__ = ['coords', 'fieldline', 'sample_data', 'tracing', 'utils', 'Input', 'Output', 'pfss']

try:
    from sunkit_magex.pfss import analytic
    __all__.append('analytic')
except ImportError:
    pass
