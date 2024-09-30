class CriresReportMixin(object):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._version = "4.0.0"

from . import crires_detmon
from . import crires_master_dark
from . import crires_wavelength
from . import crires_echelle_flatfield
from . import crires_rawdisp
from . import crires_std_star

__all__ = [
    crires_detmon,
    crires_master_dark,
    crires_wavelength,
    crires_echelle_flatfield,
    crires_rawdisp,
    crires_std_star
]
