from edps import JobParameters, List, get_parameter, ClassifiedFitsFile

from . import crires_keywords as kwd


# This function set the dynamic parameter wavelength_range with the appropriate value:
# wavelength_range="short_wavelength" if observational setup is Y,J,H, or K
# wavelength_range="long wavelength" if setup is L or M
def get_wavelength_range(files: List[ClassifiedFitsFile]):
    value = files[0].get_keyword_value(kwd.ins_wlen_id, 'x')
    if value and value[0] in ["Y", "J", "H", "K"]:
        return "short_wavelength"
    elif value and value[0] in ["L", "M"]:
        return "long_wavelength"
    else:
        return "undefined"


def is_short_wavelength(params: JobParameters) -> bool:
    return get_parameter(params, "wavelength_range") == "short_wavelength"


def is_long_wavelength(params: JobParameters) -> bool:
    return get_parameter(params, "wavelength_range") == "long_wavelength"
