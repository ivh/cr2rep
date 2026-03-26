from typing import List

from edps import JobParameters, get_parameter, ClassifiedFitsFile, Job

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


# Filter out flat calibrations for standard stars, depending on the value of the workflow parameter "only_masterflat_for_standard".
# If this parameter is set to "true", only CAL_FLAT_MASTER will be kept in the input filter, while CAL_FLAT_TW and CAL_FLAT_EXTRACT_1D will be removed.
def change_input_filter(job: Job):
    only_masterflat = str(job.parameters.get_workflow_param("only_masterflat_for_standard", "None")).lower()
    if only_masterflat == 'true':
        new_filter = [x for x in job.input_filter if x not in ['CAL_FLAT_TW', 'CAL_FLAT_EXTRACT_1D']]
        job.input_filter = new_filter
