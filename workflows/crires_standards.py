from edps import task, subworkflow, qc1calib, calchecker, IDP, ReportInput

from .crires_datasources import *
from .crires_task_functions import is_short_wavelength, is_long_wavelength, get_wavelength_range


@subworkflow("standard_stars", "")
def crires_standards(dark, flat_calibrations, wavelength_calibrations, detector_linearity):
    standard_star = (task("standard_star")
                     .with_recipe("cr2res_obs_nodding")
                     .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
                     .with_report("crires_std_star", ReportInput.RECIPE_INPUTS_OUTPUTS)
                     .with_main_input(raw_standard)
                     .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                     .with_alternatives(flat_calibrations)
                     .with_alternative_associated_inputs(wavelength_calibrations)
                     .with_associated_input(photo_flux, min_ret=0)
                     .with_associated_input(util_wave_tw, min_ret=0)
                     .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                     .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                     .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_WAVE_TW, CAL_FLAT_EXTRACT_1D,
                                        photo_flux_class, util_wave_tw_class, util_trace_tw_class)
                     .with_meta_targets([qc1calib, calchecker, IDP])
                     .build())

    spectropolarimetric_standard = (task("standard_star_pol")
                                    .with_recipe("cr2res_obs_pol")
                                    .with_main_input(raw_polarimetry_standard)
                                    .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                                    .with_alternatives(flat_calibrations)
                                    .with_alternative_associated_inputs(wavelength_calibrations)
                                    .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM],
                                                           match_rules=match_dark_for_polarimetry, min_ret=1)
                                    .with_associated_input(util_wave_tw, min_ret=0)
                                    .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                                    .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                                    .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
                                    .with_meta_targets([qc1calib, calchecker])
                                    .build())

    return standard_star, spectropolarimetric_standard
