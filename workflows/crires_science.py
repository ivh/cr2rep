from edps import task, QC0, SCIENCE, CALCHECKER, subworkflow

from .crires_datasources import *
from .crires_task_functions import *

IDP = "idp"


@subworkflow("science", "")
def process_science(dark, flat_calibrations, wavelength_calibrations, detector_linearity):
    # This workflow collects all the tasks needed to process science observations

    # - Process nodding exposures
    science_nodding = (task("science_nodding")
                       .with_recipe("cr2res_obs_nodding")
                       .with_main_input(raw_science_nodding)
                       .with_alternatives(flat_calibrations)
                       .with_alternative_associated_inputs(wavelength_calibrations)
                       .with_associated_input(photo_flux, min_ret=0)
                       .with_associated_input(util_wave_tw, min_ret=0)
                       .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                       .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                       .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM],
                                              match_rules=match_dark_for_science, min_ret=0)
                       .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
                       .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                       .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_WAVE_TW, photo_flux_class,
                                          CAL_FLAT_EXTRACT_1D, util_wave_tw_class, util_trace_tw_class,
                                          linearity_coefficients)
                       .with_meta_targets([SCIENCE, QC0, IDP, CALCHECKER])
                       .build())

    # - Process staring observations
    science_staring = (task("science_staring")
                       .with_recipe("cr2res_obs_staring")
                       .with_main_input(raw_science_staring)
                       .with_alternatives(flat_calibrations)
                       .with_alternative_associated_inputs(wavelength_calibrations)
                       .with_associated_input(photo_flux, min_ret=0)
                       .with_associated_input(util_wave_tw, min_ret=0)
                       .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                       .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                       .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM],
                                              match_rules=match_dark_for_science)
                       .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
                       .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                       .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_WAVE_TW, CAL_FLAT_EXTRACT_1D,
                                          photo_flux_class, util_wave_tw_class, util_trace_tw_class,
                                          linearity_coefficients)
                       .with_meta_targets([SCIENCE, QC0, IDP, CALCHECKER])
                       .build())

    # - Process spectro-polarimetry observations
    science_polarimetry = (task("science_polarimetry")
                           .with_recipe("cr2res_obs_pol")
                           .with_main_input(raw_sci_polarimetry)
                           .with_alternatives(flat_calibrations)
                           .with_alternative_associated_inputs(wavelength_calibrations)
                           .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM],
                                                  match_rules=match_dark_for_polarimetry)
                           .with_associated_input(photo_flux, min_ret=0)
                           .with_associated_input(util_wave_tw, min_ret=0)
                           .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                           .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                           .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
                           .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                           .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_WAVE_TW, CAL_FLAT_EXTRACT_1D,
                                              photo_flux_class, util_wave_tw_class, util_trace_tw_class,
                                              linearity_coefficients, CAL_DARK_MASTER, CAL_DARK_BPM)
                           .with_meta_targets([SCIENCE, QC0, CALCHECKER])
                           .build())

    # Process astrometric observations
    science_astrometry = (task("science_astrometry")
                          .with_recipe("cr2res_obs_nodding")
                          .with_main_input(raw_sci_astrometry)
                          .with_alternatives(flat_calibrations)
                          .with_alternative_associated_inputs(wavelength_calibrations)
                          .with_associated_input(photo_flux, min_ret=0)
                          .with_associated_input(util_wave_tw, min_ret=0)
                          .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                          .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                          .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
                          .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM],
                                                 match_rules=match_dark_for_science)
                          .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                          .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_WAVE_TW, CAL_FLAT_EXTRACT_1D,
                                             photo_flux_class, util_wave_tw_class, util_trace_tw_class,
                                             linearity_coefficients)
                          .with_meta_targets([SCIENCE, QC0, CALCHECKER])
                          .build())

    # - Process two-dimensional observations
    science_2d = (task("science_2d")
                  .with_recipe("cr2res_obs_2d")
                  .with_main_input(raw_science_2d)
                  .with_associated_input(raw_science_2d_sky, min_ret=0, max_ret=1000)
                  .with_alternatives(flat_calibrations)
                  .with_alternative_associated_inputs(wavelength_calibrations)
                  .with_associated_input(photo_flux, min_ret=0)
                  .with_associated_input(util_wave_tw, min_ret=0)
                  .with_associated_input(util_trace_tw, min_ret=0, condition=is_short_wavelength)
                  .with_associated_input(util_trace_tw, min_ret=1, condition=is_long_wavelength)
                  .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
                  .with_dynamic_parameter("wavelength_range", get_wavelength_range)
                  .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_WAVE_TW, CAL_FLAT_EXTRACT_1D, photo_flux_class,
                                     util_wave_tw_class, util_trace_tw_class, linearity_coefficients)
                  .with_meta_targets([QC0, SCIENCE, CALCHECKER])
                  .build())

    return science_nodding, science_staring, science_polarimetry, science_astrometry, science_2d
