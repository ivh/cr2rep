from edps import QC1_CALIB, CALCHECKER, ReportInput
from edps import task, alternative_associated_inputs

from .crires_datasources import *
from .crires_science import process_science
from .crires_standards import crires_standards
from .crires_task_functions import *
from .crires_wavelength import wavelength_calibration

__title__ = "CRIRES workflow"
IDP = "idp"
# --- Processing tasks -------------------------------------------------------------------

# - Task to monitor the detector linearity.
detector_linearity = (task("detector_linearity")
                      .with_recipe("cr2res_cal_detlin")
                      .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
                      .with_report("crires_detmon", ReportInput.RECIPE_INPUTS_OUTPUTS)
                      .with_main_input(raw_detlin)
                      .with_meta_targets([QC1_CALIB])
                      .build())

# - Task to process DARK frames
dark = (task("dark")
        .with_recipe("cr2res_cal_dark")
        .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
        .with_report("crires_master_dark", ReportInput.RECIPE_INPUTS_OUTPUTS)
        .with_main_input(raw_dark)
        .with_meta_targets([QC1_CALIB])
        .build())

# ----- Flat calibrations -----
# - Task to process Flat Field frames
flat = (task("flat")
        .with_recipe("cr2res_cal_flat")
        .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
        .with_report("crires_echelle_flatfield", ReportInput.RECIPE_INPUTS_OUTPUTS)
        .with_main_input(raw_flat)
        .with_associated_input(dark, [CAL_DARK_MASTER], match_rules=match_dark_for_calibs)
        .with_associated_input(util_wave_tw, min_ret=0)
        .with_associated_input(util_trace_tw, min_ret=0)
        .with_associated_input(detector_linearity, [linearity_coefficients], min_ret=0)
        .with_meta_targets([QC1_CALIB, CALCHECKER])
        .build())

# CAL_FLAT_EXTRACT_1D is not mandatory. This associates the set of calibrations
# (only mandatory or all calibrations) that is closer in time.
flat_calibrations = (alternative_associated_inputs(sort_keys=[kwd.mjd_obs])
                     .with_associated_input(flat, [CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_FLAT_EXTRACT_1D])
                     .with_associated_input(flat, [CAL_FLAT_MASTER, CAL_FLAT_TW]))

# ------------

# -- Subworkflow to run tasks for wavelength calibration
wavelength_calibrations, wavelength_fpet, wavelength_une, wavelength_gas_cell, wavelength_sky = wavelength_calibration(
    dark, flat_calibrations)

# - Task to process standard stars
standard_star, spectropolarimetric_standard = crires_standards(dark, flat_calibrations, wavelength_calibrations,
                                                               detector_linearity)

# --- Subworkflow to reduce scientific exposures --------------------------------------------------

science_nodding, science_staring, science_polarimetry, science_astrometry, science_2d = process_science(dark,
                                                                                                        flat_calibrations,
                                                                                                        wavelength_calibrations,
                                                                                                        detector_linearity)
