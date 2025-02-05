from edps import task, QC1_CALIB, CALCHECKER, subworkflow, alternative_associated_inputs, ReportInput

from .crires_datasources import *
from .crires_task_functions import *


@subworkflow("wavelength calibration", "")
def wavelength_calibration(dark, flat_calibrations):
    # This subworflow collects all the tasks associated to the processing of lamp calibrations.
    # Only the products of tasks "wavelength_fabry_perot_etalon" or "wavelength_uranium_neon" are used for subsequent steps in data reduction.
    # Gas cells wavelength calibrations obtained with gas-cells are also associated to the relevant science tasks, but
    # they are not processed by the pipeline.
    # The association rules and preferential orders on which calibration to associate to science tasks are specified
    # in the alternative association object "wavelength_calibration".
    # Wavelength calibrations obtained from sky frames are used only for instrument monitoring purposes.

    # This task reduces Fabry-Perot Etalon wavelength calibrations. It als requires a Uranium Neon reference from
    # the same template (raw_wave_une_ref)
    wavelength_fpet = (task("wave_fabry_perot_etalon")
                       .with_recipe("cr2res_cal_wave")
                       .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
                       .with_report("crires_wavelength", ReportInput.RECIPE_INPUTS_OUTPUTS)
                       .with_main_input(raw_wave_fpet)
                       .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM],
                                              match_rules=match_dark_for_calibs)
                       .with_alternatives(flat_calibrations)
                       .with_associated_input(emission_lines)
                       .with_associated_input(util_wave_tw, min_ret=0)
                       .with_associated_input(util_trace_tw, min_ret=0)
                       .with_associated_input(raw_wave_une)
                       .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_DARK_MASTER, CAL_DARK_BPM,
                                          emission_lines_class, util_trace_tw_class, util_wave_tw_class)
                       .with_meta_targets([QC1_CALIB, CALCHECKER])
                       .build())

    # This task reduces Uranium-Neon wavelength calibrations.
    wavelength_une = (task("wave_uranium_neon")
                      .with_recipe("cr2res_cal_wave")
                      .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
                      .with_report("crires_wavelength", ReportInput.RECIPE_INPUTS_OUTPUTS)
                      .with_main_input(raw_wave_une)
                      .with_associated_input(dark, [CAL_DARK_MASTER, CAL_DARK_BPM], match_rules=match_dark_for_calibs)
                      .with_alternatives(flat_calibrations)
                      .with_associated_input(emission_lines, min_ret=0)
                      .with_associated_input(util_wave_tw, min_ret=0)
                      .with_associated_input(util_trace_tw, min_ret=0)
                      .with_input_filter(CAL_FLAT_MASTER, CAL_FLAT_TW, CAL_DARK_MASTER, CAL_DARK_BPM,
                                         emission_lines_class, util_wave_tw_class)
                      .with_meta_targets([QC1_CALIB])
                      .build())

    # This collects wavelength calibrations from gas cells. It is used only for monitoring and associated to the relevant
    # science tasks, but it is not used not in data processing.
    wavelength_gas_cell = (task("wave_gas_cell")
                           .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
                           .with_main_input(raw_wave_gas_cell)
                           .with_associated_input(dark, [CAL_DARK_MASTER], match_rules=match_dark_for_calibs)
                           .with_alternatives(flat_calibrations)
                           .with_associated_input(emission_lines, min_ret=0)
                           .with_associated_input(util_wave_tw, min_ret=0)
                           .with_associated_input(util_trace_tw, min_ret=0)
                           .with_meta_targets([QC1_CALIB])
                           .build())

    # This tasks collects wavelenght calibrations from sky. It is used for monitoring but not in the data processing.
    wavelength_sky = (task("wave_sky")
                      .with_report("crires_rawdisp", ReportInput.RECIPE_INPUTS)
                      .with_main_input(raw_wave_sky)
                      # .with_alternatives(flat_calibrations)
                      # .with_associated_input(emission_lines, min_ret=0)
                      # .with_associated_input(util_wave_tw, min_ret=0)
                      # .with_associated_input(util_trace_tw, min_ret=0)
                      .with_meta_targets([QC1_CALIB])
                      .build())

    # The alternative_association "wavelength_calibrations" defines the preferences for the associations to be
    # used for wavelength calibration in the following manner:
    # A) For short-wavelength calibrations (L and M bands):
    #   - The first preference is to associate the task "wave_fabry_perot_etalon" (which combines the information of a
    #     Fabry-Perot and a Uranium-Neon observation of the same template) within 2.5 days of the science or standard
    #     star observation.
    #   - The second preference is to associate the task "wave_uranium_neon" (which contains only the information of the
    #     Uranium-Neon lamp) within 2.5 days of the science or standard star observation.
    #   - The third and fourth preferences are as the first two, but with validity ranges of 10 days.
    #   - The fifth and sixth preferences are as the first two, but without time restrictions.
    # B) For long-wavelength calibrations (Y, H, J and K bands):
    #   The task "wave_gas_cell" closest in time to the science or standard star observation is associated. The quality
    #   level of the association is determined by the time ranges.

    match_wave_7days_l0 = (match_rules()
                           .with_match_keywords(match_flat, time_range=ONE_WEEK, level=0))
    match_wave_7days_l1 = (match_rules()
                           .with_match_keywords(match_flat, time_range=ONE_WEEK, level=1))
    match_wave_10days = (match_rules()
                         .with_match_keywords(match_flat, time_range=RelativeTimeRange(-10, 10), level=2))
    match_wave_unlimited = (match_rules()
                            .with_match_keywords(match_flat, time_range=UNLIMITED, level=3))

    wavelength_calibrations = (alternative_associated_inputs()
                               .with_associated_input(wavelength_fpet, [CAL_WAVE_TW], match_rules=match_wave_7days_l0,
                                                      condition=is_short_wavelength)
                               .with_associated_input(wavelength_une, [CAL_WAVE_UNE], match_rules=match_wave_7days_l1,
                                                      condition=is_short_wavelength)
                               .with_associated_input(wavelength_fpet, [CAL_WAVE_TW], match_rules=match_wave_10days,
                                                      condition=is_short_wavelength)
                               .with_associated_input(wavelength_une, [CAL_WAVE_UNE], match_rules=match_wave_10days,
                                                      condition=is_short_wavelength)
                               .with_associated_input(wavelength_fpet, [CAL_WAVE_TW], match_rules=match_wave_unlimited,
                                                      condition=is_short_wavelength)
                               .with_associated_input(wavelength_une, [CAL_WAVE_UNE], match_rules=match_wave_unlimited,
                                                      condition=is_short_wavelength)
                               .with_associated_input(wavelength_gas_cell, [CAL_WAVE_GAS_CELL],
                                                      match_rules=match_wave_7days_l0, condition=is_long_wavelength)
                               .with_associated_input(wavelength_gas_cell, [CAL_WAVE_GAS_CELL],
                                                      match_rules=match_wave_10days, condition=is_long_wavelength)
                               .with_associated_input(wavelength_gas_cell, [CAL_WAVE_GAS_CELL],
                                                      match_rules=match_wave_unlimited, condition=is_long_wavelength))

    return wavelength_calibrations, wavelength_fpet, wavelength_une, wavelength_gas_cell, wavelength_sky
