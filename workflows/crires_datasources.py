from edps import data_source, match_rules
from edps.generator.time_range import *

from .crires_classification import *

# Convention for Data sources Association rule levels:
# Each data source can have several match function which correspond to different
# quality levels for the selected data. The level is specified as a number that
# follows this convention:
#   level < 0: more restrictive than the calibration plan
#   level = 0: follows the calibration plan
#   level = 1: quality sufficient for science reduction
#   level = 2: probably still acceptable quality
#   level = 3: significant risk of bad quality results


arm = [kwd.ins_path]
binning = [kwd.det_binx, kwd.det_biny]
arm_binning = arm + binning
instrument = [kwd.instrume]
group_cal_common = [kwd.tpl_start, kwd.ins_wlen_id, kwd.ins_slit1_id, kwd.det_seq1_dit]
group_detlin = instrument + [kwd.det_read_curname, kwd.tpl_start]

# match_flat and matching keywords:
match_flat = [kwd.ins_wlen_id, kwd.ins_slit1_id]
setup = match_flat + [kwd.det_seq1_dit]
setup_detlin = instrument + [kwd.det_read_curname]
match_grating_detector = arm_binning + [kwd.instrume, kwd.grat1_wlen, kwd.grat2_wlen]

# --- Datasources for raw exposures ------------------------------------------------------

# Datasource for detector monitoring
raw_detlin = (data_source("DETECTOR_LINEARITY")
              .with_classification_rule(detlin_lamp_class)
              .with_classification_rule(detlin_dark_class)
              .with_setup_keywords(setup_detlin)
              .with_grouping_keywords(group_detlin)
              .with_match_keywords([kwd.instrume], time_range=RelativeTimeRange(-100, 100), level=0)
              .with_match_keywords([kwd.instrume], time_range=RelativeTimeRange(-180, 180), level=1)
              .build())

# -- darks --
# To measure dark current
raw_dark = (data_source('DARK')
            .with_classification_rule(dark_class)
            .with_grouping_keywords(group_cal_common)
            .with_setup_keywords(setup)
            .build())

# Calibration tasks and science tasks have different requirements for matching darks
match_dark_for_calibs = (match_rules()
                         .with_match_keywords(setup, time_range=RelativeTimeRange(-1.5, 1.5), level=0)
                         .with_match_keywords(setup, time_range=TWO_DAYS, level=1)
                         .with_match_keywords(setup, time_range=UNLIMITED, level=3))

match_dark_for_science = (match_rules()
                          .with_match_function(rules.assoc_dark, time_range=ONE_DAY, level=0)
                          .with_match_function(rules.assoc_dark, time_range=RelativeTimeRange(-1.5, 1.5), level=1)
                          .with_match_function(rules.assoc_dark, time_range=UNLIMITED, level=3))

match_dark_for_polarimetry = (match_rules()
                              .with_match_keywords(match_flat, time_range=ONE_DAY, level=0)
                              .with_match_keywords(match_flat, time_range=RelativeTimeRange(-1.5, 1.5), level=1)
                              .with_match_keywords(match_flat, time_range=UNLIMITED, level=3))

raw_flat = (data_source('FLAT')
            .with_classification_rule(flat_class)
            .with_grouping_keywords(group_cal_common)
            .with_setup_keywords(setup + [kwd.det_ndit])
            .with_match_function(rules.assoc_flat, time_range=ONE_WEEK, level=0)
            .with_match_function(rules.assoc_flat, time_range=RelativeTimeRange(-10, 10), level=2)
            .with_match_function(rules.assoc_flat, time_range=UNLIMITED, level=3)
            .build())

raw_wave_une = (data_source('WAVE_Uranium_Neon')
                .with_classification_rule(wave_une_class)
                .with_grouping_keywords(group_cal_common)
                .with_setup_keywords(setup + [kwd.object, kwd.ocs_mtrlgy_st])
                .with_match_keywords(match_flat + [kwd.tpl_start], time_range=ONE_DAY, level=0)
                .build())

raw_wave_fpet = (data_source('WAVE_FabryPerot_Etalon')
                 .with_classification_rule(wave_fpet_class)
                 .with_grouping_keywords(group_cal_common + [kwd.dpr_type])
                 .with_setup_keywords(setup + [kwd.object, kwd.ocs_mtrlgy_st])
                 .build())

# These wavelength calbirations are not processed by any recipe, but they are used for monitoring purposes
raw_wave_gas_cell = (data_source('WAVE_GAS_CELL')
                     .with_classification_rule(wave_sgc_class)
                     .with_classification_rule(wave_n20_class)
                     .with_classification_rule(wave_other_class)
                     .with_grouping_keywords(group_cal_common + [kwd.dpr_type])
                     .with_setup_keywords(setup)
                     .build())

raw_wave_sky = (data_source('WAVE_SKY')
                .with_classification_rule(wave_sky)
                .with_grouping_keywords(group_cal_common + [kwd.dpr_type])
                .with_setup_keywords(setup)
                .build())

raw_wave_une_ref = (data_source('REFERENCE_Uranium_Neon')
                    .with_classification_rule(wave_une_rassoc_class)
                    .with_grouping_keywords([kwd.mjd_obs])
                    .with_match_keywords(match_flat + [kwd.tpl_start], time_range=ONE_DAY, level=0)
                    .build())

raw_standard = (data_source('STANDARD_STAR')
                .with_classification_rule(std_nod_other_class)
                .with_classification_rule(std_nod_jitter_class)
                .with_min_group_size(2)
                .with_grouping_keywords(group_cal_common + [kwd.dpr_tech])
                .with_setup_keywords(match_flat)
                .with_match_keywords(match_flat, time_range=RelativeTimeRange(-90, 90), level=0)
                .with_match_keywords(match_flat, time_range=UNLIMITED, level=3)
                .build())

raw_polarimetry_standard = (data_source('STANDARD_STAR_POL')
                            .with_classification_rule(std_polarimetry_class)
                            .with_min_group_size(8)
                            .with_grouping_keywords(group_cal_common)
                            .with_setup_keywords(match_flat)
                            .build())

raw_science_nodding = (data_source('SCIENCE_NODDING')
                       .with_classification_rule(sci_nod_other_class)
                       .with_classification_rule(sci_nod_jitter_class)
                       .with_min_group_size(2)
                       .with_grouping_keywords(group_cal_common + [kwd.dpr_tech])
                       .with_setup_keywords(match_flat)
                       .build())

raw_science_staring = (data_source('SCIENCE_STARING')
                       .with_classification_rule(sci_staring_jitter_class)
                       .with_classification_rule(sci_staring_other_class)
                       .with_grouping_keywords(group_cal_common + [kwd.dpr_tech])
                       .with_setup_keywords(match_flat)
                       .build())

raw_sci_polarimetry = (data_source('SCIENCE_POLARIMETRY')
                       .with_classification_rule(sci_polarimetry_class)
                       .with_grouping_keywords(group_cal_common)
                       .with_setup_keywords(match_flat)
                       .build())

raw_sci_astrometry = (data_source('SCIENCE_ASTROMETRY')
                      .with_classification_rule(sci_astrometry_other_class)
                      .with_classification_rule(sci_astrometry_jitter_class)
                      .with_grouping_keywords(group_cal_common + [kwd.dpr_tech])
                      .with_setup_keywords(match_flat)
                      .build())

raw_science_2d = (data_source('SCIENCE_2D')
                  .with_classification_rule(sci_2d_object_class)
                  .with_grouping_keywords(group_cal_common + [kwd.dpr_type])
                  .with_setup_keywords(match_flat)
                  .build())

raw_science_2d_sky = (data_source('SKY_2D')
                      .with_classification_rule(sci_2d_sky_class)
                      .with_grouping_keywords(group_cal_common + [kwd.dpr_type])
                      .with_match_keywords([kwd.tpl_start], level=0)
                      .build())

# --- Data sources for static calibrations

util_wave_tw = (data_source()
                .with_classification_rule(util_wave_tw_class)
                .with_match_keywords([kwd.ins_wlen_id], time_range=UNLIMITED, level=0)
                .build())

util_trace_tw = (data_source()
                 .with_classification_rule(util_trace_tw_class)
                 .with_match_keywords([kwd.ins_wlen_id], time_range=UNLIMITED, level=0)
                 .build())

emission_lines = (data_source()
                  .with_classification_rule(emission_lines_class)
                  .with_match_keywords([kwd.ins_wlen_id], time_range=UNLIMITED, level=0)
                  .build())

photo_flux = (data_source()
              .with_classification_rule(photo_flux_class)
              .with_match_keywords(instrument, time_range=UNLIMITED, level=0)
              .build())
