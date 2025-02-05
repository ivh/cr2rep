from edps import classification_rule

from . import crires_keywords as kwd
from . import crires_rules as rules

# Dictionaries containing the values of header keywords that define calibrations and science data
# in the red and blue arms.
crires = {kwd.instrume: "CRIRES"}
calib_keywords = {**crires, kwd.pro_catg: None, kwd.dpr_catg: "CALIB"}
science_keywords = {**crires, kwd.pro_catg: None, kwd.dpr_catg: "SCIENCE"}
calib_ima_keywords = {**calib_keywords, kwd.dpr_tech: "IMAGE", kwd.det_read_curname: "New_RR_UpTheRamp"}
calib_spc_keywords = {**calib_keywords, kwd.dpr_tech: "SPECTRUM", kwd.det_read_curname: "New_RR_UpTheRamp"}
calib_other_keywords = {**calib_keywords, kwd.dpr_tech: "SPECTRUM,NODDING,OTHER",
                        kwd.det_read_curname: "New_RR_UpTheRamp"}
calib_jitter_keywords = {**calib_keywords, kwd.dpr_tech: "SPECTRUM,NODDING,JITTER",
                         kwd.det_read_curname: "New_RR_UpTheRamp"}
science_nod_other_keywords = {**science_keywords, kwd.dpr_tech: "SPECTRUM,NODDING,OTHER",
                              kwd.det_read_curname: "New_RR_UpTheRamp"}
science_nod_jitter_keywords = {**science_keywords, kwd.dpr_tech: "SPECTRUM,NODDING,JITTER",
                               kwd.det_read_curname: "New_RR_UpTheRamp"}
science_stare_other_keywords = {**science_keywords, kwd.dpr_tech: "SPECTRUM,DIRECT,OTHER",
                                kwd.det_read_curname: "New_RR_UpTheRamp"}
science_stare_jitter_keywords = {**science_keywords, kwd.dpr_tech: "SPECTRUM,DIRECT,JITTER",
                                 kwd.det_read_curname: "New_RR_UpTheRamp"}

# ---- Classification rules  -----------------
# Raw types
dark_class = classification_rule('DARK', {**calib_ima_keywords, kwd.dpr_type: "DARK"})
flat_class = classification_rule('FLAT', rules.is_flat)
wave_une_class = classification_rule('WAVE_UNE', {**calib_spc_keywords, kwd.dpr_type: "WAVE,UNE"})
wave_fpet_class = classification_rule('WAVE_FPET', {**calib_spc_keywords, kwd.dpr_type: "WAVE,FPET"})
wave_sgc_class = classification_rule('WAVE_SGC', {**calib_spc_keywords, kwd.dpr_type: "WAVE,ABSORPTION_SGC"})
wave_n20_class = classification_rule('WAVE_N20', {**calib_spc_keywords, kwd.dpr_type: "WAVE,ABSORPTION_N2O"})
wave_sky = classification_rule('WAVE_SKY', {**calib_keywords, kwd.dpr_type: "WAVE,SKY"})
wave_other_class = classification_rule('WAVE_OTHER', {**calib_spc_keywords, kwd.dpr_type: "WAVE,ABSORPTION_OTHER"})
detlin_lamp_class = classification_rule("DETLIN_LAMP", {**calib_keywords, kwd.dpr_type: "FLAT,LAMP,DETCHECK"})
detlin_dark_class = classification_rule("DETLIN_DARK", {**calib_keywords, kwd.dpr_type: "DARK,DETCHECK"})
wave_une_rassoc_class = classification_rule('WAVE_UNE_RASSOC', {**crires, kwd.dpr_type: "WAVE,UNE"})

std_nod_other_class = classification_rule('CAL_NODDING_OTHER', {**calib_other_keywords,
                                                                kwd.dpr_type: "STD"})

std_nod_jitter_class = classification_rule('CAL_NODDING_JITTER', {**calib_jitter_keywords,
                                                                  kwd.dpr_type: "STD"})
std_polarimetry_class = classification_rule("OBS_POLARIMETRY_OTHER", {**crires, kwd.dpr_catg: "CALIB",
                                                                      kwd.dpr_type: "STD,POL"})

sci_nod_other_class = classification_rule('OBS_NODDING_OTHER', {**science_nod_other_keywords,
                                                                kwd.dpr_type: "OBJECT"})
sci_nod_jitter_class = classification_rule('OBS_NODDING_JITTER', {**science_nod_jitter_keywords,
                                                                  kwd.dpr_type: "OBJECT"})
sci_staring_other_class = classification_rule('OBS_STARING_OTHER', {**science_stare_other_keywords,
                                                                    kwd.dpr_type: "OBJECT"})
sci_staring_jitter_class = classification_rule('OBS_STARING_JITTER', {**science_stare_jitter_keywords,
                                                                      kwd.dpr_type: "OBJECT"})
sci_astrometry_other_class = classification_rule('OBS_ASTROMETRY_OTHER', rules.is_science_astrometry_other)
sci_astrometry_jitter_class = classification_rule('OBS_ASTROMETRY_JITTER', rules.is_science_astrometry_jitter)
sci_polarimetry_class = classification_rule('OBS_POLARIMETRY_OTHER', rules.is_science_polarimetry_other)
sci_2d_object_class = classification_rule('OBS_2D_OBJECT', rules.is_science_2d_object)
sci_2d_sky_class = classification_rule('OBS_2D_SKY', rules.is_science_2d_sky)

lamp_metrology_class = classification_rule('METROLOGY', {**calib_spc_keywords, kwd.dpr_type: "LAMP,METROLOGY"})
acquisition_class = classification_rule("ACQUISITION", {**crires, kwd.dpr_catg: "ACQUISITION"})

# master calibrations
CAL_DARK_MASTER = classification_rule('CAL_DARK_MASTER', {**crires, kwd.pro_catg: 'CAL_DARK_MASTER'})
CAL_DARK_BPM = classification_rule('CAL_DARK_BPM', {**crires, kwd.pro_catg: 'CAL_DARK_BPM'})
CAL_FLAT_MASTER = classification_rule('CAL_FLAT_MASTER', {**crires, kwd.pro_catg: 'CAL_FLAT_MASTER'})
CAL_FLAT_BPM = classification_rule('CAL_FLAT_BPM', {**crires, kwd.pro_catg: 'CAL_FLAT_BPM'})
CAL_FLAT_TW = classification_rule('CAL_FLAT_TW', {**crires, kwd.pro_catg: 'CAL_FLAT_TW'})
CAL_FLAT_EXTRACT_1D = classification_rule('CAL_FLAT_EXTRACT_1D', {**crires, kwd.pro_catg: 'CAL_FLAT_EXTRACT_1D'})

CAL_WAVE_TW = classification_rule('CAL_WAVE_TW', {**crires, kwd.pro_catg: 'CAL_WAVE_TW', kwd.object: "WAVE,FPET"})
CAL_WAVE_UNE = classification_rule('CAL_WAVE_TW', {**crires, kwd.pro_catg: 'CAL_WAVE_TW', kwd.object: "WAVE,UNE"})
CAL_WAVE_GAS_CELL = classification_rule('CAL_WAVE_TW', rules.is_gas_cell)

# Static calibrations
util_wave_tw_class = classification_rule('UTIL_WAVE_TW', {**crires, kwd.pro_catg: 'UTIL_WAVE_TW'})
util_trace_tw_class = classification_rule('UTIL_TRACE_TW', {**crires, kwd.pro_catg: 'UTIL_TRACE_TW'})
emission_lines_class = classification_rule('EMISSION_LINES', {**crires, kwd.pro_catg: 'EMISSION_LINES'})
photo_flux_class = classification_rule('PHOTO_FLUX', {**crires, kwd.pro_catg: 'PHOTO_FLUX'})
linearity_coefficients = classification_rule('CAL_DETLIN_COEFFS', {**crires, kwd.pro_catg: 'CAL_DETLIN_COEFFS'})
