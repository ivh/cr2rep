from edps import match

from . import crires_keywords as kwd


def is_crires(f):
    return f[kwd.instrume] == "CRIRES"


def is_raw(f):
    return f[kwd.pro_catg] is None


def is_calib(f):
    return is_crires(f) and is_raw(f) and f[kwd.dpr_catg] == "CALIB"


def is_science(f):
    return is_crires(f) and is_raw(f) and f[kwd.dpr_catg] == "SCIENCE"


def is_technical(f):
    return is_crires(f) and is_raw(f) and f[kwd.dpr_catg] == "TECHNICAL"


def is_test(f):
    return is_crires(f) and is_raw(f) and f[kwd.dpr_catg] == "TEST"


def is_flat(f):
    return is_calib(f) and \
           (f[kwd.dpr_type] == "FLAT" and f[kwd.dpr_tech] == 'SPECTRUM' and
            f[kwd.det_read_curname] == "New_RR_UpTheRamp" and f[kwd.ins_opti8_id] == "Open"
            and kwd.det_ndit != 50)


def is_science_astrometry_other(f):
    return is_science(f) and \
           (f[kwd.dpr_type] == "OBJECT" and (f[kwd.dpr_tech] == 'SPECTRUM,NODDING,OTHER,ASTROMETRY' or
                                             f[kwd.dpr_tech] == 'SPECTRUM,NODDING,OTHER,ASTROM') and
            f[kwd.det_read_curname] == "New_RR_UpTheRamp")


def is_science_astrometry_jitter(f):
    return is_science(f) and \
           (f[kwd.dpr_type] == "OBJECT" and (f[kwd.dpr_tech] == 'SPECTRUM,NODDING,JITTER,ASTROMETRY' or
                                             f[kwd.dpr_tech] == 'SPECTRUM,NODDING,JITTER,ASTROM') and
            f[kwd.det_read_curname] == "New_RR_UpTheRamp")


def is_science_polarimetry_other(f):
    return is_science(f) and \
           (f[kwd.dpr_type] == "OBJECT" and (f[kwd.dpr_tech] == 'SPECTRUM,NODDING,OTHER,POLARIMETRY' or
                                             f[kwd.dpr_tech] == 'SPECTRUM,NODDING,OTHER,POLARI') and
            f[kwd.det_read_curname] == "New_RR_UpTheRamp")


def is_science_2d_object(f):
    return is_science(f) and \
           (f[kwd.dpr_type] == "OBJECT" and (f[kwd.dpr_tech] == 'SPECTRUM,GENERIC') and
            f[kwd.det_read_curname] == "New_RR_UpTheRamp")


def is_science_2d_sky(f):
    return is_science(f) and \
           (f[kwd.dpr_type] == "SKY" and (f[kwd.dpr_tech] == 'SPECTRUM,GENERIC') and
            f[kwd.det_read_curname] == "New_RR_UpTheRamp")


def is_gas_cell(f):
    return f[kwd.pro_catg] == "CAL_WAVE_TW" and f[kwd.object] not in ["WAVE,FPET", "WAVE,UNE", "WAVE,SKY"]


# ASSOCIATION RULES
# first:  ref=trigger (e.g. science)
# second: f  =file to associate (e.g. calibration)

def assoc_dark(ref, f):
    return match(ref, f, [kwd.ins_wlen_id, kwd.det_seq1_dit]) and f[kwd.ins_lamp11_st] == "closed"
