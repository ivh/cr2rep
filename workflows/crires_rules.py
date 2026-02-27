from edps import match

from . import crires_keywords as kwd


def is_gas_cell(f):
    return f[kwd.pro_catg] == "CAL_WAVE_TW" and f[kwd.object] not in ["WAVE,FPET", "WAVE,UNE", "WAVE,SKY"]


# ASSOCIATION RULES
# first:  ref=trigger (e.g. science)
# second: f  =file to associate (e.g. calibration)

def assoc_dark(ref, f):
    return match(ref, f, [kwd.ins_wlen_id, kwd.det_seq1_dit]) and f[kwd.ins_slit1_id] == "closed"


def assoc_flat(ref, f):
    # Deep-flats with det.ndit=50 are not considered for association
    return match(ref, f, [kwd.ins_wlen_id, kwd.ins_slit1_id]) and f[kwd.det_ndit] != 50
