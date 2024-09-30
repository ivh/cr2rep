from adari_core.utils.utils import fetch_kw_or_default

extensions = ["CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1"]

def md_1(hdul):
    metadata = [
        "INS.WLEN.ID: " + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO INS WLEN ID", default="N/A")),
        "INS.SLIT1.ID: " + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO INS SLIT1 ID", default="N/A")),
        "DET.SEQ1.DIT: " + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO DET SEQ1 DIT", default="N/A")),
    ]
    return metadata

def md_2(hdul):
    metadata = md_1(hdul)
    metadata.append("PRO.DATANCOM: " + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO PRO DATANCOM", default= "N/A")))
    return metadata

class CriresSetupInfo:

    def detmon(hdul):
        metadata = [
            "DET.READ.CURNAM: "
            + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO DET READ CURNAME",
 default="N/A")),
        ]
        return metadata


    def master_dark(hdul):
        metadata = md_2(hdul)
        return metadata

    def master_flat(hdul):
        metadata = md_2(hdul)
        return metadata


    def wavelength(hdul):
        metadata = [
            "INS.WLEN.ID: "
            + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO INS WLEN ID",
 default="N/A")),
            "INS.SLIT1.ID: "
            + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO INS SLIT1 ID",
 default="N/A")),
            "DET.SEQ1.DIT: "
            + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO DET SEQ1 DIT",
 default="N/A")),
        ]
        return metadata

    def flat(hdul):
        metadata = md_1(hdul)
        return metadata

    def dark(hdul):
        metadata = md_1(hdul)
        return metadata

    def detector_linearity(hdul):
        metadata = [
            "DET.READ.CURNAME: " + str(fetch_kw_or_default(hdul["PRIMARY"], "HIERARCH ESO DET READ CURNAME", default= "N/A")),
        ]
        return metadata

    def wave_uranium_neon(hdul):
        metadata = md_1(hdul)
        return metadata

    def wave_fabry_perot_etalon(hdul):
        metadata = md_1(hdul)
        return metadata

    def wave_gas_cell(hdul):
        metadata = md_1(hdul)
        return metadata

    def wave_sky(hdul):
        metadata = md_1(hdul)
        return metadata

    def standard_star(hdul):
        metadata = md_1(hdul)
        return metadata




