import math
import numpy as np
import os

from adari_core.data_libs.master_spec_science import MasterSpecScienceReport
from adari_core.plots.images import ImagePlot
from adari_core.plots.panel import Panel
from adari_core.plots.points import LinePlot
from adari_core.plots.text import TextPlot
from adari_core.utils.clipping import clipping_percentile, clipping_sigma
from adari_core.utils.utils import fetch_kw_or_default

from .crires_utils import CriresReportMixin

class CriresScienceReport(CriresReportMixin, MasterSpecScienceReport):
    def __init__(self):
        super().__init__("crires_science")
        self.nodding = False

    def parse_sof(self):
        sci_spec = None
        det_img = []
        tw = []
        caltw = None
        comb = None

        for filename, catg in self.inputs:
            if catg == "OBS_NODDING_EXTRACTC_IDP" or catg == "OBS_STARING_IDP":
                sci_spec = filename
            if catg == "OBS_NODDING_COMBINEDA" or catg == "OBS_NODDING_COMBINEDB" or catg == "OBS_STARING_COMBINED":
                det_img.append(filename) 
            if catg == "OBS_NODDING_TWA" or catg == "OBS_NODDING_TWB":
                tw.append(filename)
            if catg == "CAL_WAVE_TW":
                caltw = filename
            if catg == "OBS_NODDING_EXTRACT_COMB" or catg == "OBS_STARING_COMBINED":
                comb = filename

       
        file_dict = {}
        if sci_spec is not None:
            file_dict["science"] = sci_spec
        if len(det_img) == 1:
            file_dict["images"] = det_img[0]
        else:
            self.nodding = True
            file_dict["imagesA"] = sorted(det_img)[0]
            file_dict["imagesB"] = sorted(det_img)[1]
        if self.nodding:
            file_dict["twA"] = sorted(tw)[0]
            file_dict["twB"] = sorted(tw)[1]
        else:
            file_dict["tw"] = caltw
        if comb is not None:
            file_dict["comb"] = comb

        return [
            file_dict,
        ]

    def generate_panels(self, **kwargs):
        panels = {}
        
        science = self.hdus[0]["science"]
        science_fname = os.path.basename(str(science.filename()))

        data = science["SPECTRUM"].data
        orderidx = np.unique(data["ORDER"])[::-1]
        detecidx = np.unique(data["DETEC"])

        p = Panel(6, 7, height_ratios=[1.5, 3, 3, 3, 2, 2, 1], y_stretch=0.6)
   
        # Extracted spectra plot
        px = 0
        py = 1
        for i in orderidx:
            specplot = LinePlot(
                title="Extracted spectrum (order=%d)"%i,
            )
            for j in detecidx:
                m = (data["ORDER"]==i)&(data["DETEC"]==j)
                w = data["WAVE"][m]
                f = data["FLUX"][m]
                specplot.add_data([w, f], label="Detector %d"%j, linewidth=0.4)
            specplot.set_ylim(0, clipping_percentile(data["FLUX"][data["ORDER"]==i], 98)[1] * 1.2)
            specplot.x_label = "Wavelength (nm)"
            specplot.y_label = "ADU"
            p.assign_plot(specplot, px, py, xext=2)
            px += 2
            if px % 6 == 0:
                px = 0
                py += 1

        # Combined detector images
        px = 0
        py = 4
        if self.nodding:
            nod = ["A", "B"]
            aspect = 4.
        else:
            nod = [""]
            aspect = 2.
        nodtitle = ""
        for k in nod:
            tw = self.hdus[0]["tw"+k]
            img = self.hdus[0]["images"+k]
            if k != "":
                nodtitle = "Nod " + k + ", "
            for j in detecidx:
                detplot = ImagePlot(
                    title="Combined image: "+nodtitle+"Det "+str(j)+", Ord 5",
                    aspect=aspect,
                )
                ext = "CHIP"+str(j)+".INT1"
                imgdata = img[ext].data
                xsize = imgdata[:, 0].size
                twdata = tw[ext].data
                m = (twdata["Order"]==5)
                lower = twdata[m]["Lower"][0]
                upper = twdata[m]["Upper"][0]
                l0 = lower[0]
                l1 = lower[0] + lower[1] * xsize + lower[2] * xsize * xsize
                u0 = upper[0]
                u1 = upper[0] + upper[1] * xsize + upper[2] * xsize * xsize
                y0 = math.floor(min(l0, l1, u0, u1))
                y1 = math.ceil(max(l0, l1, u0, u1))
                detplot.add_data(imgdata)
                detplot.y_min = y0
                detplot.y_max = y1
                detplot.set_vlim(clipping_sigma(imgdata[y0:y1, :], nsigma=1))
                p.assign_plot(detplot, px, py, xext=2)
                px += 2
                if px % 6 == 0:
                    px = 0
                    py += 1

        # Upper Text Plot
        vspace = 0.3
        t0 = TextPlot(columns=1, v_space=vspace)
        col0 = (
            str(fetch_kw_or_default(science["PRIMARY"], "INSTRUME", default="N/A"))
            + " science product preview",
            "Product: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO PRO CATG", default="N/A")
            ),
            "Raw file: "
            + str(
                fetch_kw_or_default(
                    science["PRIMARY"], "ESO PRO REC1 RAW1 NAME", default="N/A"
                )
            ),
            "MJD-OBS: "
            + str(fetch_kw_or_default(science["PRIMARY"], "MJD-OBS", default="N/A")),
            "TPL ID: "
            + str(fetch_kw_or_default(science["PRIMARY"], "ESO TPL ID", default="N/A")),
            "RUN ID: "
            + str(
                fetch_kw_or_default(
                    science["PRIMARY"], "ESO OBS PROG ID", default="N/A"
                )
            ),
        )
        t0.add_data(col0, fontsize=13)
        p.assign_plot(t0, 0, 0, xext=1)

        t1 = TextPlot(columns=1, v_space=vspace)
        col1 = (
            "Target: "
            + str(fetch_kw_or_default(science["PRIMARY"], "OBJECT", default="N/A")),
            "OB ID: "
            + str(fetch_kw_or_default(science["PRIMARY"], "ESO OBS ID", default="N/A")),
            "OB NAME: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO OBS NAME", default="N/A")
            ),
            "WLEN ID: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO INS WLEN ID", default="N/A")
            ),
            "Slit width [arcsec]: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO INS SLIT1 WID", default="N/A")
            ),
            "AO loop: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO AOS RTC LOOP STATE", default="N/A")
            ),
        )
        t1.add_data(col1, fontsize=13)
        p.assign_plot(t1, 2, 0, xext=1)
        
        t2 = TextPlot(columns=1, v_space=vspace)
        col2 = (
            "Rotator mode: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO INS1 DROT MODE", default="N/A")
            ),
            "PRO TECH: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO PRO TECH", default="N/A")
            ),
            "Carrier position: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO INS1 OPTI1 NAME", default="N/A")
            ),
            "Metrology status: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO OCS MTRLGY ST", default="N/A")
            ),
            "Nod throw [arcsec]: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO SEQ NODTHROW", default="N/A")
            ),
            "Jitter width [arcsec]: "
            + str(
                fetch_kw_or_default(science["PRIMARY"], "ESO SEQ JITTER WIDTH", default="N/A")
            ),
        )
        t2.add_data(col2, fontsize=13)
        p.assign_plot(t2, 4, 0, xext=1)
        


        # Bottom Text Plot
        vspace = 0.5
        t4 = TextPlot(columns=1, v_space=vspace)
        col4 = (
            "Exp. time [s]: "
            + "%.1f"
            % fetch_kw_or_default(science["PRIMARY"], "TEXPTIME", default="N/A"),
            "N exposures: "
            + "%i"
            % fetch_kw_or_default(science["PRIMARY"], "NCOMBINE", default="N/A"),
            "S/N average: "
            + "%.2f"
            % fetch_kw_or_default(science["PRIMARY"], "SNR", default="N/A"),
            "Seeing: "
            + "%.2f"
            % fetch_kw_or_default(
                science["PRIMARY"], "ESO TEL IA FWHMLINOBS", default="N/A"
            ),
            "Airmass: "
            + "%.2f"
            % fetch_kw_or_default(
                science["PRIMARY"], "ESO TEL AIRM START", default="N/A"
            ),
        )
        t4.add_data(col4, fontsize=13)
        p.assign_plot(t4, 0, 6, xext=1)

        t5 = TextPlot(columns=1, v_space=vspace)
        col5 = (
            "Resolving power: "
            + "%.2f"
            % fetch_kw_or_default(science["PRIMARY"], "SPEC_RES", default="N/A"),
            "Lambda start [nm]: "
            + "%.2f"
            % fetch_kw_or_default(
                science["PRIMARY"], "WAVELMIN", default="N/A"
            ),
            "Lambda end [nm]: "
            + "%.2f"
            % fetch_kw_or_default(
                science["PRIMARY"], "WAVELMAX", default="N/A"
            ),
            "N pix sat: "
            + str(fetch_kw_or_default(
                science["PRIMARY"], "ESO QC NUMSAT", default="N/A")
            ),
        )
        t5.add_data(col5, fontsize=13)
        p.assign_plot(t5, 2, 6, xext=1)


        # Get FWHM
        comb = self.hdus[0]["comb"]
        fwhm1px = fetch_kw_or_default(comb["CHIP1.INT1"], "ESO QC SLITFWHM MED", default="N/A")
        fwhm1arc = 0.056 * fwhm1px
        fwhm2px = fetch_kw_or_default(comb["CHIP2.INT1"], "ESO QC SLITFWHM MED", default="N/A")
        fwhm2arc = 0.056 * fwhm2px
        fwhm3px = fetch_kw_or_default(comb["CHIP3.INT1"], "ESO QC SLITFWHM MED", default="N/A")
        fwhm3arc = 0.056 * fwhm3px

        t6 = TextPlot(columns=1, v_space=vspace)
        col6 = (
            "FWHM [px] det 1: " + "%.2f" % (fwhm1px),
            "FWHM [arcsec] det 1: " + "%.2f" % (fwhm1arc),
            "FWHM [px] det 2: " + "%.2f" % (fwhm2px),
            "FWHM [arcsec] det 2: " + "%.2f" % (fwhm2arc),
            "FWHM [px] det 3: " + "%.2f" % (fwhm3px),
            "FWHM [arcsec] det 3: " + "%.2f" % (fwhm3arc),
        )
        t6.add_data(col6, fontsize=13)
        p.assign_plot(t6, 4, 6, xext=1)

        input_files = [science.filename()]
        if self.nodding:
            input_files.append(self.hdus[0]["twA"].filename())
            input_files.append(self.hdus[0]["imagesA"].filename())
            input_files.append(self.hdus[0]["twB"].filename())
            input_files.append(self.hdus[0]["imagesB"].filename())
            input_files.append(self.hdus[0]["comb"].filename())
        else:
            input_files.append(self.hdus[0]["tw"].filename())
            input_files.append(self.hdus[0]["images"].filename())
            input_files.append(self.hdus[0]["comb"].filename())

        addme = {
            "report_name": f"CRIRES_{str(science_fname).removesuffix('.fits').lower()}",
            "report_description": "Science panel",
            "report_tags": [],
            "report_prodcatg": "ANCILLARY.PREVIEW",
            "input_files": input_files,
        }

        panels[p] = addme

        return panels


rep = CriresScienceReport()


