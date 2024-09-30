from adari_core.plots.panel import Panel
from adari_core.plots.points import LinePlot
from adari_core.plots.text import TextPlot
from adari_core.plots.cut import CutPlot
from adari_core.plots.histogram import HistogramPlot
from adari_core.plots.images import ImagePlot, CentralImagePlot
from adari_core.plots.points import ScatterPlot
from adari_core.data_libs.master_std_star_ifu import MasterSpecphotStdReport
from adari_core.utils.utils import fetch_kw_or_error, fetch_kw_or_default

from .crires_utils import CriresSetupInfo
import astropy.io.fits as fits
import os
import numpy as np

from . import CriresReportMixin

class CriresSpecphotStdReport(CriresReportMixin, MasterSpecphotStdReport):
    def __init__(self):
        super().__init__("crires_std_star")
        self.data_extensions = {
            "CHIP1.INT1": "CHIP1",
            "CHIP2.INT1": "CHIP2",
            "CHIP3.INT1": "CHIP3"
        }

    def parse_sof(self):
        raw_flat = None
        combine_a = None
        combine_b = None
        extract_comb = None
        
        for filename, catg in self.inputs:
            if catg == f"CAL_NODDING_OTHER" or catg == f"CAL_NODDING_JITTER":
                raw_flat = filename
                            
            if catg == f"OBS_NODDING_COMBINEDA":
                combine_a = filename
        
            elif catg == f"OBS_NODDING_COMBINEDB":
                combine_b = filename
        
            elif catg == f"OBS_NODDING_EXTRACT_COMB":
                extract_comb = filename
        
        file_lists = []
        raw_extension_default = []
        if combine_a is not None and combine_b is not None and extract_comb is not None and raw_flat is not None:
            file_lists.append(
                {
                    "combine_a": combine_a,
                    "combine_b": combine_b,
                    "extract_comb": extract_comb,
                    "raw_flat": raw_flat
                }
            )
        return file_lists

    def generate_panels(self, **kwargs):
        panels = {}
        combine_a = self.hdus[0]["combine_a"]
        combine_b = self.hdus[0]["combine_b"]
        ext_comb = self.hdus[0]["extract_comb"]
        raw_flat = self.hdus[0]["raw_flat"]
        # Generate indidual raw-cuts
        for ext in self.data_extensions:
            data_a = combine_a[ext].data
            data_b = combine_b[ext].data
            ext_comb_cols = ext_comb[ext].columns.names

            spec_arr = []
            wl_arr = []
            col_names = self.hdus[0]['extract_comb'][ext].columns.names
            for i in col_names:
                if "SPEC" in i:
                    spec_arr.append(i)
                    #spec = hdus[0]['extract_comb']["CHIP1.INT1"].data.field(i)
                if "WL" in i:
                    wl_arr.append(i)
                    #wl = hdus[0]['extract_comb']["CHIP1.INT1"].data.field(i)
            hr=[1]
            hr.extend([4]*len(spec_arr))
            p = Panel(4, len(spec_arr)+1, height_ratios=hr)
            
            # Text Plot
            vspace = 0.5
            t1 = TextPlot(columns=1, v_space=vspace)
            fname = os.path.basename(str(ext_comb.filename()))
            rname = os.path.basename(str(raw_flat.filename()))
            pext = "PRIMARY"
            hdr = ext_comb[pext].header
            col1 = (
                str(hdr.get("INSTRUME")),
                "EXTNAME: " + ext,
                "PRO CATG: " + str(hdr.get("HIERARCH ESO PRO CATG")),
                "FILE NAME: " + fname,
                "RAW1.NAME: " + str(hdr.get("HIERARCH ESO PRO REC1 RAW1 NAME"))
            )
            t1.add_data(col1)
            p.assign_plot(t1, 0, 0, xext=2)
    
            t2 = TextPlot(columns=1, v_space=vspace, xext=1)
            self.metadata = CriresSetupInfo.standard_star(ext_comb)
            col2 = self.metadata
            t2.add_data(col2)
            p.assign_plot(t2, 2, 0, xext=1)

            for i in range(len(spec_arr)):
                spec = self.hdus[0]['extract_comb'][ext].data.field(spec_arr[i])
                wl = self.hdus[0]['extract_comb'][ext].data.field(wl_arr[i])
                if np.all(np.isnan(spec)==True):
                    print("No data to plot.")
                else:
                    self.wavelength = wl
                    self.wavelength_unit = "nm"
                    line_plot = self.plot_vs_wavelength(data=spec, label=spec_arr[i], title="OBS_NODDING_EXTRACT_COMB")
                    p.assign_plot(line_plot, 0, i+1, xext=2)

            image_plot_a = ImagePlot(
                combine_a[ext].data,
                v_clip="percentile",
                v_clip_kwargs={"percentile": 95.0},                
                title=f"OBS_NODDING_COMBINEDA",
            )
            image_plot_a.interp = "nearest"
            p.assign_plot(image_plot_a, 2, 1, xext=2, yext=2)

            image_plot_b = ImagePlot(
                combine_b[ext].data,
                v_clip="percentile",
                v_clip_kwargs={"percentile": 95.0},
                title=f"OBS_NODDING_COMBINEDB",
            )
            image_plot_b.interp = "nearest"
            p.assign_plot(image_plot_b, 2, 3, xext=2, yext=2)

            raw_hist = HistogramPlot(
                raw_data=raw_flat[ext].data,
                title="Raw value counts",
                bins=50,
                v_min=-5000,
                v_max=70000,
                legend=False,
            )
            p.assign_plot(raw_hist, 2, 5, xext=2, yext=2)
            
            panels[p] = {
                "combine_a": combine_a.filename(),
                "combine_b": combine_b.filename(),
                "extract_comb": ext_comb.filename(),
                "ext": ext,
                "report_name": f"CRIRES_{ext}",
                "report_description": f"Specphotometric Std panel - ({combine_a.filename()}, "f"{combine_b.filename()}, "f"{ext})",
                "report_tags": []
            }
            
        return panels

rep = CriresSpecphotStdReport()
