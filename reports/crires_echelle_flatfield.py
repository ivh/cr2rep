from adari_core.plots.text import TextPlot
from adari_core.plots.panel import Panel
from adari_core.data_libs.echelle_flatfield import MasterEchelleFlatfieldReport
from adari_core.data_libs.master_dark_bias import MasterDarkBiasReport

from adari_core.plots.images import ImagePlot
from astropy.io import fits
import numpy as np
from .crires_utils import CriresSetupInfo
from adari_core.utils.utils import fetch_kw_or_default
import logging

import os

from . import CriresReportMixin

logger = logging.getLogger(__name__)

center_size = 200

class CriresEchelleFlatfieldReport(CriresReportMixin, MasterEchelleFlatfieldReport, MasterDarkBiasReport):
    image_category = "master_product"
    raw_extension_default = 0
    im_clipping = "percentile"
    im_clipping_kwargs = {"percentile": 99}
    center_size = 200
    
    def __init__(self):
        super().__init__("crires_echelle_flatfield")
        self.center_size = center_size
        self.hist_bins_max = 20
        self.data_extensions = [
            "CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1"
        ]
        
    def parse_sof(self):
        raw_flat = None
        master_flat = None
        
        for filename, catg in self.inputs:
            if catg == f"CAL_FLAT_MASTER":
                master_flat = filename
        
            if catg == f"FLAT":
                raw_flat = filename
        
        file_lists = []
        raw_extension_default = []
        if master_flat is not None and raw_flat is not None:
            file_lists.append(
                {
                    "master_product": master_flat,
                    "raw": raw_flat,
                }
            )        
        return file_lists


    def generate_first_panel(self):
        panels = {}
        self.panel_kwargs = dict(x=3, y=3, height_ratios=[1, 4, 4])
        # Generate indidual raw-cuts
        for ext in self.data_extensions:
            new_panels = super().generate_panels(
                master_product_ext=ext, 
                raw_ext=self.raw_extension_default,
                direction="x",
                panel_kwargs = self.panel_kwargs,
                interpolation = "nearest",
            )

            for i, (panel, panel_descr) in enumerate(new_panels.items()):
                panel_descr["report_description"] = (
                    f"CRIRES echelle flatfield panel - "
                    f"{os.path.basename(panel_descr['master_product'])}, "
                    f"{panel_descr['master_product_ext']}"
                )
            
                master = self.hdus[i]["master_product"]
                raw = self.hdus[i]["raw"]
                hdr = master[ext].header

                plot = panel.retrieve(2, 1)
                plot.legend = True
                plot.add_data(np.nan_to_num(raw[ext].data))
    
                px = 0
                py = 0
                # which hdul and ext to use
                vspace = 0.3
                fname = os.path.basename(str(master.filename()))
                t1 = TextPlot(columns=1, v_space=vspace)
                col1 = (
                    str(master["PRIMARY"].header.get("INSTRUME")),
                    "EXTNAME: " + str(master[ext].header.get("EXTNAME", "N/A")),
                    "PRO CATG: "
                    + str(master["PRIMARY"].header.get("HIERARCH ESO PRO CATG")),
                    "FILE NAME: " + fname,
                    "RAW1 NAME: "
                    + str(
                        master["PRIMARY"].header.get(
                            "HIERARCH ESO PRO REC1 RAW1 NAME"
                        )
                    ),
                )
                t1.add_data(col1)
                panel.assign_plot(t1, px, py, xext=2)
    
                px = px + 2
                t2 = TextPlot(columns=1, v_space=vspace, xext=1)
                self.metadata = CriresSetupInfo.master_flat(master)
                col2 = self.metadata
                t2.add_data(col2)
                panel.assign_plot(t2, px, py, xext=1)

            panels = {**panels, **new_panels}

        return panels

    def generate_second_panel(self):

        panels = {}
        p = Panel(x=3, y=3, height_ratios=[1, 4, 4])

        # Generate the multi-extension panels
        channel = ["CHIP1", "CHIP2", "CHIP3"]

                # Metadata in Text Plot
        px, py = 0, 0
        vspace = 0.3
        t1 = TextPlot(columns=1, v_space=vspace)
        fname = os.path.basename(str(self.hdus[0]["master_product"].filename()))
        master_im = self.hdus[0]["master_product"]


        col1 = (
            str(master_im["PRIMARY"].header.get("INSTRUME")),
            "PRO CATG: "
            + str(master_im["PRIMARY"].header.get("HIERARCH ESO PRO CATG")),
            "FILE NAME: " + fname,
            "RAW1 NAME: "
            + str(
                master_im["PRIMARY"].header.get(
                    "HIERARCH ESO PRO REC1 RAW1 NAME"
                )
            ),
        )
        t1.add_data(col1)
        p.assign_plot(t1, 0, 0, xext=2)

        t2 = TextPlot(columns=1, v_space=vspace, xext=1)
        self.metadata = CriresSetupInfo.master_flat(master_im)
        col2 = self.metadata
        t2.add_data(col2)
        p.assign_plot(t2, 2, 0, xext=1)
        self.interpolation = "nearest"

        for i in range(len(self.data_extensions)):
            extname = self.data_extensions[i]

            full_plot, zoom_plot = super().image_plot(
                self.hdus[0][self.image_category][extname],
                zoom_in=True,
                zoom_in_extent=self.center_size,
                img_kwargs={
                    "title": extname,
                    "v_clip": self.im_clipping,
                    "v_clip_kwargs": self.im_clipping_kwargs,
                },
                zoom_img_kwargs={"v_clip": self.im_clipping,
                    "v_clip_kwargs": self.im_clipping_kwargs,
                                },
            )
            p.assign_plot(full_plot, i, 1, xext=1)
            p.assign_plot(zoom_plot, i, 2, xext=1)

            plt = p.retrieve(i, 1)
            plt.interp = "nearest"
            
            hdul = self.hdus[0][self.image_category]
            setup = str(hdul["PRIMARY"].header.get("HIERARCH ESO INS MODE", "N/A"))

            addme = {
                "report_name": "crires_echelle_flatfield_multi",
                "report_description": "CRIRES echelle flatfield multi panel",
                "report_tags": [],
            }

            panels[p] = addme

        return panels

    def generate_panels(self, **kwargs):
        """Create both single and multiple extension panels."""

        panels = {
            **self.generate_first_panel(),
            **self.generate_second_panel(),
        }
        return panels



rep = CriresEchelleFlatfieldReport()
