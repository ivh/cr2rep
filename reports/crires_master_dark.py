from adari_core.plots.text import TextPlot
from adari_core.data_libs.master_dark_bias import MasterDarkBiasReport
from adari_core.plots.panel import Panel

from astropy.io import fits
import numpy as np
from .crires_utils import CriresSetupInfo
from adari_core.utils.utils import fetch_kw_or_default
import logging

import os

from . import CriresReportMixin

logger = logging.getLogger(__name__)


class CriresMasterDarkReport(CriresReportMixin, MasterDarkBiasReport):
    image_category = "master_im"
    im_clipping = "percentile"
    im_clipping_kwargs = {"percentile": 99}
    zoom_in = True
    zoom_in_extent=200
    
    def __init__(self):
        super().__init__("crires_master_dark")
        self.center_size = 200
        self.hist_bins_max = 20
        self.data_extensions = [
            "CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1"
        ]
        
    def parse_sof(self):
        # Need to generate two report sets:
        # CAL_DARK_MASTER
        master_dark = None
        
        for filename, catg in self.inputs:
            if catg == "CAL_DARK_MASTER" and master_dark is None:
                master_dark = filename       
        
        # Build and return the file name list
        file_lists = []
        if master_dark is not None:
            file_lists.append(
                {
                    "master_im": master_dark,
                }
            )
        return file_lists

    def generate_first_panel(self):
        panels = {}
        # Generate indidual raw-cuts
        for ext in self.data_extensions:
            new_panels = super().generate_raw_cuts_panels(
                master_im_ext=ext,
                master_title="Master dark",
                master_im_clipping="percentile",
                master_im_n_clipping=95,
                hist_clipping="sigma",
                hist_n_clipping=5,
                cut_clipping="percentile",
                cut_n_clipping=95,
                cut_cent_clipping="percentile",
                cut_cent_n_clipping=95.0,
                collapse_clipping="sigma",
                collapse_n_clipping=5,
                interpolation = "nearest",
                )
            
            for i, (panel, panel_descr) in enumerate(new_panels.items()):
                panel_descr["report_description"] = (
                    f"CRIRES dark panel - "
                    f"{os.path.basename(panel_descr['master_im'])}, "
                    f"{panel_descr['master_im_ext']}"
                )
                master_im = self.hdus[i]["master_im"]
                

                # Text Plot
                px = 0
                py = 0
                # which hdul and ext to use
                vspace = 0.3
                fname = os.path.basename(str(master_im.filename()))
                t1 = TextPlot(columns=1, v_space=vspace)
                col1 = (
                    str(master_im["PRIMARY"].header.get("INSTRUME")),
                    "EXTNAME: " + str(master_im[ext].header.get("EXTNAME", "N/A")),
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
                panel.assign_plot(t1, px, py, xext=2)

                px = px + 2
                t2 = TextPlot(columns=1, v_space=vspace, xext=1)
                self.metadata = CriresSetupInfo.master_dark(master_im)
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
        fname = os.path.basename(str(self.hdus[0]["master_im"].filename()))
        master_im = self.hdus[0]["master_im"]


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
        self.metadata = CriresSetupInfo.master_dark(master_im)
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

            #plt = p.retrieve(i, 1)
            #plt.interp = "nearest"
            
            hdul = self.hdus[0][self.image_category]
            setup = str(hdul["PRIMARY"].header.get("HIERARCH ESO INS MODE", "N/A"))

            addme = {
                "report_name": "crires_master_dark_multi",
                "report_description": "CRIRES dark multi panel",
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

rep = CriresMasterDarkReport()
