from . import CriresReportMixin
from .crires_utils import CriresSetupInfo, extensions
from adari_core.data_libs.echelle_flatfield import MasterEchelleFlatfieldReport
from adari_core.data_libs.master_dark_bias import MasterDarkBiasReport
from adari_core.utils.utils import fetch_kw_or_default
from adari_core.plots.panel import Panel
from adari_core.plots.text import TextPlot
import os
import numpy as np

class CriresWavelengthReport(CriresReportMixin, MasterEchelleFlatfieldReport, MasterDarkBiasReport):

   def __init__(self):
        super().__init__("crires_wavelength")
        

   def parse_sof(self):
        map = None
        raw = None
 
        for filename, catg in self.inputs:
            if catg == "CAL_WAVE_MAP" and map is None:
               map = filename
            if (catg == "WAVE_FPET" or catg == "WAVE_UNE") and raw is None:
               raw = filename
        # Build and return the file name list
        file_lists = []
        if map is not None and raw is not None:
            file_lists.append(
                {
                    "master_product": map,
                    "raw": raw,
                    "master_im": map,
                }
            )
        return file_lists
 
   def generate_first_panel(self):
        panels = {}
        self.panel_kwargs = dict(x=3, y=3, height_ratios=[1, 4, 4])

        for ext1 in extensions:
            new_panels = super().generate_panels(
                master_product_ext=ext1, 
                raw_ext=ext1,
                direction="x",
                ylabel="wavelength / nm",
                panel_kwargs = self.panel_kwargs,
                interpolation = "nearest",
            )
            for i, (panel, panel_descr) in enumerate(new_panels.items()):
               master = self.hdus[i]["master_product"]
               rawname = str(master["PRIMARY"].header.get("HIERARCH ESO PRO REC1 RAW1 NAME"))

               panel_descr["report_description"] = (
                       f"CRIRES - wavelength map"
                       f"{os.path.basename(panel_descr['master_product'])}, "
                       f"{panel_descr['master_product_ext']}"
               )
   
               # Text Plot
               px = 0
               py = 0
               vspace = 0.3
               fname = os.path.basename(str(master.filename()))
               procatg = str(master["PRIMARY"].header.get("HIERARCH ESO PRO CATG"))
               t1 = TextPlot(columns=1, v_space=vspace)
               col1 = (
                   str(master["PRIMARY"].header.get("INSTRUME")),
                   "EXTNAME: " + ext1,
                   "PRO CATG: " + procatg,
                   "FILE NAME: " + fname,
                   "RAW1 NAME: " + rawname,
               )
               t1.add_data(col1)
               panel.assign_plot(t1, px, py, xext=2)
               
               px = px + 2
               t2 = TextPlot(columns=1, v_space=vspace, xext=1)
               col2 = self.metadata
               t2.add_data(col2)
               panel.assign_plot(t2, px, py, xext=1)

               # Adjust plots
               p = panel.pop(2,2)

               panels = {**panels, **new_panels}
   
        return panels

   def generate_second_panel(self):

        panels = {}
        p = Panel(x=3, y=3, height_ratios=[1, 4, 4])

        master = self.hdus[0]["master_im"]
        rawname = str(master["PRIMARY"].header.get("HIERARCH ESO PRO REC1 RAW1 NAME"))

        # Text Plot
        px = 0
        py = 0
        vspace = 0.25
        fname = os.path.basename(str(master.filename()))
        procatg = str(master["PRIMARY"].header.get("HIERARCH ESO PRO CATG"))
        t1 = TextPlot(columns=1, v_space=vspace)
        col1 = (
                   str(master["PRIMARY"].header.get("INSTRUME")),
                   "EXTNAME: N/A",
                   "PRO CATG: " + procatg,
                   "FILE NAME: " + fname,
                   "RAW1 NAME: " + rawname,
               )

        t1.add_data(col1)
        p.assign_plot(t1, 0, 0, xext=2)

        t2 = TextPlot(columns=1, v_space=vspace, xext=1)
        col2 = self.metadata
        t2.add_data(col2)
        p.assign_plot(t2, 2, 0, xext=1)


        for i,ext1 in enumerate(extensions):

            full_plot, zoom_plot = super().image_plot(
                self.hdus[0]["master_im"][ext1],
                zoom_in=True,
                zoom_in_extent=200,
                img_kwargs={
                    "title": ext1,
                },
            )
            p.assign_plot(full_plot, i, 1, xext=1)
            p.assign_plot(zoom_plot, i, 2, xext=1)

            plt = p.retrieve(i, 1)
            plt.interp = "nearest"

            addme = {
                "report_name": "crires_wavelength_CAL_WAVE_MAP_multi",
                "report_description": "CRIRES - wavelength map - multi",
                "report_tags": [],
            }

            panels[p] = addme

        return panels

   def generate_panels(self, **kwargs):
        """Create both single and multiple extension panels."""

        self.metadata = CriresSetupInfo.wavelength(list(self.hdus[0].values())[0])

        panels = {
            **self.generate_first_panel(),
            **self.generate_second_panel(),
        }
        return panels

rep = CriresWavelengthReport()

