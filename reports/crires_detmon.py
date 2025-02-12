from .crires_utils import CriresSetupInfo, extensions
from adari_core.data_libs.master_dark_bias import MasterDarkBiasReport
from adari_core.utils.utils import fetch_kw_or_default
from adari_core.plots.points import LinePlot
from adari_core.plots.text import TextPlot
import os
import numpy as np
import astropy.io.fits as fits

class CriresDetmonReport(MasterDarkBiasReport):

   def __init__(self):
        super().__init__("crires_detmon")
        self.data_readers["master_im"] = self.master_im_data_reader

   def master_im_data_reader(self, filename):
        hdu = fits.open(filename, mode="readonly")
        for ext1 in extensions:
            if len(hdu[ext1].data.shape) > 2:
                hdu[ext1].data = hdu[ext1].data[0]
        return hdu

   def parse_sof(self):
        bpm = None
        coeffs = None

        for filename, catg in self.inputs:
            if catg == "CAL_DETLIN_BPM" and bpm is None:
                bpm = filename
            if catg == "CAL_DETLIN_COEFFS" and coeffs is None:
                coeffs = filename
        # Build and return the file name list
        file_lists = []
        if bpm is not None:
            file_lists.append(
                {
                    "master_im": bpm,
                }
            )
        if coeffs is not None:
            file_lists.append(
                {
                    "master_im": coeffs,
                }
            )

        return file_lists

   def generate_panels(self, **kwargs):
        panels = {}
        self.metadata = CriresSetupInfo.detmon(list(self.hdus[0].values())[0])

        for ext1 in extensions:
            new_panels = super().generate_raw_cuts_panels(
               master_im_ext=ext1,
               master_title="detlin",
               master_im_clipping=None,
               master_im_n_clipping=None,
               master_im_zoom_clipping=None,
               master_im_zoom_n_clipping=None,
               interpolation="nearest",
            )

            for i, (panel, panel_descr) in enumerate(new_panels.items()):
               panel_descr["report_description"] = (
                       f"CRIRES - "
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
               procatg = str(master_im["PRIMARY"].header.get("HIERARCH ESO PRO CATG"))
               t1 = TextPlot(columns=1, v_space=vspace)
               col1 = (
                   str(master_im["PRIMARY"].header.get("INSTRUME")),
                   "EXTNAME: " + ext1,
                   "PRO CATG: " + procatg,
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
               col2 = self.metadata
               t2.add_data(col2)
               panel.assign_plot(t2, px, py, xext=1)
   
               # Histogram - replace default plot
               if procatg == "CAL_DETLIN_BPM":
                   bins = np.arange(0,18,1)
               elif procatg == "CAL_DETLIN_COEFFS":
                   bins = np.arange(0.98,1.02,0.0025)
               h,b = np.histogram(master_im[ext1].data, bins=bins)
               hist = LinePlot(title="Histogram", x_label="y")
               
               hist.x_label = fetch_kw_or_default(master_im["PRIMARY"], "BUNIT", default="ADU")
               hist.y_label = "Frequency"
               hist.y_scale = "log"
               hist.add_data(
                   (    
                       np.reshape([b[:-1],b[1:]],2*b[:-1].size,order='F'),
                       np.reshape([h,h],2*b[:-1].size,order='F'),
                   ),
                   color="red",
                   label=procatg,
               )
               panel.assign_plot(hist, px+1, py+1) 

               # Adjust panels
               # change of colorbar and y-axis ranges 
               # this approach relies on the internal layout specified by the master report
               # TO DO: use another mechanism to retrieve plots
               p = panel.retrieve(0,1)
               if procatg == "CAL_DETLIN_BPM":
                   p.set_vlim(0,16)
               if procatg == "CAL_DETLIN_COEFFS":
                   p.set_vlim(0.98,1.02)
               p = panel.retrieve(0,2)
               if procatg == "CAL_DETLIN_BPM":
                   p.set_vlim(0,16)
               if procatg == "CAL_DETLIN_COEFFS":
                   p.set_vlim(0.98,1.02)

               if procatg == "CAL_DETLIN_COEFFS":
                  p = panel.retrieve(1,1)
                  p.set_ylim(0.98,1.02)
                  p = panel.retrieve(1,2)
                  p.set_ylim(0.98,1.02)
                  p = panel.retrieve(2,1)
                  p.set_ylim(0.98,1.02)
                  p = panel.retrieve(2,2)
                  p.set_ylim(0.98,1.02)                
               panels = {**panels, **new_panels}
   
        return panels


rep = CriresDetmonReport()
