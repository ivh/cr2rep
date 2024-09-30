from adari_core.plots.text import TextPlot
from adari_core.data_libs.master_rawdisp import MasterRawdispReport
from astropy.io import fits
from .crires_utils import CriresSetupInfo
from adari_core.utils.utils import fetch_kw_or_default
from . import CriresReportMixin
import logging

import os

class CriresRawdispReport(CriresReportMixin, MasterRawdispReport):
	def __init__(self):
		super().__init__("crires_rawdisp")
		#rest here
		#channels
		self.raw_extensions = {
            "CHIP1.INT1": "CHIP1",
            "CHIP2.INT1": "CHIP2",
            "CHIP3.INT1": "CHIP3"
            }

		self.extensions = []
		self.tasks = {
            #
            "DARK": "dark",
            #
            "DETLIN_DARK": "detector_linearity",
            "DETLIN_LAMP": "detector_linearity",
            #
            "FLAT": "flat",
            #
            "WAVE_UNE": "wave_uranium_neon",
            #
            "WAVE_FPET": "wave_fabry_perot_etalon",
            #
            "WAVE_N20": "wave_gas_cell",
            "WAVE_SGC": "wave_gas_cell",
            "WAVE_OTHER": "wave_gas_cell",
            #
            "WAVE_SKY": "wave_sky",
            #
            "CAL_NODDING_JITTER": "standard_star",
            "CAL_NODDING_OTHER": "standard_star"
        }
        
		self.select_raw_files = {}

		self.task_scaling = {}

		self.setup_info = CriresSetupInfo

	def parse_sof(self):
		# we building multiple report sets, so we append multiple reports to file_lists
		# get a list of tags
		tags = list(self.tasks.keys())
		added = {}
		file_lists = []
		for filename, catg in self.inputs:
			if catg in tags:
				if filename is not None and catg not in added:
					file_lists.append({"filename": filename})
					added[catg] = self.tasks[catg]
					self.sof_tag.append(catg)        
		return file_lists

	def get_extensions(self):
		"""Find the data extensions required for each file.

		Description
		-----------
		After the SOF has been parsed, this method iterates over the different
		HDUS files to find which extension(s) contains the data.
		"""
		new_hdus_list = []
		new_sof_tag = []
		for i, filedict in enumerate(self.hdus):
			hdul = filedict["filename"]
			all_hdu_names = [hdu.name for hdu in hdul]
			for chip_name in self.raw_extensions.keys():
				if chip_name in all_hdu_names:
					new_hdus_list.append(filedict)
					new_sof_tag.append(self.sof_tag[i])
					self.extensions.append(chip_name)

		self.hdus = new_hdus_list
		self.sof_tag = new_sof_tag

	def generate_panels(self, **kwargs):
		panels = {}
		self.get_extensions()
		new_panels = super().generate_panels(ext=self.extensions, **kwargs)
		for i, (panel, panel_descr) in enumerate(new_panels.items()):
			# Alter the cut pos, or remove CutPlot(s) completely,
			# depending on task name
			try:
				task_name = panel_descr["task_name"]
			except KeyError as e:
				raise RuntimeError(
				    "A report has been created by "
				    "MasterRawdispReport that did "
				    "not come back with a task name "
				    "attached!"
				)
			panel_descr["report_name"] = "crires_rawdisp_{}_{}_{}".format(
				task_name,
				self.raw_extensions[self.extensions[i]],
				self.sof_tag[i].lower(),
			)
			panel_descr["report_description"] = (
				f"CRIRES rawdisp panel - "
				f"{panel_descr['task_name']}, "
				f"{panel_descr['tag']}, "
				f"{os.path.basename(panel_descr['filename'])}, "
				f"{panel_descr['ext']}"
			)
		panels = {**panels, **new_panels}

		return panels


rep = CriresRawdispReport()
