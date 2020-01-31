from __future__ import with_statement
from __future__ import print_function
import sys


try:
    import numpy
    import reflex
    from pipeline_product import PipelineProduct
    from reflex_plot_widgets import *
    import pipeline_display
    import_success = True

except ImportError:
    import_success = False
    print("Error importing modules astropy, wx, matplotlib, numpy")


def paragraph(text, width=None):
    """ wrap text string into paragraph
       text:  text to format, removes leading space and newlines
       width: if not None, wraps text, not recommended for tooltips as
              they are wrapped by wxWidgets by default
    """
    import textwrap
    if width is None:
        return textwrap.dedent(text).replace('\n', ' ').strip()
    else:
        return textwrap.fill(textwrap.dedent(text), width=width)


class DataPlotterManager(object):
    """
    This class must be added to the PipelineInteractiveApp with setPlotManager
    It must have following member functions which will be called by the app:
     - setInteractiveParameters(self)
     - readFitsData(self, fitsFiles):
     - addSubplots(self, figure):
     - plotProductsGraphics(self)
    Following members are optional:
     - setWindowHelp(self)
     - setWindowTitle(self)
     - setCurrentParameterHelper(self, helper)
    """

    # static members
    recipe_name = "rrrecipe"
    img_cat = "RRRECIPE_DOCATG_RESULT"

    def setInteractiveParameters(self):
        """
        This function specifies which are the parameters that should be presented
        in the window to be edited.  Note that the parameter has to also be in the
        in_sop port (otherwise it won't appear in the window). The descriptions are
        used to show a tooltip. They should match one to one with the parameter
        list.
        """
        return [
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="boolopt",
                                   group="group1", description="Desc2"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="intopt",
                                   group="group1", description="Desc2"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="floatopt",
                                   group="group1", description="Desc2"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="rangeopt",
                                   group="group1", description="Valid range is 0-10"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="enumopt",
                                   group="group1", description="Desc2"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="floatrangeopt",
                                   group="group1", description="Valid range is -5.5 - 5.5"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="stropt",
                                   group="group2", description="Desc1"),
            reflex.RecipeParameter(recipe=self.recipe_name, displayName="fileopt",
                                   group="group2", description="Desc1"),
        ]

    def readFitsData(self, fitsFiles):
        """
        This function should be used to read and organize the raw fits files
        produced by the recipes.
        It receives as input a list of reflex.FitsFiles
        """
        # organize the files into a dictionary, here we assume we only have 
        # one file per category if there are more, one must use a
        # dictionary of lists
        self.frames = dict()
        for f in fitsFiles:
            self.frames[f.category] = PipelineProduct(f)

        # we only have two states, we have data or we don't
        # define the plotting functions we want to use for each
        if self.img_cat in self.frames:
            self._add_subplots = self._add_subplots
            self._plot = self._data_plot
        else:
            self._add_subplots = self._add_nodata_subplots
            self._plot = self._nodata_plot

    def addSubplots(self, figure):
        """
        This function should be used to setup the subplots of the gui.  The the
        matplotlib documentation for a description of subplots.
        """
        self._add_subplots(figure)

    def plotProductsGraphics(self):
        """
        This function should be used to plot the data onto the subplots.
        """
        self._plot()

    def setWindowHelp(self):
        return 'Help for rrrecipe interactive window'

    def setWindowTitle(self):
        return 'rrrecipe interactive window'

    def _add_nodata_subplots(self, figure):
        self.txt_plot = figure.add_subplot(111)

    def _add_subplots(self, figure):
        self.img_plot = figure.add_subplot(211)
        self.spec_plot = figure.add_subplot(212)

    def _data_plot(self):
        # get the right category file from our dictionary
        p = self.frames[self.img_cat]
        p.readImage()
        # setup the image display
        imgdisp = pipeline_display.ImageDisplay()
        imgdisp.setLabels('X', 'Y')
        tooltip = paragraph("""\
        Bias (and dark) corrected arc lamp pinhole
        image. Blue points are the predicted positions of the arc lines
        in the reference arc line list using the initial physical
        model. Yellow points are the positions of those predicted
        lines that are successfully found and fitted by a two-dimensional
        Gaussian. Green points are the line positions which remain after
        thresholding the Gaussian fits for signal-to-noise and
        refining the physical model fit.
        """)

        imgdisp.setCmap('hot')
        imgdisp.display(self.img_plot, "Image title", tooltip, p.image)

        # plot a spectrum, we just us the sum of the rows here
        spec = pipeline_display.SpectrumDisplay()
        spec.setLabels('Lambda', 'Flux')
        flux = p.image.sum(axis=0)
        wave = numpy.linspace(1, 20, num=flux.size)
        spec.display(self.spec_plot, "title", "tooltip", wave, flux)

    def _nodata_plot(self):
        # could be moved to reflex library?
        self.txt_plot.set_axis_off()
        text_nodata = "Data not found. Input files should contain these" \
                       " types:\n%s" % self.img_cat
        self.txt_plot.text(0.1, 0.6, text_nodata, color='#11557c',
                      fontsize=18, ha='left', va='center', alpha=1.0)
        self.txt_plot.tooltip = 'No data found'

    def plotWidgets(self) :
        widgets = list()
        # Updating number
        self.clickablespec = InteractiveClickableSubplot(
            self.spec_plot, self.increaseFloatNumber)
        widgets.append(self.clickablespec)
        return widgets

    def increaseFloatNumber(self, point) :
        floatpoint = self.getCurrentParameterHelper('floatopt') + 1
        new_params = list()
        new_params.append(reflex.RecipeParameter(self.recipe_name,'floatopt',
                                                 value=str(floatpoint)))
        return new_params

    def setCurrentParameterHelper(self, helper) :
        self.getCurrentParameterHelper = helper

#This is the 'main' function
if __name__ == '__main__':
    from reflex_interactive_app import PipelineInteractiveApp

    # Create interactive application
    interactive_app = PipelineInteractiveApp(enable_init_sop=True)

    # get inputs from the command line
    interactive_app.parse_args()

    #Check if import failed or not
    if not import_success:
        interactive_app.setEnableGUI(False)

    #Open the interactive window if enabled
    if interactive_app.isGUIEnabled():
        #Get the specific functions for this window
        dataPlotManager = DataPlotterManager()

        interactive_app.setPlotManager(dataPlotManager)
        interactive_app.showGUI()
    else:
        interactive_app.set_continue_mode()

    #Print outputs. This is parsed by the Reflex python actor to
    #get the results. Do not remove
    interactive_app.print_outputs()
    sys.exit()
