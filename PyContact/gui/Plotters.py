import warnings

from PyQt5.QtWidgets import QSizePolicy
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
    as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import cm

# from ..settings import *
# from ..core.Biochemistry import *
from ..core.ContactFilters import *
from ..core.Biochemistry import AccumulationMapIndex


class MplPlotter(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class ContactPlotter(MplPlotter):

    def plot_contact_figure(self, contact):
        self.axes.plot(contact.scoreArray)
        self.axes.set_xlabel("frame")
        self.axes.set_ylabel("score")

    def plot_all_contacts_figure(self, contacts):
        values = []

        for frame in range(len(contacts[0].scoreArray)):
            current = 0
            for c in contacts:
                current += c.scoreArray[frame]
            values.append(current)
        self.axes.plot(values)
        self.axes.set_xlabel("frame")
        self.axes.set_ylabel("score")


class ContactPlotParameters():
    mean, median, lifetime, median_life_time, hbond_percentage = range(5)
    mapping = ["Mean Score", "Median Score", "Mean Lifetime", "Median Lifetime", "Hbond percentage"]


class HistPlotter(MplPlotter):
    """Simple canvas with an histogram plot."""

    def plotGeneralHist(self, currentContacts, attribute, threshold, nsPerFrame):
        values = []

        if attribute == "Mean Score":
            for c in currentContacts:
                values.append(c.mean_score())
        elif attribute == "Median Score":
            for c in currentContacts:
                values.append(c.median_score())
        elif attribute == "Mean Lifetime":
            for c in currentContacts:
                values.append(c.mean_life_time(nsPerFrame, threshold))
        elif attribute == "Median Lifetime":
            for c in currentContacts:
                values.append(c.median_life_time(nsPerFrame, threshold))

        valuesNp = np.array(values, dtype = float)
        self.axes.hist(valuesNp, bins=20)
        self.axes.set_ylabel("N")
        self.axes.set_xlabel(attribute + " bins")
        self.fig.subplots_adjust(bottom=0.2, top=0.95, left=0.15, right=0.85)

    def plotContactHist(self, currentContacts, attribute, threshold, nsPerFrame,xticksfontsize):
        values = []
        titles = []

        if attribute == "Mean Score":
            for c in currentContacts:
                values.append(c.mean_score())
                titles.append(c.title)
        elif attribute == "Median Score":
            for c in currentContacts:
                values.append(c.median_score())
                titles.append(c.title)
        elif attribute == "Mean Lifetime":
            for c in currentContacts:
                values.append(c.mean_life_time(nsPerFrame, threshold))
                titles.append(c.title)
        elif attribute == "Median Lifetime":
            for c in currentContacts:
                values.append(c.median_life_time(nsPerFrame, threshold))
                titles.append(c.title)
        elif attribute == "Hbond percentage":
            for c in currentContacts:
                values.append(c.hbond_percentage())
                titles.append(c.title)

        valuesNp = np.array(values, dtype = float)
        titlesNp = np.array(titles, dtype = str)

        x = range(len(currentContacts))
        h = self.axes.bar(x, valuesNp, color="red")
        xticks_pos = [0.7071 * patch.get_width() + patch.get_xy()[0] for patch in h]
        self.axes.set_xticklabels(titlesNp, ha='right', size=8, rotation=45)
        self.axes.set_xticks(xticks_pos)
        for tick in self.axes.xaxis.get_major_ticks():
            tick.label.set_fontsize(xticksfontsize)
        self.axes.set_ylabel(attribute)
        self.fig.subplots_adjust(bottom=0.2, top=0.95, left=0.1, right=0.9)

    def saveFigure(self, path, outputFormat):
        self.fig.savefig(path + "." + outputFormat, format=outputFormat)


class MapPlotter(MplPlotter):
    """Simple canvas with an 2d heatmap plot."""
    def plotMap(self, contacts, map1, map2, label1, label2, attribute, threshold, nsPerFrame):
        minmaxresids1 = []
        minmaxresids2 = []
        if not map1[AccumulationMapIndex.resid] or not map2[AccumulationMapIndex.resid]:
            return -1

        for cont in contacts:
            minmaxresids1.append(int(cont.key1[AccumulationMapIndex.resid]))
            minmaxresids2.append(int(cont.key2[AccumulationMapIndex.resid]))

        x = np.arange(np.min(minmaxresids1), np.max(minmaxresids1)+1)
        y = np.arange(np.min(minmaxresids2), np.max(minmaxresids2)+1)
        minx = np.min(minmaxresids1)
        miny = np.min(minmaxresids2)
        data = np.zeros((len(y), len(x)))

        if attribute == "Mean Score":
            for c in contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid])-minx
                r2 = int(c.key2[AccumulationMapIndex.resid])-miny
                data[r2, r1] = c.mean_score()
        elif attribute == "Median Score":
            for c in contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid])-minx
                r2 = int(c.key2[AccumulationMapIndex.resid])-miny
                data[r2, r1] = c.median_score()
        elif attribute == "Mean Lifetime":
            for c in contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid])-minx
                r2 = int(c.key2[AccumulationMapIndex.resid])-miny
                data[r2, r1] = c.mean_life_time(nsPerFrame, threshold)
        elif attribute == "Median Lifetime":
            for c in contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid])-minx
                r2 = int(c.key2[AccumulationMapIndex.resid])-miny
                data[r2, r1] = c.median_life_time(nsPerFrame, threshold)
        elif attribute == "Hbond percentage":
            for c in contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid])-minx
                r2 = int(c.key2[AccumulationMapIndex.resid])-miny
                data[r2, r1] = c.hbond_percentage()

        cax = self.axes.matshow(data, cmap=cm.Greys, label=attribute)

        # TODO: do this automatically
        stridex = 5
        stridey = 5
        self.axes.set_xticks(np.arange(0, x.size, stridex))
        self.axes.set_xticklabels(np.arange(minx, x.size+minx, stridex))
        # self.axes.set_title("Contact Map")
        self.axes.set_xlabel(label1)
        self.axes.set_ylabel(label2)

        self.axes.set_yticks(np.arange(0, y.size, stridey))
        self.axes.set_yticklabels(np.arange(miny, y.size+miny, stridey))
        cb = self.fig.colorbar(cax)
        cb.set_label(attribute)
        self.fig.tight_layout()

    def saveFigure(self, path, outputFormat):
        self.fig.savefig(path + "." + outputFormat, format=outputFormat)


class SimplePlotter(MplPlotter):
    def plot(self, x, y):
        self.axes.plot(x, y)
        self.axes.set_xlabel("x")
        self.axes.set_ylabel("f(x)")
        self.axes.xaxis.set_label_position('top')
