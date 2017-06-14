from PyQt5.QtWidgets import QSizePolicy, QApplication
import numpy as np
from math import factorial
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg \
    as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import cm
from matplotlib import animation as ani
from matplotlib.widgets import Slider, Button
import matplotlib

from ..core.ContactFilters import *
from ..core.Biochemistry import AccumulationMapIndex
from matplotlib import pyplot as plt
plt.style.use('ggplot')

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)


class MplPlotter(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class ContactPlotter(MplPlotter):
    """Plots a frame-score plot with lines."""

    def plot_contact_figure(self, contact, nsPerFrame):
        frames = np.arange(0, len(contact.scoreArray), 1) * nsPerFrame
        self.axes.plot(frames, contact.scoreArray)
        self.axes.set_xlabel("time [ns]")
        self.axes.set_ylabel("score")

    def plot_all_contacts_figure(self, contacts, smooth, nsPerFrame):
        values = []

        for frame in range(len(contacts[0].scoreArray)):
            current = 0
            for c in contacts:
                current += c.scoreArray[frame]
            values.append(current)
        frames = np.arange(0, len(values), 1) * nsPerFrame
        if smooth:
            val = self.savitzky_golay(np.array(values), smooth, 2)
            smaller = np.where(val < 0)
            val[smaller] = 0
            self.axes.plot(val)
        else:
            self.axes.plot(values)
        self.axes.set_xlabel("time [ns]")
        self.axes.set_ylabel("score")

    def plot_hbondNumber(self, contacts, smooth, nsPerFrame):
        values = []
        for c in contacts:
            c.hbondFramesScan()

        for frame in range(len(contacts[0].scoreArray)):
            current = 0
            for c in contacts:
                current += c.hbondFrames[frame]
            values.append(current)
        frames = np.arange(0, len(values), 1) * nsPerFrame
        if smooth:
            val = self.savitzky_golay(np.array(values), smooth, 2)
            smaller = np.where(val < 0)
            val[smaller] = 0
            self.axes.plot(frames, val)
        else:
            self.axes.plot(frames, values)
        self.axes.set_xlabel("time [ns]")
        self.axes.set_ylabel("hbond number")

    # from http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    def savitzky_golay(self, y, window_size, order, deriv=0, rate=1):

        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError, msg:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve( m[::-1], y, mode='valid')

    def saveFigure(self, path, outputFormat):
        self.fig.savefig(path + "." + outputFormat, format=outputFormat)


class ContactPlotParameters:
    """Parameter for the contact plotter."""
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

        valuesNp = np.array(values, dtype=float)
        self.axes.hist(valuesNp, bins=20)
        self.axes.set_ylabel("N")
        self.axes.set_xlabel(attribute + " bins")
        self.fig.subplots_adjust(bottom=0.2, top=0.95, left=0.15, right=0.85)

    def plotContactHist(self, currentContacts, attribute, threshold, nsPerFrame, xticksfontsize):
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

        valuesNp = np.array(values, dtype=float)
        titlesNp = np.array(titles, dtype=str)

        x = range(len(currentContacts))
        # h = self.axes.bar(x, valuesNp, color="red")
        h = self.axes.bar(x, valuesNp)
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
    """Simple plotter, used in the SASA view."""
    def plot(self, x, y):
        self.axes.plot(x, y)
        self.axes.set_xlabel("x")
        self.axes.set_ylabel("f(x)")
        # self.axes.xaxis.set_label_position('top')

    def saveFigure(self, path, outputFormat):
        self.fig.savefig(path + "." + outputFormat, format=outputFormat)

    def clearFigure(self):
        self.fig.clf()



class AnimateMapPlotter(MplPlotter):
    """Animated canvas with an 2d heatmap plot."""

    def plotMap(self, contacts, map1, map2, label1, label2, attribute, threshold, nsPerFrame):
        self.contacts = contacts
        self.map1 = map1
        self.map2 = map2
        self.label1 = label1
        self.label2 = label2
        self.attribute = attribute
        self.threshold = threshold
        self.nsPerFrame = nsPerFrame

        self.minmaxresids1 = []
        self.minmaxresids2 = []

        if not self.map1[AccumulationMapIndex.resid] or not self.map2[AccumulationMapIndex.resid]:
            return -1

        for cont in self.contacts:
            self.minmaxresids1.append(int(cont.key1[AccumulationMapIndex.resid]))
            self.minmaxresids2.append(int(cont.key2[AccumulationMapIndex.resid]))

        self.x = np.arange(np.min(self.minmaxresids1), np.max(self.minmaxresids1) + 1)
        self.y = np.arange(np.min(self.minmaxresids2), np.max(self.minmaxresids2) + 1)
        self.minx = np.min(self.minmaxresids1)
        self.miny = np.min(self.minmaxresids2)
        rng = np.arange(len(contacts[0].scoreArray))
        self.cax = None
        self.ttl = self.axes.set_title("Map")
        # TODO: do this automatically
        stridex = 5
        stridey = 5
        self.axes.set_xticks(np.arange(0, self.x.size, stridex))
        self.axes.set_xticklabels(np.arange(self.minx, self.x.size + self.minx, stridex))
        self.axes.set_xlabel(self.label1)
        self.axes.set_ylabel(self.label2)
        self.axes.set_yticks(np.arange(0, self.y.size, stridey))
        self.axes.set_yticklabels(np.arange(self.miny, self.y.size + self.miny, stridey))
        self.fig.tight_layout()
        self.cb = None
        animation = ani.FuncAnimation(self.fig, self.updateMap, rng, init_func=self.initFig, interval=200, blit=True, repeat=False)
        # animation = ani.FuncAnimation(self.fig, self.updateMap, 50, interval=200, blit=True)
        QApplication.processEvents()
        self.show()

        # writer = ani.writers['ffmpeg'](fps=30)
        # animation.save('demo.mp4',writer=writer,dpi=100)

    def initFig(self):
        data = np.zeros((len(self.y), len(self.x)))
        frame = 0
        print("frame: ", frame)
        attribute = "Mean Score"
        if attribute == "Mean Score":
            for c in self.contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid]) - self.minx
                r2 = int(c.key2[AccumulationMapIndex.resid]) - self.miny
                data[r2, r1] = c.scoreArray[frame]
        self.cax = self.axes.matshow(data, cmap=cm.Greys, label=self.attribute)
        self.ttl.set_text("Frame 0")
        if self.cb is None:
            self.cb = self.fig.colorbar(self.cax)
            self.cb.set_label(self.attribute)
            # mp_slider_ax  = self.fig.add_axes([1, 1, 0.65, 0.03])
            # sfreq = Slider(mp_slider_ax, 'Freq', 0.1, 30.0, valinit=10)
        self.draw()
        return self.cax,

    def updateMap(self, frame):
        print(frame)
        data = np.zeros((len(self.y), len(self.x)))
        attribute = "Mean Score"
        if attribute == "Mean Score":
            for c in self.contacts:
                r1 = int(c.key1[AccumulationMapIndex.resid]) - self.minx
                r2 = int(c.key2[AccumulationMapIndex.resid]) - self.miny
                data[r2, r1] = c.scoreArray[frame]
        self.cax.set_data(data)  # update the data
        self.ttl.set_text("Frame " + str(frame))

        # elif attribute == "Median Score":
        #     for c in contacts:
        #         r1 = int(c.key1[AccumulationMapIndex.resid])-minx
        #         r2 = int(c.key2[AccumulationMapIndex.resid])-miny
        #         data[r2, r1] = c.median_score()
        # elif attribute == "Mean Lifetime":
        #     for c in contacts:
        #         r1 = int(c.key1[AccumulationMapIndex.resid])-minx
        #         r2 = int(c.key2[AccumulationMapIndex.resid])-miny
        #         data[r2, r1] = c.mean_life_time(nsPerFrame, threshold)
        # elif attribute == "Median Lifetime":
        #     for c in contacts:
        #         r1 = int(c.key1[AccumulationMapIndex.resid])-minx
        #         r2 = int(c.key2[AccumulationMapIndex.resid])-miny
        #         data[r2, r1] = c.median_life_time(nsPerFrame, threshold)
        # elif attribute == "Hbond percentage":
        #     for c in contacts:
        #         r1 = int(c.key1[AccumulationMapIndex.resid])-minx
        #         r2 = int(c.key2[AccumulationMapIndex.resid])-miny
        #         data[r2, r1] = c.hbond_percentage()

        self.draw()

        return self.cax,

    def saveFigure(self, path, outputFormat):
        self.fig.savefig(path + "." + outputFormat, format=outputFormat)
