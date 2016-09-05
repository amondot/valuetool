"""
/***************************************************************************
         Value Tool       - A QGIS plugin to get values at the mouse pointer
                             -------------------
    begin                : 2008-08-26
    copyright            : (C) 2008 by G. Picard
                           (C) 2016 by CS Systemes d'information (CS SI)
    contributors         : G. Picard
                           A. Mondot
    email                :
                           alexia.mondot@c-s.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""

import logging
logger = logging.getLogger('ValueWidget')
logger.setLevel(logging.DEBUG)

# change the level back to logging.WARNING(the default) before releasing
logging.basicConfig(level=logging.DEBUG)

from PyQt4 import QtCore, QtGui

# Import the PyQt and QGIS libraries
from PyQt4.QtCore import Qt, QObject, SIGNAL, QSize, QSettings
from PyQt4.QtGui import QApplication, QWidget, QTableWidgetItem, QBrush, QPen, QToolButton, QActionGroup, QMenu, QAction

# import GDAL and QGIS libraries
from qgis.core import (QGis,
                       QgsMapLayerRegistry,
                       QgsApplication,
                       QgsPoint,
                       QgsRasterBlock,
                       QgsRaster,
                       QgsMapLayer,
                       QgsRasterDataProvider,
                       QgsCoordinateReferenceSystem,
                       QgsCoordinateTransform,
                       QgsRectangle,
                       QgsCsException,
                       QgsRasterBandStats
                       )

from osgeo import osr, gdal
import gdalconst

from ui_valuewidgetbase import Ui_ValueWidgetBase as Ui_Widget

hasqwt = True
try:
    from PyQt4.Qwt5 import QwtPlot, QwtPlotCurve, QwtScaleDiv, QwtSymbol, QwtText
except:
    hasqwt = False

# test if matplotlib >= 1.0
hasmpl = True
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
except:
    hasmpl = False
if hasmpl:
    if int(matplotlib.__version__[0]) < 1:
        hasmpl = False

debug = 0


class ValueWidget(QWidget, Ui_Widget):
    def __init__(self, iface):

        self.hasqwt = hasqwt
        self.hasmpl = hasmpl
        self.layerMap = dict()
        self.statsChecked = False
        self.ymin = 0
        self.ymax = 250
        self.isActive = False

        # Statistics (>=1.9)
        self.statsSampleSize = 2500000
        self.stats = {}  # stats per layer

        self.layersSelected = []
        self.layerBands = dict()

        self.iface = iface
        self.canvas = self.iface.mapCanvas()
        self.legend = self.iface.legendInterface()
        self.logger = logging.getLogger('.'.join((__name__,
                                                  self.__class__.__name__)))


        # defines the size of the neighborhood to work on
        self.neighboroodSize = 1
        # defines if we display the neighborhood (if we calculate it)
        self.neighborood = False

        QWidget.__init__(self)
        self.setupUi(self)
        self.tabWidget.setEnabled(False)
        self.cbxClick.setChecked(QSettings().value('plugins/valuetool/mouseClick', False, type=bool))

        # self.setupUi_plot()
        # don't setup plot until Plot tab is clicked - workaround for bug #7450
        # qgis will still crash in some cases, but at least the tool can be used in Table mode
        self.qwtPlot = None
        self.mplPlot = None
        self.mplLine = None

        QObject.connect(self.plotSelector, SIGNAL("currentIndexChanged ( int )"), self.changePlot)
        QObject.connect(self.tabWidget, SIGNAL("currentChanged ( int )"), self.tabWidgetChanged)
        QObject.connect(self.cbxLayers, SIGNAL("currentIndexChanged ( int )"), self.updateLayers)
        QObject.connect(self.cbxBands, SIGNAL("currentIndexChanged ( int )"), self.updateLayers)
        QObject.connect(self.tableWidget2, SIGNAL("cellChanged ( int , int )"), self.layerSelected)

    def setupUi_plot(self):

        # plot
        self.plotSelector.setVisible(False)
        self.cbxStats.setVisible(False)
        # stats by default because estimated are fast
        self.cbxStats.setChecked(True)
        self.plotSelector.addItem('Qwt')
        self.plotSelector.addItem('mpl')

        # Page 2 - qwt
        if self.hasqwt:
            self.qwtPlot = QwtPlot(self.stackedWidget)
            self.qwtPlot.setAutoFillBackground(False)
            self.qwtPlot.setObjectName("qwtPlot")
            self.curve = QwtPlotCurve()
            self.curve.setSymbol(
                QwtSymbol(QwtSymbol.Ellipse,
                          QBrush(Qt.white),
                          QPen(Qt.red, 2),
                          QSize(9, 9)))
            self.curve.attach(self.qwtPlot)
        else:
            self.qwtPlot = QtGui.QLabel("Need Qwt >= 5.0 or matplotlib >= 1.0 !")

        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.qwtPlot.sizePolicy().hasHeightForWidth())
        self.qwtPlot.setSizePolicy(sizePolicy)
        self.qwtPlot.updateGeometry()
        self.stackedWidget.addWidget(self.qwtPlot)

        # Page 3 - matplotlib
        self.mplLine = None  # make sure to invalidate when layers change
        if self.hasmpl:
            # mpl stuff
            # should make figure light gray
            self.mplBackground = None  # http://www.scipy.org/Cookbook/Matplotlib/Animations
            self.mplFig = plt.Figure(facecolor='w', edgecolor='w')
            self.mplFig.subplots_adjust(left=0.1, right=0.975, bottom=0.13, top=0.95)
            self.mplPlt = self.mplFig.add_subplot(111)
            self.mplPlt.tick_params(axis='both', which='major', labelsize=12)
            self.mplPlt.tick_params(axis='both', which='minor', labelsize=10)
            # qt stuff
            self.pltCanvas = FigureCanvasQTAgg(self.mplFig)
            self.pltCanvas.setParent(self.stackedWidget)
            self.pltCanvas.setAutoFillBackground(False)
            self.pltCanvas.setObjectName("mplPlot")
            self.mplPlot = self.pltCanvas
        else:
            self.mplPlot = QtGui.QLabel("Need Qwt >= 5.0 or matplotlib >= 1.0 !")

        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mplPlot.sizePolicy().hasHeightForWidth())
        self.mplPlot.setSizePolicy(sizePolicy)
        self.mplPlot.updateGeometry()
        self.stackedWidget.addWidget(self.mplPlot)

        if (self.hasqwt and self.hasmpl):
            self.plotSelector.setEnabled(True)
            self.plotSelector.setVisible(True)
            self.plotSelector.setCurrentIndex(0);
        else:
            if self.hasqwt:
                self.plotSelector.setCurrentIndex(0);
            else:
                self.plotSelector.setCurrentIndex(1);
        self.changePlot()

    def keyPressEvent(self, e):
        if (e.modifiers() == Qt.ControlModifier or e.modifiers() == Qt.MetaModifier) and e.key() == Qt.Key_C:
            items = ''
            for rec in range(self.tableWidget.rowCount()):
                items += '"' + self.tableWidget.item(rec, 0).text() + '",' + self.tableWidget.item(rec, 1).text() + "\n"
            if not items == '':
                clipboard = QApplication.clipboard()
                clipboard.setText(items)
        else:
            QWidget.keyPressEvent(self, e)

    def changePlot(self):
        if (self.plotSelector.currentText() == 'mpl'):
            self.stackedWidget.setCurrentIndex(1)
        else:
            self.stackedWidget.setCurrentIndex(0)

    def changeActive(self, active, gui=True):
        self.isActive = active

        if (active):
            self.cbxEnable.setCheckState(Qt.Checked)
            QObject.connect(self.canvas, SIGNAL("layersChanged ()"), self.invalidatePlot)
            if not self.cbxClick.isChecked():
                QObject.connect(self.canvas, SIGNAL("xyCoordinates(const QgsPoint &)"), self.printValue)
                #connection computing neighborhood
                QObject.connect(self.cbxVoisinage, SIGNAL("stateChanged(int)"), self.changeNeighborActive)
                QObject.connect(self.spinBoxVoisinage, SIGNAL("valueChanged(int)"), self.defNeighborood)
        else:
            self.cbxEnable.setCheckState(Qt.Unchecked)
            QObject.disconnect(self.canvas, SIGNAL("layersChanged ()"), self.invalidatePlot)
            QObject.disconnect(self.canvas, SIGNAL("xyCoordinates(const QgsPoint &)"), self.printValue)

        if gui:
            self.tabWidget.setEnabled(active)
            if active:
                self.labelStatus.setText(self.tr("Value tool is enabled"))
                if self.tabWidget.currentIndex() == 2:
                    self.updateLayers()
            else:
                self.labelStatus.setText(self.tr(""))
                # use this to clear plot when deactivated
                # self.values=[]
                # self.showValues()

    def changeNeighborActive(self, state):
        """
        Updates self.neighborood according to the user choice.
        Allows or not to run the calculation of the average

        Keyword arguments:
            state -- state of the checkbox

        """
        if (state == Qt.Checked):
            self.neighborood = True
        else:
            self.neighborood = False

    def activeRasterLayers(self, index=None):
        layers = []
        allLayers = []

        if not index:
            index = self.cbxLayers.currentIndex()
        if index == 0:
            allLayers = self.canvas.layers()
        elif index == 1:
            allLayers = self.legend.layers()
        elif index == 2:
            for layer in self.legend.layers():
                if layer.id() in self.layersSelected:
                    allLayers.append(layer)

        for layer in allLayers:
            if layer != None and layer.isValid() and \
                            layer.type() == QgsMapLayer.RasterLayer and \
                    layer.dataProvider() and \
                    (layer.dataProvider().capabilities() & QgsRasterDataProvider.IdentifyValue):
                layers.append(layer)

        return layers

    def activeBandsForRaster(self, layer):
        activeBands = []

        if self.cbxBands.currentIndex() == 1 and layer.renderer():
            activeBands = layer.renderer().usesBands()
        elif self.cbxBands.currentIndex() == 2:
            if layer.bandCount() == 1:
                activeBands = [1]
            else:
                activeBands = self.layerBands[layer.id()] if (layer.id() in self.layerBands) else []
        else:
            activeBands = range(1, layer.bandCount() + 1)

        return activeBands

    def printValue(self, position):
        # coordinate pixel/line to display
        coordx = "0"
        coordy = "0"

        if debug > 0:
            print(position)

        if not position:
            return
        if self.tabWidget.currentIndex() == 2:
            return

        if debug > 0:
            print("%d active rasters, %d canvas layers" % (len(self.activeRasterLayers()), self.canvas.layerCount()))
        layers = self.activeRasterLayers()
        if len(layers) == 0:
            if self.canvas.layerCount() > 0:
                self.labelStatus.setText(self.tr("No valid layers to display - change layers in options"))
            else:
                self.labelStatus.setText(self.tr("No valid layers to display"))
            self.values = []
            self.showValues()
            return

        self.labelStatus.setText(self.tr('Coordinate:') + ' (%f, %f)' % (position.x(), position.y()))

        needextremum = (self.tabWidget.currentIndex() == 1)  # if plot is shown

        # count the number of requires rows and remember the raster layers
        nrow = 0
        rasterlayers = []
        layersWOStatistics = []

        for layer in layers:

            nrow += layer.bandCount()
            rasterlayers.append(layer)

            # check statistics for each band
            if needextremum:
                for i in range(1, layer.bandCount() + 1):
                    has_stats = self.getStats(layer, i) is not None
                    if not layer.id() in self.layerMap and not has_stats \
                            and not layer in layersWOStatistics:
                        layersWOStatistics.append(layer)

        if layersWOStatistics and not self.statsChecked:
            self.calculateStatistics(layersWOStatistics)

        irow = 0
        self.values = []
        self.ymin = 1e38
        self.ymax = -1e38

        mapCanvasSrs = self.iface.mapCanvas().mapRenderer().destinationCrs()

        # TODO - calculate the min/max values only once, instead of every time!!!
        # keep them in a dict() with key=layer.id()


        # pull out wavelength if it exists in metadata
        # piece to pull out wavelength information if present in metadata
        rasterMeta = rasterlayers[0].metadata()
        self.wavelengths = {}
        self.wavelength_units = ''
        if ('wavelength' in rasterMeta):
            mdSplit = rasterMeta.split('</p>')
            for d in mdSplit:
                if ('Band_' in d and 'glossy' not in d and '=' in d):
                    variableName, valueWavelength = d.split('=')
                    bandNumber = int(variableName.split('_')[1])
                    self.wavelengths[bandNumber] = float(valueWavelength.split(' ')[-2].replace('(', ''))
                elif ('wavelength_units' in d):
                    variableName, v = d.split('=')
                    self.wavelength_units = v
                    ####
        for layer in rasterlayers:

            layername = unicode(layer.name())
            layerSrs = layer.crs()

            pos = position

            # if given no position, get dummy values
            if position is None:
                pos = QgsPoint(0, 0)
            # transform points if needed
            elif not mapCanvasSrs == layerSrs and self.iface.mapCanvas().hasCrsTransformEnabled():
                srsTransform = QgsCoordinateTransform(mapCanvasSrs, layerSrs)
                try:
                    pos = srsTransform.transform(position)
                except QgsCsException, err:
                    # ignore transformation errors
                    continue

            if True:  # for QGIS >= 1.9
                if not layer.dataProvider():
                    continue

                ident = None
                identNeighborhood = None
                if position is not None:
                    canvas = self.iface.mapCanvas()

                    # first test if point is within map layer extent
                    # maintain same behaviour as in 1.8 and print out of extent
                    if not layer.dataProvider().extent().contains(pos):
                        ident = dict()
                        for iband in range(1, layer.bandCount() + 1):
                            ident[iband] = str(self.tr('out of extent'))
                    # we can only use context if layer is not projected
                    elif canvas.hasCrsTransformEnabled() and layer.dataProvider().crs() != canvas.mapRenderer().destinationCrs():
                        ident = layer.dataProvider().identify(pos, QgsRaster.IdentifyFormatValue).results()
                        # if computing of neighborhood is checked
                        if self.neighborood:
                            identNeighborhood = self.calculateNeighborood(pos, layer)
                    else:
                        extent = canvas.extent()
                        width = round(extent.width() / canvas.mapUnitsPerPixel());
                        height = round(extent.height() / canvas.mapUnitsPerPixel());

                        extent = canvas.mapRenderer().mapToLayerCoordinates(layer, extent);
                        # if computing of neighborhood is checked
                        if self.neighborood:
                            identNeighborhood = self.calculateNeighborood(pos, layer, canvas.extent(), width, height)

                        ident = layer.dataProvider().identify(pos, QgsRaster.IdentifyFormatValue, canvas.extent(),
                                                              width, height).results()
                    if not len(ident) > 0:
                        continue

                # if given no position, set values to 0
                if position is None and ident is not None and ident.iterkeys() is not None:
                    for key in ident.iterkeys():
                        ident[key] = layer.dataProvider().noDataValue(key)

                # bands displayed depends on cbxBands (all / active / selected)
                activeBands = self.activeBandsForRaster(layer)


                # get the pixel line of the current position
                coordx, coordy = self.calculatePixelLine( layer, pos )
                # self.printLatlongInStatusBar()


                for indexBand, iband in enumerate(activeBands):  # loop over the active bands
                    average = "#"
                    layernamewithband = layername
                    if ident is not None and len(ident) > 1:
                        layernamewithband += ' ' + layer.bandName(iband)

                    if not ident or not ident.has_key(iband):  # should not happen
                        bandvalue = "?"
                    else:
                        bandvalue = ident[iband]
                        if bandvalue is None:
                            bandvalue = "no data"

                    # if we need to display neighborhood
                    if self.neighborood and identNeighborhood is not None:
                        average = QgsRasterBlock.printValue( identNeighborhood[indexBand-1] )

                    # add the pixel line to the end of "values "
                    self.values.append( (layernamewithband, str(bandvalue), coordx, coordy, average) )

                    if needextremum:
                        # estimated statistics
                        stats = self.getStats(layer, iband)
                        if stats:
                            self.ymin = min(self.ymin, stats.minimumValue)
                            self.ymax = max(self.ymax, stats.maximumValue)

        if len(self.values) == 0:
            self.labelStatus.setText(self.tr("No valid bands to display"))

        self.showValues()

    def showValues(self):
        if self.tabWidget.currentIndex() == 1:
            # TODO don't plot if there is no data to plot...
            self.plot()
        else:
            self.printInTable()

    def calculateStatistics(self, layersWOStatistics):

        self.invalidatePlot(False)

        self.statsChecked = True

        layerNames = []
        for layer in layersWOStatistics:
            if not layer.id() in self.layerMap:
                layerNames.append(layer.name())

        if (len(layerNames) != 0):
            if not self.cbxStats.isChecked():
                for layer in layersWOStatistics:
                    self.layerMap[layer.id()] = True
                return
        else:
            print('ERROR, no layers to get stats for')

        save_state = self.isActive
        self.changeActive(False, False)  # deactivate

        # calculate statistics
        for layer in layersWOStatistics:
            if not layer.id() in self.layerMap:
                self.layerMap[layer.id()] = True
                for i in range(1, layer.bandCount() + 1):
                    self.getStats(layer, i, True)

        if save_state:
            self.changeActive(True, False)  # activate if necessary

    def getStats(self, layer, bandNo, force=False):
        """
        Get cached statistics for layer and band or None if not calculated
        :param layer:
        :param bandNo:
        :param force:
        :return:
        """
        if self.stats.has_key(layer):
            if self.stats[layer].has_key(bandNo):
                return self.stats[layer][bandNo]
        else:
            self.stats[layer] = {}

        if force or layer.dataProvider().hasStatistics(bandNo, QgsRasterBandStats.Min | QgsRasterBandStats.Min,
                                                       QgsRectangle(), self.statsSampleSize):
            self.stats[layer][bandNo] = layer.dataProvider().bandStatistics(bandNo,
                                                                            QgsRasterBandStats.Min | QgsRasterBandStats.Min,
                                                                            QgsRectangle(), self.statsSampleSize)
            return self.stats[layer][bandNo]

        return None

    def printInTable(self):
        """
        print self.values in the table
        """
        # set table widget row count
        self.tableWidget.setRowCount(len(self.values))

        # current line
        irow = 0
        for row in self.values:
            layername, value, x, y, average = row
            average = str( average )

            # limit number of decimal places if requested
            if self.cbxDigits.isChecked():
                try:
                    value = str("{0:." + str(self.spinDigits.value()) + "f}").format(float(value))
                except ValueError:
                    pass

            if (self.tableWidget.item(irow, 0) == None):
                # create the item
                self.tableWidget.setItem(irow, 0, QTableWidgetItem())
                self.tableWidget.setItem(irow, 1, QTableWidgetItem())
                self.tableWidget.setItem(irow, 2, QTableWidgetItem())
                self.tableWidget.setItem(irow, 3, QTableWidgetItem())
                #if self.neighborood :
                self.tableWidget.setItem(irow, 4, QTableWidgetItem())

            self.tableWidget.item(irow, 0).setText(layername)
            self.tableWidget.item(irow, 1).setText(value)
            self.tableWidget.item(irow, 2).setText( str( x ) )
            self.tableWidget.item(irow, 3).setText( str( y ) )
            if self.neighborood :
                self.tableWidget.item(irow, 4).setText(average)
            irow += 1

    def plot(self):
        numvalues = []
        wlvalues = []
        bandsUsed = []

        if (self.hasqwt or self.hasmpl):
            for row in self.values:
                layername, value = row
                bandsUsed.append(int(layername.split(' ')[-1]))
                try:
                    numvalues.append(float(value))
                except:
                    numvalues.append(0)
            if (len(self.wavelengths) != 0):
                for i in bandsUsed:
                    wlvalues.append(self.wavelengths[i])
            else:
                wlvalues = range(1, len(numvalues) + 1)

        ymin = self.ymin
        ymax = self.ymax
        xmin = float(min(wlvalues))
        xmax = float(max(wlvalues))
        if self.leYMin.text() != '' and self.leYMax.text() != '':
            ymin = float(self.leYMin.text())
            ymax = float(self.leYMax.text())

        if (self.hasqwt and (self.plotSelector.currentText() == 'Qwt')):

            self.qwtPlot.setAxisMaxMinor(QwtPlot.xBottom, 0)
            # self.qwtPlot.setAxisMaxMajor(QwtPlot.xBottom,0)
            self.qwtPlot.setAxisScale(QwtPlot.xBottom, xmin, xmax)
            # self.qwtPlot.setAxisScale(QwtPlot.yLeft,self.ymin,self.ymax)
            self.qwtPlot.setAxisScale(QwtPlot.yLeft, ymin, ymax)

            self.curve.setData(wlvalues, numvalues)
            if (self.wavelength_units.lower() == 'micrometers' or self.wavelength_units.lower() == 'microns'):
                xlabel = QwtText('Wavelength (um)')
            elif (self.wavelength_units.lower() == 'nanometers'):
                xlabel = QwtText('Wavelength (nm)')
            else:
                xlabel = QwtText('Bands')
            self.qwtPlot.setAxisTitle(QwtPlot.xBottom, xlabel)
            self.qwtPlot.replot()
            self.qwtPlot.setVisible(len(numvalues) > 0)

        elif (self.hasmpl and (self.plotSelector.currentText() == 'mpl')):

            self.mplPlt.clear()
            self.mplPlt.plot(wlvalues, numvalues, marker='o', color='k', mfc='b', mec='b')
            # self.mplPlt.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            self.mplPlt.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            self.mplPlt.set_xlim((xmin, xmax))
            self.mplPlt.set_ylim((ymin, ymax))

            if (self.wavelength_units.lower() == 'micrometers' or self.wavelength_units.lower() == 'microns'):
                xlabel = 'Wavelength ($\mu$m)'
            elif (self.wavelength_units.lower() == 'nanometers'):
                xlabel = 'Wavelength (nm)'
            else:
                xlabel = 'Bands'
            self.mplPlt.set_xlabel(xlabel)
            self.mplFig.canvas.draw()

    def invalidatePlot(self, replot=True):
        if self.tabWidget.currentIndex() == 2:
            self.updateLayers()
        if not self.isActive:
            return
        self.statsChecked = False
        if self.mplLine is not None:
            del self.mplLine
            self.mplLine = None
        # update empty plot
        if replot and self.tabWidget.currentIndex() == 1:
            # self.values=[]
            self.printValue(None)

    def resizeEvent(self, event):
        self.invalidatePlot()

    def tabWidgetChanged(self):
        if self.tabWidget.currentIndex() == 1 and not self.qwtPlot:
            self.setupUi_plot()
        if self.tabWidget.currentIndex() == 2:
            self.updateLayers()

    def updateLayers(self):
        """
        update active layers in table
        :return:
        """
        if self.tabWidget.currentIndex() != 2:
            return

        if self.cbxLayers.currentIndex() == 0:
            layers = self.activeRasterLayers(0)
        else:
            layers = self.activeRasterLayers(1)

        self.tableWidget2.blockSignals(True)
        self.tableWidget2.clearContents()
        self.tableWidget2.setRowCount(len(layers))
        self.tableWidget2.horizontalHeader().resizeSection(0, 20)
        self.tableWidget2.horizontalHeader().resizeSection(2, 20)

        j = 0
        for layer in layers:

            item = QTableWidgetItem()
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            if self.cbxLayers.currentIndex() != 2:
                item.setFlags(item.flags() & ~ Qt.ItemIsEnabled)
                item.setCheckState(Qt.Checked)
            else:
                if layer.id() in self.layersSelected:
                    item.setCheckState(Qt.Checked)
                else:
                    item.setCheckState(Qt.Unchecked)
            self.tableWidget2.setItem(j, 0, item)
            item = QTableWidgetItem(layer.name())
            item.setData(Qt.UserRole, layer.id())
            self.tableWidget2.setItem(j, 1, item)
            activeBands = self.activeBandsForRaster(layer)
            button = QToolButton()
            button.setText("#")  # TODO add edit? icon
            button.setPopupMode(QToolButton.InstantPopup)
            group = QActionGroup(button)
            group.setExclusive(False)
            QObject.connect(group, SIGNAL("triggered(QAction*)"), self.bandSelected)
            if self.cbxBands.currentIndex() == 2 and layer.bandCount() > 1:
                menu = QMenu()
                menu.installEventFilter(self)

                for iband in range(1, layer.bandCount() + 1):
                    action = QAction(str(layer.bandName(iband)), group)
                    action.setData([layer.id(), iband, j, False])
                    action.setCheckable(True)
                    action.setChecked(iband in activeBands)
                    menu.addAction(action)
                if layer.bandCount() > 1:
                    action = QAction(str(self.tr("All")), group)
                    action.setData([layer.id(), -1, j, True])
                    action.setCheckable(False)
                    menu.addAction(action)
                    action = QAction(str(self.tr("None")), group)
                    action.setData([layer.id(), -1, j, False])
                    action.setCheckable(False)
                    menu.addAction(action)

                button.setMenu(menu)
            else:
                button.setEnabled(False)
            self.tableWidget2.setCellWidget(j, 2, button)
            item = QTableWidgetItem(str(activeBands))
            item.setToolTip(str(activeBands))
            self.tableWidget2.setItem(j, 3, item)
            j = j + 1

        self.tableWidget2.blockSignals(False)

    def layerSelected(self, row, column):
        """
        slot for when active layer selection has changed
        :param row:
        :param column:
        :return:
        """
        if column != 0:
            return

        self.layersSelected = []
        for i in range(0, self.tableWidget2.rowCount()):
            item = self.tableWidget2.item(i, 0)
            layerID = self.tableWidget2.item(i, 1).data(Qt.UserRole)
            if item and item.checkState() == Qt.Checked:
                self.layersSelected.append(layerID)
            elif layerID in self.layersSelected:
                self.layersSelected.remove(layerID)

    def bandSelected(self, action):
        """
        slot for when active band selection has changed
        :param action:
        :return:
        """
        layerID = action.data()[0]
        layerBand = action.data()[1]
        j = action.data()[2]
        toggleAll = action.data()[3]
        activeBands = self.layerBands[layerID] if (layerID in self.layerBands) else []

        # special actions All/None
        if layerBand == -1:
            for layer in self.legend.layers():
                if layer.id() == layerID:
                    if toggleAll:
                        activeBands = range(1, layer.bandCount() + 1)
                    else:
                        activeBands = []
                    # toggle all band# actions
                    group = action.parent()
                    if group and not isinstance(group, QtGui.QActionGroup):
                        group = None
                    if group:
                        group.blockSignals(True)
                        for a in group.actions():
                            if a.isCheckable():
                                a.setChecked(toggleAll)
                        group.blockSignals(False)

        # any Band# action
        else:
            if action.isChecked():
                activeBands.append(layerBand)
            else:
                if layerBand in activeBands:
                    activeBands.remove(layerBand)
            activeBands.sort()

        self.layerBands[layerID] = activeBands

        # update UI
        item = QTableWidgetItem(str(activeBands))
        item.setToolTip(str(activeBands))
        self.tableWidget2.setItem(j, 3, item)

    def eventFilter(self, obj, event):
        """
        event filter for band selection menu, do not close after toggling each band

        :param obj:
        :param event:
        :return:
        """
        if event.type() in [QtCore.QEvent.MouseButtonRelease]:
            if isinstance(obj, QtGui.QMenu):
                if obj.activeAction():
                    if not obj.activeAction().menu():  # if the selected action does not have a submenu
                        # eat the event, but trigger the function
                        obj.activeAction().trigger()
                        return True
        return super(ValueWidget, self).eventFilter(obj, event)

    def shouldPrintValues(self):
        return self.isVisible() and not self.visibleRegion().isEmpty() \
               and self.isActive and self.tabWidget.currentIndex() != 2

    def toolMoved(self, position):
        if self.shouldPrintValues() and not self.cbxClick.isChecked():
            self.printValue(self.canvas.getCoordinateTransform().toMapCoordinates(position))

    def toolPressed(self, position):
        if self.shouldPrintValues() and self.cbxClick.isChecked():
            self.printValue(self.canvas.getCoordinateTransform().toMapCoordinates(position))

    def defNeighborood(self, size):
        """
        Updates self.neighboroodSize
            
        Keyword arguments:
            size -- size entered by the user

        """
        self.neighboroodSize = size
        logger.debug(self.neighboroodSize)

    def calculatePixelLine(self, layer, pos):
        """
        Computes the pixel line of the mouse position

        Keyword arguments:
            layer  -- The current layer the caller function is working on
            pos    -- The position in a given coordinate system

        Returns : nothing
        """

        self.posx = pos.x()
        self.posy = pos.y()
        self.posQgsPoint = pos

        # using GDAL
        try:
            dataset = gdal.Open(str(layer.source()), gdalconst.GA_ReadOnly)
        except RuntimeError:
            message = "Failed to open " + layer.source() + ". You may have to many data opened in QGIS."
            logger.error(message)
        else:
            spatialReference = osr.SpatialReference()
            spatialReference.ImportFromWkt(dataset.GetProjectionRef())

            # if the layer is not georeferenced
            if len(str(spatialReference)) == 0:
                # never go in this loop because qgis always georeference a file at opening
                coordx = pos.x()
                coordy = pos.y()
            else:
                # getting the extent and the spacing
                geotransform = dataset.GetGeoTransform()
                if not geotransform is None:
                    origineX = geotransform[0]
                    origineY = geotransform[3]
                    spacingX = geotransform[1]
                    spacingY = geotransform[5]

                    coordx = (pos.x() - origineX) / spacingX
                    coordy = (pos.y() - origineY) / spacingY
                    if coordy < 0:
                        coordy = - coordy

            if not layer.dataProvider().extent().contains(pos):
                #        if coordx < 0 or coordx > layer.width() or coordy < 0 or coordy > layer.height() :
                coordx = "Out of image"
                coordy = "Out of image"
            else:
                #            coordx = QgsRasterBlock.printValue( coordx )#QString(pos.x())
                coordx = str(int(coordx))  # QString(pos.x())
                #            coordy = QgsRasterBlock.printValue( coordy )#QString(pos.y())
                coordy = str(int(coordy))  # QString(pos.y())
            return coordx, coordy

    def changeNeighborActive(self, state):
        """
        Updates self.neighborood according to the user choice.
        Allows or not to run the calculation of the average

        Keyword arguments:
            state -- state of the checkbox

        """
        if (state == Qt.Checked):
            self.neighborood = True
        else:
            self.neighborood = False

    def defNeighborood(self, size):
        """
        Updates self.neighboroodSize

        Keyword arguments:
            size -- size entered by the user

        """
        self.neighboroodSize = size
        logger.debug(self.neighboroodSize)

    def calculateNeighborood(self, pos, layer, extent=None, width=0, height=0):
        """
        Computes the average in a given neighborhood of the mouse position

        Keyword arguments:
            pos    -- The position in a given coordinate system
            layer  -- The current layer the caller function is working on
            extent -- Current extent (default None)
            witdh  -- Computed width with pixel size (default 0)
            height -- Computed height with pixel size (default 0)

        Returns : nothing
        """
        logger.debug("calculate neighborood")
        # Values of the neighborhood for all layers

        identNeighborhood = []
        for iband in range(0, layer.bandCount()):
            identNeighborhood.append(0)

        #        logger.debug( identNeighborhood )

        nbElement = 0

        for i in range(0 - self.neighboroodSize, self.neighboroodSize):
            for j in range(0 - self.neighboroodSize, self.neighboroodSize):
                # going throw neighbors

                #                logger.debug( "pos x : " + str( pos.x() ) + ", pos y : " + str( pos.y() ) )

                point = QgsPoint(pos.x() + i, pos.y() + j)

                #                logger.debug( "pixel : i :" + str( i ) + " j : " + str( j ) )

                #                logger.debug( "point x : " + str( point.x() ) + ", point y : " + str( point.y() ) )

                if extent == None:
                    ident = layer.dataProvider().identify(point, QgsRaster.IdentifyFormatValue).results()
                else:
                    ident = layer.dataProvider().identify(point, QgsRaster.IdentifyFormatValue, extent,
                                                          width, height).results()
                nbElement += 1

                '''
                elif canvas.hasCrsTransformEnabled() and layer.dataProvider().crs() != canvas.mapRenderer().destinationCrs():
                    ident = layer.dataProvider().identify(pos, QgsRaster.IdentifyFormatValue).results()
                    # if computing of neighborhood is checked
                    if self.neighborood:
                        identNeighborhood = self.calculateNeighborood(pos, layer)
                else:
                    extent = canvas.extent()
                    width = round(extent.width() / canvas.mapUnitsPerPixel());
                    height = round(extent.height() / canvas.mapUnitsPerPixel());

                    extent = canvas.mapRenderer().mapToLayerCoordinates(layer, extent);
                    # if computing of neighborhood is checked
                    if self.neighborood:
                        identNeighborhood = self.calculateNeighborood(pos, layer, canvas.extent(), width, height)

                    ident = layer.dataProvider().identify(pos, QgsRaster.IdentifyFormatValue, canvas.extent(),
                                                          width, height).results()
                '''




                # add for each band the value of ident in identneighborhood
                for key in ident.iterkeys():
                    #                    logger.debug( "index in dictionary : " + str( key ) )
                    if ident is not None and len(ident) > 1:
                        #                        logger.debug( "dictionary is not empty" )
                        pass
                    if not ident or not ident.has_key(key):  # should not happen
                        #                        logger.debug( "dictionary is empty or had a missing key" )
                        pass
                    else:
                        if not ident[key]:
                            bandvalue="#"
                        # test if value is str (out of extent)
                        # this is kind of contrived, but trying to minimize changes
                        elif isinstance(ident[key], str):
                            #                            logger.debug( "value is a string" )
                            bandvalue = ident[key]
                            identNeighborhood[key - 1] += bandvalue
                            logger.debug("identNeighborhood {}, key: {}".format(identNeighborhood[key - 1], key))
                        else:
                            logger.debug("ident[key] {}".format(ident[key]))
                            doubleValue = float(ident[key])
                            # TODO
                            # if not layer.dataProvider().isNoDataValue(key, doubleValue):
                                #                                logger.debug( "value en i, j : " + str( doubleValue ) )
                                # add doublevalue to identvoisin
                            identNeighborhood[key - 1] += doubleValue
                            #                                logger.debug( "added value : " + str( identNeighborhood[key-1] ) )

        # averaging each band of identneighborhood
        for iband in range(0, layer.bandCount()):  # loop over the bands
            identNeighborhood[iband] /= nbElement

        #        stringLabel = QString( "Average of pixel (" + str( pos.x() ) + ", " + str( pos.y() )
        #                               + ") on the neighborhood of " + str( self.neighboroodSize ) + " is : " +  str( identNeighborhood ) )

        return identNeighborhood
