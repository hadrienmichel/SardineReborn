## TODO/IDEAS
# Add slider for zoom control on zoom panel
# Add option to merge similar sources/receivers at loading of geometry file
# Add posibilities for auto-picking
# Add robustness for reading/writing operations
# Debug Modelling with multiple datasets

## Imports for the inner functions
import sys
import os
import re
from copy import deepcopy
import numpy as np
import time
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
NavigationToolbar2QT.toolitems = [('Home', 'Reset original view', 'home', 'home'),
                                  (None, None, None, None), 
                                  ('Pan', 'Left button pans, Ri...xes aspect', 'move', 'pan'), 
                                  ('Zoom', 'Zoom to rectangle\nx/...xes aspect', 'zoom_to_rect', 'zoom'),
                                  (None, None, None, None), 
                                  ('Save', 'Save the figure', 'filesave', 'save_figure')]
from matplotlib import animation
from matplotlib.figure import Figure
from matplotlib.backend_bases import MouseButton
## Imports for the seismic data input
from obspy import read
from obspy.signal.trigger import aic_simple
## Imports for the data inversion
import pygimli as pg
from pygimli.physics import TravelTimeManager as TTMgr
from pygimli.physics.traveltime import ratools, drawFirstPicks
from pygimli.viewer import show as pgshow
## Imports for saving/reading states:
import pickle
## Imports for the GUI
from PyQt5.QtWidgets import (
    QApplication,  
    QWidget,
    QTabWidget,
    QVBoxLayout,
    QGridLayout,
    QAction,
    QMainWindow,
    QStatusBar,
    QMessageBox,
    QFileDialog,
    QInputDialog,
    QProgressBar,
    QComboBox,
    QPushButton,
    QSpinBox,
    QLabel,
    QGroupBox,
    QLineEdit)
from PyQt5.QtCore import Qt
from PyQt5 import QtGui

defaultStatus = "Idle."

def buildModel(sourceX, receiversX, times, nbLayers=2, orientation=1):
    '''BUILDMODEL is a function that computes a simple optimal model for 
    a given set of source, receivers and traveltimes.
    It computes the curvature of the TT and detects breaks in the hodochrones.
    It serves only as a first guess and should not be used as an optimal model!

    Inputs:
        - sourceX:
        - receiversX:
        - times:
        - nbLayers (default=2):
        - orientation (default=1):

    Outputs:
        - interceptTimes:
        - apparentVelocities:

    '''
    # Change the coordinates if requiered and normalize to sourceX = 0:
    receiversXSave = receiversX
    if orientation < 0:
        receiversX = np.flip(np.abs(receiversX - sourceX))
        times = np.flip(times)
    else:
        receiversX = receiversX - sourceX
    curvature = np.abs(np.gradient(np.gradient(times, receiversX), receiversX))
    # Find the highest curvatures in the set:
    maxCurvIndex = np.sort(np.argpartition(curvature, -(nbLayers-1))[-(nbLayers-1):])
    maxCurvIndex = np.concatenate(([0], maxCurvIndex, [len(times)-1]))
    v = np.zeros((nbLayers,))
    inter = np.zeros((nbLayers,))
    points = np.zeros((nbLayers+1,2)) # Points to store the intersections of interest
    for i in range(nbLayers):
        x = receiversX[maxCurvIndex[i] : maxCurvIndex[i+1]]
        y = times[maxCurvIndex[i] : maxCurvIndex[i+1]]
        if i == 0:
            x = x[:, np.newaxis]
            p = np.linalg.lstsq(x, y, rcond=None)
            v[i] = 1/p[0][0]
            inter[i] = 0
        else:
            p = np.polyfit(x, y, 1)
            v[i] = 1/p[0]
            inter[i] = p[1]
    
    points[0,0] = sourceX
    for i in range(nbLayers):
        if i < nbLayers-1:
            xTemp = (inter[i]-inter[i+1])/((1/v[i+1])-(1/v[i]))
            tTemp = inter[i] + xTemp*(1/v[i])
            points[i+1,0] = sourceX + orientation*xTemp
            points[i+1,1] = tTemp
        else:
            if orientation > 0:
                points[i+1,0] = max(receiversXSave)
            else:
                points[i+1,0] = min(receiversXSave)
            points[i+1,1] = inter[-1] + np.abs(sourceX - points[i+1,0])*(1/v[-1])

    return inter, v, points

def model1D(inter, v):
    h = np.ones((len(v)-1,))
    pos = np.zeros((len(v), 2))
    for i in np.arange(1, len(v)):
        iCr = np.arcsin(v[:i]/v[i])
        tInterp = inter[i]
        for j in np.arange(i-1):
            tInterp -= 2*h[j]*np.cos(iCr[j])/v[j]
            if j < i-1:
                pos[i, 0] += np.tan(iCr[j]) * h[j]
        h[i-1] = v[i-1] * tInterp / (2*np.cos(iCr[i-1]))
        pos[i, 1] = pos[i-1, 1] + h[i-1]
        pos[i, 0] = pos[i-1, 0] + np.tan(iCr[i-1]) * h[i-1]
    return h, v, pos

def modelWithSlope(interS, vS):
    '''Model from MOTA L. (1954) entiteld "Determination of dips and depths of geological layers by the seismic refraction method"
    '''
    nbLayers = np.shape(vS)[0]
    v0 = np.sum(vS[0,:])/2 # We take the mean velocity between the two possibilities 
    # Find alpha and v2
    ipA = np.arcsin(v0/vS[1,0])
    imA = np.arcsin(v0/vS[1,1])
    iCr = (ipA + imA)/2
    theta1 = (imA - ipA)/2
    v1 = v0/np.sin(iCr)
    zL = [v0*interS[1,0]/(2*np.cos(iCr))]
    zR = [v0*interS[1,1]/(2*np.cos(iCr))]
    hL = [zL[0]/np.cos(theta1)] # Depth of the interface below the source A                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    hR = [zR[0]/np.cos(theta1)] # Depth of the interface below the source B
    v = [v0, v1]
    if nbLayers > 2:
        alpha = np.arcsin(v0/vS[2,1]) - theta1
        beta = np.arcsin(v0/vS[2,0]) + theta1
        gamma = np.arcsin(v1/v0 * np.sin(alpha))
        delta = np.arcsin(v1/v0 * np.sin(beta))
        iCr = (gamma+delta)/2
        theta2 = (gamma-delta)/2 + theta1
        v2 = v1/np.sin(iCr)
        zL.append(v1*(interS[2,0] - zL[0]/v0 * (np.cos(alpha+beta)+1)/np.cos(alpha))/(2*np.cos(iCr)))
        zR.append(v1*(interS[2,1] - zR[0]/v0 * (np.cos(alpha+beta)+1)/np.cos(beta))/(2*np.cos(iCr)))
        hL.append(1/np.cos(theta2) * (zL[0]*np.cos(alpha-theta2+theta1)/np.cos(alpha) + zL[1]))
        hR.append(1/np.cos(theta2) * (zR[0]*np.cos(beta+theta2-theta1)/np.cos(beta) + zR[1]))
        v = [v0, v1, v2]
    return v, hL, hR


## Need to take a closer look at this: https://programmerall.com/article/10751929193/
# https://matplotlib.org/devdocs/gallery/widgets/polygon_selector_demo.html#polygon-selector (for the line selection)
# https://build-system.fman.io/ (for building into executable)
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=75):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.cBar = None # For the (eventual) colorbar
        super(MplCanvas, self).__init__(self.fig)
        self.setParent(parent)

# class VerticalNavigationToolbar2QT(NavigationToolbar2QT):
#     def __init__(self, canvas, parent, coordinates=True):
#         super().__init__(canvas, parent, coordinates)
class paths():
    def __init__(self) -> None:
        self.directory = []
        self.geometryFile = []
        self.nbFiles = []
        self.seg2Files = []
class geometry():
    def __init__(self) -> None:
        self.sensors = []
        self.sourcesId = []
        self.receivers = []
class model():
    def __init__(self) -> None:
        self.nbLayers = 2
        self.thickLeft = np.ones((self.nbLayers-1,))*5
        self.dipAngles = np.zeros((self.nbLayers-1,))
        self.Vp = np.linspace(1000,3000,num=self.nbLayers)
    def changeLayers(self, nbLayers):
        '''Initialize the model for a given number of layers'''
        self.nbLayers = nbLayers
        self.thickLeft = np.ones((self.nbLayers-1,))*5
        self.dipAngles = np.zeros((self.nbLayers-1,))
        self.Vp = np.linspace(1000,3000,num=self.nbLayers)
class animationPicking():
    def __init__(self) -> None:
        self.timeOnClick = 0        # To know if the click is in fixed position of to pan/zoom
        self.mousePosition = 0      # Storing the current position of the mouse
        self.currSelect = 0         # Storing the current trace number beiing analyzed
        self.changedSelect = True   # Variable to tell if the trace selected is different
        self.first = True           # Variable to tell if plotting for the first time
        self.maxClickLength = 0.5   # Maximum time (in sec) to consider a click on place
class inversionData():
    def __init__(self) -> None:
        self.lam = 20.0
        self.zWeight = 0.5
        self.vTop = 500.0
        self.vBottom = 3000.0
        self.vMin = 10.0
        self.vMax = 5000.0
        self.startModel = None
        self.meshMaxCellSize = 5.0
        self.meshDepthMax = 50.0
        # Pygimli inversion features:
        self.mesh = None
        self.data = None
        self.manager = None
    def setStartModelGradient(self, data, mesh):
        self.startModel = pg.Vector(ratools.createGradientModel2D(data, mesh, self.vTop, self.vBottom))
class modellingAnimation():
    def __init__(self) -> None:
        self.currPosition = [0, 0]
        self.maxClickLength = 0.5
        self.timeOnClick = 0
        self.offsetVelocity = 0.5
        self.tol = [1, 0.002] # x in meteres, y in seconds
        self.currPointId = 1
        self.changingPts = False
        self.namesSources = []
        self.namesOrientations = []
class modellingData():
    def __init__(self) -> None:
        self.sensors = []
        self.measurements = []
        self.hodoPoints = []
        self.combinationSR = []
        self.model = []
        self.nbLayers = 2
        self.sourcesX = []
        self.orientations = []
        self.appVelocities = []
        self.interceptTime = []
class dataStorage():
    def __init__(self) -> None:
        ## Data variables:
        self.paths = paths()
        self.geometry = geometry()
        self.sisDataOriginal = []
        self.sisData = []
        self.picking = []
        self.pickingError = []
        self.model = model()
        self.invData = inversionData()
        self.modellingAnimation = modellingAnimation()
        self.modellingData = modellingData()
        ## Status variables:
        self.dataLoaded = False
        self.pickingDone = False
        self.inversionDone = False
        self.meshLoaded = False
        self.sisFileId = 0
        self.beginTime = []
        ## Graphical animation variables:
        self.animationPicking = animationPicking()

class picklingStatus():
    def __init__(self) -> None:
        self.paths = paths()
        self.geometry = geometry()
        self.sisDataOriginal = []
        self.sisData = []
        self.beginTime = []
        self.picking = []
        self.pickingError = []
        self.sisFileId = 0
class Window(QMainWindow):
    def __init__(self) -> None:
        super().__init__()

        ## Initializing the data structure
        self.dataUI = dataStorage()

        ## Initialize the UI
        self.setWindowTitle('Sardine Reborn')
        self.setWindowIcon(QtGui.QIcon('./images/SardineRebornLogo_100ppp.png'))
        self.resize(1100,600)

        ## Adding tabs to the layout
        self.tabs = QTabWidget(self)
        self.tabs.addTab(self._pickTracesTabUI(), 'Traces picking')
        self.tabs.addTab(self._inversionTabUI(), 'Inversion')
        self.tabs.addTab(self._modelTabUI(), 'Model')

        self.setCentralWidget(self.tabs) 

        ## Defining menu actions:
        ### File openning / saving picking
        openFile = QAction("&Open Geometry File",self)
        openFile.setShortcut('Ctrl+O') # Setting ctrl+o as the shortcut
        openFile.setStatusTip('Open the *.geometry file')
        openFile.triggered.connect(self._openGeometry)

        savePicking = QAction("Save Current &Picking",self)
        savePicking.setShortcut('Ctrl+P')
        savePicking.setStatusTip('Save the picking state in a *.sgt file')
        savePicking.triggered.connect(self._savePicking)

        loadPicking = QAction("&Load Picking File", self)
        loadPicking.setShortcut('Ctrl+L')
        loadPicking.setStatusTip('Load an existing *.sgt file for inversion')
        loadPicking.triggered.connect(self._loadPicking)

        savePickingStatus = QAction("Save the picking state", self)
        savePickingStatus.setStatusTip('Save the status into a pickle')
        savePickingStatus.triggered.connect(self.saveStatePicking)

        loadPickingStatus = QAction("Load a picking state", self)
        loadPickingStatus.setStatusTip('Load a pickle status')
        loadPickingStatus.triggered.connect(self.loadStatePicking)

        ### Filtering of the signal:
        dcFilter = QAction("Apply a DC filter", self)
        dcFilter.setStatusTip('Apply a DC filter to the signal')
        dcFilter.triggered.connect(self.dcFilter)

        trimFilter = QAction("Trim the datasets", self)
        trimFilter.setStatusTip('Apply trimming to the signal')
        trimFilter.triggered.connect(self.trimFilter)

        highpassFilter = QAction("Apply a high-pass filter", self)
        highpassFilter.setStatusTip('Apply a high pass filter to the signal')
        highpassFilter.triggered.connect(self.highpassFilter)

        lowpassFilter = QAction("Apply a low-pass filter", self)
        lowpassFilter.setStatusTip('Apply a low pass filter to the signal')
        lowpassFilter.triggered.connect(self.lowpassFilter)

        resetFilters = QAction("Reset the filters", self)
        resetFilters.setStatusTip('Reset all previously applied filters')
        resetFilters.triggered.connect(self.resetFilters)

        autoPicking = QAction("Automated picking", self)
        autoPicking.setStatusTip('Automatically pick the fisrt arrival for each signals ().')
        autoPicking.triggered.connect(self.autoPicking)

        ### Inversion through pigimli:
        loadInvMesh = QAction('Load Inversion Mesh', self)
        loadInvMesh.setStatusTip('Load an inversion mesh that was already created (*.poly)')
        loadInvMesh.triggered.connect(self._loadInvMesh) 

        loadInitialModel = QAction('Load Initial Model', self)
        loadInitialModel.setStatusTip('Load an existing model as the starting model (*.vector) for the inversion')
        loadInitialModel.triggered.connect(self._loadInitModel)
        
        saveInvMesh = QAction('Save Inversion Mesh', self)
        saveInvMesh.setStatusTip('Save the inversion mesh (*.poly)')
        saveInvMesh.triggered.connect(self._saveInvMesh) 
        
        saveInvAsVTK = QAction('Save Inversion as VTK', self)
        saveInvAsVTK.setStatusTip('Save the inversion results into a VTK file (for Paraview)')
        saveInvAsVTK.triggered.connect(self._saveInvAsVTK)

        saveInvResponse = QAction('Save the Inverse Response', self)
        saveInvResponse.setStatusTip('Save the model response for the last iteration (*.vector)')
        saveInvResponse.triggered.connect(self._saveInvResponse) 

        saveInvResult = QAction('Save the Inverse Results', self)
        saveInvResult.setStatusTip('Save the model for the last iteration (*.vector)')
        saveInvResult.triggered.connect(self._saveInvResult)

        # saveModel = QAction("Save Current &Model",self)
        # saveModel.setShortcut('Ctrl+M')
        # saveModel.setStatusTip('Save the current model in a *.txt file')
        # saveModel.triggered.connect(self._saveModel)

        # Adding the menu bar atop
        menuBarInternal = self.menuBar()
        menuBarInternal.setNativeMenuBar(False)
        fileMenu = menuBarInternal.addMenu("Picking")
        fileMenu.addAction(openFile)
        fileMenu.addAction(savePicking)
        fileMenu.addAction(loadPicking)
        fileMenu.addAction(savePickingStatus)
        fileMenu.addAction(loadPickingStatus)
        fileMenu = menuBarInternal.addMenu("Filters")
        fileMenu.addAction(dcFilter)
        fileMenu.addAction(trimFilter)
        fileMenu.addAction(highpassFilter)
        fileMenu.addAction(lowpassFilter)
        fileMenu.addAction(resetFilters)
        fileMenu.addAction(autoPicking)
        # fileMenu.addAction(saveModel)
        fileMenu = menuBarInternal.addMenu("Inversion")
        fileMenu.addAction(loadInvMesh)
        fileMenu.addAction(loadInitialModel)
        fileMenu.addAction(saveInvMesh)
        fileMenu.addAction(saveInvResult)
        fileMenu.addAction(saveInvResponse)
        fileMenu.addAction(saveInvAsVTK)

        ## Defining the status bar:
        self.statusBar = QStatusBar(self)
        self.statusBar.showMessage(defaultStatus)
        self.setStatusBar(self.statusBar)

        pass

    ## Oppening and closing message boxes:
    def showEvent(self, event):
        msgBox = QMessageBox(self)
        msgBox.setIconPixmap(QtGui.QPixmap('./images/SardineRebornLogo_100ppp.png').scaled(200,100,aspectRatioMode=Qt.KeepAspectRatio))
        msgBox.setText('Welcome to Sardine Reborn!')
        msgBox.setInformativeText('by Hadrien Michel (2022)')
        msgBox.setWindowTitle('Welcome!')
        msgBox.show()
        event.accept()
    
    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Closing ...', 'Are you sure you want to quit?', QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Ok)
        if reply == QMessageBox.Ok:
            event.accept()
            matplotlib.pyplot.close('all') # To get rid of the figure openning after the closing.
        else:
            event.ignore()
        
    ## Animations definitions:
    # PyQt5 animations
    
    # Matplotlib animations:
    def changeMouse(self, event):
        if event.inaxes is not None:
            self.dataUI.animationPicking.mousePosition = event.xdata
        return 0

    def onPress(self, event): # For both windows
        if event.button == MouseButton.LEFT or event.button == MouseButton.RIGHT:
            self.dataUI.animationPicking.timeOnClick = time.time()
        return 0

    def onRelease(self, event):
        if event.button == MouseButton.LEFT and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength): # If left click and not dragging accross the pannel
            if self.buttonTabPickingSetT0.isChecked():
                self.dataUI.beginTime[self.dataUI.sisFileId] -= self.dataUI.animationPicking.mousePosition
                # Change all picked signals for this stream:
                self.dataUI.picking[self.dataUI.sisFileId,:] -= self.dataUI.animationPicking.mousePosition
                self.buttonTabPickingSetT0.setChecked(False)
            else:
                if self.dataUI.animationPicking.mousePosition < 0: # To remove a picked trace, click on times below 0
                    self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.nan
                    self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.nan
                else:
                    self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = self.dataUI.animationPicking.mousePosition
                    self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = max(self.dataUI.animationPicking.mousePosition*0.03, 0.000001) # Default error is 3%
            self.dataUI.animationPicking.changedSelect = True
        if event.button == MouseButton.RIGHT and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength): # If right click and not dragging accross the pannel
            if not(np.isnan(self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect])):
                if self.dataUI.animationPicking.mousePosition > 0:
                    self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.abs(self.dataUI.animationPicking.mousePosition-self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect])/self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] # Default error is 3%
                    self.dataUI.animationPicking.changedSelect = True
        return 0
    
    def onPressModelling(self, event):
        if event.button == MouseButton.LEFT:
            sourceId = self.sourceSelector.currentIndex()
            receiversOrientation = self.receiversSelector.currentIndex()
            points = np.asarray(self.dataUI.modellingData.hodoPoints[sourceId][receiversOrientation])
            nbPoints = points.shape[0]
            mousePosition = np.repeat(np.asarray([event.xdata, event.ydata]).reshape((1,2)), nbPoints, axis=0)
            diff = np.abs(points-mousePosition)
            diffXBool = diff[:,0] < self.dataUI.modellingAnimation.tol[0]
            diffYBool = diff[:,1] < self.dataUI.modellingAnimation.tol[1]
            ptsId = np.logical_and(diffXBool, diffYBool)
            ptsId = np.where(ptsId)[0]
            if len(ptsId) > 1:
                ptsId = ptsId[-1] # We take the last from both possibilities
            if len(ptsId) > 0:
                if ptsId > 0:
                    self.dataUI.modellingAnimation.currPointId = ptsId
                    self.dataUI.modellingAnimation.changingPts = True
                else:
                    self.dataUI.modellingAnimation.changingPts = False
        return

    def onReleaseModelling(self, event):
        self.dataUI.modellingAnimation.changingPts = False    

    def changeMouseModelling(self, event): # For the modelling window
        if (event.inaxes is not None) and self.dataUI.modellingAnimation.changingPts:
            sourceId = self.sourceSelector.currentIndex()
            receiversOrientation = self.receiversSelector.currentIndex()
            ptsId = self.dataUI.modellingAnimation.currPointId
            # print(f'Mouse at x={event.xdata}, y={event.ydata} changed')
            self.dataUI.modellingData.hodoPoints[sourceId][receiversOrientation][ptsId,:] = [event.xdata, event.ydata]
        return 0
    
    def animationZoom(self, i):
        axZoom = self.zoomGraph.axes
        # Get axis variables
        deltaT = float(self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        timeSEG2 = np.arange(self.dataUI.beginTime[self.dataUI.sisFileId], self.dataUI.beginTime[self.dataUI.sisFileId]+nbPoints*deltaT, deltaT)
        # Change plot to go at the correct position:
        axZoom.clear()
        idx = np.greater_equal(timeSEG2,self.dataUI.animationPicking.mousePosition-150*deltaT) & np.less_equal(timeSEG2,self.dataUI.animationPicking.mousePosition+150*deltaT)
        timeZoom = timeSEG2[idx]
        axZoom.plot(timeZoom,self.dataUI.sisData[self.dataUI.sisFileId][self.dataUI.animationPicking.currSelect].data[idx],color='k')
        axZoom.set_xlim(left=self.dataUI.animationPicking.mousePosition-150*deltaT,right=self.dataUI.animationPicking.mousePosition+150*deltaT)
        axZoom.autoscale(axis='y')
        z = axZoom.get_ylim()
        axZoom.plot([self.dataUI.animationPicking.mousePosition, self.dataUI.animationPicking.mousePosition],z,color='r')
        currPicking = self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect]
        if not(np.isnan(currPicking)):
            currError = self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] * currPicking
            axZoom.plot([currPicking, currPicking], z, 'g')
            axZoom.plot([currPicking - currError, currPicking - currError], z, ':g')
            axZoom.plot([currPicking + currError, currPicking + currError], z, ':g')
        axZoom.set_frame_on(False)
        axZoom.tick_params(
            axis='both',       # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left=False,
            right=False,
            labelleft=False,
            labelbottom=False) # labels along the bottom edge are off
        self.zoomGraph.draw()
        return 0
    
    def animationMain(self, i):
        axMain = self.mainGraph.axes
        # Get axis variables
        deltaT = float(self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        timeSEG2 = np.arange(self.dataUI.beginTime[self.dataUI.sisFileId], self.dataUI.beginTime[self.dataUI.sisFileId]+nbPoints*deltaT, deltaT)
        if self.dataUI.animationPicking.changedSelect:
            # Change red graph + pick
            if not(self.dataUI.animationPicking.first):
                limitsY = axMain.get_ylim()
                limitsX = axMain.get_xlim()
            axMain.clear()
            i = 0
            for tr in self.dataUI.sisData[self.dataUI.sisFileId]:
                data = tr.data 
                data = data/(max(data)-min(data))+i
                if i == self.dataUI.animationPicking.currSelect:
                    axMain.plot(timeSEG2,data,color='r')
                else:
                    axMain.plot(timeSEG2,data,color='k')
                currPicking = self.dataUI.picking[self.dataUI.sisFileId, i]
                if not(np.isnan(currPicking)):
                    currError = self.dataUI.pickingError[self.dataUI.sisFileId, i] * currPicking
                    axMain.plot([currPicking, currPicking], [i-0.5, i+0.5],color='g')
                    axMain.plot([currPicking - currError, currPicking - currError], [i-0.25, i+0.25], ':g')
                    axMain.plot([currPicking + currError, currPicking + currError], [i-0.25, i+0.25], ':g')
                i += 1
            if not(self.dataUI.animationPicking.first):
                axMain.set_xlim(left=limitsX[0],right=limitsX[1])
                axMain.set_ylim(bottom=limitsY[0],top=limitsY[1])
            self.dataUI.animationPicking.first=False
            self.dataUI.animationPicking.changedSelect = False
        self.mainGraph.draw()
        return 0
    
    ## UI objects definition
    def _openGeometry(self):
        self.statusBar.showMessage('Openning Geometry file . . .')
        ## Opening the geometry file:
        fname, _ = QFileDialog.getOpenFileName(self,'Open geometry file',filter='Geometry file (*.geometry)')# The first argument returned is the filename and path
        if fname != "":
            headTail = os.path.split(fname)
            path = headTail[0]
            file = headTail[1]
            # nameSave = file[:-9] # remove the *.geometry
            # Retreive the datafiles names:
            SEG2Files = []
            SourcePosition = []
            ReceiversPosition = []
            sources = False
            receivers = False
            with open(fname) as f:
                Lines = f.read().splitlines()
            for line in Lines:
                if len(line.strip('\t')) != 0:
                    if line.startswith("SOURCES"):
                        sources = True
                        receivers = False
                    elif line.startswith("RECEIVERS"):
                        sources = False
                        receivers = True
                    else:
                        if sources:
                            tmp = re.split(r'\t+', line.strip('\t'))
                            name = tmp[0]
                            CurrSource = [float(i) for i in tmp[1:]]
                            SEG2Files.append(name)
                            SourcePosition.append(CurrSource)
                        elif receivers:
                            CurrReceiver = [float(i) for i in re.split(r'\t+', line.strip('\t'))]
                            ReceiversPosition.append(CurrReceiver)
            # Check if Sources in List of Receivers --> Constitute Sensors array for output file:
            sensors = deepcopy(ReceiversPosition)
            sourcesId = np.zeros((len(SourcePosition),))
            i = 0
            for source in SourcePosition:
                if not(ReceiversPosition.count(source) == 1): # If the source is not in the receivers array
                    sensors.append(source)
                sourcesId[i] = sensors.index(source) # We store the position of the Id position of the current source in the sensor array
                i += 1
        
            if self.dataUI.dataLoaded == True:
                reply = QMessageBox.question(self, 'Overwritting data . . .', 'Are you sure you want to overwrite the current dataset?', QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
                if reply == QMessageBox.Yes:
                    ## Setting up the paths:
                    self.saveDataUI(path, file, SEG2Files, ReceiversPosition, sensors, sourcesId)

                    ## Updating status bar
                    self.dataUI.dataLoaded = True
                    self.statusBar.showMessage(f'{len(SEG2Files)} data files retreived from the geometry file with {len(sensors)} sensors.')
                else:
                    self.statusBar.showMessage(f'No data loaded')
            else:
                self.saveDataUI(path, file, SEG2Files, ReceiversPosition, sensors, sourcesId)

                ## Updating status bar
                self.dataUI.dataLoaded = True
                self.statusBar.showMessage(f'{len(SEG2Files)} data files retreived from the geometry file with {len(sensors)} sensors.')
            
            ## Return to the picking tab
            self.updateTab0()
            self.tabs.setCurrentIndex(0)
        else:
            self.statusBar.showMessage('No file loaded!')

    def saveDataUI(self, path, file, SEG2Files, ReceiversPosition, sensors, sourcesId):
        ## Setting up the paths:
        self.dataUI.paths.directory = path
        self.dataUI.paths.geometryFile = file
        self.dataUI.paths.nbFiles = len(sourcesId)
        self.dataUI.paths.seg2Files = SEG2Files

            ## Setting up the data:
        self.dataUI.geometry.sensors = sensors
        self.dataUI.geometry.sourcesId = sourcesId
        self.dataUI.geometry.receivers = ReceiversPosition

            ## Reading the datasets:
        for name  in SEG2Files:
            _, ext = os.path.splitext(name)
            if ext == '.segy' or ext == '.sgy':
                st = read(os.path.join(path,name), 'SEGY')
                self.dataUI.beginTime.append(0)
            elif ext == '.seg2' or ext == '.sg2':
                st = read(os.path.join(path,name), 'SEG2')
                self.dataUI.beginTime.append(float(st[0].stats.seg2["DELAY"].replace(',', '.')))
            else:
                st = read(os.path.join(path,name), 'SEG2')
                self.dataUI.beginTime.append(float(st[0].stats.seg2["DELAY"].replace(',', '.')))
            if len(st) != len(ReceiversPosition): # If the number of geophones does not match between the loaded array and the gemoetry
                raise Exception('The file referenced in the geometry file does not match the geometry of the array!')
            self.dataUI.sisData.append(st)
        
        self.dataUI.sisDataOriginal = deepcopy(self.dataUI.sisData)
        # self.filterSignal()
        self.dataUI.picking = np.empty((len(self.dataUI.paths.seg2Files),len(self.dataUI.sisData[0])))
        self.dataUI.picking[:] = np.nan
        self.dataUI.pickingError = np.empty((len(self.dataUI.paths.seg2Files),len(self.dataUI.sisData[0])))
        self.dataUI.pickingError[:] = np.nan
        self.spinBoxCurrSelect.setMinimum(0)
        self.spinBoxCurrSelect.setMaximum(len(self.dataUI.sisData[0])-1)
        self.spinBoxCurrSelect.setPrefix('Trace number ')
        self.spinBoxCurrSelect.valueChanged.connect(self.traceNumberChanged)
        self.spinBoxCurrSelect.setEnabled(True)

    def dcFilter(self):
        for st in self.dataUI.sisData:
            for tr in st:
                tr.detrend('constant')
        self.dataUI.animationPicking.changedSelect = True

    def highpassFilter(self):
        band, ok = QInputDialog.getDouble(self, "High-pass filter", "Frequency [Hz]", 5.0, 0.0, 10000.0)
        if ok:
            for st in self.dataUI.sisData:
                for tr in st:
                    tr.filter('highpass', freq=band)
        self.dataUI.animationPicking.changedSelect = True

    def lowpassFilter(self):
        band, ok = QInputDialog.getDouble(self, "Low-pass filter", "Frequency [Hz]", 100.0, 0.0, 10000.0)
        if ok:
            for st in self.dataUI.sisData:
                for tr in st:
                    tr.filter('lowpass', freq=band)
        self.dataUI.animationPicking.changedSelect = True

    def trimFilter(self):
        deltaT = float(self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        tEnd = self.dataUI.beginTime[self.dataUI.sisFileId]+nbPoints*deltaT
        trim, ok = QInputDialog.getDouble(self, "Trimming", "Time: ", tEnd/2, 0.0, tEnd, 2)
        if ok:
            for st in self.dataUI.sisData:
                st.trim(starttime=st[0].stats.starttime, endtime=st[0].stats.starttime + trim)
        self.dataUI.animationPicking.changedSelect = True

    def resetFilters(self):
        self.dataUI.sisData = deepcopy(self.dataUI.sisDataOriginal)
        self.dataUI.animationPicking.changedSelect = True

    def filterSignal(self, dcFilter:bool=True, highpass:float=5.0, lowpass:float=100.0):
        self.dataUI.sisData = deepcopy(self.dataUI.sisDataOriginal)
        if not(dcFilter) and (highpass > 0.0 or lowpass < np.Inf):
            QMessageBox.warning(self, 'Warning !', 'Using the lowpass and/or highpass filter\n without the DC filter is not advised!')
        for st in self.dataUI.sisData:
            for tr in st:
                if dcFilter:
                    tr.detrend('constant')
                if highpass > 0.0:
                    tr.filter('highpass', freq=highpass)
                if lowpass < np.Inf:
                    tr.filter('lowpass', freq=lowpass)
        self.dataUI.animationPicking.changedSelect = True
    
    def autoPicking(self):
        pBar = QProgressBar(self)
        nbTraces = len(self.dataUI.sisData)*len(self.dataUI.sisData[0])
        pBar.setMaximum(nbTraces)
        for i, st in enumerate(self.dataUI.sisData):
            for j, tr in enumerate(st):
                aic_f = aic_simple(tr.data)
                p_idx = aic_f.argmin()
                self.dataUI.picking[i, j] = p_idx/tr.stats.sampling_rate
                self.dataUI.pickingError[i, j] = max(self.dataUI.picking[i, j] * 0.05, 0.000001)
                pBar.setValue((i*len(self.dataUI.sisData[0]) + (j)))
        pBar.close()
        self.dataUI.animationPicking.changedSelect = True
    
    def traceNumberChanged(self, value):
        self.dataUI.animationPicking.currSelect = value
        self.dataUI.animationPicking.changedSelect = True

    def updateTab0(self):
        self.comboBoxFilesPicking.clear()
        for name in self.dataUI.paths.seg2Files:
            self.comboBoxFilesPicking.addItem(name)
        self.dataUI.sisFileId = 0
        self.comboBoxFilesPicking.currentIndexChanged.connect(self.comboBoxChange)
        self.connectMouse = self.mainGraph.mpl_connect('motion_notify_event', self.changeMouse)
        self.connectPress = self.mainGraph.mpl_connect('button_press_event', self.onPress)
        self.connectRelease = self.mainGraph.mpl_connect('button_release_event', self.onRelease)
        self.aniMain = animation.FuncAnimation(self.mainGraph.fig, self.animationMain, interval=16.7)
        self.aniZoom = animation.FuncAnimation(self.zoomGraph.fig, self.animationZoom, interval=16.7)
        self.mainGraph.draw()
        self.zoomGraph.draw()

    def comboBoxChange(self, newId):
        self.dataUI.sisFileId = newId
        self.dataUI.animationPicking.changedSelect = True

    def _savePicking(self):
        self.statusBar.showMessage('Save current picking . . .')
        ## Building the sgt array:
        sensors = self.dataUI.geometry.sensors
        sourcesId = self.dataUI.geometry.sourcesId
        receivers = self.dataUI.geometry.receivers
        picksSave = []
        for nbFile in range(len(self.dataUI.paths.seg2Files)):
            for i in range(len(self.dataUI.sisData[0])):
                sId = int(sourcesId[nbFile])
                rId = int(sensors.index(receivers[i]))
                if sId != rId: # The traveltime for source = receiever is 0 and not usefull for inversion!
                    t = self.dataUI.picking[nbFile, i]
                    err = self.dataUI.pickingError[nbFile, i] * t
                    if not(np.isnan(t)):
                        picksSave.append([sId, rId, t, err])
        # Remove unused sensors from the list:
        usedSensors = [False]*len(sensors)
        for pick in picksSave:
            usedSensors[pick[0]] = True
            usedSensors[pick[1]] = True
        oldId = range(len(sensors))
        oldId = [i for i in range(len(sensors)) if usedSensors[i]]
        sensors = [sensors[i] for i in range(len(sensors)) if usedSensors[i]]
        newId = range(len(sensors))
        for pick in picksSave:
            pick[0] = newId[oldId.index(pick[0])]
            pick[1] = newId[oldId.index(pick[1])]
        ## Saving the file:
        fname, _ = QFileDialog.getSaveFileName(self,'Select file to save',filter='Source-Receiver-Time file (*.sgt)')# The first argument returned is the filename and path
        if fname != "":
            f = open(os.path.join(fname),'w')# Create a new file
            nbSensors = len(sensors)
            f.write('%d # shot/geophone points\n' % nbSensors)
            f.write('#x\ty\n')
            for i in range(nbSensors):
                f.write('%.2f\t%.2f\n' % (sensors[i][0], sensors[i][1]))
            nbMeas = len(picksSave)
            f.write('%d # measurements\n' % nbMeas)
            f.write('#s\tg\tt\terr\n')
            for i in range(nbMeas):
                f.write('%d\t%d\t%f\t%f\n' % (picksSave[i][0]+1, picksSave[i][1]+1, picksSave[i][2], picksSave[i][3]))
            f.close()
            self.statusBar.showMessage(defaultStatus)
            self.filePicksPath.setText(fname)
            self._initPygimli(fname)
        else:
            self.statusBar.showMessage('No file saved!')


    def _loadPicking(self):
        self.statusBar.showMessage('Loading picking file . . .')
        # 1) Load a file with the first arrival:
        fName, _ = QFileDialog.getOpenFileName(self,'Select file to load',filter='Source-Receiver-Time file (*.sgt)')
        # We retreived a first-arrival file --> geometry of the sensors + first arrivals 
        if fName != "":
            with open(fName) as f:
                lines = f.read().splitlines()
            markerNbSensors = "# shot/geophone points"
            markerMeasurements = "# measurements"
            for line in lines:
                if line.endswith(markerNbSensors):
                    nbSensors = int(line[:-len(markerNbSensors)])
                    sensors = np.zeros((nbSensors,2))
                    idxSensor = 0
                elif line.endswith("#x\ty"):
                    pass
                elif idxSensor < nbSensors:
                    sensors[idxSensor,:] = re.split(r'\t+', line)
                    idxSensor += 1
                elif line.endswith(markerMeasurements):
                    nbMeasurements = int(line[:-len(markerMeasurements)])
                    measurements = np.zeros((nbMeasurements,4)) # s g t err
                    idxMeas = 0
                elif line.endswith('#s\tg\tt\terr'):
                    measurements = np.zeros((nbMeasurements,4)) # s g t err
                elif line.endswith('#s\tg\tt'):
                    measurements = np.zeros((nbMeasurements,3)) # s g t
                elif idxMeas < nbMeasurements:
                    measurements[idxMeas,:] = re.split(r'\t+', line)
                    idxMeas += 1
            self.dataUI.modellingData.sensors = sensors
            self.dataUI.modellingData.measurements = measurements
            loadPicking = False
            if self.dataUI.dataLoaded:
                # Check if the sgt file corresponds to the already loaded dataset (similar array)
                sensorsInList = True
                for i in range(nbSensors):
                    currSensor = sensors[i,:]
                    equals = np.all(self.dataUI.geometry.sensors == currSensor, axis=1)
                    if np.sum(equals) < 1:
                        # The sensor is not in the list
                        sensorsInList = False
                        break
                if sensorsInList:
                    if not(np.all(np.isnan(self.dataUI.picking))):
                        # Load the picking to add atop the traces:
                        reply = QMessageBox.question(self, 'Overwritting picking . . .', 'Are you sure you want to overwrite the current picking?', QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
                        if reply == QMessageBox.Yes:
                            ## loading the picking:
                            loadPicking = True
                        else:
                            ## Do nothing for the picking window:
                            loadPicking = False
                            self.statusBar.showMessage(f'Data loaded but picking not presented')
                    else:
                        loadPicking = True
            if loadPicking:
                sources = np.asarray(self.dataUI.geometry.sensors)
                sources = sources[list(self.dataUI.geometry.sourcesId.astype(int)), :]
                receivers = np.asarray(self.dataUI.geometry.receivers)
                self.dataUI.pickingError[:] = np.nan
                self.dataUI.picking[:] = np.nan
                for i in range(nbMeasurements):
                    sCurr = sensors[int(measurements[i,0])-1,:]
                    rCurr = sensors[int(measurements[i,1])-1,:]
                    pickCurr = measurements[i,2]
                    if len(measurements[i,:]) > 3:
                        errCurr = measurements[i,3]
                    sId = np.where(np.all(sources == sCurr, axis=1))[0]
                    rId = np.where(np.all(receivers == rCurr, axis=1))[0]
                    self.dataUI.picking[sId, rId] = pickCurr
                    if len(measurements[i,:]) > 3:
                        self.dataUI.pickingError[sId, rId] = errCurr/pickCurr
                    else:
                        self.dataUI.pickingError[sId, rId] = 0.03
                self.statusBar.showMessage(f'Data loaded with picking on graphs')
                self.dataUI.animationPicking.changedSelect = True
            self.filePicksPath.setText(fName)
            self._initPygimli(fName)
            self._initModelling()
    
    def _initPygimli(self, fname):
        ## Preparing inversion of data (pygimli)
        self.dataUI.invData.data = pg.DataContainer(fname, sensorTokens='s g')
        self.dataUI.invData.manager = TTMgr(self.dataUI.invData.data)
        self.dataUI.invData.mesh = self.dataUI.invData.manager.createMesh(data=self.dataUI.invData.data, paraMaxCellSize=self.dataUI.invData.meshMaxCellSize, paraDepth=self.dataUI.invData.meshDepthMax)
        self.invModelGraph.axes.clear()
        self.dataGraph.axes.clear()
        self.fitGraph.axes.clear()
        pgshow(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
        drawFirstPicks(ax=self.dataGraph.axes, data=self.dataUI.invData.data)
        self.dataGraph.fig.tight_layout()
        self.dataGraph.draw()
        self.invModelGraph.fig.tight_layout()
        self.invModelGraph.draw()

    def _initModelling(self):
        if len(self.dataUI.modellingData.sensors) == 0:
            return
        sensors, measurements = self.dataUI.modellingData.sensors, self.dataUI.modellingData.measurements
        if np.all(np.abs(sensors[:,-1]) < 1e-6):
            # Plotting the hodochrones:
            self.plotHodochrones(sensors, measurements)
        else:
            QMessageBox.warning(self, 'Warning !', 'Impossible to model using the intercept time method!\nNo topography authorized.')
            return
        # Initialize the modelling:
        self.dataUI.modellingData.nbLayers = self.nbLayersSelector.value()
        self.dataUI.modellingData.hodoPoints = []
        self.dataUI.modellingData.combinationSR = []
        self.dataUI.modellingData.sourcesX = []
        self.dataUI.modellingData.appVelocities = []
        self.dataUI.modellingData.interceptTime = []
        self.dataUI.modellingData.orientations = []
        # Gather the possible sources and orientations:
        self.dataUI.modellingAnimation.namesSources = []
        self.dataUI.modellingAnimation.namesOrientations = []
        sources = np.unique(measurements[:,0]).astype(int)
        for sId in sources:
            sourceX = sensors[sId-1,0]
            self.dataUI.modellingAnimation.namesSources.append(f'Source at {sourceX} m.')
            self.dataUI.modellingData.sourcesX.append(sourceX)
            index = measurements[:,0].astype(int) == sId
            sourceX = sensors[sId-1, 0]
            receiversX = sensors[measurements[index, 1].astype(int) - 1, 0]
            times = measurements[index, 2]
            receiversLeft = receiversX[receiversX < sourceX]
            receiversRight = receiversX[receiversX > sourceX]
            orientationsText = []
            orientations = []
            appVel = []
            intercept = []
            hodoPts = []
            if receiversLeft.size != 0:
                orientationsText.append('Left')
                orientations.append(-1)
                self.dataUI.modellingData.combinationSR.append([sourceX, -1])
                inter, v, points = buildModel(sourceX, receiversLeft, times[receiversX < sourceX], self.dataUI.modellingData.nbLayers, -1)
                appVel.append(v)
                intercept.append(inter)
                hodoPts.append(points)
            if receiversRight.size != 0:
                orientationsText.append('Right')
                orientations.append(1)
                self.dataUI.modellingData.combinationSR.append([sourceX, 1])
                inter, v, points = buildModel(sourceX, receiversRight, times[receiversX > sourceX], self.dataUI.modellingData.nbLayers, 1)
                appVel.append(v)
                intercept.append(inter)
                hodoPts.append(points)
            self.dataUI.modellingAnimation.namesOrientations.append(orientationsText)
            self.dataUI.modellingData.appVelocities.append(appVel)
            self.dataUI.modellingData.interceptTime.append(intercept)
            self.dataUI.modellingData.orientations.append(orientations)
            self.dataUI.modellingData.hodoPoints.append(hodoPts)
        self._updateHodoGraph()
        self._updateModelGraph()
        # Implement the source/receiver selector:
        self.sourceSelector.clear()
        self.receiversSelector.clear()
        self.sourceSelector.addItems(self.dataUI.modellingAnimation.namesSources)
        self.receiversSelector.addItems(self.dataUI.modellingAnimation.namesOrientations[0])
        self.sourceSelector.currentIndexChanged.connect(self.sourceSelectorChanged)
    
    def _updateHodoGraph(self):
        axHod = self.hodochronesGraph.axes
        axHod.cla()
        self.plotHodochrones(self.dataUI.modellingData.sensors, self.dataUI.modellingData.measurements)
        colors = matplotlib.pyplot.cm.tab10(np.arange(10))
        xAxisShow = np.linspace(np.min(self.dataUI.modellingData.sensors[:,0]), np.max(self.dataUI.modellingData.sensors[:,0]),1000)
        for i, sourceX in enumerate(self.dataUI.modellingData.sourcesX):
            for j, orientation in enumerate(self.dataUI.modellingData.orientations[i]):
                vel = self.dataUI.modellingData.appVelocities[i][j]
                inter = self.dataUI.modellingData.interceptTime[i][j]
                points = self.dataUI.modellingData.hodoPoints[i][j]
                xShow = xAxisShow[(xAxisShow - sourceX)*orientation >= 0]
                for k in range(self.dataUI.modellingData.nbLayers):
                    times = inter[k] + (1/vel[k])*np.abs(xShow-sourceX)
                    axHod.plot(xShow, times, color=colors[i %10])
                axHod.plot(points[:,0], points[:,1], linestyle='none', marker='o', color='k', markersize=5)
        self.hodochronesGraph.draw()

    def _updateModelGraph(self):    
        axMod = self.modelGraph.axes
        axMod.cla()
        # Show the sensors array at the surface:
        axMod.plot(self.dataUI.modellingData.sensors[:,0], self.dataUI.modellingData.sensors[:,1], marker='v', color='k')
        axMod.grid()
        offsetVelocity = self.dataUI.modellingAnimation.offsetVelocity
        sourcesOrientation = []
        for i, sourceX in enumerate(self.dataUI.modellingData.sourcesX):
            for j, orientation in enumerate(self.dataUI.modellingData.orientations[i]):
                sourcesOrientation.append([i, j, sourceX, orientation])
                if orientation > 0:
                    ha='left'
                else:
                    ha='right'
                vel = self.dataUI.modellingData.appVelocities[i][j]
                inter = self.dataUI.modellingData.interceptTime[i][j]
                _, v, pos = model1D(inter, vel)
                axMod.plot(sourceX + pos[:,0]*orientation, pos[:,1], linestyle='none', color='k', marker='x')
                for k in range(self.dataUI.modellingData.nbLayers):
                    axMod.text(sourceX+pos[j,0] + offsetVelocity, pos[k,1] + 2*offsetVelocity, f'$v_{k}$ = {round(v[k],2)} m/s', ha=ha)
        sourcesOrientation = np.asarray(sourcesOrientation)
        sourcesOrientation[np.lexsort((sourcesOrientation[:,3], sourcesOrientation[:,2]))]
        sourcesLeft = sourcesOrientation[sourcesOrientation[:,3]<0, :]
        sourcesRight = sourcesOrientation[sourcesOrientation[:,3]>0, :]
        for s1 in sourcesRight.tolist():
            sourceX = s1[2]
            for s2 in sourcesLeft[sourcesLeft[:,2] > sourceX].tolist():
                vel1 = self.dataUI.modellingData.appVelocities[int(s1[0])][int(s1[1])]
                inter1 = self.dataUI.modellingData.interceptTime[int(s1[0])][int(s1[1])]
                vel2 = self.dataUI.modellingData.appVelocities[int(s2[0])][int(s2[1])]
                inter2 = self.dataUI.modellingData.interceptTime[int(s2[0])][int(s2[1])]
                vS = np.asarray([vel1, vel2]).T
                interpS = np.asarray([inter1, inter2]).T
                vSlope, hL, hR = modelWithSlope(interpS, vS)
                for i in range(self.dataUI.modellingData.nbLayers):
                    if i != self.dataUI.modellingData.nbLayers-1:
                        axMod.plot([s1[2], s2[2]], [hL[i], hR[i]], color='k')
                    if i == 0:
                        axMod.text(np.mean([s1[2], s2[2]]), 2*offsetVelocity, f'$v_{i}$ = {round(vSlope[i],2)} m/s', ha='center')
                    else:
                        axMod.text(np.mean([s1[2], s2[2]]), (hL[i-1] + hR[i-1])/2 + 2*offsetVelocity, f'$v_{i}$ = {round(vSlope[i],2)} m/s', ha='center')
        axMod.set_xlabel('X [m]')
        axMod.set_ylabel('Depth [m]')
        xRange = axMod.get_xlim()
        xRange = xRange[1] - xRange[0]
        axMod.set_ylim((-1, np.ceil(xRange/10)))
        axMod.invert_yaxis()
        self.modelGraph.draw()

    def _updateModelling(self, frame):
        nbLayers = self.dataUI.modellingData.nbLayers # Gather the number of layers in the model
        ## Gathering the info about intercept times and apparent velocities:
        sourceId = self.sourceSelector.currentIndex()
        receiversOrientation = self.receiversSelector.currentIndex()
        orientation = self.dataUI.modellingData.orientations[sourceId][receiversOrientation]
        points = np.asarray(self.dataUI.modellingData.hodoPoints[sourceId][receiversOrientation])
        nbLayers = self.dataUI.modellingData.nbLayers
        vel = np.zeros((nbLayers,))
        inter = np.zeros((nbLayers,))
        sourceX = points[0, 0]
        for i in range(self.dataUI.modellingData.nbLayers):
            vel[i] = np.abs((points[i+1,0] - points[i,0])/(points[i+1,1] - points[i,1]))
            inter[i] = 1/vel[i] * (sourceX-points[i,0]) * orientation + points[i,1]
        self.dataUI.modellingData.appVelocities[sourceId][receiversOrientation] = vel
        self.dataUI.modellingData.interceptTime[sourceId][receiversOrientation] = inter
        ## Updating the graphs:
        self._updateHodoGraph()
        self._updateModelGraph()

    def animateModelling(self):
        self.connectMouse = self.hodochronesGraph.mpl_connect('motion_notify_event', self.changeMouseModelling)
        self.connectPress = self.hodochronesGraph.mpl_connect('button_press_event', self.onPressModelling)
        self.connectRelease = self.hodochronesGraph.mpl_connect('button_release_event', self.onReleaseModelling)
        self.aniModelling = animation.FuncAnimation(self.hodochronesGraph.fig, self._updateModelling)# , interval=16.7)
        self.modelGraph.draw()
        self.hodochronesGraph.draw()

    def setAnimationModelling(self):
        if self.aniModelling is None:
            if self.movePoints.isChecked():
                self.animateModelling()
        else:
            if self.movePoints.isChecked():
                self.aniModelling.event_source.start()
                self.connectMouse = self.hodochronesGraph.mpl_connect('motion_notify_event', self.changeMouseModelling)
                self.connectPress = self.hodochronesGraph.mpl_connect('button_press_event', self.onPressModelling)
                self.connectRelease = self.hodochronesGraph.mpl_connect('button_release_event', self.onReleaseModelling)
            else:
                self.aniModelling.event_source.stop()
                self.hodochronesGraph.mpl_disconnect(self.connectMouse)
                self.hodochronesGraph.mpl_disconnect(self.connectPress)
                self.hodochronesGraph.mpl_disconnect(self.connectRelease)

    def sourceSelectorChanged(self, i):
        self.receiversSelector.clear()
        self.receiversSelector.addItems(self.dataUI.modellingAnimation.namesOrientations[i])
    
    def plotHodochrones(self, sensors, measurements):
        ax=self.hodochronesGraph.axes
        sources = np.unique(measurements[:,0]).astype(int)
        colors = matplotlib.pyplot.cm.tab10(np.arange(10))
        hasError = len(measurements[0,:])==4
        for i, sId in enumerate(sources):
            index = measurements[:,0].astype(int) == sId
            ti = measurements[index, 2]
            ri = sensors[measurements[index, 1].astype(int) - 1, 0]
            if hasError:
                erri = measurements[index, 3]
                ax.errorbar(ri, ti, yerr=erri, linestyle='none', color=colors[i % 10], marker='s', markersize=5)
            else:
                ax.plot(ri, ti, linestyle='none', color=colors[i % 10], marker='s', markersize=5)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Traveltime (s)')
        ax.set_ylim(bottom=0.0)
        ax.grid(True)
        self.hodochronesGraph.fig.tight_layout()
        self.hodochronesGraph.draw()

    # def _saveModel(self):
    #     # TODO
    #     self.statusBar.showMessage('Saving current model . . .')
    #     time.sleep(10)
    #     self.statusBar.showMessage(defaultStatus)

    def _loadInvMesh(self):
        # Load an inversion mesh that was already created (*.poly)
        fName, _ = QFileDialog.getOpenFileName(self,'Select file to load',filter='GIMLi mesh file (*.bms)')
        if fName != "":
            self.dataUI.invData.mesh = pg.Mesh(fName)
            self.invModelGraph.axes.cla()
            pg.show(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
            self.invModelGraph.draw()
            self.dataUI.meshLoaded = True
            self.statusBar.showMessage(f'Mesh loaded from {fName}')
        else:
            self.statusBar.showMessage('No mesh could be loaded!')

    def _loadInitModel(self):
        # Load an existing model as the starting model (*.vector) for the inversion
        fName, _ = QFileDialog.getOpenFileName(self,'Select file to load',filter='Result Vector file (*.vector)')
        if fName != "":
            self.dataUI.invData.startModel = pg.Vector(np.loadtxt(fName))
            self.invModelGraph.axes.cla()
            try:
                pg.show(self.dataUI.invData.mesh, self.dataUI.invData.startModel, ax=self.invModelGraph.axes)
                self.invModelGraph.draw()
                self.statusBar.showMessage(f'Initial model loaded from {fName}')
            except:
                QMessageBox.warning(self, 'Warning !', 'The initial model and the current mesh do not match in size!')
                self.dataUI.invData.startModel = None
                self.statusBar.showMessage('No initial model loaded!')
        else:
            self.statusBar.showMessage('No initial model loaded!')

    def _saveInvMesh(self):
        # Save the inversion mesh (*.poly)
        fName, _ = QFileDialog.getSaveFileName(self,'Select file to save',filter='GIMLi mesh file (*.bms)')
        if fName != "":
            self.dataUI.invData.mesh.save(fName)
            self.statusBar.showMessage(f'Mesh saved to {fName}')
        else:
            self.statusBar.showMessage('No mesh saved!')

    def _saveInvAsVTK(self):
        # Save the inversion results into a VTK file (for Paraview)
        fName, _ = QFileDialog.getSaveFileName(self,'Select file to save',filter='Paraview mesh file (*.vtk)')
        if fName != "":
            mgr = self.dataUI.invData.manager
            m = self.dataUI.invData.mesh
            m.addData("Velocity [m/s]", mgr.model)
            coverage = mgr.fop.jacobian().transMult(np.ones(mgr.fop.data.size()))
            m.addData("Coverage [/]", coverage)
            C = mgr.fop.constraintsRef()
            m.addData("Standardized Coverage [/]", np.sign(np.absolute(C.transMult(C * coverage))))
            m.exportVTK(fName)
            self.statusBar.showMessage(f'Results saved to {fName}')
        else:
            self.statusBar.showMessage('Results where NOT saved!')

    def _saveInvResponse(self):
        # Save the model response for the last iteration (*.vector)
        fName, _ = QFileDialog.getSaveFileName(self,'Select file to save',filter='Response Vector file (*.vector)')
        if fName != "":
            np.savetxt(fName, self.dataUI.invData.manager.inv.response)
            self.statusBar.showMessage(f'Response saved to {fName}')
        else:
            self.statusBar.showMessage('Response was NOT saved!')

    def _saveInvResult(self):
        # Save the model for the last iteration (*.vector)
        fName, _ = QFileDialog.getSaveFileName(self,'Select file to save',filter='Result Vector file (*.vector)')
        if fName != "":
            np.savetxt(fName, self.dataUI.invData.manager.model)
            self.statusBar.showMessage(f'Result saved to {fName}')
        else:
            self.statusBar.showMessage('Result was NOT saved!')
    
    def _pickTracesTabUI(self):
        importTab = QWidget(self.tabs)
        layout = QGridLayout(self.tabs) # Grid of 10-by-15

        ## Comb Box for choosing the correct file to pick.
        self.comboBoxFilesPicking = QComboBox()
        for name in self.dataUI.paths.seg2Files:
            self.comboBoxFilesPicking.addItem(name)
        layout.addWidget(self.comboBoxFilesPicking,0,0,1,15)

        ## Main graph with all the traces
        self.mainGraph = MplCanvas(importTab, width=8, height=7, dpi=100)
        self.mainGraph.axes.plot([0,1,2,3,4], [10,1,20,3,40], animated=True)
        self.mainGraphToolbar = NavigationToolbar2QT(self.mainGraph, importTab)
        mainGraphLayout = QVBoxLayout()
        mainGraphLayout.addWidget(self.mainGraphToolbar)
        mainGraphLayout.addWidget(self.mainGraph)
        layout.addLayout(mainGraphLayout,2,0,9,10)
        ## Zoom graph with the current trace
        self.zoomGraph = MplCanvas(importTab, width=5, height=4, dpi=100)
        self.zoomGraph.axes.plot([0,1,2,3,4], [10,1,20,3,40], animated=True)
        self.zoomGraph.axes.tick_params(
                axis='both',       # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                left=False,
                right=False,
                labelleft=False,
                labelbottom=False) # labels along the bottom edge are off
        self.aniMain = None
        self.aniZoom = None
        layout.addWidget(self.zoomGraph,2,10,4,5)
        self.buttonTabPickingSet = QPushButton('Set picking')
        self.buttonTabPickingReset = QPushButton('Reset picking')
        self.spinBoxCurrSelect = QSpinBox(importTab)
        self.spinBoxCurrSelect.setEnabled(False)
        self.buttonTabPickingSetT0 = QPushButton('Set t=0')
        textBox = QLabel(importTab)
        textBox.setText('Select the trace number:')
        layout.addWidget(textBox,6,10,1,5)
        layout.addWidget(self.spinBoxCurrSelect,7,10,1,5)
        layout.addWidget(self.buttonTabPickingSetT0,8,10,1,5)
        layout.addWidget(self.buttonTabPickingSet,9,10,1,5)
        layout.addWidget(self.buttonTabPickingReset,10,10,1,5)
        self.buttonTabPickingSet.setCheckable(True)
        self.buttonTabPickingSet.clicked.connect(self.setPicking)
        self.buttonTabPickingReset.clicked.connect(self.resetPicking)
        self.buttonTabPickingSetT0.setCheckable(True)
        importTab.setLayout(layout)
        return importTab
    
    def setPicking(self):
        if self.dataUI.dataLoaded:
            if self.buttonTabPickingSet.isChecked():
                # We stop the animation:
                self.aniZoom.event_source.stop()
                self.aniMain.event_source.stop()
                self.mainGraph.mpl_disconnect(self.connectMouse)
                self.mainGraph.mpl_disconnect(self.connectPress)
                self.mainGraph.mpl_disconnect(self.connectRelease)
                self.spinBoxCurrSelect.setEnabled(False)
            else:
                # We begin back the animation:
                self.aniZoom.event_source.start()
                self.aniMain.event_source.start()
                self.connectMouse = self.mainGraph.mpl_connect('motion_notify_event', self.changeMouse)
                self.connectPress = self.mainGraph.mpl_connect('button_press_event', self.onPress)
                self.connectRelease = self.mainGraph.mpl_connect('button_release_event', self.onRelease)
                self.spinBoxCurrSelect.setEnabled(True)

    def resetPicking(self):
        if self.dataUI.dataLoaded:
            self.dataUI.picking = np.empty((len(self.dataUI.paths.seg2Files),len(self.dataUI.sisData[0])))
            self.dataUI.picking[:] = np.nan 
            self.dataUI.animationPicking.changedSelect = True

    def _inversionTabUI(self):
        ## Tab for the inversion of the data using pygimli api.
        inversionTab = QWidget(self.tabs)
        layout = QGridLayout(self.tabs)
        ## Selecting the picking file:
        self.filePicksPath = QLabel('File path', inversionTab)
        self.filePicksPath.setAlignment(Qt.AlignCenter)
        self.buttonGetSgtFile = QPushButton('...', inversionTab)
        self.buttonGetSgtFile.clicked.connect(self._loadPicking)
        layout.addWidget(self.filePicksPath, 0, 0, 1, 9)
        layout.addWidget(self.buttonGetSgtFile,0, 9, 1, 1)
        self.invModelGraph = MplCanvas(inversionTab)
        invModelGraphToolbar = NavigationToolbar2QT(self.invModelGraph, inversionTab)
        invModelGraphLayout = QVBoxLayout()
        invModelGraphLayout.addWidget(invModelGraphToolbar)
        invModelGraphLayout.addWidget(self.invModelGraph)
        layout.addLayout(invModelGraphLayout, 1, 0, 5, 5)
        self.dataGraph = MplCanvas(inversionTab)
        dataGraphToolbar = NavigationToolbar2QT(self.dataGraph, inversionTab)
        dataGraphLayout = QVBoxLayout()
        dataGraphLayout.addWidget(dataGraphToolbar, alignment=Qt.AlignRight)
        dataGraphLayout.addWidget(self.dataGraph)
        layout.addLayout(dataGraphLayout, 1, 5, 5, 5)
        self.fitGraph = MplCanvas(inversionTab)
        fitGraphToolbar = NavigationToolbar2QT(self.fitGraph, inversionTab)
        fitGraphLayout = QVBoxLayout()
        fitGraphLayout.addWidget(self.fitGraph)
        fitGraphLayout.addWidget(fitGraphToolbar, alignment=Qt.AlignRight)
        layout.addLayout(fitGraphLayout, 6, 5, 5, 5)
        self.groupeOption = QGroupBox(inversionTab)
        self.groupeOption.setTitle('Inversion options')
        ## List of inversion options to enable:
        groupOptionLayout = QGridLayout(self.groupeOption) # Grid of 5 by 4
        self.setLambda = QLineEdit(str(self.dataUI.invData.lam),self.groupeOption)
        self.setLambda.setValidator(QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setLambda))
        lambdaText = QLabel('Lambda :')
        self.setZWeight = QLineEdit(str(self.dataUI.invData.zWeight),self.groupeOption)
        self.setZWeight.setValidator(QtGui.QDoubleValidator(0.01, 100.0, 3, self.setZWeight))
        zWeightText = QLabel('Z-weight (/) :')
        self.setVTop = QLineEdit(str(self.dataUI.invData.vTop),self.groupeOption)
        self.setVTop.setValidator(QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVTop))
        vTopText = QLabel('V<sub>top</sub> (m/s) :')
        self.setVBottom = QLineEdit(str(self.dataUI.invData.vBottom),self.groupeOption)
        self.setVBottom.setValidator(QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVBottom))
        vBottomText = QLabel('V<sub>bottom</sub> (m/s) :')
        # self.loadInitModel = QPushButton('Load Initial Model', self.groupeOption)
        minVText = QLabel('Min. velocity (m/s) :')
        self.setVMin = QLineEdit(str(self.dataUI.invData.vMin),self.groupeOption)
        self.setVMin.setValidator(QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVMin))
        maxVText = QLabel('Max. velocity (m/s) :')
        self.setVMax = QLineEdit(str(self.dataUI.invData.vMax),self.groupeOption)
        self.setVMax.setValidator(QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVMax))
        maxCellText = QLabel('Mesh max. cell size (m??) :')
        self.setMaxCell = QLineEdit(str(self.dataUI.invData.meshMaxCellSize),self.groupeOption)
        self.setMaxCell.setValidator(QtGui.QDoubleValidator(0.1, 100.0, 2, self.setMaxCell))
        maxDepthText = QLabel('Mesh max. depth (m) :')
        self.setMaxDepth = QLineEdit(str(self.dataUI.invData.meshDepthMax),self.groupeOption)
        self.setMaxDepth.setValidator(QtGui.QDoubleValidator(5.0, 1000.0, 2, self.setMaxDepth))
        self.runInversion = QPushButton('Run inversion', self.groupeOption)
        groupOptionLayout.addWidget(lambdaText, 0, 0, 1, 1)
        groupOptionLayout.addWidget(self.setLambda,0, 1, 1, 1)
        groupOptionLayout.addWidget(zWeightText, 0, 2, 1, 1)
        groupOptionLayout.addWidget(self.setZWeight, 0, 3, 1, 1)
        groupOptionLayout.addWidget(vTopText, 1, 0, 1, 1)
        groupOptionLayout.addWidget(self.setVTop, 1, 1, 1, 1)
        groupOptionLayout.addWidget(vBottomText, 1, 2, 1, 1)
        groupOptionLayout.addWidget(self.setVBottom, 1, 3, 1, 1)
        # groupOptionLayout.addWidget(self.loadInitModel, 2, 0, 1, 4)
        # self.loadInitModel.clicked.connect(self._setStartModel)
        groupOptionLayout.addWidget(minVText, 2, 0, 1, 1)
        groupOptionLayout.addWidget(self.setVMin, 2, 1, 1, 1)
        groupOptionLayout.addWidget(maxVText, 2, 2, 1, 1)
        groupOptionLayout.addWidget(self.setVMax, 2, 3, 1, 1)
        groupOptionLayout.addWidget(maxCellText, 3, 0, 1, 1)
        groupOptionLayout.addWidget(self.setMaxCell, 3, 1, 1, 1)
        groupOptionLayout.addWidget(maxDepthText, 3, 2, 1, 1)
        groupOptionLayout.addWidget(self.setMaxDepth, 3, 3, 1, 1)
        groupOptionLayout.addWidget(self.runInversion, 4, 0, 1, 4)
        self.runInversion.clicked.connect(self._runInversion)
        self.groupeOption.setLayout(groupOptionLayout)
        # - Lambda ('lam')
        # - InitialModel
        # - ErrorModel
        # Button for running the inversion
        layout.addWidget(self.groupeOption, 6, 0, 5, 5)
        # self.inversionText = QTextEdit(inversionTab)
        # layout.addWidget(self.inversionText, 10, 0, 1, 5)
        inversionTab.setLayout(layout)
        return inversionTab
    
    # def _setStartModel(self):
    #     pass

    def _runInversion(self):
        # Parameters for inversion:
        self.dataUI.invData.lam = float(self.setLambda.text())
        self.dataUI.invData.zWeight = float(self.setZWeight.text())
        self.dataUI.invData.vTop = float(self.setVTop.text())
        self.dataUI.invData.vBottom = float(self.setVBottom.text())
        self.dataUI.invData.vMin = float(self.setVMin.text())
        self.dataUI.invData.vMax = float(self.setVMax.text())
        # Mesh model and start model:
        self.dataUI.invData.meshMaxCellSize = float(self.setMaxCell.text())
        self.dataUI.invData.meshDepthMax = float(self.setMaxDepth.text())
        if self.dataUI.invData.data is not None:
            if self.dataUI.inversionDone:
                self.dataUI.invData.manager = TTMgr(data = self.dataUI.invData.data)
            # Creating mesh with mesh parameters:
            # self.dataUI.invData.data = pg.DataContainer(self.filePicksPath.text())
            # self.dataUI.invData.manager = TTMgr(data=self.dataUI.invData.data)
            if self.dataUI.meshLoaded == False or self.dataUI.invData.mesh == None:
                self.dataUI.invData.mesh = self.dataUI.invData.manager.createMesh(data=self.dataUI.invData.data, paraMaxCellSize=self.dataUI.invData.meshMaxCellSize, paraDepth=self.dataUI.invData.meshDepthMax)
                mesh = self.dataUI.invData.mesh
            pgshow(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
            self.invModelGraph.draw()
            if (self.dataUI.invData.startModel is None): # and (self.dataUI.invData.startModel.shape()[0]):
                self.dataUI.invData.setStartModelGradient(data=self.dataUI.invData.data, mesh=self.dataUI.invData.mesh)
            # Running the inversion
            self.dataUI.invData.manager.invert(data = self.dataUI.invData.data,
                                               mesh = self.dataUI.invData.mesh,
                                               zWeight = self.dataUI.invData.zWeight,
                                               lam = self.dataUI.invData.lam,
                                               startModel = self.dataUI.invData.startModel,
                                               limits = [self.dataUI.invData.vMin, self.dataUI.invData.vMax],
                                               verbose = False)
            self.invModelGraph.fig.clear()
            self.invModelGraph.axes = self.invModelGraph.fig.add_subplot(111)
            self.fitGraph.axes.cla()
            drawFirstPicks(ax=self.fitGraph.axes, data=self.dataUI.invData.data, tt=(np.abs(np.asarray(self.dataUI.invData.data('t')-np.asarray(self.dataUI.invData.manager.inv.response)))/np.asarray(self.dataUI.invData.data('t')))*100)
            self.fitGraph.axes.set_xlabel('X (m)')
            self.fitGraph.axes.set_ylabel('Data misfit (%)')
            self.fitGraph.fig.tight_layout()
            self.fitGraph.draw()
            _, cBar = self.dataUI.invData.manager.showResult(ax=self.invModelGraph.axes, cMap='cividis')
            self.invModelGraph.axes.set_title('Inversion result (chi?? = {: .2f}, RMS = {: .2f} %, RRMS = {: .2f} %)'.format(self.dataUI.invData.manager.inv.chi2(), self.dataUI.invData.manager.inv.absrms()*100, self.dataUI.invData.manager.inv.relrms()*100))
            self.invModelGraph.cBar = cBar
            self.dataUI.invData.manager.drawRayPaths(ax=self.invModelGraph.axes, color='w', lw=0.3, alpha=0.5)
            self.invModelGraph.fig.tight_layout()
            self.invModelGraph.draw()
            self.dataUI.invData.startModel = None
            self.dataUI.inversionDone = True
            self.dataUI.meshLoaded = False

    def _modelTabUI(self):
        '''In this tab , we will propose to draw the hodochrones on top
         of the picking and build the corresponding layered model.
        Those models are build using the intercept time-method.
        '''
        modellingTab = QWidget()
        # Pick loading
        self.filePicksPath = QLabel('File path', modellingTab)
        self.filePicksPath.setAlignment(Qt.AlignCenter)
        self.buttonGetSgtFile = QPushButton('...', modellingTab)
        self.buttonGetSgtFile.clicked.connect(self._loadPicking)
        # Graph with the hodochrones:
        self.hodochronesGraph = MplCanvas(modellingTab)
        hodochronesToolbar = NavigationToolbar2QT(self.hodochronesGraph, modellingTab)
        hodochronesWidget = QVBoxLayout()
        hodochronesWidget.addWidget(hodochronesToolbar, alignment=Qt.AlignLeft)
        hodochronesWidget.addWidget(self.hodochronesGraph)
        # Graph with the model
        self.modelGraph = MplCanvas(modellingTab)
        modelToolbar = NavigationToolbar2QT(self.modelGraph, modellingTab)
        modelWidget = QVBoxLayout()
        modelWidget.addWidget(modelToolbar, alignment=Qt.AlignRight)
        modelWidget.addWidget(self.modelGraph)
        # Options for the tab:
        self.groupOptionModelling = QGroupBox()
        self.groupOptionModelling.setTitle('Options')
        self.sourceSelector = QComboBox(self.groupOptionModelling)
        self.receiversSelector = QComboBox(self.groupOptionModelling)
        labelSourceReceiver = QLabel('Select source and orientation :',self.groupOptionModelling)
        layoutOptions = QGridLayout()
        layoutOptions.addWidget(labelSourceReceiver, 1, 0, 1, 4, alignment=Qt.AlignRight)
        layoutOptions.addWidget(self.sourceSelector, 1, 4, 1, 4)
        layoutOptions.addWidget(self.receiversSelector, 1, 8, 1, 2)
        self.movePoints = QPushButton('Move Points', self.groupOptionModelling)
        self.movePoints.setCheckable(True)
        self.movePoints.clicked.connect(self.setAnimationModelling)
        self.nbLayersSelector = QSpinBox(self.groupOptionModelling)
        self.nbLayersSelector.setMaximum(3)
        self.nbLayersSelector.setMinimum(2)
        self.nbLayersSelector.setValue(self.dataUI.modellingData.nbLayers)
        self.nbLayersSelector.valueChanged.connect(self._initModelling)
        nbLayersLabel = QLabel('Select the number of layers :')
        layoutOptions.addWidget(nbLayersLabel, 0, 0, 1, 5, alignment=Qt.AlignRight)
        layoutOptions.addWidget(self.nbLayersSelector, 0, 5, 1, 5)
        layoutOptions.addWidget(self.movePoints, 2, 0, 1, 10)
        self.groupOptionModelling.setLayout(layoutOptions)
        # Setup of the layout:
        layout = QGridLayout(modellingTab)
        layout.addWidget(self.filePicksPath, 0, 0, 1, 9)
        layout.addWidget(self.buttonGetSgtFile,0, 9, 1, 1)
        layout.addLayout(hodochronesWidget, 1, 0, 7, 5)
        layout.addLayout(modelWidget, 1, 5, 7, 5)
        layout.addWidget(self.groupOptionModelling, 9, 0, 2, 5)
        modellingTab.setLayout(layout)
        self.aniModelling = None
        return modellingTab
    
    def saveStatePicking(self):
        fName, _ = QFileDialog.getSaveFileName(self,'Select file to save',filter='Pickled structure (*.pkl)')
        if fName != "":
            dataLoaded = picklingStatus()
            dataLoaded.paths = self.dataUI.paths
            dataLoaded.geometry = self.dataUI.geometry
            dataLoaded.sisData = self.dataUI.sisData
            dataLoaded.beginTime = self.dataUI.beginTime
            dataLoaded.picking = self.dataUI.picking
            dataLoaded.pickingError = self.dataUI.pickingError
            dataLoaded.sisFileId = self.dataUI.sisFileId
            file = open(fName, 'wb')
            pickle.dump(dataLoaded,file)
            file.close()
            self.statusBar.showMessage(f'Result saved to {fName}')
        else:
            self.statusBar.showMessage('Result was NOT saved!')

    def loadStatePicking(self):
        fName, _ = QFileDialog.getOpenFileName(self,'Select file to load',filter='Pickled structure (*.pkl)')
        if fName != "":
            file = open(fName, 'rb')
            dataLoaded = pickle.load(file)
            file.close()
            # Assign the different elements :
            if len(dataLoaded.paths.directory) != 0:
                if self.dataUI.dataLoaded:
                    self.aniZoom.event_source.stop()
                    self.aniMain.event_source.stop()
                    self.mainGraph.mpl_disconnect(self.connectMouse)
                    self.mainGraph.mpl_disconnect(self.connectPress)
                    self.mainGraph.mpl_disconnect(self.connectRelease)
                self.dataUI.paths = dataLoaded.paths
                self.dataUI.geometry = dataLoaded.geometry
                self.dataUI.sisData = dataLoaded.sisData
                self.dataUI.beginTime = dataLoaded.beginTime
                self.dataUI.picking = dataLoaded.picking
                self.dataUI.pickingError = dataLoaded.pickingError
                self.dataUI.sisFileId = dataLoaded.sisFileId
                self.dataUI.dataLoaded = True
                # Change the interface:
                self.spinBoxCurrSelect.setMinimum(0)
                self.spinBoxCurrSelect.setMaximum(len(self.dataUI.sisData[0])-1)
                self.spinBoxCurrSelect.setPrefix('Trace number ')
                self.spinBoxCurrSelect.valueChanged.connect(self.traceNumberChanged)
                self.spinBoxCurrSelect.setEnabled(True)
                self.dataUI.animationPicking.changedSelect = True
                self.updateTab0()
            else:
                self.statusBar.showMessage('Empty strate loaded')
        else:
            self.statusBar.showMessage('No status loaded!')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = Window()
    window.show()
    sys.exit(app.exec_())