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
## Imports for the data inversion
import pygimli as pg
from pygimli.physics import TravelTimeManager as TTMgr
from pygimli.physics.traveltime import ratools, drawFirstPicks
from pygimli.viewer import show as pgshow
## Imports for the GUI
from PyQt5.QtWidgets import (
    QApplication,  
    QWidget,
    QTabWidget,
    QHBoxLayout,
    QVBoxLayout,
    QGridLayout,
    QCheckBox,
    QAction,
    QMainWindow,
    QStatusBar,
    QMessageBox,
    QFileDialog,
    QComboBox,
    QPushButton,
    QSpinBox,
    QLabel,
    QGroupBox,
    QLineEdit)
from PyQt5.QtCore import Qt
from PyQt5 import QtGui

defaultStatus = "Idle."

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
        self.drawingLine = False
        self.movingPoint = False
        self.currPosition = [0, 0]
        self.maxClickLength = 0.5
class modellingData():
    def __init__(self) -> None:
        self.hodo = []
        self.model = []
        self.nbLayers = []
        self.appVelocities = []
        self.interceptTime = []
class dataStorage():
    def __init__(self) -> None:
        ## Data variables:
        self.paths = paths()
        self.geometry = geometry()
        self.sisData = []
        self.picking = []
        self.pickingError = []
        self.model = model()
        self.invData = inversionData()
        ## Status variables:
        self.dataLoaded = False
        self.pickingDone = False
        self.inversionDone = False
        self.sisFileId = 0
        self.beginTime = 0
        ## Graphical animation variables:
        self.animationPicking = animationPicking()
class Window(QMainWindow):
    def __init__(self) -> None:
        super().__init__()

        ## Initializing the data structure
        self.dataUI = dataStorage()

        ## Initialize the UI
        self.setWindowTitle('Sardine Reborn')
        self.setWindowIcon(QtGui.QIcon('SardineRebornLogo.png'))
        self.resize(1100,600)

        ## Adding tabs to the layout
        self.tabs = QTabWidget(self)
        self.tabs.addTab(self._pickTracesTabUI(), 'Traces picking')
        self.tabs.addTab(self._inversionTabUI(), 'Inversion')
        self.tabs.addTab(self._modelTabUI(), 'Model')

        self.setCentralWidget(self.tabs) 

        ## Defining menu actions:
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

        saveModel = QAction("Save Current &Model",self)
        saveModel.setShortcut('Ctrl+M')
        saveModel.setStatusTip('Save the current model in a *.txt file')
        saveModel.triggered.connect(self._saveModel)
        # Adding the menu bar atop
        menuBarInternal = self.menuBar()
        menuBarInternal.setNativeMenuBar(False)
        fileMenu = menuBarInternal.addMenu("File")
        fileMenu.addAction(openFile)
        fileMenu.addAction(savePicking)
        fileMenu.addAction(loadPicking)
        fileMenu.addAction(saveModel)

        ## Defining the status bar:
        self.statusBar = QStatusBar(self)
        self.statusBar.showMessage(defaultStatus)
        self.setStatusBar(self.statusBar)

        pass

    ## Oppening and closing message boxes:
    def showEvent(self, event):
        msgBox = QMessageBox(self)
        msgBox.setIconPixmap(QtGui.QPixmap('SardineRebornLogo.png').scaled(200,100,aspectRatioMode=Qt.KeepAspectRatio))
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

    def onPress(self, event):
        if event.button == MouseButton.LEFT or event.button == MouseButton.RIGHT:
            self.dataUI.animationPicking.timeOnClick = time.time()
        return 0

    def onRelease(self, event):
        if event.button == MouseButton.LEFT and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength): # If left click and not dragging accross the pannel
            if self.dataUI.animationPicking.mousePosition < 0: # To remove a picked trace, click on times below 0
                self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.nan
                self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.nan
            else:
                self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = self.dataUI.animationPicking.mousePosition
                self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = self.dataUI.animationPicking.mousePosition*0.03 # Default error is 3%
            self.dataUI.animationPicking.changedSelect = True
        if event.button == MouseButton.RIGHT and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength): # If right click and not dragging accross the pannel
            if not(np.isnan(self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect])):
                if self.dataUI.animationPicking.mousePosition > 0:
                    self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.abs(self.dataUI.animationPicking.mousePosition-self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect])/self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] # Default error is 3%
                    self.dataUI.animationPicking.changedSelect = True
        return 0
    
    def animationZoom(self, i):
        axZoom = self.zoomGraph.axes
        # Get axis variables
        deltaT = float(self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        timeSEG2 = np.arange(self.dataUI.beginTime, self.dataUI.beginTime+nbPoints*deltaT, deltaT)
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
        timeSEG2 = np.arange(self.dataUI.beginTime, self.dataUI.beginTime+nbPoints*deltaT, deltaT)
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
            if line.startswith("SOURCES"):
                sources = True
                receivers = False
            elif line.startswith("RECEIVERS"):
                sources = False
                receivers = True
            else:
                if sources:
                    tmp = re.split(r'\t+', line)
                    name = tmp[0]
                    CurrSource = [float(i) for i in tmp[1:]]
                    SEG2Files.append(name)
                    SourcePosition.append(CurrSource)
                elif receivers:
                    CurrReceiver = [float(i) for i in re.split(r'\t+', line)]
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
                self.dataUI.beginTime = 0
            elif ext == '.seg2' or ext == '.sg2':
                st = read(os.path.join(path,name), 'SEG2')
                self.dataUI.beginTime = float(st[0].stats.seg2["DELAY"])
            else:
                st = read(os.path.join(path,name), 'SEG2')
                self.dataUI.beginTime = float(st[0].stats.seg2["DELAY"])
            if len(st) != len(ReceiversPosition): # If the number of geophones does not match between the loaded array and the gemoetry
                raise Exception('The file referenced in the geometry file does not match the geometry of the array!')
            self.dataUI.sisData.append(st)
        
        self.dataUI.picking = np.empty((len(self.dataUI.paths.seg2Files),len(self.dataUI.sisData[0])))
        self.dataUI.picking[:] = np.nan
        self.dataUI.pickingError = np.empty((len(self.dataUI.paths.seg2Files),len(self.dataUI.sisData[0])))
        self.dataUI.pickingError[:] = np.nan
        self.spinBoxCurrSelect.setMinimum(0)
        self.spinBoxCurrSelect.setMaximum(len(self.dataUI.sisData[0])-1)
        self.spinBoxCurrSelect.setPrefix('Trace number ')
        self.spinBoxCurrSelect.valueChanged.connect(self.traceNumberChanged)
        self.spinBoxCurrSelect.setEnabled(True)
    
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

    def _loadPicking(self):
        self.statusBar.showMessage('Loading picking file . . .')
        # 1) Load a file with the first arrival:
        fname, _ = QFileDialog.getOpenFileName(self,'Select file to load',filter='Source-Receiver-Time file (*.sgt)')
        # We retreived a first-arrival file --> geometry of the sensors + first arrivals 
        with open(fname) as f:
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
        self.filePicksPath.setText(fname)
        self._initPygimli(fname)
        self._initModelling(sensors, measurements)
    
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

    def _initModelling(self, sensors, measurements):
        
        pass

    def _saveModel(self):
        # TODO
        self.statusBar.showMessage('Saving current model . . .')
        time.sleep(10)
        self.statusBar.showMessage(defaultStatus)
    
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
        layout.addWidget(self.zoomGraph,2,10,5,5)
        self.buttonTabPickingSet = QPushButton('Set picking')
        self.buttonTabPickingReset = QPushButton('Reset picking')
        self.spinBoxCurrSelect = QSpinBox(importTab)
        self.spinBoxCurrSelect.setEnabled(False)
        textBox = QLabel(importTab)
        textBox.setText('Select the trace number:')
        layout.addWidget(textBox,7,10,1,5)
        layout.addWidget(self.spinBoxCurrSelect,8,10,1,5)
        layout.addWidget(self.buttonTabPickingSet,9,10,1,5)
        layout.addWidget(self.buttonTabPickingReset,10,10,1,5)
        self.buttonTabPickingSet.setCheckable(True)
        self.buttonTabPickingSet.clicked.connect(self.setPicking)
        self.buttonTabPickingReset.clicked.connect(self.resetPicking)
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
        maxCellText = QLabel('Mesh max. cell size (m²) :')
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
        # - ReferenceModel --> impossible with the invert method used for TT --> line 422 of inversion.py --> self.inv.setReferenceModel(self.startModel)!
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
            self.dataUI.invData.mesh = self.dataUI.invData.manager.createMesh(data=self.dataUI.invData.data, paraMaxCellSize=self.dataUI.invData.meshMaxCellSize, paraDepth=self.dataUI.invData.meshDepthMax)
            pgshow(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
            self.invModelGraph.draw()
            # if (self.dataUI.invData.startModel is None): # and (self.dataUI.invData.startModel.shape()[0]):
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
            drawFirstPicks(ax=self.fitGraph.axes, data=self.dataUI.invData.data, tt=np.abs(np.asarray(self.dataUI.invData.data('t')-np.asarray(self.dataUI.invData.manager.inv.response)))/np.asarray(self.dataUI.invData.data('t'))*100)
            self.fitGraph.axes.set_xlabel('x (m)')
            self.fitGraph.axes.set_ylabel('Data misfit (%)')
            self.fitGraph.fig.tight_layout()
            self.fitGraph.draw()
            _, cBar = self.dataUI.invData.manager.showResult(ax=self.invModelGraph.axes, cMap='cividis')
            self.invModelGraph.cBar = cBar
            self.dataUI.invData.manager.drawRayPaths(ax=self.invModelGraph.axes, color='w', lw=0.3, alpha=0.5)
            self.invModelGraph.fig.tight_layout()
            self.invModelGraph.draw()
            self.dataUI.inversionDone = True

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
        # Setup of the layout:
        layout = QGridLayout(modellingTab)
        layout.addWidget(self.filePicksPath, 0, 0, 1, 9)
        layout.addWidget(self.buttonGetSgtFile,0, 9, 1, 1)
        layout.addLayout(hodochronesWidget, 1, 0, 7, 5)
        layout.addLayout(modelWidget, 1, 5, 7, 5)
        layout.addWidget(self.groupOptionModelling, 9, 0, 2, 10)
        modellingTab.setLayout(layout)
        return modellingTab

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = Window()
    window.show()
    sys.exit(app.exec_())