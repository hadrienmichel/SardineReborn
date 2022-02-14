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
## Imports for the seismic data input
from obspy import read
## Imports for the data inversion
import pygimli as pg
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
    QTextEdit)
from PyQt5.QtCore import Qt

defaultStatus = "Idle."

## Need to take a closer look at this: https://programmerall.com/article/10751929193/
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.fig.tight_layout()
        super(MplCanvas, self).__init__(self.fig)
        self.setParent(parent)
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
class dataStorage():
    def __init__(self) -> None:
        ## Data variables:
        self.paths = paths()
        self.geometry = geometry()
        self.sisData = []
        self.picking = []
        self.model = model()
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
        self.resize(600,350)

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
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText('Welcome to Sardine Reborn!')
        msgBox.setInformativeText('by Hadrien Michel (2022)')
        msgBox.setWindowTitle('Welcome!')
        msgBox.show()
        event.accept()
    
    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Closing ...', 'Are you sure you want to quit?', QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Ok)
        if reply == QMessageBox.Ok:
            event.accept()
        else:
            event.ignore()
        
    ## Animations definitions:
    # PyQt5 animations
    def keyPressEvent(self, event):
        if self.tabs.currentIndex() == 0: # Only used if the key is pressed in the first tab
            if event.key() == Qt.Key_Up: # Change i += 1
                self.dataUI.animationPicking.currSelect += 1
                if self.dataUI.animationPicking.currSelect >= len(self.dataUI.sisData[0]):
                    self.dataUI.animationPicking.currSelect = 0
                self.dataUI.animationPicking.changedSelect = True
            elif event.key() == Qt.Key_Down: # Change i -= 1
                self.dataUI.animationPicking.currSelect -= 1
                if self.dataUI.animationPicking.currSelect < 0:
                    self.dataUI.animationPicking.currSelect = len(self.dataUI.sisData[0])-1
                self.dataUI.animationPicking.changedSelect = True    
            event.accept()
        else:
            event.ignore()
    
    # Matplotlib animations:
    def changeMouse(self, event):
        if event.inaxes is not None:
            self.dataUI.animationPicking.mousePosition = event.xdata
        return 0

    def onPress(self, event):
        if event.button == 1:
            self.dataUI.animationPicking.timeOnClick = time.time()
        return 0

    def onRelease(self, event):
        if event.button == 1 and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength): # If left click and not dragging accross the pannel
            self.dataUI.picking[self.dataUI.sisFileId][self.dataUI.animationPicking.currSelect] = self.dataUI.animationPicking.mousePosition
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
        if self.dataUI.picking[self.dataUI.sisFileId][self.dataUI.animationPicking.currSelect] is not None:
            axZoom.plot([self.dataUI.picking[self.dataUI.sisFileId][self.dataUI.animationPicking.currSelect], self.dataUI.picking[self.dataUI.sisFileId][self.dataUI.animationPicking.currSelect]], z,color='g')
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
                if self.dataUI.picking[self.dataUI.sisFileId][i] is not None:
                    axMain.plot([self.dataUI.picking[self.dataUI.sisFileId][i], self.dataUI.picking[self.dataUI.sisFileId][i]], [i-0.5, i+0.5],color='g')
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
        
        self.dataUI.picking = [[None]*len(self.dataUI.sisData[0])]*len(SEG2Files)       

    def updateTab0(self):
        self.comboBoxFilesPicking.clear()
        for name in self.dataUI.paths.seg2Files:
            self.comboBoxFilesPicking.addItem(name)
        self.dataUI.sisFileId = 0
        self.comboBoxFilesPicking.currentIndexChanged.connect(self.comboBoxChange)
        self.mainGraph.mpl_connect('motion_notify_event', self.changeMouse)
        # self.mainGraph.mpl_connect('key_press_event', self.onKey)
        self.mainGraph.mpl_connect('button_press_event', self.onPress)
        self.mainGraph.mpl_connect('button_release_event', self.onRelease)
        self.aniMain = animation.FuncAnimation(self.mainGraph.fig, self.animationMain, interval=16.7)
        self.aniZoom = animation.FuncAnimation(self.zoomGraph.fig, self.animationZoom, interval=16.7)
        self.mainGraph.draw()
        self.zoomGraph.draw()

    def comboBoxChange(self, newId):
        self.dataUI.sisFileId = newId
        self.dataUI.animationPicking.changedSelect = True

    def _savePicking(self):
        self.statusBar.showMessage('Save current picking . . .')
        time.sleep(10)
        self.statusBar.showMessage(defaultStatus)

    def _loadPicking(self):
        self.statusBar.showMessage('Loading picking file . . .')
        time.sleep(10)
        self.statusBar.showMessage(defaultStatus)

    def _saveModel(self):
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
        # plt.tight_layout()
        self.mainGraphToolbar = NavigationToolbar2QT(self.mainGraph, importTab)
        mainGraphLayout = QVBoxLayout()
        mainGraphLayout.addWidget(self.mainGraphToolbar)
        mainGraphLayout.addWidget(self.mainGraph)
        layout.addLayout(mainGraphLayout,2,0,9,10)
        
        # self.mainGraph.connect

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
        # plt.tight_layout()
        self.aniMain = None
        self.aniZoom = None
        layout.addWidget(self.zoomGraph,2,10,5,5)
        layout.addWidget(QTextEdit(),7,10,2,5)
        layout.addWidget(QPushButton('Option 1'),9,10,1,5)
        layout.addWidget(QPushButton('Option 2'),10,10,1,5)
        importTab.setLayout(layout)
        return importTab

    def _inversionTabUI(self):
        importTab = QWidget()
        layout = QHBoxLayout()
        layout.addWidget(QCheckBox('Option 1'))
        layout.addWidget(QCheckBox('Option 2'))
        importTab.setLayout(layout)
        return importTab

    def _modelTabUI(self):
        importTab = QWidget()
        layout = QHBoxLayout()
        layout.addWidget(QCheckBox('Option 1'))
        layout.addWidget(QCheckBox('Option 2'))
        importTab.setLayout(layout)
        return importTab
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = Window()
    window.show()
    sys.exit(app.exec_())