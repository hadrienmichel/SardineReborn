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

defaultStatus = "Idle."

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        fig.tight_layout()
        super(MplCanvas, self).__init__(fig)

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
    
    def _openGeometry(self):
        self.statusBar.showMessage('Openning Geometry file . . .')
        ## Opening the geometry file:
        fname, _ = QFileDialog.getOpenFileName(self,'Open geometry file',filter='Geometry file (*.geometry)')# The first argument returned is the filename and path
        headTail = os.path.split(fname)
        path = headTail[0]
        file = headTail[1]
        nameSave = file[:-9] # remove the *.geometry
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
        self.updateTab()
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
            elif ext == '.seg2' or ext == '.sg2':
                st = read(os.path.join(path,name), 'SEG2')
            else:
                st = read(os.path.join(path,name), 'SEG2')
            if len(st) != len(ReceiversPosition): # If the number of geophones does not match between the loaded array and the gemoetry
                raise Exception('The file referenced in the geometry file does not match the geometry of the array!')
            self.dataUI.sisData.append(st)

    def updateTab(self):
        self.comboBoxFilesPicking.clear()
        for name in self.dataUI.paths.seg2Files:
            self.comboBoxFilesPicking.addItem(name)

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
        self.mainGraph.axes.plot([0,1,2,3,4], [10,1,20,3,40])
        # plt.tight_layout()
        self.mainGraphToolbar = NavigationToolbar2QT(self.mainGraph, importTab)
        mainGraphLayout = QVBoxLayout()
        mainGraphLayout.addWidget(self.mainGraphToolbar)
        mainGraphLayout.addWidget(self.mainGraph)
        layout.addLayout(mainGraphLayout,2,0,9,10)
        # self.mainGraph.

        ## Zoom graph with the current trace
        self.zoomGraph = MplCanvas(importTab, width=5, height=4, dpi=100)
        self.zoomGraph.axes.plot([0,1,2,3,4], [10,1,20,3,40])
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