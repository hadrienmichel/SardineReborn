import sys
from PyQt5.QtWidgets import (
    QApplication,  
    QWidget,
    QTabWidget,
    QHBoxLayout,
    QCheckBox,
    QAction,
    QMenu,
    QMainWindow)

class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Sardine Reborn')
        self.resize(600,350)

        # Adding tabs to the layout
        tabs = QTabWidget(self)
        tabs.addTab(self._pickTracesTabUI(), 'Traces picking')
        tabs.addTab(self._inversionTabUI(), 'Inversion')
        tabs.addTab(self._modelTabUI(), 'Model')

        self.setCentralWidget(tabs) 

        # Defining menu actions:
        openFile = QAction("&Open Geometry File",self)
        openFile.setShortcut('Ctrl+O') # Setting ctrl+o as the shortcut
        openFile.setStatusTip('Open the *.geometry file')
        openFile.triggered.connect(self._openGeometry)

        savePicking = QAction("Save Current &Picking",self)
        savePicking.setShortcut('Ctrl+P')
        savePicking.setStatusTip('Save the picking state in a *.sgt file')
        savePicking.triggered.connect(self._savePicking)

        saveModel = QAction("Save Current &Model",self)
        saveModel.setShortcut('Ctrl+M')
        saveModel.setStatusTip('Save the current model in a *.txt file')
        saveModel.triggered.connect(self._saveModel)
        # Adding the menu bar atop
        menuBarInternal = self.menuBar()
        menuBarInternal.setNativeMenuBar(False)
        fileMenu = menuBarInternal.addMenu("&File")
        fileMenu.addAction(openFile)
        fileMenu.addAction(savePicking)
        fileMenu.addAction(saveModel)

    
    def _openGeometry(self):
        print('Open')

    def _savePicking(self):
        print('Save Picking')

    def _saveModel(self):
        print('Save Model')
    
    def _pickTracesTabUI(self):
        importTab = QWidget()
        layout = QHBoxLayout()
        layout.addWidget(QCheckBox('Option 1'))
        layout.addWidget(QCheckBox('Option 2'))
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