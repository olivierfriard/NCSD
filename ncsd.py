#!/usr/bin/env python

"""
Neurolucida Companion for Spatial Distribution Analysis (NCSD)

a GUI interface (with PyQt5)
Copyright Olivier Friard 2010-2011
Universita' di Torino

This file is part of "Neurolucida Companion for Spatial Distribution Analysis".

"Neurolucida Companion for Spatial Distribution Analysis" is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or any later version.

"Neurolucida Companion for Spatial Distribution Analysis" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "Neurolucida Companion for Spatial Distribution Analysis"; see the file COPYING.TXT.  If not, write to
the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
"""

__version__ = '7'

"""
- 2011-05-10: improved functions error
- 2011-04-28: added ncsd_distance module 
- 2011-03-09: added "Open directory" function
- 2011-02-23: changed output file extension to .txt (from .tsv)
"""

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtSvg import *

from ncsd_ui import Ui_MainWindow
import sys
import os



class MainWindow(QMainWindow, Ui_MainWindow):

    filename = ''
    directoryName = ''

    def __init__(self, parent=None):
        '''
        QWidget.__init__(self, parent)
        self.ui = ncsd_ui.Ui_MainWindow()
        self.ui.setupUi(self)
        '''
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        
        #QObject.connect(self.pbMarkersAnalysis, SIGNAL("clicked()"), self.pbMarkersAnalysis_clicked)
        self.pbMarkersAnalysis.clicked.connect(self.pbMarkersAnalysis_clicked)

        #QObject.connect(self.pbCellBodiesAnalysis, SIGNAL("clicked()"), self.pbCellBodiesAnalysis_clicked)
        self.pbCellBodiesAnalysis.clicked.connect(self.pbCellBodiesAnalysis_clicked)

        #QObject.connect(self.pbDistanceAnalysis, SIGNAL("clicked()"), self.pbDistanceAnalysis_clicked)
        self.pbDistanceAnalysis.clicked.connect(self.pbDistanceAnalysis_clicked)

        #QObject.connect(self.label, SIGNAL("clicked()"), self.openXML)
        #self.label.clicked.connect(self.openXML)
        

        # menu items
        # QObject.connect(self.actionOpen_Neurolucida_XML_file, SIGNAL("activated()"), self.openXML)
        self.actionOpen_Neurolucida_XML_file.triggered.connect(self.openXML)
        # QObject.connect(self.actionOpenAll, SIGNAL("activated()"), self.openAll)
        self.actionOpenAll.triggered.connect(self.openAll)
        # QObject.connect(self.actionQuit, SIGNAL("activated()"), self.quit)
        self.actionQuit.triggered.connect(self.quit)
        # QObject.connect(self.actionAbout, SIGNAL("activated()"), self.about)
        self.actionAbout.triggered.connect(self.about)

    def pbMarkersAnalysis_clicked(self):
        '''Analyze markers from XML file'''

        try:
            import ncsd_markers
        except ImportError:
            QMessageBox.warning(self, 'NCSD', '"ncsd_markers" module not found!')
            return

        if (self.filename == '') and (self.directoryName == ''):
            QMessageBox.warning(self, 'NCSD', 'Select a XML file or a directory!')
            return

        self.statusbar.showMessage("Markers analysis", 0)

        angle = self.sbAngle.value()

        if self.filename:
            files_list = [str(self.filename)]
        if self.directoryName:
            import glob
            files_list = glob.glob(str(self.directoryName) + '/*.xml')

            if not files_list:
                QMessageBox.warning(self, 'NCSD', 'The directory %s does not contain any XML file' % str(self.directoryName))
                return

        flagError = False
        for fileName in files_list:

            outputFileName = fileName.replace('.xml', '.txt')

            print(fileName, outputFileName)
            r, txt = ncsd_markers.main(angle, False, (self.cbSVG.checkState() == Qt.Checked), fileName, outputFileName)

            if r:
                flagError = True
                if txt:
                    QMessageBox.warning(self, "Error", txt)
                else:
                    QMessageBox.warning(self, "Error", 'Error #%d' % r)

        if not flagError:
            QMessageBox.about(self, "Markers analysis", 'Analysis done successfully')

        self.statusbar.showMessage("Done", 2000)

    def pbCellBodiesAnalysis_clicked(self):
        '''Analyze cell bodies from XML file'''

        try:
            import ncsd_cellbodies
        except ImportError:
            QMessageBox.warning(self, 'NCSD', '"ncsd_cellbodies" module not found!')
            return

        if self.filename == '':
            QMessageBox.warning(self, 'NCSD', 'You do not select any XML file!')
            return


        self.statusbar.showMessage("Cell bodies analysis", 0)

        angle = self.sbAngle.value()

        ### reference contour
        ref_contour = str(self.leRefContour_2.text())


        ### center contour
        center_contour = str(self.le_centerContour_2.text())


        if self.filename:
            files_list = [str(self.filename)]
        if self.directoryName:
            import glob
            files_list = glob.glob(str(self.directoryName) + '/*.xml')

            if not files_list:
                QMessageBox.warning(self, 'NCSD', 'The directory %s does not contain any XML file' % str(self.directoryName))
                return

        flagError = False
        for fileName in files_list:

            outputFileName = fileName.replace('.xml', '.txt')

            print(fileName, outputFileName)

            r, txt = ncsd_cellbodies.main(angle, ref_contour, center_contour, False, (self.cbSVG.checkState() == Qt.Checked), fileName, outputFileName)

            if r:
                flagError = True
                if txt:
                    QMessageBox.warning(self, 'Error', txt)
                else:
                    QMessageBox.warning(self, 'Error', 'Error #%d' % r)

        if not flagError:
            QMessageBox.about(self, 'Cell bodies analysis', 'Analysis done successfully')

        self.statusbar.showMessage('Done', 2000)

    def pbDistanceAnalysis_clicked(self):
        '''Analyze markers in sub contours of neurolucida XML file'''

        try:
            import ncsd_distance
        except ImportError:
            QMessageBox.warning(self, 'NCSD', '"ncsd_distance" module not found!')
            return

        if (self.filename == '') and (self.directoryName == ''):
            QMessageBox.warning(self, 'NCSD', 'Select a XML file or a directory!')
            return

        self.statusbar.showMessage("Distance analysis", 0)

        angle = self.sbAngle.value()


        # reference contour
        ref_contour = str(self.leRefContour.text())

        # center contour
        center_contour = str(self.le_centerContour.text())

        # number of sub contours
        n_subcontours = self.sbSubContoursNb.value()


        if self.filename:
            files_list = [str(self.filename)]
        if self.directoryName:
            import glob
            files_list = glob.glob(str(self.directoryName) + '/*.xml')

            if not files_list:
                QMessageBox.warning(self, 'NCSD', 'The directory %s does not contain any XML file' % str(self.directoryName))
                return

        flagError = False
        for fileName in files_list:

            outputFileName = fileName.replace('.xml', '.txt')

            print(fileName, outputFileName)
            r, txt = ncsd_distance.main(angle, ref_contour, center_contour, n_subcontours, False, (self.cbSVG.checkState() == Qt.Checked), fileName, outputFileName)

            if r:
                flagError = True
                if txt:
                    QMessageBox.warning(self, "Error", txt)
                else:
                    QMessageBox.warning(self, "Error", 'Error #%d' % r)

        if not flagError:
            QMessageBox.about(self, "Distance analysis", 'Analysis done successfully')

        self.statusbar.showMessage("Done", 2000)



    def openXML(self):
        """Open XML"""
        self.statusbar.showMessage('Open Neurolucida XML file', 0)
        fd = QFileDialog(self)
        self.filename = fd.getOpenFileName(self, 'Open a Neurolucida XML file', '', 'XML file (*.xml);;All files (*)')
        if type(self.filename) == type(()):  ### a tuple?
            self.filename = self.filename[0]
        
        if self.filename:
            self.lbFileName.setText(self.filename)
            # cancel directory
            self.directoryName = ''
            self.lbDirectoryName.setText('')

            output_filename = str(self.filename)
            output_filename = output_filename.replace('.xml', '.txt')

        self.statusbar.showMessage(self.filename, 0)


    def openAll(self):
        """Open all XML from directory"""
        self.statusbar.showMessage('Open all XML from directory', 0)
        fd = QFileDialog(self)

        self.directoryName = QFileDialog.getExistingDirectory(self, 'Open Directory', '.', QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        print(self.directoryName)
        if self.directoryName:
            self.lbDirectoryName.setText(self.directoryName)
            # cancel file name
            self.filename = ''
            self.lbFileName.setText('')
        self.statusbar.showMessage(self.directoryName, 0)


    def about(self):
        import platform

        try:
            import ncsd_cellbodies
            cellbodies_version = ncsd_cellbodies.version
        except ImportError:
            cellbodies_version = 'not found'

        try:
            import ncsd_markers
            markers_version = ncsd_markers.version
        except ImportError:
            markers_version = 'not found'

        try:
            import ncsd_distance
            distance_version = ncsd_distance.version
        except ImportError:
            distance_version = 'not found'


        QMessageBox.about(self, "About Neurolucida Companion for Spatial Distribution Analysis",
        """<b>Neurolucida Companion for Spatial Distribution Analysis</b><br>version %s<br><br>
        Markers analysis version v. %s<br>
        Cell bodies analysis version v. %s<br>
        Distance analysis version v. %s
        <p>Copyright &copy; 2010-2011 Olivier Friard - Universit&agrave; di Torino.<br>
        All rights reserved.
        <p>Python %s - Qt %s - PySide %s on %s""" % (__version__, markers_version, cellbodies_version, distance_version, platform.python_version(), PySide.QtCore.__version__, PySide.__version__, platform.system()))


    def quit(self):
        self.close()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    neurolucida = MainWindow()
    neurolucida.show()
    neurolucida.raise_()
    sys.exit(app.exec_())


