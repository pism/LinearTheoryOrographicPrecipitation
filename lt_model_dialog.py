# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LinearTheoryOrographicPrecipitationDialog
                                 A QGIS plugin
 Implements the Smith & Barstad (2004) LT model
                             -------------------
        begin                : 2016-09-15
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Andy Aschwanden
        email                : andy.aschwanden@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os

from PyQt4 import QtGui, uic
from PyQt4.QtCore import QFileInfo
from qgis.utils import iface

import numpy as np
from osgeo import gdal
from linear_orog_precip import OrographicPrecipitation, saveRaster, ReadRaster

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'lt_model_dialog_base.ui'))


class LinearTheoryOrographicPrecipitationDialog(QtGui.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(LinearTheoryOrographicPrecipitationDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.configUI()
        self.physical_constants = dict()
        self.connectSignals()
        self.inFileName = None
        self.outFileName = None


    def configUI(self):
        self.inputLineEdit.setReadOnly(True)
        self.outputLineEdit.setReadOnly(True)
        self.runButton.setDisabled(True)


    def connectSignals(self):
        self.closeButton.clicked.connect(self.close)
        self.runButton.clicked.connect(self.run)
        self.inputButton.clicked.connect(self.showOpenDialog)
        self.outputButton.clicked.connect(self.showSaveDialog)
        self.gaussianCheckBox.clicked.connect(self.change_input)


    def get_Cw(self):
        Cw = float(str(self.Cw.text()))
        return Cw


    def get_Hw(self):
        Hw = float(str(self.Hw.text()))
        return Hw


    def get_latitude(self):
        latitude = float(str(self.latitude.text()))
        return latitude


    def get_Nm(self):
        Nm = float(str(self.Nm.text()))
        return Nm

    
    def get_tau_c(self):
        tau_c = float(str(self.tau_c.text()))
        return tau_c

    
    def get_tau_f(self):
        tau_f = float(str(self.tau_f.text()))
        return tau_f

    
    def get_truncate(self):
        # Implement: True if box is checked
        return True

    
    def get_P0(self):
        P0 = float(str(self.P0.text()))
        return P0

    
    def get_P_scale(self):
        P_scale = float(str(self.P_scale.text()))
        return P_scale

    
    def get_wind_direction(self):
        wind_direction = float(str(self.wind_direction.text()))
        return wind_direction


    def get_wind_speed(self):
        wind_speed = float(str(self.wind_speed.text()))
        return wind_speed

    
    def gaussian_bump(self, xmin, xmax, ymin, ymax, dx, dy):
        # Reproduce Fig 4c in SB2004
        x, y = np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dy)
        h_max = 500.
        x0 = -25e3
        y0 = 0
        sigma_x = sigma_y = 15e3
        X, Y = np.meshgrid(x, y)
        Orography = h_max * np.exp(-(((X-x0)**2/(2*sigma_x**2))+((Y-y0)**2/(2*sigma_y**2))))
        return X, Y, Orography

    def run(self):

        Cw = self.get_Cw()
        Hw = self.get_Hw()
        latitude = self.get_latitude()
        Nm = self.get_Nm()
        P0 = self.get_P0()
        P_scale = self.get_P_scale()
        tau_c = self.get_tau_c()
        tau_f = self.get_tau_f()
        truncate = self.get_truncate()
        wind_direction = self.get_wind_direction()
        wind_speed = self.get_wind_speed()
        physical_constants = self.physical_constants
        physical_constants['tau_c'] = tau_c
        physical_constants['tau_f'] = tau_f
        physical_constants['f'] = 2 * 7.2921e-5 * np.sin(latitude * np.pi / 180)
        physical_constants['Nm'] = Nm
        physical_constants['Cw'] = Cw
        physical_constants['Hw'] = Hw
        physical_constants['u'] = -np.sin(wind_direction * 2 * np.pi / 360) * wind_speed
        physical_constants['v'] = np.cos(wind_direction * 2 * np.pi / 360) * wind_speed
        physical_constants['P0'] = P0
        physical_constants['P_scale'] = P_scale
        
        inFileName = self.inFileName
        if self.gaussianCheckBox.isChecked():
            dx = dy = 750
            xmin = -100e3
            xmax = 200e3
            ymin = -150e3
            ymax = 150e3
            X, Y, Orography = self.gaussian_bump(xmin, xmax, ymin, ymax, dx, dx)
            geoTrans = [0., dx, 0., 0., 0., -dy]
            proj4 = ''
        else:
            gd = ReadRaster(inFileName)
            X = gd.X
            Y = gd.Y
            Orography = gd.RasterArray
            geoTrans = gd.geoTrans
            proj4 = gd.proj4
        OP = OrographicPrecipitation(X, Y, Orography, physical_constants, truncate=self.truncateCheckBox.isChecked())
        outFileName = self.outFileName
        saveRaster(outFileName, geoTrans, proj4, OP.P)
        if self.addResultCheckBox.isChecked():
            fileInfo = QFileInfo(outFileName)
            baseName = fileInfo.baseName()
            iface.addRasterLayer(outFileName, baseName)
        self.close()

    def change_input(self):
        self.inputLineEdit.setText('non-georeferenced Gaussian bump')
        
    def showOpenDialog(self):
        fileName = str(QtGui.QFileDialog.getOpenFileName(self,
                                                        "Input Raster File:"))

        if len(fileName) is 0:
            return
        else:
            self.inFileName = fileName;

        gdal.AllRegister()
        dataset = gdal.Open(str(self.inFileName))
        self.rasterBands = dataset.RasterCount
        dataset = None


        if self.inFileName is not None and self.outFileName is not None:
            self.runButton.setDisabled(False)

        self.inputLineEdit.clear()
        self.inputLineEdit.setText(self.inFileName)


    def showSaveDialog(self):
        #self.outFileName = str(QtGui.QFileDialog.getSaveFileName(self,
        #        'Output Raster File:', '', '*.tif'))

        # Declare the filetype in which to save the output file
        # Currently the plugin only supports tif files
        fileTypes = 'GeoTIFF Files (*.tif *.tiff)'
        fileName, filter = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self, 'Output Raster File:', '', fileTypes)

        if len(fileName) is 0:
            return
        else:
            # Extract the base filename without the suffix if it exists
            # Convert the fileName from QString to python string
            fileNameStr = str(fileName)

            # Split the fileNameStr where/if a '.' exists
            splittedFileName = fileNameStr.split('.')

            # Finally extract the base filename from the splitted filename
            baseFileName = splittedFileName[0]

            # Initialize the suffix string
            suffixStr = ''

            # Check if the user entered a suffix
            suffixExists = False
            existingSuffix = ''
            if len(splittedFileName) != 1:
                existingSuffix = splittedFileName[len(splittedFileName) - 1]
                if existingSuffix is not None:
                    suffixExists = True


            # Extract the suffix from the selected filetype filter
            # Convert the selected filter from QString to python string
            filterStr = str(filter)

            # Split the filter string where/if an asterisk (*) exists
            # I do this to find where the first suffix of the selected filetype
            # occurs
            splittedFilter = filterStr.split('*')

            # If a suffix is not supplied by the user it will be automatically
            # added to the filename. The default suffix will be the first
            # available suffix for the chosen filetype
            if not suffixExists:
                # Extract the 'dirty' suffix string where the first suffix is located
                dirtySuffixStr = splittedFilter[1]

                # Find out the number of the available suffixes
                suffixNum = len(splittedFilter) - 1

                if suffixNum == 1:
                    # Split the dirty suffix string where a ')' occurs
                    # which indicates where the selected filetype ends
                    splittedDirtySuffixStr = dirtySuffixStr.split(')')
                else:
                    # Split the dirty suffix string where a space occurs which
                    # indicates where the selected filetype suffix ends
                    splittedDirtySuffixStr = dirtySuffixStr.split(' ')
                suffixStr = splittedDirtySuffixStr[0]
            else:
                # WE NEED TO CHECK IF THE SUPPLIED SUFFIX CORRESPONDS TO THE
                # SELECTED FILETYPE

                # Extract all the suffixes available for the selected filetype
                # First find out the number of the available suffixes
                suffixNum = len(splittedFilter) - 1

                if suffixNum == 1:
                    # Extract the 'dirty' suffix string where the suffix is located
                    dirtySuffixStr = splittedFilter[1]

                    # Split the dirty suffix string where a space occurs which
                    # indicates where the selected filetype suffix ends
                    splittedDirtySuffixStr = dirtySuffixStr.split(' ')
                    suffixStr = splittedDirtySuffixStr[0]


                else:
                    suffixList = []
                    if suffixNum == 2:
                        # Extract the first suffix and put it in the list
                        dirtySuffixStr = splittedFilter[1]
                        splittedDirtySuffixStr = dirtySuffixStr.split(' ')
                        suffixList.append(splittedDirtySuffixStr[0])

                        # Extract the second suffix and put it in the list
                        dirtySuffixStr = splittedFilter[2]
                        splittedDirtySuffixStr = dirtySuffixStr.split(')')
                        suffixList.append(splittedDirtySuffixStr[0])

                    else:
                        # Extract the first suffix and put it in the list
                        dirtySuffixStr = splittedFilter[1]
                        splittedDirtySuffixStr = dirtySuffixStr.split(' ')
                        suffixList.append(splittedDirtySuffixStr[0])

                        # Extract the last suffix and put it in the list
                        dirtySuffixStr = splittedFilter[suffixNum]
                        splittedDirtySuffixStr = dirtySuffixStr.split(')')
                        suffixList.append(splittedDirtySuffixStr[0])

                        # Extract the rest of the suffixes and put them in the list
                        for i in xrange(2, suffixNum):
                            dirtySuffixStr = splittedFilter[i]
                            splittedDirtySuffixStr = dirtySuffixStr.split(' ')
                            suffixList.append(splittedDirtySuffixStr[0])

                    # Find if the user supplied suffix is valid for the
                    # chosen filetype and set it as the filename suffix
                    isValidSuffix = False
                    userSuffix = '.' + existingSuffix
                    for i in xrange(suffixNum + 1):
                        if userSuffix == suffixList[i]:
                            isValidSuffix = True
                            suffixStr = userSuffix
                            break

                    # If the supplied suffix is not valid replace it
                    # with the default suffix for the chosen filetype
                    if not isValidSuffix:
                        suffixStr = suffixList[0]

            self.outFileName = baseFileName + suffixStr

        if (self.inFileName is not None or self.gaussianCheckBox.isChecked()) and self.outFileName is not None:
            self.runButton.setDisabled(False)
        self.outputLineEdit.clear()
        self.outputLineEdit.setText(self.outFileName)
