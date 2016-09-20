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

import logging
import re, math
import numpy as np
from osgeo import gdal, osr
from netcdftime import utime

from linear_orog_precip import OrographicPrecipitation, saveRaster

debug = True
_units = ['days', 'hours', 'minutes', 'seconds', 'day', 'hour', 'minute', 'second']

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'lt_model_dialog_base.ui'))

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


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
        self.prefix = ''
        self.variables = None
        self.uri = None
        self.isNetCDF = False
        self.rasterBand = None
        self.netCDFVariable = None

    def configUI(self):
        self.inputLineEdit.setReadOnly(True)
        self.outputLineEdit.setReadOnly(True)
        self.runButton.setDisabled(True)
        self.timesComboBox.setDisabled(True)


    def connectSignals(self):
        self.closeButton.clicked.connect(self.close)
        self.runButton.clicked.connect(self.run)
        self.inputButton.clicked.connect(self.showOpenDialog)
        self.outputButton.clicked.connect(self.showSaveDialog)
        self.gaussianCheckBox.clicked.connect(self.change_input)
        self.varsComboBox.currentIndexChanged.connect(self.update_variable)
        self.timesComboBox.currentIndexChanged.connect(self.update_time)


    def update_variable(self):
        if self.isNetCDF:
            self.updateNetCDFVariable()
        else:
            self.updateRasterBand()


    def update_time(self):
        self.rasterBand = self.timesComboBox.currentIndex()
        print 'current index {}'.format(self.rasterBand)


    def clear(self):
        self.isNetCDF = False
        self.rasterBand = None
        self.varsComboBox.clear()
        self.timesComboBox.clear()


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
            self.readRaster()
            X = self.X
            Y = self.Y
            Orography = self.RasterArray
            geoTrans = self.geoTrans
            proj4 = self.proj4
        OP = OrographicPrecipitation(X, Y, Orography, physical_constants, truncate=self.truncateCheckBox.isChecked())
        outFileName = self.outFileName
        saveRaster(outFileName, geoTrans, proj4, OP.P)
        if self.addResultCheckBox.isChecked():
            fileInfo = QFileInfo(outFileName)
            baseName = fileInfo.baseName()
            iface.addRasterLayer(outFileName, baseName)
        self.close()


    def updateURI(self):
        if debug>0:
            print('updateURI')
        # update URI
        fileInfo = QFileInfo(self.inFileName)
        uri = 'NETCDF:"%s":%s' % (fileInfo.fileName(), self.varsComboBox.currentText())
        self.uri = uri


    def change_input(self):
        if self.gaussianCheckBox.isChecked():
            self.inputLineEdit.setText('<non-georeferenced Gaussian bump>')
            self.varsComboBox.setDisabled(True)
            self.timesComboBox.setDisabled(True)
        else:
            self.inputLineEdit.clear()
            self.varsComboBox.setDisabled(True)
            self.timesComboBox.setDisabled(True)

    def updateFile(self):
        # self.clear()
        fileName = self.inFileName
        if debug>0:
            print('updateFile ' + fileName)
        if fileName == '':
            return

        self.prefix = ''
        self.variables = []

        gdal.PushErrorHandler('CPLQuietErrorHandler')
        ds = gdal.Open(fileName)
        gdal.PopErrorHandler()
        if ds is None:
            return
        if self.isNetCDF:
            md = ds.GetMetadata("SUBDATASETS")
            for key in sorted(md.iterkeys()):
                #SUBDATASET_1_NAME=NETCDF:"file.nc":var
                if re.match('^SUBDATASET_[0-9]+_NAME$', key) is None:
                    continue
                m = re.search('^(NETCDF:".+"):(.+)', md[key])
                if m is None:
                    continue
                self.prefix = m.group(1)
                self.variables.append(m.group(2))
        else:
            for band in range(ds.RasterCount):
                band += 1
                srcband = ds.GetRasterBand(band)
                if srcband is None:
                    continue
                self.variables.append('Band {}'.format(band))

        self.varsComboBox.blockSignals(True)
        self.varsComboBox.clear()
        for var in self.variables:
            self.varsComboBox.addItem(var)
        if self.isNetCDF:
            if debug:
                print 'update netCDFVariable'
            self.updateNetCDFVariable()
        else:
            self.updateRasterBand()    
        self.varsComboBox.blockSignals(False)

        if debug > 0:
            print('done updateFile ' + fileName)


    def showOpenDialog(self):
        self.timesComboBox.setDisabled(True)
        fileName = str(QtGui.QFileDialog.getOpenFileName(self,
                                                        "Input Raster File:"))

        if len(fileName) is not 0:
            self.inFileName = fileName;

        gdal.AllRegister()
        dataset = gdal.Open(str(self.inFileName))
        self.rasterBands = dataset.RasterCount
        dataset = None

        if self.inFileName is not None and self.outFileName is not None:
            self.runButton.setDisabled(False)

        self.inputLineEdit.clear()
        self.inputLineEdit.setText(self.inFileName)
        inFileNameSuffix  = QFileInfo(self.inFileName).suffix()
        if inFileNameSuffix in ('nc', 'nc3', 'nc4'):
            self.isNetCDF = True
        self.updateFile()
        self.timesComboBox.setDisabled(False)


    def updateRasterBand(self):
        fullBandName = self.varsComboBox.currentText()
        rasterBand = num(fullBandName.split('Band ')[1])
        self.rasterBand = rasterBand
        self.uri = self.inFileName


    def updateNetCDFVariable(self):
        self.netCDFVariable =  self.varsComboBox.currentText()
        dim_map = dict()
        self.dim_names = []
        self.dim_values = dict()
        self.dim_values2 = dict()
        self.dim_def = dict()
        self.dim_band = dict()
        # self.clear()
        uri = 'NETCDF:"%s":%s' % (self.inFileName, self.netCDFVariable)

        if debug>0:
            print('updateVariable ' + str(uri))

        #look for extra dim definitions
        #  NETCDF_DIM_EXTRA={time,tile}
        #  NETCDF_DIM_tile_DEF={3,6}
        #  NETCDF_DIM_tile_VALUES={1,2,3}
        #  NETCDF_DIM_time_DEF={12,6}
        #  NETCDF_DIM_time_VALUES={1,32,60,91,121,152,182,213,244,274,305,335}

        gdal.PushErrorHandler('CPLQuietErrorHandler')
        ds = gdal.Open(uri)
        gdal.PopErrorHandler()
        if ds is None:
            return
        md = ds.GetMetadata()
        if md is None:
            return

        for key in sorted(md.iterkeys()):
            if key.startswith('NETCDF_DIM_'):
                line="%s=%s" % (key, md[key])
                m = re.search('^(NETCDF_DIM_.+)={(.+)}', line)
                if m is not None:
                    dim_map[m.group(1)] = m.group(2)

        if not 'NETCDF_DIM_EXTRA' in dim_map:
            self.warning()
            return
        
        tok = dim_map['NETCDF_DIM_EXTRA']
        if debug:
            print tok
        if tok is not None:
            for dim in tok.split(','):
                self.dim_names.append( dim )
                tok2 = dim_map.get('NETCDF_DIM_' + dim + '_VALUES')
                self.dim_values[dim] = []
                if tok2 is not None:
                    for s in tok2.split(','):
                        self.dim_values[dim].append(num(s))
                tok2 = dim_map.get('NETCDF_DIM_' + dim + '_DEF')
                self.dim_def[ dim ] = []
                if tok2 is not None:
                    for s in tok2.split(','):
                        self.dim_def[ dim ].append(num(s))

        # remove any dims which have only 1 element
        dim_names = self.dim_names
        self.dim_names = []
        for dim in dim_names:
            if self.dim_def[dim][0] <= 1:
                del self.dim_values[dim]
                del self.dim_def[dim]
            else:
                self.dim_names.append(dim)

        for dim in dim_names:
            if dim in self.dim_values:
                if (dim + "#units") in md:
                    timestr = md[dim + "#units"]
                    units = timestr.split()[0].lower()
                    if (dim + "#calendar") in md:
                        calendar = md[dim + "#calendar"]
                        cdftime = utime(timestr, calendar=calendar)
                    if units in _units:
                        try:
                            dates = cdftime.num2date(self.dim_values[dim])
                        except ValueError:
                            continue
                        self.dim_values2[dim] = []
                        only_days = True
                        for date in dates:
                            val = date.strftime("%Y-%m-%d %H:%M:%S")
                            if not val.endswith(" 00:00:00"):
                                only_days = False
                            self.dim_values2[dim].append(val)
                        if only_days:
                            for i in range(0,len(self.dim_values2[ dim ])):
                                self.dim_values2[dim][i] = self.dim_values2[dim][i][0:10]
        # print "dim_values2 {}".format(self.dim_values2)
        # if debug>0:
        #     print(str(dim_map))
        #     print(str(self.dim_names))
        #     print(str(self.dim_def))
        #     print(str(self.dim_values))
        #     print(str(self.dim_values2))
        self.updateNetCDFTime()


    def updateNetCDFTime(self):
        self.timesComboBox.setDisabled(False)
        self.timesComboBox.blockSignals(True)
        self.timesComboBox.clear()
        dim = 'time'
        print "dim_values2 {}".format(self.dim_values2)
        
        if dim in self.dim_values2:
            for k, item  in enumerate(self.dim_values2[dim]):
                self.timesComboBox.addItem('time {}: {}'.format(k, item))
        else:
            for k, item  in enumerate(self.dim_values['time']):
                self.timesComboBox.addItem('time {}: {}'.format(k, item))
        self.timesComboBox.blockSignals(False)
        self.update_time()
        if debug:
            print('done updateNetCDFTime ' + str(item))


    def showSaveDialog(self):
        # Declare the filetype in which to save the output file
        # Currently the plugin only supports tif files
        fileTypes = 'All Supported Raster Files (*.tif *.tiff)'
        fileName, filter = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self, 'Output Raster File:', '', fileTypes)

        if len(fileName) is 0:
            return
        else:
            # Extract the base filename without the suffix if it exists
            # Convert the fileName from QString to python string
            fileNameStr = str(self.outFileName)

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


        
    def readRaster(self):
        uri = self.uri
        ds = gdal.Open(uri)
        rasterBand = self.rasterBand
        rb = ds.GetRasterBand(rasterBand)
        self.RasterArray = rb.ReadAsArray()
        self.projection = ds.GetProjection()

        geoTrans = ds.GetGeoTransform()
        pxwidth = ds.RasterXSize
        pxheight = ds.RasterYSize
        ulx = geoTrans[0]
        uly = geoTrans[3]
        rezX = geoTrans[1]
        rezY = geoTrans[5]
        self.geoTrans = geoTrans
        rx = ulx + pxwidth * rezX
        ly = uly + pxheight * rezY
        osr_ref = osr.SpatialReference()
        osr_ref.ImportFromWkt(self.projection)
        self.proj4 = osr_ref.ExportToProj4()

        easting = np.arange(ulx, rx + rezX, rezX)
        northing = np.arange(ly, uly - rezY, -rezY)
        self.X, self.Y = np.meshgrid(easting, northing)


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
