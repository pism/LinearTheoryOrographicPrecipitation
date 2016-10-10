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

NOTE:

Several functions, including, but not limited to, updateFile and updateVariable
are adaptions of the functions found in the netcdfbrowser plugin
(https://github.com/etiennesky/netcdfbrowser). Credit for these funtions is due
to Etienne Tourigny (etourigny.dev@gmail.com).
"""

import os
# QWdiget bug has issues with uic importing, the except should be removed
# once this is fixed.
try:
    # QGIS VERSION >= 2.14
    from qgis.PyQt import QtGui, QtCore, uic
    from qgis.PyQt.QtCore import QFileInfo
except:
    from PyQt4 import QtGui, QtCore, uic
    from PyQt4.QtCore import QFileInfo

from qgis.utils import iface

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

import logging
import re
import math
import numpy as np
from osgeo import gdal, osr
from netcdftime import utime
import json

hasCFunits = False
try:
    import cf_units
    hasCFunits = True
except:
    pass

from linear_orog_precip import OrographicPrecipitation

_units = [
    'days',
    'hours',
    'minutes',
    'seconds',
    'day',
    'hour',
    'minute',
    'second']


# create logger
logger = logging.getLogger('LTOP')
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
# fh = logging.RotatingFileHandler('/Users/andy/ltop.log')
fh = logging.StreamHandler()
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s')

# add formatter to ch and fh
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
#logger.addHandler(ch)
logger.addHandler(fh)

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'lt_model_dialog_base.ui'))


def uniquify_list(seq, idfun=None):
    '''
    Remove duplicates from a list, order preserving.
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    '''

    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def get_all_raster_drivers():
    gdal.AllRegister()
    raster_driver_dict = dict()
    for no in range(gdal.GetDriverCount()):
        driver = gdal.GetDriver(no)
        if driver.GetMetadataItem(gdal.DCAP_RASTER):
            raster_driver_dict[
                driver.GetMetadataItem(
                    gdal.DMD_LONGNAME)] = driver.GetMetadataItem(
                gdal.DMD_EXTENSIONS)
    return raster_driver_dict


def readRaster(uri, rasterBand):
    proj4 = None
    gdal.AllRegister()
    logger.debug('readRaster: Open uri {} and band {}'.format(uri, rasterBand))
    ds = gdal.Open(uri)
    rb = ds.GetRasterBand(rasterBand)
    RasterArray = rb.ReadAsArray()
    projection = ds.GetProjection()

    geoTrans = ds.GetGeoTransform()
    pxwidth = ds.RasterXSize
    pxheight = ds.RasterYSize
    ulx = geoTrans[0]
    uly = geoTrans[3]
    rezX = geoTrans[1]
    rezY = geoTrans[5]
    geoTrans = geoTrans
    rx = ulx + pxwidth * rezX
    ly = uly + pxheight * rezY
    osr_ref = osr.SpatialReference()
    osr_ref.ImportFromWkt(projection)
    proj4 = osr_ref.ExportToProj4()

    easting = np.arange(ulx, rx + rezX, rezX)
    northing = np.arange(ly, uly - rezY, -rezY)
    X, Y = np.meshgrid(easting, northing)
    logger.debug('readRaster: geoTrans {} and proj4 {}'.format(geoTrans, proj4))

    return X, Y, RasterArray, geoTrans, proj4


def saveRaster(newRasterfn, geoTrans, proj4, array):
    '''
    Function to export geo-coded raster

    Parameters
    ----------

    '''

    cols = array.shape[1]
    rows = array.shape[0]
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(geoTrans)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(proj4)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


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
        logger.info(
            'Initializing instance of LinearOrographicPrecipitationDialog')
        self.inFileName = None
        self.outFileName = None
        self.prefix = ''
        self.variables = None
        self.uri = None
        self.isNetCDF = False
        self.rasterBand = None
        self.netCDFVariable = None
        self.plugin_dir = os.path.dirname(__file__)
        self.setupUi(self)
        self.configUI()
        self.physical_constants = dict()
        self.connectSignals()
        logger.info(
            'Finished initializing instance of LinearOrographicPrecipitationDialog')

    def configUI(self):
        logger.info('Running configUI')
        self.restoreDefaults()
        self.inputLineEdit.setReadOnly(True)
        self.outputLineEdit.setReadOnly(True)
        self.runButton.setDisabled(True)
        self.timesComboBox.setDisabled(True)
        self.setupUnitsComboBox()

    def connectSignals(self):
        logger.info('Running connectSignals')
        self.closeButton.clicked.connect(self.close)
        self.closeButton_help.clicked.connect(self.close)
        self.runButton.clicked.connect(self.run)
        self.inputButton.clicked.connect(self.showOpenDialog)
        self.outputButton.clicked.connect(self.showSaveDialog)
        self.gaussianCheckBox.clicked.connect(self.change_input)
        self.restoreButton.clicked.connect(self.restoreDefaults)
        self.resetButton.clicked.connect(self.reset)
        self.varsComboBox.currentIndexChanged.connect(self.update_variable)
        self.timesComboBox.currentIndexChanged.connect(self.update_time)
        self.unitsComboBox.currentIndexChanged.connect(self.update_units)

    def setupUnitsComboBox(self):
        logger.info('Running setupUnitsComboBox')
        self.unitsComboBox.addItem('mm hr-1')
        self.out_units = self.unitsComboBox.currentText()
        if hasCFunits:
            self.unitsComboBox.addItem('mm s-1')
            self.unitsComboBox.addItem('mm yr-1')
            self.unitsComboBox.addItem('m hr-1')
            self.unitsComboBox.addItem('m s-1')
            self.unitsComboBox.addItem('m yr-1')
            self.unitsComboBox.setDisabled(False)
            self.cf_units_text.clear()

    def restoreDefaults(self):
        logger.info('Running restoreDefaults')
        defaults_file = os.path.join(self.plugin_dir, 'default_constants.json')
        with open(defaults_file, 'r') as jsonfile:
            defaults = json.load(jsonfile)
        self.Cw.setText(str(defaults['Cw']))
        self.Hw.setText(str(defaults['Hw']))
        self.latitude.setText(str(defaults['latitude']))
        self.tau_c.setText(str(defaults['tau_c']))
        self.tau_f.setText(str(defaults['tau_f']))
        self.P0.setText(str(defaults['P0']))
        self.P_scale.setText(str(defaults['P_scale']))

    def reset(self):
        logger.info('Running reset')
        self.restoreDefaults
        self.inFileName = None
        self.outFileName = None
        self.prefix = ''
        self.variables = None
        self.uri = None
        self.isNetCDF = False
        self.rasterBand = None
        self.netCDFVariable = None
        self.inputLineEdit.clear()
        self.outputLineEdit.clear()
        self.varsComboBox.blockSignals(True)
        self.varsComboBox.clear()
        self.varsComboBox.blockSignals(False)
        self.timesComboBox.blockSignals(True)
        self.timesComboBox.clear()
        self.timesComboBox.blockSignals(True)
        self.unitsComboBox.setCurrentIndex(0)

    def update_units(self):
        logger.info('Running update_units')
        self.out_units = self.unitsComboBox.currentText()

    def update_variable(self):
        logger.info('Running update_variable')
        logger.debug('isnetCDF = '.format(self.isNetCDF))
        if self.isNetCDF:
            self.updateNetCDFVariable()
        else:
            self.updateRasterBand()

    def update_time(self):
        logger.info('Running update_time')
        # GDAL rasterBands are 1 indexed
        self.rasterBand = self.timesComboBox.currentIndex() + 1
        logger.debug('update_time: current rasterBand is {}'.format(self.rasterBand))

    def run(self):
        logger.info('Running run')
        logger.debug('Preparing physical_constants dict')
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
        physical_constants['f'] = 2 * 7.2921e-5 * \
            np.sin(latitude * np.pi / 180)
        physical_constants['Nm'] = Nm
        physical_constants['Cw'] = Cw
        physical_constants['Hw'] = Hw
        physical_constants['u'] = - \
            np.sin(wind_direction * 2 * np.pi / 360) * wind_speed
        physical_constants['v'] = np.cos(
            wind_direction * 2 * np.pi / 360) * wind_speed
        physical_constants['P0'] = P0
        physical_constants['P_scale'] = P_scale

        inFileName = self.inFileName
        if self.gaussianCheckBox.isChecked():
            dx = dy = 750
            xmin = -100e3
            xmax = 200e3
            ymin = -150e3
            ymax = 150e3
            X, Y, Orography = self.gaussian_bump(
                xmin, xmax, ymin, ymax, dx, dx)
            geoTrans = [0., dx, 0., 0., 0., -dy]
            proj4 = ''
        else:
            X, Y, Orography, geoTrans, proj4 = readRaster(
                self.uri, self.rasterBand)
        logger.info('Calling OrographicPrecipitation')
        logger.debug('hasCFunits = {}'.format(hasCFunits))
        if hasCFunits:
            OP = OrographicPrecipitation(
                X,
                Y,
                Orography,
                physical_constants,
                truncate=self.truncateCheckBox.isChecked(),
                ounits=self.unitsComboBox.currentText())
        else:
            OP = OrographicPrecipitation(
                X,
                Y,
                Orography,
                physical_constants,
                truncate=self.truncateCheckBox.isChecked())
        logger.debug('run: geoTrans {} and proj4 {}'.format(geoTrans, proj4))
        outFileName = self.outFileName
        saveRaster(outFileName, geoTrans, proj4, OP.P)
        if self.addResultCheckBox.isChecked():
            fileInfo = QFileInfo(outFileName)
            baseName = fileInfo.baseName()
            iface.addRasterLayer(outFileName, baseName)
        self.close()

    def updateURI(self):
        logger.info('updateURI')
        # update URI
        fileInfo = QFileInfo(self.inFileName)
        if self.isNetCDF:
            uri = 'NETCDF:{}:{}'.format(
                fileInfo.fileName(),
                self.varsComboBox.currentText())
        else:
            uri = self.inFileName
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
        logger.info('updateFile')
        fileName = self.inFileName
        logger.info('updateFile ' + fileName)
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
            try:
                md = ds.GetMetadata("SUBDATASETS")
                for key in sorted(md.iterkeys()):
                    # SUBDATASET_1_NAME=NETCDF:"file.nc":var
                    if re.match('^SUBDATASET_[0-9]+_NAME$', key) is None:
                        continue
                    m = re.search('^(NETCDF:".+"):(.+)', md[key])
                    if m is None:
                        continue
                    self.prefix = m.group(1)
                    self.variables.append(m.group(2))
            except:
                md = ds.GetMetadata()
                my_vars = []
                for key in sorted(md.iterkeys()):
                    ncvar = key.split('#')[0]
                    if ncvar not in ('NC_GLOBAL'):
                        my_vars.append(ncvar)
                my_vars = uniquify_list(my_vars)
                for myvar in my_vars:
                    self.variables.append(myvar)
        else:
            for band in range(ds.RasterCount):
                band += 1
                srcband = ds.GetRasterBand(band)
                if srcband is None:
                    continue
                self.variables.append('Band {}'.format(band))

        self.varsComboBox.blockSignals(True)
        self.varsComboBox.clear()
        for var in sorted(self.variables):
            self.varsComboBox.addItem(var)
        if self.isNetCDF:
            self.updateNetCDFVariable()
        else:
            self.updateRasterBand()
        self.varsComboBox.blockSignals(False)


    def showOpenDialog(self):
        self.timesComboBox.setDisabled(True)
        fileName = _fromUtf8(
            QtGui.QFileDialog.getOpenFileName(
                self, "Input Raster File:"))

        if len(fileName) is not 0:
            self.inFileName = fileName

        gdal.AllRegister()
        dataset = gdal.Open(_fromUtf8(self.inFileName))
        self.rasterBands = dataset.RasterCount
        dataset = None

        if self.inFileName is not None and self.outFileName is not None:
            self.runButton.setDisabled(False)

        self.inputLineEdit.clear()
        self.inputLineEdit.setText(self.inFileName)
        inFileNameSuffix = QFileInfo(self.inFileName).suffix()
        # FIXME: use gdal to test if netcdffile
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
        self.netCDFVariable = self.varsComboBox.currentText()
        dim_map = dict()
        self.dim_names = []
        self.dim_values = dict()
        self.dim_values2 = dict()
        self.dim_def = dict()
        self.dim_band = dict()
        # self.clear()
        uri = 'NETCDF:"%s":%s' % (self.inFileName, self.netCDFVariable)
        self.uri = uri

        logger.info('updateVariable ' + _fromUtf8(uri))

        # look for extra dim definitions
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
                line = "{}={}".format(key, md[key])
                m = re.search('^(NETCDF_DIM_.+)={(.+)}', line)
                if m is not None:
                    dim_map[m.group(1)] = m.group(2)
                else:
                    # netCDF file has time dimension but only 1 entry
                    # 'NETCDF_DIM_time_VALUES=-155088000'
                    m = re.search('^(NETCDF_DIM_.+)=(.+)', line)
                    if m is not None:
                        dim_map[m.group(1)] = m.group(2)

        if 'NETCDF_DIM_EXTRA' in dim_map:

            tok = dim_map['NETCDF_DIM_EXTRA']
            if tok is not None:
                for dim in tok.split(','):
                    self.dim_names.append(dim)
                    tok2 = dim_map.get('NETCDF_DIM_' + dim + '_VALUES')
                    self.dim_values[dim] = []
                    if tok2 is not None:
                        for s in tok2.split(','):
                            self.dim_values[dim].append(num(s))
                    tok2 = dim_map.get('NETCDF_DIM_' + dim + '_DEF')
                    self.dim_def[dim] = []
                    if tok2 is not None:
                        for s in tok2.split(','):
                            self.dim_def[dim].append(num(s))

            dim_names = self.dim_names
            self.dim_names = []
            for dim in dim_names:
                self.dim_names.append(dim)
            for dim in dim_names:
                if dim in self.dim_values:
                    if (dim + "#units") in md:
                        timestr = md[dim + "#units"]
                        units = timestr.split()[0].lower()
                        if (dim + "#calendar") in md:
                            calendar = md[dim + "#calendar"]
                            if calendar in ('none'):
                                # PISMS writes 'none' as calendar
                                calendar = '365_day'
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
                                for i in range(0, len(self.dim_values2[dim])):
                                    self.dim_values2[dim][
                                        i] = self.dim_values2[dim][i][0:10]
            self.updateNetCDFTime()
        else:
            # No EXTRA dim found, has only 1 raster band
            band = 1
            self.rasterBand = band

    def updateNetCDFTime(self):
        self.timesComboBox.setDisabled(False)
        self.timesComboBox.blockSignals(True)
        self.timesComboBox.clear()
        dim = 'time'

        if dim in self.dim_values2:
            for k, item in enumerate(self.dim_values2[dim]):
                self.timesComboBox.addItem('time {}: {}'.format(k, item))
        else:
            for k, item in enumerate(self.dim_values['time']):
                self.timesComboBox.addItem('time {}: {}'.format(k, item))
        self.timesComboBox.blockSignals(False)
        self.update_time()
        if debug:
            print('done updateNetCDFTime ' + _fromUtf8(item))

    def showSaveDialog(self):
        # Declare the filetype in which to save the output file
        # Currently the plugin only supports tif files
        fileTypes = 'All Supported Raster Files (*.vrt *.tif *.tiff *.img *.asc *.png *.jpg *.jpeg *.gif *.xpm *.bmp *.pix *.map *.mpr *.mpl *.hgt *.nc *.grb *.rst *.grd *.rda *.hdr *.dem *.blx *.sqlite *.sdat)'
        fileTypes = 'GeoTiff (*.tif *.tiff)'
        fileName, filter = QtGui.QFileDialog.getSaveFileNameAndFilter(
            self, 'Output Raster File:', '', fileTypes)
        if len(fileName) is 0:
            return
        else:
            # Extract the base filename without the suffix if it exists
            # Convert the fileName from QString to python string
            fileNameStr = _fromUtf8(fileName)

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
            filterStr = _fromUtf8(filter)

            # Split the filter string where/if an asterisk (*) exists
            # I do this to find where the first suffix of the selected filetype
            # occurs
            splittedFilter = filterStr.split('*')

            # If a suffix is not supplied by the user it will be automatically
            # added to the filename. The default suffix will be the first
            # available suffix for the chosen filetype
            if not suffixExists:
                # Extract the 'dirty' suffix string where the first suffix is
                # located
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
                    # Extract the 'dirty' suffix string where the suffix is
                    # located
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

                        # Extract the rest of the suffixes and put them in the
                        # list
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
            self.outFileNameSuffixStr = suffixStr

        if (self.inFileName is not None or self.gaussianCheckBox.isChecked()
                ) and self.outFileName is not None:
            self.runButton.setDisabled(False)
        self.outputLineEdit.clear()
        self.outputLineEdit.setText(self.outFileName)

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
        Orography = h_max * \
            np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) + ((Y - y0)**2 / (2 * sigma_y**2))))
        return X, Y, Orography
