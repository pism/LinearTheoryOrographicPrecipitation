# -*- coding: utf-8 -*-

"""
/***************************************************************************
 LTOrographicPrecipitation
                                 A QGIS plugin
 Implements the Smith & Barstad (2004) LT model
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2018-05-02
        copyright            : (C) 2018 by Andy Aschwanden and Constantine Khrulev
        email                : ckhroulev@alaska.edu
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

__author__ = 'Andy Aschwanden and Constantine Khrulev'
__date__ = '2018-05-02'
__copyright__ = '(C) 2018 by Andy Aschwanden and Constantine Khrulev'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

import os
from PyQt5.QtCore import QCoreApplication
from qgis.core import (Qgis,
                       QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsRasterFileWriter,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterPoint,
                       QgsProcessingParameterCrs,
                       QgsProcessingParameterExtent,
                       QgsProcessingParameterDefinition,
                       QgsProcessingParameterRasterDestination)

import numpy as np
from .linear_orog_precip import gaussian_bump

class LTOrographicPrecipitationTestInput(QgsProcessingAlgorithm):
    """
    Creates the gaussian bump field for testing.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OUTPUT = 'OUTPUT'
    TARGET_CRS = 'TARGET_CRS'
    EXTENT = 'EXTENT'
    DX = 'DX'
    DY = 'DY'
    SIGMA_X = 'SIGMA_X'
    SIGMA_Y = 'SIGMA_Y'
    CENTER = 'CENTER'

    def initAlgorithm(self, config):

        x_min   = -100e3
        x_max   = 200e3
        y_min   = -150e3
        y_max   = 150e3
        dx      = 750
        dy      = 750
        x0      = -25e3
        y0      = 0.0
        sigma_x = 15e3
        sigma_y = 15e3

        self.addParameter(QgsProcessingParameterCrs(self.TARGET_CRS,
                                                    self.tr('Target CRS'),
                                                    'ProjectCrs'))

        self.addParameter(QgsProcessingParameterExtent(self.EXTENT,
                                                       self.tr('Extent'),
                                                       "{}, {}, {}, {}".format(x_min, x_max,
                                                                               y_min, y_max)))

        self.addParameter(QgsProcessingParameterNumber(self.DX,
                                                       self.tr("Grid spacing (dx)"),
                                                       QgsProcessingParameterNumber.Double,
                                                       defaultValue=dx,
                                                       minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(self.DY,
                                                       self.tr("Grid spacing (dy)"),
                                                       QgsProcessingParameterNumber.Double,
                                                       defaultValue=dy,
                                                       minValue=0.0))

        self.addParameter(QgsProcessingParameterPoint(self.CENTER,
                                                      self.tr("Center"),
                                                      "{}, {}".format(x0, y0)))

        s_x = QgsProcessingParameterNumber(self.SIGMA_X,
                                           self.tr("Spread in the X direction (sigma_x)"),
                                           QgsProcessingParameterNumber.Double,
                                           defaultValue=sigma_x,
                                           minValue=0.0)
        s_x.setFlags(s_x.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(s_x)

        s_y = QgsProcessingParameterNumber(self.SIGMA_Y,
                                           self.tr("Spread in the Y direction (sigma_y)"),
                                           QgsProcessingParameterNumber.Double,
                                           defaultValue=sigma_y,
                                           minValue=0.0)
        s_y.setFlags(s_y.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(s_y)

        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT,
                                                                  self.tr('Output')))

    def prepareAlgorithm(self, parameters, context, feedback):
        self.crs    = self.parameterAsCrs(parameters, self.TARGET_CRS, context)
        self.extent = self.parameterAsExtent(parameters, self.EXTENT, context, self.crs)
        self.center = self.parameterAsPoint(parameters, self.CENTER, context, self.crs)

        self.dx = self.parameterAsDouble(parameters, self.DX, context)
        self.dy = self.parameterAsDouble(parameters, self.DY, context)

        self.sigma_x = self.parameterAsDouble(parameters, self.SIGMA_X, context)
        self.sigma_y = self.parameterAsDouble(parameters, self.SIGMA_Y, context)

        if self.dx > self.extent.width():
            feedback.reportError(self.tr("Grid spacing cannot exceed grid extent (x direction)."))
            return False
        if self.dy > self.extent.height():
            feedback.reportError(self.tr("Grid spacing cannot exceed grid extent (y direction)."))
            return False

        self.outputFile = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)

        return True

    def processAlgorithm(self, parameters, context, feedback):

        extent = self.extent

        x_min = extent.xMinimum()
        x_max = extent.xMaximum()
        y_min = extent.yMinimum()
        y_max = extent.yMaximum()

        x0 = self.center.x()
        y0 = self.center.y()

        outputFormat = QgsRasterFileWriter.driverForExtension(os.path.splitext(self.outputFile)[1])

        rows = max([np.ceil(extent.height() / self.dy), 1.0])
        cols = max([np.ceil(extent.width() / self.dx), 1.0])

        writer = QgsRasterFileWriter(self.outputFile)
        writer.setOutputProviderKey('gdal')
        writer.setOutputFormat(outputFormat)
        provider = writer.createOneBandRaster(Qgis.Float64, cols, rows, extent, self.crs)
        provider.setNoDataValue(1, -9999)

        _, _, data = gaussian_bump(x_min, x_max, y_min, y_max, self.dx, self.dy,
                                   x0=x0, y0=y0, sigma_x=self.sigma_x, sigma_y=self.sigma_y)

        provider.write(bytes(data.data),
                       1,       # band
                       cols,    # width
                       rows,    # height
                       0, 0)    # offset

        provider.setEditable(False)

        return {self.OUTPUT: self.outputFile}

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'bump'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr("Gaussian bump")

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr("Testing")

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'testing'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return LTOrographicPrecipitationTestInput()

    def shortHelpString(self):
        return """
Create a raster layer containing a Gaussian "bump" orography that can be used to test the model.

Default parameter values produce the field needed to reproduce Figure 4c in Smith and Barstad (2004)."""
