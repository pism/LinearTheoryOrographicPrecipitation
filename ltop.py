# -*- coding: utf-8 -*-

"""
/***************************************************************************
 LTOrographicPrecipitation
                                 A QGIS plugin
 Implements the Smith & Barstad (2004) LT model
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2018-05-02
        copyright            : (C) 2018-2020 by Andy Aschwanden and Constantine Khrulev
        email                : uaf-pism@alaska.edu
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

__author__ = "Andy Aschwanden and Constantine Khrulev"
__date__ = "2018-05-02"
__copyright__ = "(C) 2018-2020 by Andy Aschwanden and Constantine Khrulev"

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = "$Format:%H$"

import inspect
import os
import sys

from qgis.core import QgsApplication

from .ltop_provider import LTOrographicPrecipitationProvider

cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]

if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)


class LTOrographicPrecipitationPlugin(object):
    def __init__(self):
        self.provider = LTOrographicPrecipitationProvider()

    def initGui(self):
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
