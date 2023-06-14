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
 This script initializes the plugin, making it known to QGIS.
"""

__author__ = "Andy Aschwanden and Constantine Khrulev"
__date__ = "2018-05-02"
__copyright__ = "(C) 2018-2020 by Andy Aschwanden and Constantine Khrulev"


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load LTOrographicPrecipitation class from file LTOrographicPrecipitation.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .ltop import LTOrographicPrecipitationPlugin

    return LTOrographicPrecipitationPlugin()
