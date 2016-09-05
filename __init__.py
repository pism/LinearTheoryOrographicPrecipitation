# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LinearTheoryOrographicPrecipitation
                                 A QGIS plugin
 Implements the Smith and Barstad (2004) Model
                             -------------------
        begin                : 2016-09-05
        copyright            : (C) 2016 by Andy Aschwanden
        email                : andy.aschwanden@gmail.com
        git sha              : $Format:%H$
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


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load LinearTheoryOrographicPrecipitation class from file LinearTheoryOrographicPrecipitation.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .lt_model import LinearTheoryOrographicPrecipitation
    return LinearTheoryOrographicPrecipitation(iface)
