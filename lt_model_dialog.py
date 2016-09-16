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
import numpy as np
from linear_orog_precip import OrographicPrecipitation, array2raster, GdalFile

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
        self.in_file = '/Users/andy/base/LinearOrographicPrecipitation/olympics_500m.tif'
        self.out_file = 'foo.nc'
        self.setupUi(self)
        self.physical_constants = dict()
        self.connectSignals()

    def connectSignals(self):
        self.button_close.clicked.connect(self.close)
        self.button_run.clicked.connect(self.run)
        
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

    def gaussian_bump(self):
        # Reproduce Fig 4c in SB2004
        dx = dy = 750.
        x, y = np.arange(-100e3, 200e3, dx), np.arange(-150e3, 150e3, dy)
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
        physical_constants['tau_c'] = tau_c      # conversion time [s]
        physical_constants['tau_f'] = tau_f      # fallout time [s]
        physical_constants['f'] = 2 * 7.2921e-5 * np.sin(latitude * np.pi / 180) # Coriolis force
        physical_constants['Nm'] = Nm   # moist stability frequency [s-1]
        physical_constants['Cw'] = Cw  # uplift sensitivity factor [kg m-3]
        physical_constants['Hw'] = Hw  # vapor scale height (m)
        physical_constants['u'] = -np.sin(wind_direction*2*np.pi/360) * wind_speed  # x-component of wind vector [m s-1]
        physical_constants['v'] = np.cos(wind_direction*2*np.pi/360) * wind_speed   # y-component of wind vector [m s-1]
        physical_constants['P0'] = P0   # background precip [mm hr-1]
        physical_constants['P_scale'] = P_scale   # precip scale factor [1]
        # FIXME: do only if no raster is select, this reproduces Fig 4c
        X, Y, Orography = self.gaussian_bump()
        # X = self.gd.X
        # Y = self.gd.Y
        # Orography = self.gd.RasterArray
        OP = OrographicPrecipitation(X, Y, Orography, physical_constants, truncate=truncate)
        P = OP.P
        units = OP.P_units
        # get this from dialog
        if self.in_file is not None:
            array2raster(self.out_file, self.gd.geoTrans, self.gd.proj4, units, P)
        else:
            geoTrans = [0., OP.dx, 0., 0., 0., -OP.dy]
            array2raster(self.out_file, geoTrans, '', units, P)
        self.close()
