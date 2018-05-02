import numpy as np
from osgeo import gdal, osr

import logging
logger = logging.getLogger('LTOP')

np.seterr(divide='ignore', invalid='ignore')

class Constants(object):
    "Spatially-constant inputs of the model"

    tau_c = 1000.0
    "conversion time [s]"

    tau_f = 1000.0
    "fallout time [s]"

    P0 = 0.0
    'Background precipitation rate [mm hr-1]'

    P_scale = 1.0
    'Precipitation scale factor'

    Nm = 0.005
    'moist stability frequency [s-1]'

    Hw = 2500
    'Water vapor scale height [m]'

    latitude = 45.0
    "Latitude used to compute the Coriolis force"

    direction = 270.0
    "Wind direction, 0 is north, 270 is west"

    speed = 15.0
    "Wind speed"

    f = None
    "Coriolis force"

    u = None
    "u component of the wind velocity"

    v = None
    "v component of the wind velocity"

    Cw = None
    "Uplift sensitivity factor [kg m-3]"

    def update(self):
        "Update derived constants"

        Theta_m  = -6.5     # K / km
        rho_Sref = 7.4e-3   # kg m-3
        gamma    = -5.8     # K / km

        self.f = 2 * 7.2921e-5 * np.sin(self.latitude * np.pi / 180.0)

        self.u = -np.sin(self.direction * 2 * np.pi / 360) * self.speed
        self.v = np.cos(self.direction * 2 * np.pi / 360) * self.speed

        self.Cw = rho_Sref * Theta_m / gamma

    def __init__(self):
        self.update()

def orographic_precipitation(dx, dy, orography, constants, truncate):
    """
    Calculates orographic precipitation following Smith & Barstad (2004).

    """
    logger.info('Computing orographic precipitation')

    eps = 1e-18
    pad = 250

    ny, nx = orography.shape
    logger.debug('Raster shape before padding ({},{})'.format(nx, ny))

    padded_orography = np.pad(orography, pad, 'constant')
    ny, nx = padded_orography.shape
    logger.debug('Raster shape after padding ({},{})'.format(ny, nx))

    logger.info('Fourier transform orography')
    padded_orography_fft = np.fft.fft2(padded_orography)

    x_n_value = np.fft.fftfreq(ny, (1.0 / ny))
    y_n_value = np.fft.fftfreq(nx, (1.0 / nx))

    x_len = nx * dx
    y_len = ny * dy
    kx_line = np.divide(np.multiply(2.0 * np.pi, x_n_value), x_len)
    ky_line = np.divide(
        np.multiply(
            2.0 * np.pi,
            y_n_value),
        y_len)[
        np.newaxis].T

    kx = np.tile(kx_line, (nx, 1))
    ky = np.tile(ky_line, (1, ny))

    # Intrinsic frequency sigma = U*k + V*l
    u0 = constants.u
    v0 = constants.v

    logger.info('Calculate sigma')
    sigma = np.add(np.multiply(kx, u0), np.multiply(ky, v0))
    sigma_sqr_reg = sigma ** 2
    m_denom = np.power(sigma, 2.) - constants.f**2

    sigma_sqr_reg[
        np.logical_and(
            np.fabs(sigma_sqr_reg) < eps,
            np.fabs(
                sigma_sqr_reg >= 0))] = eps
    sigma_sqr_reg[
        np.logical_and(
            np.fabs(sigma_sqr_reg) < eps,
            np.fabs(
                sigma_sqr_reg < 0))] = -eps

    # The vertical wave number
    # Eqn. 12
    # Regularization
    m_denom[
        np.logical_and(
            np.fabs(m_denom) < eps,
            np.fabs(m_denom) >= 0)] = eps
    m_denom[
        np.logical_and(
            np.fabs(m_denom) < eps,
            np.fabs(m_denom) < 0)] = -eps

    m1 = np.divide(
        np.subtract(
            constants.Nm**2,
            np.power(
                sigma,
                2.)),
        m_denom)
    m2 = np.add(np.power(kx, 2.), np.power(ky, 2.))
    m_sqr = np.multiply(m1, m2)
    logger.info('Calculating m')
    m = np.sqrt(-1 * m_sqr)
    # Regularization
    m[np.logical_and(m_sqr >= 0, sigma == 0)] = np.sqrt(
        m_sqr[np.logical_and(m_sqr >= 0, sigma == 0)])
    m[np.logical_and(m_sqr >= 0, sigma != 0)] = np.sqrt(m_sqr[np.logical_and(
        m_sqr >= 0, sigma != 0)]) * np.sign(sigma[np.logical_and(m_sqr >= 0, sigma != 0)])
    # Numerator in Eqn. 49
    P_karot_num = np.multiply(np.multiply(np.multiply(
        constants.Cw, 1j), sigma), padded_orography_fft)
    P_karot_denom_Hw = np.subtract(1, np.multiply(
        np.multiply(constants.Hw, m), 1j))
    P_karot_denom_tauc = np.add(1, np.multiply(np.multiply(
        sigma, constants.tau_c), 1j))
    P_karot_denom_tauf = np.add(1, np.multiply(np.multiply(
        sigma, constants.tau_f), 1j))
    # Denominator in Eqn. 49
    P_karot_denom = np.multiply(
        P_karot_denom_Hw, np.multiply(
            P_karot_denom_tauc, P_karot_denom_tauf))
    P_karot = np.divide(P_karot_num, P_karot_denom)

    # Converting from wave domain back to space domain
    logger.info('Performing inverse Fourier transform')
    P = np.fft.ifft2(P_karot)
    spy = 31556925.9747
    logger.info('De-pad array')
    P = P[pad:-pad, pad:-pad]
    P = np.multiply(np.real(P), 3600)   # mm hr-1
    # Add background precip
    P0 = constants.P0
    logger.info('Adding background precpipitation {} mm hr-1'.format(P0))
    P += P0
    # Truncation

    if truncate:
        logger.info('Truncate precipitation')
        P[P < 0] = 0
    P_scale = constants.P_scale
    logger.info('Scale precipitation P = P * {}'.format(P_scale))
    P *= P_scale

    return P


class ReadRaster(object):

    '''
    A class to read a GDAL File

    Parameters
    ----------

    filename: a valid gdal file
    '''

    def __init__(self, file_name):
        self.file_name = file_name
        try:
            print("\n  opening file %s" % file_name)
            ds = gdal.Open(file_name)
        except:
            print("  could not open file %s" % file_name)

        self.RasterArray = ds.ReadAsArray()
        self.projection = ds.GetProjection()

        geoT = ds.GetGeoTransform()
        pxwidth = ds.RasterXSize
        pxheight = ds.RasterYSize
        ulx = geoT[0]
        uly = geoT[3]
        rezX = geoT[1]
        rezY = geoT[5]
        rx = ulx + pxwidth * rezX
        ly = uly + pxheight * rezY
        osr_ref = osr.SpatialReference()
        osr_ref.ImportFromWkt(self.projection)
        self.proj4 = osr_ref.ExportToProj4()

        self.geoTrans = geoT
        self.width = np.abs(pxwidth * rezX)
        self.height = np.abs(pxheight * rezY)
        self.center_x = ulx + pxwidth * rezX / 2
        self.center_y = uly + pxheight * rezY / 2
        self.easting = np.arange(ulx, rx + rezX, rezX)
        self.northing = np.arange(ly, uly - rezY, -rezY)
        self.X, self.Y = np.meshgrid(self.easting, self.northing)


def array2raster(newRasterfn, geoTrans, proj4, units, array):
    '''
    Function to export geo-coded raster

    Parameters
    ----------

    '''

    cols = array.shape[1]
    rows = array.shape[0]

    driver = gdal.GetDriverByName('netCDF')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((geoTrans))
    outband = outRaster.GetRasterBand(1)
    outband.SetMetadata('units', units)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromProj4(proj4)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def gaussian_bump(xmin, xmax, ymin, ymax, dx, dy, h_max=500.0,
                  x0=-25e3, y0=0.0, sigma_x=15e3, sigma_y=15e3):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    # Reproduce Fig 4c in SB2004
    x = np.arange(xmin, xmax, dx)
    y = np.arange(ymin, ymax, dy)
    X, Y = np.meshgrid(x, y)
    Orography = h_max * np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                                 ((Y - y0)**2 / (2 * sigma_y**2))))
    return X, Y, Orography

if __name__ == "__main__":
    print('Linear Theory Orographic Precipitation Model by Smith & Barstad (2004)')

    import itertools
    from argparse import ArgumentParser

    import logging
    import logging.handlers

    # create logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    fh = logging.handlers.RotatingFileHandler('ltop.log')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(lineno)d - %(message)s')
    # add formatter to ch and fh
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)
    logger.addHandler(fh)

    constants = Constants()

    parser = ArgumentParser()
    parser.add_argument('-i', dest='in_file',
                        help='Gdal-readable DEM', default=None)
    parser.add_argument('-o', dest='out_file',
                        help='Output file', default='foo.nc')
    parser.add_argument(
        '--background_precip',
        dest='P0',
        type=float,
        help='Background precipitation rate [mm hr-1].',
        default=0.)
    parser.add_argument('--precip_scale_factor', dest='P_scale', type=float,
                        help='Precipitation scale factor.', default=constants.P_scale)
    parser.add_argument('--no_trunc', dest='truncate', action='store_false',
                        help='Do not truncate precipitation.', default=True)
    parser.add_argument('--latitude', dest='lat', type=float,
                        help='Latitude to compute Coriolis term.', default=constants.latitude)
    parser.add_argument('--tau_c', dest='tau_c', type=float,
                        help='conversion time [s].', default=constants.tau_c)
    parser.add_argument('--tau_f', dest='tau_f', type=float,
                        help='fallout time [s].', default=constants.tau_f)
    parser.add_argument('--moist_stability', dest='Nm', type=float,
                        help='moist stability frequency [s-1].', default=constants.Nm)
    parser.add_argument('--vapor_scale_height', dest='Hw', type=float,
                        help='Water vapor scale height [m].', default=constants.Hw)
    parser.add_argument('--wind_direction', dest='direction', type=float,
                        help='Direction from which the wind is coming.',
                        default=constants.direction)
    parser.add_argument('--wind_magnitude', dest='magnitude', type=float,
                        help='Magnitude of wind velocity [m/s].', default=constants.speed)
    parser.add_argument('--units', dest='units', type=str,
                        help='Output precipitation units.', default="mm hr-1")

    options = parser.parse_args()
    in_file = options.in_file
    out_file = options.out_file

    constants.tau_c = options.tau_c
    constants.tau_f = options.tau_f
    constants.Nm = options.Nm
    constants.Hw = options.Hw
    constants.P0 = options.P0
    constants.P_scale = options.P_scale
    constants.update()

    if in_file is not None:
        gd = ReadRaster(in_file)
        X = gd.X
        Y = gd.Y
        orography = gd.RasterArray
    else:
        # Reproduce Fig 4c in SB2004
        dx = 750
        dy = 750
        xmin = -100e3
        xmax = 200e3
        ymin = -150e3
        ymax = 150e3
        X, Y, orography = gaussian_bump(xmin, xmax, ymin, ymax, dx, dy)

    dx = X[0,1] - X[0,0]
    dy = Y[1,0] - Y[0,0]

    P = orographic_precipitation(dx, dy, orography, constants,
                                 truncate=options.truncate)

    units = options.units
    if units != "mm hr-1":
        from cf_units import Unit
        internal = Unit("mm hr-1")
        P = Unit("mm hr-1").convert(P, Unit(units))

    if in_file is not None:
        array2raster(out_file, gd.geoTrans, gd.proj4, units, P)
    else:
        geoTrans = [0., dx, 0., 0., 0., -dy]
        array2raster(out_file, geoTrans, '', units, P)
