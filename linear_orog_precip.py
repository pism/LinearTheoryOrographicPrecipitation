from scipy.fftpack import fft2, fftfreq
import numpy as np
import cmath
from osgeo import gdal, osr
import logging
logger = logging.getLogger('LTOP')

np.seterr(divide='ignore', invalid='ignore')


class OrographicPrecipitation(object):

    """
    Calculates orographic precipitation following Smith & Barstad (2004).

    """

    def __init__(
            self,
            X,
            Y,
            Orography,
            physical_constants,
            truncate=True,
            ounits=None):
        self.X = X
        self.Y = Y
        self.Orography = Orography
        self.physical_constants = physical_constants
        self.dx = np.diff(X)[0, 0]
        self.dy = np.diff(Y, axis=0)[0, 0]
        self.nx = len(Orography[1, :])
        self.ny = len(Orography)
        self.truncate = truncate

        self.P = self._compute_precip(ounits)
        if ounits is not None:
            self.P_units = ounits
        else:
            self.P_units = 'mm hr-1'

    def _compute_precip(self, ounits):
        logger.info('Running _compute_precip')
        physical_constants = self.physical_constants
        eps = 1e-18
        pad_max = 200
        pad = int(np.ceil(((self.nx + self.ny) / 2) / 100)) * 100
        if pad > pad_max:
            pad = pad_max
        logger.debug(
            'Raster shape before padding ({},{})'.format(
                self.nx, self.ny))
        Orography = np.pad(self.Orography, pad, 'constant')
        nx, ny = Orography.shape
        logger.debug('Raster shape after padding ({},{})'.format(ny, nx))
        logger.info('Fourier transform orography')
        Orography_fft = np.fft.fft2(Orography)
        dx = self.dx
        dy = self.dy
        x_n_value = np.fft.fftfreq(
            ny, (1.0 / ny))
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
        u0 = physical_constants['u']
        v0 = physical_constants['v']
        logger.info('Calculate sigma')
        sigma = np.add(np.multiply(kx, u0), np.multiply(ky, v0))
        sigma_sqr_reg = sigma ** 2
        m_denom = np.power(sigma, 2.) - physical_constants['f']**2

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
                physical_constants['Nm']**2,
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
            physical_constants['Cw'], cmath.sqrt(-1)), sigma), Orography_fft)
        P_karot_denom_Hw = np.subtract(1, np.multiply(
            np.multiply(physical_constants['Hw'], m), cmath.sqrt(-1)))
        P_karot_denom_tauc = np.add(1, np.multiply(np.multiply(
            sigma, physical_constants['tau_c']), cmath.sqrt(-1)))
        P_karot_denom_tauf = np.add(1, np.multiply(np.multiply(
            sigma, physical_constants['tau_f']), cmath.sqrt(-1)))
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
        P0 = physical_constants['P0']
        logger.info('Adding background precpipitation {} mm hr-1'.format(P0))
        P += P0
        # Truncation
        truncate = self.truncate
        if truncate is True:
            logger.info('Truncate precipitation')
            P[P < 0] = 0
        P_scale = physical_constants['P_scale']
        logger.info('Scale precipitation P = P * {}'.format(P_scale))
        P *= P_scale

        if ounits is not None:
            import cf_units
            in_units = cf_units.Unit('mm hr-1')
            out_units = cf_units.Unit(ounits)
            logger.info('Converting mm hr-1 to {}'.format(ounits))
            P = in_units.convert(P, out_units)
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
                        help='Precipitation scale factor.', default=1.)
    parser.add_argument('--no_trunc', dest='truncate', action='store_false',
                        help='Do not truncate precipitation.', default=True)
    parser.add_argument('--latitude', dest='lat', type=float,
                        help='Latitude to compute Coriolis term.', default=45.)
    parser.add_argument('--tau_c', dest='tau_c', type=float,
                        help='conversion time [s].', default=1000)
    parser.add_argument('--tau_f', dest='tau_f', type=float,
                        help='fallout time [s].', default=1000)
    parser.add_argument('--moist_stability', dest='Nm', type=float,
                        help='moist stability frequency [s-1].', default=0.005)
    parser.add_argument('--vapor_scale_height', dest='Hw', type=float,
                        help='Water vapor scale height [m].', default=2500)
    parser.add_argument(
        '--wind_direction',
        dest='direction',
        type=float,
        help='Direction from which the wind is coming.',
        default=270)
    parser.add_argument('--wind_magnitude', dest='magnitude', type=float,
                        help='Magnitude of wind velocity [m/s].', default=15)

    options = parser.parse_args()
    in_file = options.in_file
    out_file = options.out_file
    direction = options.direction
    lat = options.lat
    magnitude = options.magnitude
    tau_c = options.tau_c
    tau_f = options.tau_f
    truncate = options.truncate
    Nm = options.Nm
    Hw = options.Hw
    P0 = options.P0
    P_scale = options.P_scale

    if in_file is not None:
        gd = ReadRaster(in_file)
        X = gd.X
        Y = gd.Y
        Orography = gd.RasterArray
    else:
        # Reproduce Fig 4c in SB2004
        dx = dy = 750.
        x, y = np.arange(-100e3, 200e3, dx), np.arange(-150e3, 150e3, dy)
        h_max = 500.
        x0 = -25e3
        y0 = 0
        sigma_x = sigma_y = 15e3
        X, Y = np.meshgrid(x, y)
        Orography = h_max * \
            np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) + ((Y - y0)**2 / (2 * sigma_y**2))))

    Theta_m = -6.5     # K / km
    rho_Sref = 7.4e-3  # kg m-3
    gamma = -5.8       # K / km

    physical_constants = dict()
    physical_constants['tau_c'] = tau_c      # conversion time [s]
    physical_constants['tau_f'] = tau_f      # fallout time [s]
    physical_constants['f'] = 2 * 7.2921e-5 * \
        np.sin(lat * np.pi / 180)  # Coriolis force
    physical_constants['Nm'] = Nm   # moist stability frequency [s-1]
    physical_constants['Cw'] = rho_Sref * Theta_m / \
        gamma  # uplift sensitivity factor [kg m-3]
    physical_constants['Hw'] = Hw         # vapor scale height
    # x-component of wind vector [m s-1]
    physical_constants['u'] = -np.sin(direction * 2 * np.pi / 360) * magnitude
    # y-component of wind vector [m s-1]
    physical_constants['v'] = np.cos(direction * 2 * np.pi / 360) * magnitude
    physical_constants['P0'] = P0   # background precip [mm hr-1]
    physical_constants['P_scale'] = P_scale   # precip scale factor [1]

    OP = OrographicPrecipitation(
        X,
        Y,
        Orography,
        physical_constants,
        truncate=truncate)
    P = OP.P
    units = OP.P_units

    if in_file is not None:
        array2raster(out_file, gd.geoTrans, gd.proj4, units, P)
    else:
        geoTrans = [0., OP.dx, 0., 0., 0., -OP.dy]
        array2raster(out_file, geoTrans, '', units, P)
