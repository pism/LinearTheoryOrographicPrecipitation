"""Linear Theory of Orographic Precipitation (LTOP) model"""

import logging

import numpy as np
import pylab as plt
logger = logging.getLogger("LTOP")

np.seterr(divide="ignore", invalid="ignore")


class LTOP():
    "Linear Theory of Orographic Precipitation (LTOP) model"

    def __init__(
        self,
        tau_c: float = 1000.0,
        tau_f: float = 1000.0,
        P0: float = 0.0,
        P_scale: float = 1.0,
        Nm: float = 0.005,
        Hw: float = 2500,
        latitude: float = 45.0,
        direction: float = 270.0,
        speed: float = 15.0,
        Theta_m: float = -6.5,
        rho_Sref: float = 7.4e-3,
        gamma: float = -5.8,
    ):
        super().__init__()

        self.tau_c = tau_c  # "conversion time [s]"
        self.tau_f = tau_f  # "fallout time [s]"
        self.P0 = P0  # "Background precipitation rate [mm hr-1]"
        self.P_scale = P_scale  # "Precipitation scale factor"
        self.Nm = Nm  # "moist stability frequency [s-1]"
        self.Hw = Hw  # "Water vapor scale height [m]"
        self.latitude = latitude  # "Latitude used to compute the Coriolis force"
        self.direction = direction  # "Wind direction, 0 is north, 270 is west"
        self.speed = speed  # "Wind speed [m s-1]"
        self.Theta_m = Theta_m
        self.rho_Sref = rho_Sref
        self.gamma = gamma
        self.f: float = 0.0  # "Coriolis force"
        self.u: float = 0.0  # "u component of the wind velocity"
        self.v: float = 0.0  # "v component of the wind velocity"
        self.Cw: float = 0.0

        self.update()

    @property
    def tau_c(self):
        return self._tau_c

    @tau_c.setter
    def tau_c(self, value):
        self._tau_c = value

    @property
    def tau_f(self):
        return self._tau_f

    @tau_f.setter
    def tau_f(self, value):
        self._tau_f = value

    @property
    def P0(self):
        return self._P0

    @P0.setter
    def P0(self, value):
        self._P0 = value

    @property
    def P_scale(self):
        return self._P_scale

    @P_scale.setter
    def P_scale(self, value):
        self._P_scale = value

    @property
    def Nm(self):
        return self._Nm

    @Nm.setter
    def Nm(self, value):
        self._Nm = value

    @property
    def Hw(self):
        return self._Hw

    @Hw.setter
    def Hw(self, value):
        self._Hw = value

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, value):
        self._latitude = value

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, value):
        self._direction = value

    @property
    def Theta_m(self):
        return self._Theta_m

    @Theta_m.setter
    def Theta_m(self, value):
        self._Theta_m = value

    @property
    def rho_Sref(self):
        return self._rho_Sref

    @rho_Sref.setter
    def rho_Sref(self, value):
        self._rho_Sref = value

    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._gamma = value

    def run(
        self, orography: np.ndarray, dx: float, dy: float, truncate: bool = True
    ) -> np.ndarray:
        "Compute orographic precipitation in mm/hour."
        # make sure derived constants are up to date
        self.update()

        eps = 1e-18

        nrows, ncols = orography.shape

        pad = max(nrows, ncols)

        h = np.pad(orography, pad, "constant")
        nrows, ncols = h.shape

        h_hat = np.fft.fft2(h)

        x_freq = np.fft.fftfreq(ncols, dx / (2 * np.pi))
        y_freq = np.fft.fftfreq(nrows, dy / (2 * np.pi))

        kx, ky = np.meshgrid(x_freq, y_freq)

        # Intrinsic frequency sigma = U*k + V*l
        u0 = self.u
        v0 = self.v

        # $\sigma = U k + V l$, see paragraph after eq 5.
        sigma = u0 * kx + v0 * ky

        denominator = sigma**2 - self.f**2
        denominator[np.logical_and(np.fabs(denominator) < eps, denominator >= 0)] = eps
        denominator[np.logical_and(np.fabs(denominator) < eps, denominator < 0)] = -eps

        m_squared = (self.Nm**2 - sigma**2) * (kx**2 + ky**2) / denominator

        m = np.sqrt(np.array(m_squared, dtype=np.cdouble))

        # Regularization
        nonzero = np.logical_and(m_squared >= 0, sigma != 0)
        m[nonzero] *= np.sign(sigma[nonzero])

        P_hat = h_hat * (
            self.Cw
            * 1j
            * sigma
            / (
                (1 - 1j * m * self.Hw)
                * (1 + 1j * sigma * self.tau_c)
                * (1 + 1j * sigma * self.tau_f)
            )
        )

        # Convert from wave domain back to space domain
        P = np.real(np.fft.ifft2(P_hat))

        # Remove padding
        if pad > 0:
            P = P[pad:-pad, pad:-pad]

        # convert to mm hr-1
        P *= 3600

        # Add background precipitation
        P += self.P0

        # Truncate
        if truncate:
            P[P < 0] = 0.0

        P *= self.P_scale

        return P

    def update(self):
        "Update derived constants"

        self.f = 2 * 7.2921e-5 * np.sin(self.latitude * np.pi / 180.0)

        self.u = -np.sin(self.direction * 2 * np.pi / 360) * self.speed
        self.v = -np.cos(self.direction * 2 * np.pi / 360) * self.speed

        self.Cw = self.rho_Sref * self.Theta_m / self.gamma


def triangle_ridge_grid(dx: float = 5.0e4, dy: float = 5.0e4):
    "Allocate the grid for the synthetic geometry test."

    x_min, x_max = -100e3, 100e3
    y_min, y_max = -100e3, 100e3

    Mx = int((x_max - x_min) / dx) + 1
    My = int((y_max - y_min) / dy) + 1

    x = np.linspace(x_min, x_max, Mx)
    y = np.linspace(y_min, y_max, My)

    return x, dx, y, dy


def triangle_ridge(x: np.ndarray, A: float = 500.0, d: float = 50.0e3):
    "Create the 'triangle ridge' topography"
    return np.maximum(A * (1 - np.fabs(x) / d), 0)


def triangle_ridge_exact(
    x: np.ndarray, u: float, Cw: float, tau: float, A: float = 500.0, d: float = 50.0e3
):
    """The exact precipitation corresponding to the "triangle ridge" topography."""
    assert d > 0

    C = Cw * u * A / d
    Ut = u * tau

    xc = Ut * np.log(2 - np.exp(-d / Ut))

    def P(x):
        if x < 0 and x >= -d:
            return C * (1.0 - np.exp(-(x + d) / Ut))
        elif x >= 0 and x <= xc:
            return C * (np.exp(-x / Ut) * (2 - np.exp(-d / Ut)) - 1)
        else:
            return 0

    try:
        return 3600 * np.array([P(t) for t in x])
    except TypeError:
        return 3600 * P(x)


def max_error(spacing: float, direction: float):
    """Compute the maximum precipitation error compared to the "triangle ridge" exact
    solution.

    `spacing` : grid spacing, meters
    `direction` : wind direction, degrees

    """
    model = LTOP()
    # Set conversion time to zero (we could set fallout time to zero instead: it does not
    # matter which one is zero)
    model.tau_c = 0.0
    model.Hw = 0.0
    model.direction = direction
    model.latitude = 0.0

    if direction == 90.0 or direction == 270.0:
        # east or west
        x, dx, y, dy = triangle_ridge_grid(dx=spacing)
        t = x

        h = triangle_ridge(t)
        orography = np.tile(h, (len(y), 1))

        P = model.run(orography, dx, dy)
        P = P[len(y) // 2, :]
    else:
        # north or south
        x, dx, y, dy = triangle_ridge_grid(dy=spacing)
        t = y

        h = triangle_ridge(t)
        orography = np.tile(h, (len(x), 1)).T

        P = model.run(orography, dx, dy)
        P = P[:, len(x) // 2]

    if direction == 0 or direction == 90:
        P_exact = triangle_ridge_exact(-t, model.speed, model.Cw, model.tau_f)
    else:
        P_exact = triangle_ridge_exact(t, model.speed, model.Cw, model.tau_f)

    return np.max(np.fabs(P - P_exact))


def convergence_rate(dxs, error, direction, plot):
    """Compute and plot the convergence rate given the resinement path `dxs` and errors in
    `error`.

    Specify wind direction in `direction` (in degrees).

    Set `plot` to True to plot.

    """
    errors = [error(dx, direction) for dx in dxs]

    p = np.polyfit(np.log10(dxs), np.log10(errors), 1)

    if plot:
        wind_direction = {0: "north", 90: "east", 180: "south", 270: "west"}

        plt.figure()
        plt.title(
            "Precipitation errors (wind direction: {})".format(
                wind_direction[direction]
            )
        )
        log_fit_plot(dxs, p, "polynomial fit (dx^{:1.4})".format(p[0]))
        log_plot(dxs, errors, "o", "errors")
        plt.legend()
        plt.grid()
        plt.xlabel("grid spacing (meters)")
        plt.ylabel("log10(error)")
        plt.show()

    return p[0]


def ltop_test(plot=False):
    "Comparing to the 'triangle ridge' exact solution"
    dxs = [2000.0, 1000.0, 500.0, 250.0]

    assert convergence_rate(dxs, max_error, 0, plot) > 1.9
    assert convergence_rate(dxs, max_error, 90, plot) > 1.9
    assert convergence_rate(dxs, max_error, 180, plot) > 1.9
    assert convergence_rate(dxs, max_error, 270, plot) > 1.9


def gaussian_bump(
    xmin,
    xmax,
    ymin,
    ymax,
    dx,
    dy,
    h_max=500.0,
    x0=-25e3,
    y0=0.0,
    sigma_x=15e3,
    sigma_y=15e3,
):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    # Reproduce Fig 4c in SB2004
    x = np.arange(xmin, xmax, dx)
    y = np.arange(ymin, ymax, dy)
    X, Y = np.meshgrid(x, y)
    Orography = h_max * np.exp(
        -(((X - x0) ** 2 / (2 * sigma_x**2)) + ((Y - y0) ** 2 / (2 * sigma_y**2)))
    )
    return X, Y, Orography


if __name__ == "__main__":

    def log_plot(x, y, style, label):
        plt.plot(np.log10(x), np.log10(y), style, label=label)
        plt.xticks(np.log10(x), x)

    def log_fit_plot(x, p, label):
        plt.plot(np.log10(x), np.polyval(p, np.log10(x)), label=label)

    ltop_test(plot=True)
