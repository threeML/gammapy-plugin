import logging
import astropy.units as u
import numpy as np
from astromodels.functions.function import Function1D, Function2D, FunctionMeta
from astromodels.utils.angular_distance import angular_distance
from past.utils import old_div

log = logging.getLogger(__name__)


class Log_parabola_gammapy(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A log-parabolic function, same parametrization as Gammapy
    latex :
        $K\left(\frac{x}{piv}\right)^{-\alpha-\beta\log{\left(\frac{x}{piv}\right)}}$

    parameters :

        K :
            desc : Normalization
            initial value : 1e-11
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e5

        piv :
            desc : Pivot (keep this fixed)
            initial value : 1
            fix : yes

        alpha :

            desc : index
            initial value : 2.0

        beta :

            desc : curvature
            initial value : 1.0

    """

    def _set_units(self, x_unit, y_unit):
        # K has units of y

        self.K.unit = y_unit

        # piv has the same dimension as x
        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = u.dimensionless_unscaled
        self.beta.unit = u.dimensionless_unscaled

    def evaluate(self, x, K, piv, alpha, beta):

        # print("Receiving %s" % ([K, piv, alpha, beta]))

        xx = np.divide(x, piv)
        return K * np.power(xx, (-alpha - beta * np.log(xx)))


class Exp_cutoff_powerlaw_gammapy(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A exp cutoff  function.
    latex :
        $K\left(\frac{x}{piv}\right)^{-index}\exp{-(\lambda x)^\alpha}$

    parameters :

        K :
            desc : Normalization
            initial value : 1e-11
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e5

        piv :
            desc : Pivot (keep this fixed)
            initial value : 1
            fix : yes

        index :
            desc : index
            initial value : 2.0

        lambda_ :
            desc : curvature (= 1/xc)
            initial value : 0.1

        alpha :
            desc : alpha
            initial value : 1
            fix : yes

    """

    def _set_units(self, x_unit, y_unit):
        # K has units of y

        self.K.unit = y_unit

        # piv has the same dimension as x
        self.piv.unit = x_unit

        # alpha and beta are dimensionless
        self.alpha.unit = u.dimensionless_unscaled
        self.lambda_.unit = 1 / x_unit
        self.index.unit = u.dimensionless_unscaled

    def evaluate(self, x, K, piv, index, lambda_, alpha):

        # print("Receiving %s" % ([K, piv, alpha, beta]))

        xx = np.divide(x, piv)
        xy = np.multiply(x, lambda_)
        return K * np.power(xx, -index) * np.exp(-np.power(xy, alpha))


class Gaussian_on_sphere(Function2D, metaclass=FunctionMeta):
    r"""
    description :

        A bidimensional Gaussian function on a sphere (in spherical coordinates)

    latex : $$ f(\vec{x}) = \left(\frac{180^\circ}{\pi}\right)^2 \frac{1}{2\pi \sqrt{\det{\Sigma}}} \, {\rm exp}\left( -\frac{1}{2} (\vec{x}-\vec{x}_0)^\intercal \cdot \Sigma^{-1}\cdot (\vec{x}-\vec{x}_0)\right) \\ \vec{x}_0 = ({\rm RA}_0,{\rm Dec}_0)\\ \Lambda = \left( \begin{array}{cc} \sigma^2 & 0 \\ 0 & \sigma^2 (1-e^2) \end{array}\right) \\ U = \left( \begin{array}{cc} \cos \theta & -\sin \theta \\ \sin \theta & cos \theta \end{array}\right) \\\Sigma = U\Lambda U^\intercal $$

    parameters :

        lon0 :

            desc : Longitude of the center of the source
            initial value : 0.0
            min : 0.0
            max : 360.0

        lat0 :

            desc : Latitude of the center of the source
            initial value : 0.0
            min : -90.0
            max : 90.0

        sigma :

            desc : Standard deviation of the Gaussian distribution
            initial value : 10
            min : 0
            max : 20

    """

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.sigma.unit = x_unit

    def evaluate(self, x, y, lon0, lat0, sigma):

        lon, lat = x, y

        angsep = angular_distance(lon0, lat0, lon, lat)

        s2 = sigma**2

        return (
            (old_div(180, np.pi)) ** 2
            * 1
            / (2.0 * np.pi * s2)
            * np.exp(-0.5 * angsep**2 / s2)
        )

    def get_boundaries(self, max_sigma=None):

        # Truncate the gaussian at 2 times the max of sigma allowed
        if max_sigma is None:
            max_sigma = self.sigma.max_value

        min_lat = max(-90.0, self.lat0.value - 2 * max_sigma)
        max_lat = min(90.0, self.lat0.value + 2 * max_sigma)

        max_abs_lat = max(np.absolute(min_lat), np.absolute(max_lat))

        if (
            max_abs_lat > 89.0
            or 2 * max_sigma / np.cos(max_abs_lat * np.pi / 180.0) >= 180.0
        ):

            min_lon = 0.0
            max_lon = 360.0

        else:

            min_lon = self.lon0.value - 2 * max_sigma / np.cos(
                max_abs_lat * np.pi / 180.0
            )
            max_lon = self.lon0.value + 2 * max_sigma / np.cos(
                max_abs_lat * np.pi / 180.0
            )

            if min_lon < 0.0:

                min_lon = min_lon + 360.0

            elif max_lon > 360.0:

                max_lon = max_lon - 360.0

        return (min_lon, max_lon), (min_lat, max_lat)

    def get_total_spatial_integral(self, z=None, binsz=(None, None)):
        # TODO this is an unprecise and very slow solution
        if self.sigma.value < 1e-1:
            msg = f"Your sigma value {self.sigma.value} is fairly small. "
            msg += "This may lead to a really high memory usage"
            log.warning(msg)
        if binsz == (None, None):
            binsz = (0.1 * self.sigma.value, 0.1 * self.sigma.value)
        if self.sigma.value * 20 < self.sigma.max_value:
            (min_l, max_l), (min_b, max_b) = self.get_boundaries(
                max_sigma=self.sigma.value * 20
            )
        else:
            (min_l, max_l), (min_b, max_b) = self.get_boundaries()
        if min_l > max_l:
            min_l -= 360
        lons = np.linspace(
            min_l, max_l, np.ceil((max_l - min_l) / binsz[0]).astype(int)
        )
        lons[lons > 360] -= 360
        lats = np.linspace(
            min_b, max_b, np.ceil((max_b - min_b) / binsz[1]).astype(int)
        )
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # Flatten the grid for evaluation
        lon_flat = lon_grid.ravel()
        lat_flat = lat_grid.ravel()

        # Evaluate the function at all grid points
        values = np.array(
            self.evaluate(
                lon_flat,
                lat_flat,
                self.lon0.value,
                self.lat0.value,
                self.sigma.value,
            )
        )

        # Area element in degrees^2 on the sphere
        dlon = (max_l - min_l) / (lons.shape[0] - 1)
        dlat = (max_b - min_b) / (lats.shape[0] - 1)

        cos_lat = np.cos(np.radians(lat_flat))
        area_elements = dlat * dlon * cos_lat
        if isinstance(z, u.Quantity):
            z = z.value
        total = np.sum(values * area_elements) * np.ones_like(z)
        return total
