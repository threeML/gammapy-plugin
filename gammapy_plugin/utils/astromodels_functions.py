from astromodels.functions.function import Function1D, FunctionMeta
from threeML import u
import numpy as np


class Log_parabola_gammapy(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A log-parabolic function.
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

            desc : curvature (positive is concave, negative is convex)
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


# TODO latex
class Exp_cutoff_powerlaw_gammapy(Function1D, metaclass=FunctionMeta):
    r"""
    description :

        A exp cutoff  function.
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

        index :
            desc : index
            initial value : 2.0

        lambda_ :
            desc : cur
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
