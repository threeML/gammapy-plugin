import astropy.units as u
import numpy as np
from threeML.io.plotting.data_residual_plot import ResidualPlot

from typing import TYPE_CHECK

if TYPE_CHECK:
    from gammapy.datasets import Datasets


def plot_model(
    datasets: Datasets,
    *args,
    **kwargs,
) -> ResidualPlot:
    """
    Plot the model and data for a given datasets object.

    :param datasets: Gammapy Datasets object containing the data and model
    :type datasets: Datasets
    :return: ResidualPlot object containing the plot
    """

    residual_plot = ResidualPlot(
        **kwargs,
    )
    for i in range(len(datasets)):

        y_unweighted = datasets[i].counts.get_spectrum().data.reshape(-1)
        x = datasets[i].counts.geom.axes["energy"].as_plot_center.to("keV").value
        xerr = [
            datasets[i].counts.geom.axes["energy"].as_plot_xerr[j].to("keV").value
            for j in [0, 1]
        ]

        bins = datasets[i].counts.geom.axes["energy"].as_plot_edges.to("keV").value
        widths = np.diff(bins)
        y = y_unweighted / widths
        y /= datasets[i].gti.time_sum.to(u.s).value

        residuals = (
            datasets[i].counts.get_spectrum() - datasets[i].npred().get_spectrum()
        ) / datasets[i].npred().get_spectrum()
        residuals = residuals.data.reshape(-1)

        residual_plot.add_data(
            x,
            y,
            residuals,
            xerr=xerr,
            label=datasets[i].name,
            show_data=kwargs.get("show_data", True),
        )

        residual_plot.add_model(
            x,
            datasets[i].npred().get_spectrum().data.reshape(-1)
            / (widths * datasets[i].gti.time_sum.to(u.s).value),
            label=kwargs.get("model_label", "Expected"),
        )

    return residual_plot.finalize(
        xlabel="Energy\n(keV)",
        ylabel="Counts/keV/s",
        xscale="log",
        yscale="log",
        show_legend=kwargs.get("show_legend", True),
    )
