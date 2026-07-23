import numpy as np
from threeML.io.plotting.data_residual_plot import ResidualPlot


def plot_model(
    datasets,
    *args,
    **kwargs,
) -> "ResidualPlot":

    residual_plot = ResidualPlot(
        **kwargs,
    )
    for i in range(len(datasets)):

        y_unweighted = datasets[i].counts.data.reshape(-1)
        x = datasets[i].counts.geom.axes["energy"].as_plot_center.to("keV").value
        xerr = [
            datasets[i].counts.geom.axes["energy"].as_plot_xerr[j].to("keV").value
            for j in [0, 1]
        ]

        bins = datasets[i].counts.geom.axes["energy"].as_plot_edges.to("keV").value
        widths = np.diff(bins)
        y = y_unweighted / widths

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
            datasets[i].npred().get_spectrum().data.reshape(-1) / widths,
            label=kwargs.get("model_label", "Expected"),
        )

    return residual_plot.finalize(
        xlabel="Energy\n(keV)",
        ylabel="Counts/keV",
        xscale="log",
        yscale="log",
        show_legend=kwargs.get("show_legend", True),
    )
