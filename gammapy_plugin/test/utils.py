from threeML.analysis_results import BayesianResults, MLEResults


def get_close(threeML_results, gammapy_result_dict):
    if isinstance(threeML_results, BayesianResults):
        bm = threeML_results.get_median_fit_model().free_parameters
    elif isinstance(threeML_results, MLEResults):
        bm = threeML_results.optimized_model.free_parameters
    else:
        raise NotImplementedError
    hdp = threeML_results.get_data_frame(error_type="hpd")
    for p in bm.keys():
        pn = p.split(".")[-1]
        if pn == "K":
            pn = "amplitude"
        for gp in gammapy_result_dict["spectral"]["parameters"]:
            if gp["name"] == pn:
                break
        if pn == "amplitude":
            min_v = (
                ((hdp.loc[p]["negative_error"] + bm[p].value) * bm[p].unit)
                .to("TeV-1 cm-2 s-1")
                .value
            )
            max_v = (
                ((hdp.loc[p]["positive_error"] + bm[p].value) * bm[p].unit)
                .to("TeV-1 cm-2 s-1")
                .value
            )
        elif pn == "alpha":
            min_v = hdp.loc[p]["negative_error"] - bm[p].value
            max_v = hdp.loc[p]["positive_error"] - bm[p].value
        else:
            min_v = hdp.loc[p]["negative_error"] + bm[p].value
            max_v = hdp.loc[p]["positive_error"] + bm[p].value
        return bool(gp["value"] <= max_v and min_v <= gp["value"])
