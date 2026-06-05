from pathlib import Path

import import_data
import matplotlib.pyplot as plt
import numpy as np
from numpy import typing as npt
from sklearn.decomposition import PCA

SEED: int = 600411609
DIR_SCRIPT: Path = Path(__file__).parent.resolve()

def pca_fit(features: npt.NDArray[np.float64],
            kwargs_pca: None | dict[str, int | None | float | str | bool] = None
            ) -> npt.NDArray[np.float64]:
    """Will take in a list of each frame's features and hyper-parameters,
    and project each frame onto the 2D plane

    Args:
        features (NDArray): array of each frame's features
        kwargs (dict): stores non-default pca hyperparameters being edited

    Returns:
        Each frame's N-features crunched to 2 dimensions
    """
    # ensure PCA hyperparameters are correct
    kwargs_pca_default: dict[str, int | None | float | str | bool] = {
        "n_components":2, 
        "copy":True, 
        "whiten":False, 
        "svd_solver":'auto', 
        "tol":0.0, 
        "iterated_power":'auto', 
        "n_oversamples":10, 
        "power_iteration_normalizer":'auto', 
        "random_state":SEED
    }
    
    if kwargs_pca is None:
        kwargs_pca: dict[str, int | None | float | str] = kwargs_pca_default
    else:
        kwargs_pca_default.update(kwargs_pca)
        kwargs_pca = kwargs_pca_default
    assert kwargs_pca is not None
    
    reducer = PCA(**kwargs_pca)
    embedding: npt.NDArray[np.float64] = reducer.fit_transform(features)

    return embedding


# colors = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
def pca_all_tog(embedding: dict[str, dict[str, npt.NDArray[np.float64]]], plot_output: Path):
    """Will take in a list of frames (with compressed 2D points), organized in a 2D array of phs and runs
    and plot them

    Args:
        embedding: each frame's location in 2D space organized by ph and run
        plot_output (path): where plots are stored
    """

    colors = ["#F0E442", "#E69F00", "#D55E00", "#035353", "#0072B2", "#56B4E9"]

    colorIndex = 0
    for run_type, runs in embedding.items():
        for run_num, data in runs.items():
            plt.scatter(
                data[:, 0],  # gets a list of x vars
                data[:, 1],  # gets a list of y vars
                c=colors[colorIndex],
                label=f"{run_type}-{run_num}",
            )
            colorIndex = colorIndex + 1

    plt.title("All Run PCA", fontsize=20)
    plt.legend()

    if not plot_output.parent.is_dir():
        plot_output.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(plot_output)


def pca_all_sep(embedding: dict[str, dict[str, npt.NDArray[np.float64]]], plot_output: Path):
    """Will take in a list of frames (with compressed 2D points), organized in a 2D array of phs and runs
    and plot them in seperate graphs

    Args:
        embedding (dict): each frame's location in 2D space organized by ph and run
        plot_output (Path): where PCA graph will be saved
    """
    run_types: list[str] = list(embedding.keys())
    run_nums: list[str] = list(embedding[run_types[0]].keys())

    fig, ax = plt.subplots(nrows=len(run_types), ncols=len(run_nums), figsize=(20, 8))

    colors = ["#F0E442", "#E69F00", "#D55E00", "#035353", "#0072B2", "#56B4E9"]

    colorIndex = 0
    for type_index, (run_type, runs) in enumerate(embedding.items()):
        for run_index, (run_num, data) in enumerate(runs.items()):
            ax[type_index][run_index].scatter(
                data[:, 0],  # gets a list of x vars
                data[:, 1],  # gets a list of y vars
                c=colors[colorIndex],
            )
            colorIndex = colorIndex + 1
            ax[type_index][run_index].set_title(
                f"{run_type}-{run_num} PCA", fontsize=20
            )
            ax[type_index][run_index].set_ylim(bottom=-20, top=25)
            ax[type_index][run_index].set_xlim(left=-25, right=35)
    
    if not plot_output.parent.is_dir():
        plot_output.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(plot_output)


def main():
    # get features
    run_types: list[str] = ["ph5", "ph7"]  # the different types of runs
    run_nums: list[str] = ["r1", "r2", "r3"]
    distance_path: Path = Path(DIR_SCRIPT / ".." / ".." / "002-finding-descriptors" / 
                               "data" / "interaction_distances")
    """where distance parquets are stored"""
    #ranked_path: Path = Path(DIR_SCRIPT / ".." / ".." / "002-finding-descriptors" / "data" / "ranked_interactions.csv")
    """where the rankings are held"""
    ranked_path: Path = Path(DIR_SCRIPT / ".." / ".." / "004-distance-over-time" / "data" / "marker_interactions.csv")

    plot_output_tog: Path = Path(DIR_SCRIPT / ".." / "figures" / "pca_tog_minter_graph.png")
    plot_output_sep: Path = Path(DIR_SCRIPT / ".." / "figures" / "pca_sep_minter_graph.png")

    """where PCA plot will be saved"""

    feature_dict: dict[str, dict[str, npt.NDArray[np.float64]]] = (
        import_data.import_features(distance_path, ranked_path,
                                    run_types, run_nums)
        )
    feature_list, run_splits = import_data.frame_dict_to_all(feature_dict)


    # Go to each run_type, run, and frame and find PCA coordinates
    pca_data = pca_fit(feature_list, kwargs_pca=None)
    pca_data = import_data.all_frames_to_dict(pca_data, run_splits, run_types, run_nums)

    # plot
    pca_all_tog(pca_data, plot_output_tog)
    pca_all_sep(pca_data, plot_output_sep)


if __name__ == "__main__":
    main()
    
