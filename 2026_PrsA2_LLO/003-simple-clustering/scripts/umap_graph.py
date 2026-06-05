from pathlib import Path

import import_data
import matplotlib.pyplot as plt
import numpy as np
import umap
from numpy import typing as npt
from sklearn.preprocessing import StandardScaler

SEED: int = 600411609
DIR_SCRIPT: Path = Path(__file__).parent.resolve()

def umap_fit(features: npt.NDArray[np.float64], 
            kwargs_umap: None | dict = None
            ) -> npt.NDArray[np.float64]:
    """Will take in a list of each frame's features and hyper-parameters,
    and project each frame onto the 2D plane

    Args:
        features (NDArray): array of each frame's features
        kwargs_umap (dict): dictionary holding edited UMAP hyperparameters

    Returns:
        Each frame's N-features crunched to 2 dimensions
    """
    # ensure PCA hyperparameters are correct
    kwargs_umap_default: dict = {
        "n_neighbors":12, 
        "n_components":2, 
        "metric":'euclidean', 
        "metric_kwds":None, 
        "output_metric":'euclidean', 
        "output_metric_kwds":None, 
        "n_epochs":None, 
        "learning_rate":1.0, 
        "init":'spectral', 
        "min_dist":0.1, 
        "spread":1.0, 
        "low_memory":True, 
        "n_jobs":1, 
        "set_op_mix_ratio":1.0, 
        "local_connectivity":1.0, 
        "repulsion_strength":1.0, 
        "negative_sample_rate":5, 
        "transform_queue_size":4.0, 
        "a":None, 
        "b":None, 
        "random_state":SEED, 
        "angular_rp_forest":False, 
        "target_n_neighbors":-1, 
        "target_metric":'categorical', 
        "target_metric_kwds":None, 
        "target_weight":0.5, 
        "transform_seed":42, 
        "transform_mode":'embedding', 
        "force_approximation_algorithm":False, 
        "verbose":False, 
        "tqdm_kwds":None, 
        "unique":False, 
        "densmap":False, 
        "dens_lambda":2.0, 
        "dens_frac":0.3, 
        "dens_var_shift":0.1, 
        "output_dens":False, 
        "disconnection_distance":None, 
        "precomputed_knn":(None, None, None)
    }
    
    if kwargs_umap is None:
        kwargs_umap: dict = kwargs_umap_default
    else:
        kwargs_umap_default.update(kwargs_umap)
        kwargs_umap = kwargs_umap_default

    assert kwargs_umap is not None
    
    reducer = umap.UMAP(**kwargs_umap)
    scaled_data: npt.NDArray[np.float64] = StandardScaler().fit_transform(features)
    embedding: npt.NDArray[np.float64] = reducer.fit_transform(scaled_data)

    return embedding


# colors = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
def umap_all_tog(embedding: dict[str, dict[str, npt.NDArray[np.float64]]], plot_output: Path):
    """Will take in a list of frames (with compressed 2D points), organized in a 2D array of phs and runs
    and plot them

    Args:
        embedding (dict): each frame's location in 2D space organized by ph and run
        Plot_output (Path): where the plots are stored
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

    plt.title("All Run UMAP", fontsize=20)
    plt.legend()

    if not plot_output.parent.is_dir():
        plot_output.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(plot_output)


def umap_all_sep(embedding: dict[str, dict[str, npt.NDArray[np.float64]]], plot_output: Path):
    """Will take in a list of frames (with compressed 2D points), organized in a 2D array of phs and runs
    and plot them in seperate graphs

    Args:
        embedding (dict): each frame's location in 2D space organized by ph and run
        plot_output (Path): where UMAP graph will be saved
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
                f"{run_type}-{run_num} UMAP", fontsize=20
            )
            ax[type_index][run_index].set_ylim(bottom=-15, top=20)
            ax[type_index][run_index].set_xlim(left=-10, right=25)

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
    ranked_path: Path = Path(DIR_SCRIPT / ".." / ".." / "004-distance-over-time" / "data" / "marker_interactions.csv")
    """where the rankings are held"""
    plot_output: Path = Path(DIR_SCRIPT / ".." / "figures" / "umap_minter_graph.png")
    """where PCA plot will be saved"""

    feature_dict: dict[str, dict[str, npt.NDArray[np.float64]]] = (
        import_data.import_features(distance_path, ranked_path,
                                    run_types, run_nums)
        )
    feature_list, run_splits = import_data.frame_dict_to_all(feature_dict)

    # UMAP hyperparameters
    min_distance = 0.1  # how close two different frames can be to each other
    neighbor_num = 12  # how big "neighborhoods" or similar clusters are. Want to have distinct largish neighborhoods, but not be one big blob

    # Go to each run_type, run, and frame and find UMAP coordinates
    umap_data = umap_fit(
        feature_list, 
        {"min_dist":min_distance, "n_neighbors":neighbor_num}
    ) 
    umap_data = import_data.all_frames_to_dict(
        umap_data, run_splits, run_types, run_nums
    )

    # plot
    umap_all_tog(umap_data, plot_output)


if __name__ == "__main__":
    main()
