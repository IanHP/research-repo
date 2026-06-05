from pathlib import Path

import import_data
import matplotlib.pyplot as plt
import numpy as np
from numpy import typing as npt
from sklearn.cluster import HDBSCAN

SEED: int = 600411609

DIR_SCRIPT: Path = Path(__file__).parent.resolve()
DIR_STUDY: Path = Path(DIR_SCRIPT / ".." / ".." / "..").resolve()


def hdbscan_fit(features: npt.NDArray[np.float64],
                kwargs_clustering: None | dict[str, int | None | float | str | bool] = None
                ) -> tuple[list[int], npt.NDArray[np.float64], int]:
    """Will take in a list of each frame's features and hyper-parameters,
    and determine which cluster each frame is apart of

    Args:
        features: array of each frame's features

    Returns:
        Each frame's cluster, index aligned
    """
    
    kwargs_clustering_default: dict[str, int | None | float | str | bool] = {
        "min_cluster_size": 5,
        "min_samples": None,
        "cluster_selection_epsilon": 0.0,
        "max_cluster_size": None,
        "metric": "euclidean",
        "metric_params": None,
        "alpha": 1,
        "algorithm": "auto",
        "leaf_size": 40,
        "n_jobs": None,
        "cluster_selection_method": "eom",
        "allow_single_cluster": False,
        "store_centers": "medoid",
        "copy": False,
    }
    if kwargs_clustering is None:
        kwargs_clustering: dict[str, int | None | float | str] = kwargs_clustering_default
    else:
        kwargs_clustering_default.update(kwargs_clustering)
        kwargs_clustering = kwargs_clustering_default

    assert kwargs_clustering is not None

    hdb = HDBSCAN(**kwargs_clustering)
    frame_labels: list[int] = hdb.fit_predict(features)
    cluster_num: int = max(frame_labels)

    centroids: npt.NDArray[np.float64] = hdb.medoids_
    centroid_frames: list[int] = [-1] * (cluster_num + 1)
    for group_num, centroid in enumerate(centroids):
        for frame_num, feature in enumerate(features):
            if feature[0] == centroid[0] and feature[5] == centroid[5]:
                centroid_frames[group_num] = frame_num

    return frame_labels, centroid_frames, cluster_num


def pca_all_tog(embedding: dict[str, dict[str, npt.NDArray[np.float64]]]):
    """Will take in a list of frames (with compressed 2D points), organized in a 2D array of phs and runs
    and plot them

    Args:
        embedding: each frame's location in 2D space organized by ph and run
    """

    colors: list[str] = ["#F0E442", "#E69F00", "#D55E00", "#035353", "#0072B2", "#56B4E9"]

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

    plt.show()


def print_clusters(embedding: dict[str, dict[str, npt.NDArray[np.float64]]],
                   cluster_num: int,centroids: npt.NDArray[np.float64],
                   output_file: Path):
    """Will take in a list of frames (with compressed 2D points), organized in a 2D array of phs and runs
    and plot them in seperate graphs

    Args:
        embedding: each frame's location in 2D space organized by ph and run
    """

    clusters: list[list[str]] = []
    for i in range(cluster_num + 1):
        clusters.append([])

    # Create a list of frames in each cluster
    for type_index, (run_type, runs) in enumerate(embedding.items()):
        for run_index, (run_num, data) in enumerate(runs.items()):
            for frame_num, cluster in enumerate(data):
                if cluster >= 0:
                    clusters[cluster].append(f"{run_type}-{run_num} {frame_num}")

    # Write out frames in each cluster + centroids
    with open(output_file, "w") as f:
        for cluster, frames in enumerate(clusters):
            f.write(f"Cluster {cluster}")
            for frame in frames:
                f.write(f",{frame}")
            f.write("\n")

        f.write("\n\nCentroids:")
        for group, frame in enumerate(centroids):
            f.write(f"\n{group},{frame}")

        f.write("\n\n\nComposition:")
        for group, frames in enumerate(clusters):
            type_set: set[str] = set()
            for frame in frames:
                type_set.add(frame.split(" ")[0])
            f.write(f"\n{group},{",".join(type_set)}")


def frame_to_run(run_splits: list[int], frame: int, run_types: list[str], 
                 run_nums: list[str]) -> str:
    """Will take a frame and return what run_type (pH) and run # it is from

    Args:
        run_splits: how the run types are seperated
        frame: the frame number being checked
        run_types: the different run types in a list Ex: ["pH5","pH7"]
        run_nums: list of different run numbers in a list Ex: ["r1","r2"]

    Returns:
        which run_type / run it is apart of in string form.
        Ex: ph5-r1 1000
    """
    type_index: int = 0
    for run_type in run_types:
        for run_num in run_nums:
            curr_min: int = run_splits[type_index]
            curr_max: int = run_splits[type_index + 1]
            if frame >= curr_min and frame < curr_max:
                return f"{run_type}-{run_num} {frame - curr_min}"
            type_index = type_index + 1

    return ""


def main():
    # get features
    run_types: list[str] = ["ph5", "ph7"]  # the different types of runs
    run_nums: list[str] = ["r1","r2","r3"]

    distance_path: Path = (
        DIR_STUDY / "analysis" / "002-finding-descriptors" / "data" / "all_distances"
    )
    """Outer file where the parquets are held"""
    ranked_path = (
        DIR_STUDY
        / "analysis"
        / "002-finding-descriptors"
        / "scripts"
        / "004-rank-interactions"
        / "ranked_interactions.csv"
    )
    """where the rankings are held"""
    output_file: Path = DIR_SCRIPT / ".." / "data" / "hdbscan_cluster_leaf_m25_eps_3.csv"

    feature_dict: dict[str, dict[str, npt.NDArray[np.float64]]] = (
        import_data.import_features(distance_path, ranked_path,
                                    run_types, run_nums)
        )

    feature_list, run_splits = import_data.frame_dict_to_all(feature_dict)

    # Go to each run_type, run, and frame and find PCA coordinates
    clustered_frames, centroids, cluster_num = hdbscan_fit(
        feature_list, kwargs_clustering={"cluster_selection_method":"leaf",
                                         "min_cluster_size":25,
                                         "cluster_selection_epsilon":14.20})  # each frame is now just 2 data points / 2 columns. x var, y var
    clustered_frames = import_data.all_frames_to_dict(
        clustered_frames, run_splits, run_types, run_nums)

    # get the actual location of each centroid
    for group in range(len(centroids)):
        centroid: int = centroids[group]
        centroids[group] = frame_to_run(run_splits, centroid, run_types, run_nums)

    # plot
    print_clusters(clustered_frames, cluster_num, centroids, output_file)



if __name__ == "__main__":
    main()
