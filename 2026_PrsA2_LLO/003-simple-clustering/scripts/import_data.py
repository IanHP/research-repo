from pathlib import Path

import fastparquet
import numpy as np
import pandas as pd
from numpy import typing as npt


def feature_types(feature_file: Path) -> list[list[str]]:
    """Will read in the ranked feature file to determine the interactions being used as features
    Ex: if ARG_121 and HIS_132 is an interaction used as a feature, will add [..,[121,132],...]

    Args:
        feature_file: location of ranked feature file

    Returns:
        list of each pair of residues that make up an interaction
        Ex: ["HIS_192","ASP_312"]
    """
    # read in and manipulate list of features file
    with open(feature_file, "r") as f:
        text: str = f.read()
    line_split: list[str] = text.split("\n")
    residue_list: list[list[str]] = [
        item.split(", ")[0].split(" - ") for item in line_split[1:]
    ]
    
    feature_types: list[list[str]] = []
    for line_index in range(len(residue_list)):
        if len(residue_list[line_index]) > 1:
            feature_types.append([])
            feature_types[-1].append(
                "_".join(residue_list[line_index][0].split("_")[0:2])
            )
            feature_types[-1].append(
                "_".join(residue_list[line_index][1].split("_")[0:2])
            )

    return feature_types


def read_in_parq(distance_path: Path, run_type: str, run_num: str, feature_list: list[list[str]]
                ) -> npt.NDArray[np.float64]:
    """Will read in a parquet file, and based on runType, runNum, will return each frames (row) interactions (col) in feature list

    Args:
        distance_path: where parquets are stored
        run_type: run type (ex: ph5) of interest
        run_num: run num of interest
        feature_list: list of each pair of residues that make up an interaction

    Returns:
        2D list. Rows are frames, columns each interaction, each value is distance
    """

    base_path: Path = distance_path / run_type
    features: list[list[float]] = []

    for interact_pair in feature_list:
        his: str = interact_pair[0]
        charge: str = interact_pair[1]
        file_name: Path = base_path / f"{run_num}_{his}.prq"
        df: pd.DataFrame = fastparquet.ParquetFile(file_name).to_pandas([charge])
        features.append(df[charge].tolist())

    features: npt.NDArray[np.float64] = np.array(features).T

    return features


def import_features(distance_path: Path, ranked_path: Path, 
                    run_types: list[str], run_nums: list[str]
                    ) -> dict[str, dict[str, npt.NDArray[np.float64]]]:
    """Takes in folder that holds all interactions distance + list of interaction types being analyzed.
    Will then find each run's features (interaction distances) for each frame

    Args:
        distance_path: folder that holds the run type folders (ph5 / ph7)
        ranked_path: csv file that holds interactions used as features

    Returns:
        2D dictionary for each run type (ph), then run (run #). Then there is a matrix  for that run
        where columns are interaction types and rows are frames.
    """

    # key variables

    # import list of feature types
    feature_list: list[list[str]] = feature_types(Path(ranked_path))

    # create runType / runNum double dictionary with each runs interaction data in a 2D matrix
    features: dict[str, dict[str, npt.NDArray[np.float64]]] = {}
    for run_type in run_types:
        features[run_type] = {}
        for run_num in run_nums:
            data: npt.NDArray[np.float64] = read_in_parq(
                Path(distance_path), run_type, run_num, feature_list
            )
            features[run_type][run_num] = data

    return features


def frame_dict_to_all(features_dict: dict[str, dict[str, npt.NDArray[np.float64]]],
                     ) -> tuple[npt.NDArray[np.float64], list[int]]:
    """will take the feature dictionary that holds all frames features for each ph / run number, and put it all into one
    np.array

    Args:
        features_dict: first is run type, then run num. Then holds all interaction
            distances for each frame

    Returns:
        every frame / interaction distance put into one array
        where each run_type / run_num starts in all_frame list
    """

    # run types and run nums
    run_types: list[str] = list(features_dict.keys())
    run_nums: list[str] = list(features_dict[run_types[0]].keys())

    # determine total number of frames
    frame_count: int = 0

    for run_type in run_types:
        for run_num in run_nums:
            frame_count = frame_count + features_dict[run_type][run_num].shape[0]

    interact_count: int = features_dict[run_type][run_num].shape[1]

    # create list of all frames together
    all_frames_features: npt.NDArray[np.float64] = np.empty(
        (frame_count, interact_count), dtype=np.float64
    )
    run_splits: list[
        int
    ] = []  # each number is the start of each frame. In order of runType / runNum
    logical_length: int = 0

    for run_type in run_types:
        for run_num in run_nums:
            run_splits.append(logical_length)
            frame_count: int = len(features_dict[run_type][run_num])
            all_frames_features[logical_length : logical_length + frame_count, :] = (
                features_dict[run_type][run_num]
            )

            logical_length = logical_length + frame_count

    run_splits.append(logical_length)

    return all_frames_features, run_splits


def all_frames_to_dict(
    all_frames_features: npt.NDArray[np.float64] | list[int],
    run_splits: list[int],
    run_types: list[str],
    run_nums: list[str],
) -> dict[str, dict[str, npt.NDArray[np.float64]]]:
    """Will take in matrix with every frame, and seperate into dictionary based on
    what run_type / run_num it is from

    Args:
        all_frames_features: matrix where rows are frames and columns are data
            about that frame
        run_splits: list of start of each different traj in list
        run_types: list of all run types
        run_nums: list of all run numbers

    Returns:
        2D dictionary for each run type (ph), then run (run #). Then there is a matrix  for that run
            where columns are interaction types and rows are frames.
    """

    split_index: int = 0
    features_dict: dict[str, dict[str, npt.NDArray[np.float64]]] = {}

    for run_type in run_types:
        features_dict[run_type] = {}
        for run_num in run_nums:
            features_dict[run_type][run_num] = np.array(
                all_frames_features[
                    run_splits[split_index] : run_splits[split_index + 1]
                ]
            )
            split_index += 1

    return features_dict
