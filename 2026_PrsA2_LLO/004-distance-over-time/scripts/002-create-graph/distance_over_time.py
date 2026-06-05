
import fastparquet
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from residue_conversion import *

SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()


def main(his_residue: str, charged_residue: str, run_types: list[str], runs: list[str], rolling_num: int,
        parquet_dir: Path, figure_dir: Path):
    """Will create a graph of the distances of the interaction over time

    Args:
        his_residue (str): the histidine residue in the interaction
        charged_residue (str): the other residue in the interaction
        run_types (list[str]): all the types of runs Ex: ph5
        runs (list[str]): how many runs of each type. Name in file storage. Ex: r1
        rolling_num (int): how large the rolling average is based on.
        parquet_dir (Path): Path that holds are distances of interactions
        figure_dir (Path): Name of file where figure will be placed
    """
    # get all interactions distances from all runs
    all_distances: dict[str, dict[str, list[float]]] = {}
    for run_type in run_types:
        all_distances[run_type] = {}
        for run in runs:
            base_dist: list[float] = readInParq(his_residue, charged_residue, run_type, run, parquet_dir)
            all_distances[run_type][run] = roll_dist(base_dist, rolling_num)

    # plot it
    plt.figure(figsize=(12, 6))
    colorss = [["#F0E442", "#E69F00", "#D55E00"],["#035353", "#0072B2", "#56B4E9"]]
    for run_type, colors in zip(run_types, colorss):
        for run, color in zip(runs,colors):
            final_y: list[float] = all_distances[run_type][run]
            final_x: list[float] = create_x_values(final_y)
            plt.plot(final_x, final_y, label=f'{run_type}-{run}', color=color)

    plt.title(f'{his_residue} ({prot_name_from_res_id(his_residue)}) to {charged_residue} ({prot_name_from_res_id(charged_residue)}) - Average Distance Over Time', fontsize=14)
    plt.xlabel('Time (nanoseconds)', fontsize=12)
    plt.ylabel('Average Distance (Å)', fontsize=12)
    #plt.xticks(fontsize=10)
    #plt.yticks(fontsize=10)
    plt.legend(fontsize=12)
    #plt.xlim(0, 5)
    plt.grid(True)
    plt.savefig(figure_dir, format='png')
    plt.close()
    #plt.show()

def roll_dist(base_dists: list[float], rolling_num: int) -> list[float]:
    """Will take in a list of numbers and create rolling average for each
    index

    Args:
        base_dists (list[float]): the list being rolling averaged
        rolling_num (int): how large the average is based on. How many
            frames on each side of index are included

    Returns:
        list[float]: the averaged list
    """
    avg_dist: list[float] = []
    for index, dist in enumerate(base_dists):
        min_ind: int = index - rolling_num 
        if(min_ind < 0):
            min_ind = 0
        max_ind: int = index + rolling_num 
        if(max_ind > len(base_dists)-1):
            max_ind = len(base_dists)-1
        
        avg: float = sum(base_dists[min_ind:max_ind+1]) / (max_ind - min_ind + 1)
        avg_dist.append(avg)
    return avg_dist


def create_x_values(distance_per_frame: list[float]) -> list[float]:
    """Will take in a list of interaction distances per frame and
    determine the ns of each frame

    Args:
        distance_per_frame (list[float]): list of distances per frame

    Returns:
        int: length of simulation in ns
    """
    # Each frame is 0.2ns
    x_axis: list[float] = []
    for frame in range(len(distance_per_frame)):
        x_axis.append(frame * 0.2 + 50)
    return x_axis


def readInParq(his_res: str, charge_res: str, run_type: str, run_num: str, parquet_dir: Path) -> list[float]:
    """Will read in a parquet file, and based on run_type / his / charge, will return all runs distance data

    Args:
        hisRes (String): histidine res of interest
        chargeRes (String): Charged intraction of interest
        runType (String): run type (ex: ph5) of interest

    Returns:
         list of int: list of all distances for interaction across all runs
    """
    ret: list[float] = []
    file_path: Path = Path(parquet_dir / run_type / f"{run_num}_{his_res}.prq").resolve()
    df = fastparquet.ParquetFile(file_path).to_pandas([charge_res])
    ret.extend(df[charge_res].tolist())
    return ret
    
if __name__ == "__main__":

    parquet_dir: Path = Path(DIR_SCRIPT / ".." / ".." / ".." / "002-finding-descriptors" / "data" / "interaction_distances").resolve()
    """The directory were parquet files all stored"""
    run_types: list[str] = ["ph5","ph7"]
    """All the runtypes. Ex: ph5"""
    runs: list[str] = ["r1","r2","r3"]
    """The runs of each run type. Ex" r1"""
    marker_interact_path: Path = Path(DIR_SCRIPT  / ".." / ".." / "data" / "marker_interactions.csv").resolve()
    """Interactions creating graphs for"""
    real_interact_path: Path = Path(DIR_SCRIPT  / ".." / ".." / "data" / "real_interactions.csv").resolve()
    """Interactions creating graphs for"""
    rolling_num: int = 20
    """How many frames on either side are included in each point"""

    for interactions_path in [marker_interact_path, real_interact_path]:
        with open(interactions_path, "r") as f:
            all_text:list[str] = f.read().split("\n")
            for line in all_text:
                ress: list[str] = line.split(",")[0].split(" - ")
                if(len(ress) == 2):
                    his_residue: str = "_".join(ress[0].split("_")[0:2])
                    """This is the histidine of the interaction being graphed"""
                    charged_residue: str = "_".join(ress[1].split("_")[0:2])
                    """This is the other residue of the interaction being graphed"""
                    figure_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "figures" / interactions_path.stem 
                                            / f"{his_residue}_{charged_residue}_dist_r{2*rolling_num+1}.png").resolve()
                    """Where the figure will be stored"""
                    if not figure_dir.parent.is_dir():
                        figure_dir.parent.mkdir(parents=True, exist_ok=True)

                    main(his_residue, charged_residue, run_types, runs, rolling_num, parquet_dir, figure_dir)