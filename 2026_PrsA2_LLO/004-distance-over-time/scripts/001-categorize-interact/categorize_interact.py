from matplotlib.pylab import real
from pathlib import Path
import fastparquet
import math

SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()



def main(run_types: list[str], runs: list[str], parquet_dir: Path, ranked_inter_dir: Path, interactions_path: Path):
    """Goes through every interaction and determines if it (1) might be a real interaction or (2) a good marker
    of general trends in the protein between ph5 and ph7

    Args:
        run_types (list[str]): all the types of runs Ex: ph5
        runs (list[str]): how many runs of each type. Name in file storage. Ex: r1
        parquet_dir (Path): Path that holds are distances of interactions
        ranked_inter_dir (Path): where the interactions of potential interest are listed
        interactions_path (Path): where the output files will be stored
    """

    # open and read all potential interactions (ranked_inter)
    interactions: list[str] = []
    with open(ranked_inter_dir, "r") as f:
        lines: list[str] = f.read().split("\n")
        for line in lines:
            if(len(line.split(",")) > 1):
                interactions.append(line)


    # go through each interaction and determine type. Add to correct group
    real_interact: list[str] = []
    marker: list[str] = []
    for interact in interactions[1:]:
        type: str = determine_interact_type(run_types, runs, parquet_dir, interact)
        if(type == "real_inter"):
            real_interact.append(interact)
        elif(type == "marker"):
            marker.append(interact)

    # print them out
    real_interact_path: Path = Path(interactions_path / "real_interactions.csv").resolve()
    marker_path: Path = Path(interactions_path / "marker_interactions.csv").resolve()
    if not real_interact_path.parent.is_dir():
        real_interact_path.parent.mkdir(parents=True, exist_ok=True)

    with open(real_interact_path, "w") as f:
        f.write("\n".join(real_interact))
    with open(marker_path, "w") as f:
        f.write("\n".join(marker))




def determine_interact_type(run_types: list[str], runs: list[str], parquet_dir: Path, interact: str) -> str:
    """Will take in the interaction and determine what type it is


    Args:
        run_types (list[str]): types of simulations run. Ex: ph5
        runs (list[str]): number of each simulation type run. Ex: r1
        parquet_dir (Path): directory where distance parquets are held
        interact (str): the string of the interaction

    Returns:
        str: interaction type
    """
    # real interaction requirements: at least one run has 30% beneath 8a
    # interact: the difference is greater then 1.5
    
    his_res: str = "_".join(interact.split(",")[0].split(" - ")[0].split("_")[0:2])
    charge_res: str = "_".join(interact.split(",")[0].split(" - ")[1].split("_")[0:2])
    delta_dist: float = float(interact.split(",")[1])

    for run_type in run_types:
        for run in runs:
            interact_dist: list[float] = readInParq(his_res, charge_res, run_type, run, parquet_dir)
            if_real_inter: bool = if_under_dist(interact_dist, 8, 0.3)
            if(if_real_inter):
                return "real_inter"
    if(abs(delta_dist) >= 1.5):
        return "marker"
    return "none"



def if_under_dist(interact_dist: list[float], under: float, percent: float) -> bool:
    """Will return if a certain list of numbers has a % (or more) of numbers under a number

    Args:
        interact_dist (list[float]): list of numbers
        under (float): value it needs to be under
        percent (float): % that need to be under value

    Returns:
        bool: if % of the list (or more) is under that number
    """

    length: float = float(len(interact_dist))
    count: float = 0.0
    for interact in interact_dist:
        if(interact <= under):
            count = count + 1
        if(math.isnan(interact)):
            length = length -1
    calc_percent: float = count / length
    return calc_percent >= percent





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
    ranked_inter_dir: Path = Path(DIR_SCRIPT / ".." / ".." / ".." / "002-finding-descriptors" / "data" / "ranked_interactions.csv").resolve()
    """Where the ranked interaction.csv file is placed"""
    interactions_path: Path = Path(DIR_SCRIPT  / ".." / ".." / "data").resolve()
    """Which interactions graphs are created for for"""

    if not interactions_path.is_dir():
        interactions_path.mkdir()


    main(run_types, runs, parquet_dir, ranked_inter_dir, interactions_path)