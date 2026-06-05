import fastparquet
import numpy as np
import residue_conversion
from pathlib import Path

SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()

def main(min_dist_path: Path, parquet_dir: Path, output_dir: Path):
    """Calculate the average distance of each interaction in the parquets,
    finds difference between ph5 and ph7 and saves it

    Args:
        min_dist_path (Path): holds smallest distance for each interaction
            for each run type. Used to find which interactions are being used
        parquet_dir (Path): directory that holds all the distance parquets
        output_dir (Path): where the csv with average distance of each interaction
            will be stored
    """
    # determine interactions present in parquets
    interactions = {}
    """holds what each histidine interacts with. D1: histidines L2: all charged it interacts with"""
    with open(min_dist_path) as f:
        text = f.read()
        his_blocks = text.split("\n\n")
        for his_block in his_blocks:
            if len(his_block) > 5:
                split_hit_block = his_block.split("\n")
                interactions[split_hit_block[0]] = []
                for charge_line in split_hit_block:
                    split_charge_line = charge_line.split(" - ")
                    if len(split_charge_line) > 1:
                        interactions[split_hit_block[0]].append(split_charge_line[0][2:])
    # go through every interaction in .txt and calculate average distance (mean)
    # from each run type across all runs / frames
    avg_dist_list_close: list[list[str | float]] = []
    """holds all of the average distances. [interaction, difference (7-5), ph7 distance, ph5 distance]"""
    for hisRes, chargeRess in interactions.items():
        for chargeRes in chargeRess:
            ph5_avg_dist: float = boltzAverage(readInParq(hisRes, chargeRes, "ph5", parquet_dir))
            ph7_avg_dist: float = boltzAverage(readInParq(hisRes, chargeRes, "ph7", parquet_dir))
            # if average is below 15, include it
            if ph5_avg_dist < 15 or ph7_avg_dist < 15:
                avg_dist_list_close.append(
                    [
                        f"{prot_name_from_res(hisRes)} - {prot_name_from_res(chargeRes)}",
                        ph7_avg_dist - ph5_avg_dist,
                        ph7_avg_dist,
                        ph5_avg_dist,
                    ]
                )
    # sort based on difference and output data (perhaps sort based off if on different proteins?)
    avg_dist_list_close = sorted(
        avg_dist_list_close, key=lambda interact: abs(interact[1]), reverse=True
    )
    # output
    with open(output_dir, "w") as f:
        f.write("Interaction, Difference (7-5), ph7 Average, ph5 Average\n")
        for avg_dist in avg_dist_list_close:
            nice_avg_dist = [str(round(float(item), 4)) for item in avg_dist[1:]]
            f.write(f"{avg_dist[0]}, " + ", ".join(nice_avg_dist))
            f.write("\n")


def prot_name_from_res(res_string) -> str:
    """Returns the protein it is apart of

    Args:
        res_string (str): residue that protein is being found
            for

    Returns:
        str: residue with new name that includes protein it is apart of
    """
    # 1-273, 274-546, 547-1051
    resnum = int(res_string.split("_")[1])

    prot_name = residue_conversion.prot_name_from_res_id(resnum)
    return f"{res_string}_{prot_name}"


def boltzAverage(distances) -> float:
    """Takes in the distances for a specific interaction / run type, and calculates
    average. Should be boltz, but is just normal mean.

    Args:
        distances (list[int]): list of all distances for an interaction across
            all runs and frames

    Returns:
        float: average
    """
    return np.nanmean(distances)


def readInParq(his_res: str, charge_res: str, run_type: str, parquet_dir: Path) -> list[float]:
    """Will read in a parquet file, and based on run_type / his / charge, will return all runs distance data

    Args:
        hisRes (String): histidine res of interest
        chargeRes (String): Charged intraction of interest
        runType (String): run type (ex: ph5) of interest

    Returns:
         list of int: list of all distances for interaction across all runs
    """
    ret: list[float] = []
    for runNum in ["r1", "r2", "r3"]:
        file_path: Path = Path(parquet_dir / run_type / f"{runNum}_{his_res}.prq").resolve()
        df = fastparquet.ParquetFile(file_path).to_pandas([charge_res])
        ret.extend(df[charge_res].tolist())
    return ret


if __name__ == "__main__":
    min_dist_path: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "min_dist_interactions.txt").resolve()
    parquet_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "interaction_distances").resolve()
    output_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "ranked_interactions.csv").resolve()
    main(min_dist_path, parquet_dir, output_dir)
