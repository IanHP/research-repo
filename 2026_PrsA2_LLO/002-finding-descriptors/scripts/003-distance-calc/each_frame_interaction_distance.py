import warnings

import fastparquet
import MDAnalysis as mda
import numpy as np
import pandas
from pathlib import Path
from numpy import typing as npt

warnings.filterwarnings("ignore", category=DeprecationWarning)
SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()


def main(run_dir: Path, output_dir: Path, run_types: list[str], runs: list[str]):
    all_dist_dict: dict[str, list[dict[str, dict[str, npt.NDArray[np.float64]]]]] = {}
    """holds all distances. D1: run_type, L2: each run, D3: each HIS, D4: each charged res, L5: each frame"""
    all_min_dict: dict[str, dict[str, dict[str, float]]] = {}
    """holds min distance across runs / each frame. D1: run_type, D2: each HIS, D3: each charged res"""

    # select all HIS and charged res for each one. Form of name_resnum. Charged selected by residues <25a away from the HIS in run_type 1, run 1
    potn_inter_file: Path = output_dir / "all_his_interactions.txt"
    all_his: dict[str, list[str]] = read_potential_interact_file(potn_inter_file)
    """dictionary. Keys: all his, values: all charged that are <25a away first frame"""

    # for every run type
    for run_type in run_types:
        print(run_type)
        all_dist_dict[run_type] = []
        all_min_dict[run_type] = {}
        # for every run, get all distances for all interactions
        for run_num, run in enumerate(runs):
            all_dist_dict[run_type].append({})
            dcd_path: Path = run_dir / run_type / f"prsa2-llo-{run_type}-{run}_op.dcd"
            psf_path: Path = run_dir / run_type / f"prsa2-llo-{run_type}-{run}_op.psf"
            all_dist_dict[run_type][run_num] = calculate_1run_all_distances(
                dcd_path, psf_path, all_his)
        # find minimum distance across every run of this run type
        for his_res, charges in all_his.items():
            all_min_dict[run_type][his_res] = {}
            for charged_res in charges:
                dist_list = []
                for run_num, run in enumerate(runs):
                    for dist in all_dist_dict[run_type][run_num][his_res][charged_res]:
                        dist_list.append(dist)

                minVal = np.nanmin(dist_list)
                all_min_dict[run_type][his_res][charged_res] = minVal

    # export all potential interactions across both types
    min_dir: Path = Path(output_dir / "min_dist_interactions.txt").resolve()
    with open(min_dir, "w") as f:
        for his_res in all_min_dict[run_types[0]].keys():
            f.write(f"{his_res}\n")
            for charged_res in all_min_dict[run_types[0]][his_res].keys():
                runType1Dist = all_min_dict[run_types[0]][his_res][charged_res]
                runType2Dist = all_min_dict[run_types[1]][his_res][charged_res]
                if runType1Dist < 10 or runType2Dist < 10:
                    f.write(
                        f"  {charged_res} - 5: {round(runType1Dist, 3)} - 7: {round(runType2Dist, 3)}\n"
                    )
            f.write("\n")

    # export all distances for later use. Save in parquet
    output_path: Path = Path(output_dir / "interaction_distances").resolve()
    output_dist_parquet(output_path, all_dist_dict, all_min_dict)



def output_dist_parquet(output_path: Path, all_dist_dict: dict[str, list[dict[str, dict[str, npt.NDArray[np.float64]]]]], 
                        all_min_dict: dict[str, dict[str, dict[str, float]]]): 
    """Will take every distance at output into parquets. Each run type has a folder, each run number / histidine has a file,
    in the file columns are interaction types, rows are frames.

    Args:
        output_path (Path): Where everything will be output
        all_dist_dict (dict[str, list[dict[str, dict[str, npt.NDArray[np.float64]]]]]): all the distances
        all_min_dict (dict[str, dict[str, dict[str, float]]]): min distances for each interaction / run_type
    """
    for run_type, runs in all_dist_dict.items():
        for run, all_his in enumerate(runs):
            for his, allCharged in all_his.items():
                # find which interactions have 1 instance where dist <10
                cor_charged_res: list[str] = []
                """Which interactions are valid"""
                for charged_res in allCharged.keys():
                    for run_type_ in all_dist_dict.keys():
                        run_type_dist: float = all_min_dict[run_type_][his][charged_res]
                        if run_type_dist < 10:
                            cor_charged_res.append(charged_res)
                            break

                # column (different interaction types) / row (frame #) names
                col: npt.NDArray = np.array(cor_charged_res)
                row: npt.NDArray[np.integer] = np.array(range(len(allCharged[cor_charged_res[0]])))

                # create np array of data
                data = np.ndarray([len(row), len(col)])
                for frame in row:
                    for chargeInd, charge in enumerate(col):
                        data[frame][chargeInd] = allCharged[charge][frame]
                # make dirs and output
                if not output_path.is_dir():
                    output_path.mkdir()
                run_output_path = Path(output_path / run_type).resolve()
                if not run_output_path.is_dir():
                    run_output_path.mkdir()
                #csv_op: Path = Path(run_output_path / f"r{run+1}_{his}.csv")
                prq_op: Path = Path(run_output_path / f"r{run + 1}_{his}.prq")
                df = pandas.DataFrame(data, row, col)
                #df.to_csv(csv_op)
                fastparquet.write(prq_op, df)



def calculate_1run_all_distances(dcd_path: Path, psf_path: Path, 
                                 all_his: dict[str, list[str]]) -> dict[str, dict[str, npt.NDArray[np.float64]]]:
    """Will calculate distances for every interaction across every frame for this specific run

    Args:
        dcd_path (Path): Location of run dcd
        psf_path (Path): Location of run psf
        all_his (dict[str, list[str]]): List of all interaction, organized by histidine

    Returns:
        dict[str, dict[str, npt.NDArray[np.float64]]]: Each histidine, 
            each charged it interacts with, list of distances
    """
    u: mda.Universe = mda.Universe(psf_path, dcd_path)

    n_frames: int = len(u.trajectory)
    run_dist_dict: dict[str, dict[str, npt.NDArray[np.float64]]] = {}
    # for every frame and interaction
    for i, ts in enumerate(u.trajectory):
        if i % 50 == 0:
            print(f"{str(dcd_path.name)} - {i}")
        for his, charges in all_his.items():
            if his not in run_dist_dict.keys():
                run_dist_dict[his] = {}
            for charge in charges:
                if charge not in run_dist_dict[his].keys():
                    run_dist_dict[his][charge] = np.full(
                        (n_frames,), np.nan, dtype=np.float64
                    )
                # find distance of interaction
                dist: float = select_dist_his_charge_atoms(u, his, charge)
                run_dist_dict[his][charge][i] = dist

    return run_dist_dict


def select_dist_his_charge_atoms(u: mda.Universe, his: str, charge: str) -> float:
    """Will determine distance between the histidine and charged residue

    Args:
        u (mda.Universe): the trajectory being analyzed. Put on correct frame.
        his (str): the histidine. In format "HIS_[resnum]"
        charge (str): the charged residue. In format "[resname]_[resnum]"

    Returns:
        float: the distance between charged atom of HIS and charged atom for this traj / frame
    """

    atom_1  = u.select_atoms(f"resid {his.split('_')[1]} and name NE2")

    charge_name = charge.split("_")[0]
    if charge_name[0:-1] == "AS":
        charged_atom = "OD2"
    if charge_name[0:-1] == "GL":
        charged_atom = "OE2"
    if charge_name[0:-1] == "LY":
        charged_atom = "NZ"
    if charge_name[0:-1] == "AR":
        charged_atom = "NH1"
    if charge_name[0:-1] == "HI":
        charged_atom = "NE2"

    atom_2 = u.select_atoms(f"resid {charge.split('_')[1]} and name {charged_atom}")
    dist:float = float(np.linalg.norm(atom_1.positions - atom_2.positions))

    return dist


def read_potential_interact_file(interact_file: Path) -> dict[str, list[str]]:
    """Reads in the potential interaction file. This defines which interaction distances will be 
    calculated and stored

    Args:
        interact_file (Path): Location of potential interaction file. Outut of 002-get-potential

    Returns:
        dict[str, list[str]]: each histidine (key) with list of each all residues it interacts with (value)
    """
    with open(interact_file, "r") as f:
        full_text: str = f.read()

    lines: list[str] = full_text.split("\n")[:-1]
    all_his: dict[str, list[str]] = {}
    for line in lines:
        his: str = line.split(": ")[0]
        charges: list[str] = line.split(": ")[1].split(", ")
        all_his[his] = charges

    return all_his


if __name__ == "__main__":
    run_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "output_trajectories").resolve()  
    """The directory with all the combined trajectories"""
    output_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "data").resolve()  
    """Directory where min distance and parquets will be put"""
    run_types = ["ph5", "ph7"]
    """The different run types. Ex: ph5, ph7"""
    runs = ["r1", "r2", "r3"]
    """Number of runs. Follows file naming. Ex: r1, r2"""

    print("starting")
    main(run_dir, output_dir, run_types, runs)
