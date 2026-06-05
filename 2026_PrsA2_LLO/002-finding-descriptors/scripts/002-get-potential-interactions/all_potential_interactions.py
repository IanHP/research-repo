import random
import warnings
import MDAnalysis as mda
import numpy as np
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning)
SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()


def main(run_dir: Path, run_types: list[str], runs: list[str], frame_num: int, op_file: Path) -> dict[str, list[str]]:
    """Looking at specified # of frames from each trajectory, will find all HIS-charged interactions
    that are less then <20a at least once

    Args:
        run_dir (Path): Path of the combined trajectories
        run_types (list[str]): Each run type (ex: ph5, ph7)
        runs (list[str]): Each run number's name. Shared in each run type (ex: r1)
        frame_num (int): number of frames from each trajectory to look at

    Returns:
        dict[str, list[str]]: each his is key, value is list of residues it may interact with
    """

    # define the runs
    random.seed(SEED)

    # select all HIS and charged res for each one. Form of name_resnum. Charged selected by residues <25a away from the HIS in run_type 1, run 1
    all_his: dict[str, list[str]] = {}  # dictionary. Keys: all his, values: all charged that are <25a away first frame
    all_charged: list[str] = []  # list of all charged residues in string form
    for run_type in run_types:
        for run in runs:
            dcd_path: Path = Path(run_dir / run_type / f"prsa2-llo-{run_type}-{run}_op.dcd")
            psf_path: Path = Path(run_dir / run_type / f"prsa2-llo-{run_type}-{run}_op.psf")
            all_his, all_charged = get_all_potential_interact(dcd_path, psf_path, frame_num, all_his, all_charged)
        print(f"{run_type} complete")
    
    # write it out
    write_each_his_interacts(op_file, all_his)
    return all_his


def get_all_potential_interact(dcd_path: Path, psf_path: Path, frame_num: int, all_his: dict[str, list[str]], all_charged: list[str]) -> tuple[dict[str, list[str]], list[str]]:
    """Will get all potential his charged interactions from this trajectory

    Args:
        dcd_path (String): trajectory DCD
        psf_path (String): trajectory PSF
        frame_num (int): how many frames to look at
        all_his (dictionary of list of strings): each his is a key, value is list of residues it may interact with
        all_charged (list of strings): all charged residues in the trajectory

    Returns:
        allHis, allCharged: dict / list updated with data from this trajectory
    """

    u: mda.Universe = mda.Universe(psf_path, dcd_path)

    # get all charged resids (including HIS) and put in string format
    if len(all_charged) < 1:
        charged_residues = u.select_atoms(
            "resname ASP or resname ASH or resname GLU or resname GLH or resname LYS or resname LYH or resname ARG or resname HIE or resname HID or resname HIP or resname HIS"
        ).residues
        for charge in charged_residues:
            all_charged.append(f"{clean_res_name(charge.resname)}_{charge.resid}")

    # get all HIS and put into list
    if len(all_his.keys()) < 1:
        his_residues = u.select_atoms(
            "resname HIE or resname HID or resname HIP or resname HIS"
        ).residues
        for his in his_residues:
            all_his[f"HIS_{his.resid}"] = []

    # go through random frames, check each interaction. If less then 25 store
    all_frames = random.sample(range(0, len(u.trajectory)), frame_num)
    for frame in all_frames:
        u.trajectory[frame]
        for his in all_his.keys():
            for charge in all_charged:
                # stops duplicate histidines or histidines "interacting" with itself
                if(charge.startswith("HIS")):
                    if(int(charge.split("_")[1]) <= int(his.split("_")[1])):
                        continue
                # find the distance and if valid, add to list
                dist = select_dist_his_charge_atoms(u, his, charge)
                if dist < 20 and charge not in all_his[his]:
                    all_his[his].append(charge)

    return all_his, all_charged


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


def write_each_his_interacts(op_file: Path, all_his: dict[str, list[str]]):
    """Takes in all interaction pairs and writes it to a .txt

    Args:
        allHis (dictionary of list of strings): each histidine has a list of charged residues it may interact with
    """
    # create directory to hold output files
    if not op_file.parent.is_dir():
        op_file.parent.mkdir(parents=True, exist_ok=True)

    with open(op_file, "w") as f:
        for his, charges in all_his.items():
            f.write(f"{his}: {', '.join(charges)}\n")


def clean_res_name(res_name: str) -> str:
    """Will take in a residue's name and return it's default name

    Args:
        res_name (str): residue name

    Returns:
        str: Cleaned residue name
    """
    if res_name.startswith("AS"):
        return "ASP"
    if res_name.startswith("GL"):
        return "GLU"
    if res_name.startswith("LY"):
        return "LYS"
    if res_name.startswith("HI"):
        return "HIS"
    return res_name


if __name__ == "__main__":
    # goal is to determine the HIS - charged interactions distances should be found for

    run_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "output_trajectories").resolve()
    op_file: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "all_his_interactions.txt").resolve()
    run_types: list[str] = ["ph5", "ph7"]
    runs: list[str] = ["r1", "r2", "r3"]

    if_test: bool = False  # this means a bunch of frame numbers will be tested instead of just one input below
    frame_num: int = 25  # how many frames to check distances in
    frame_nums: list[int] = [15, 20, 25]  # if testing, list of frames to chefck

    if not if_test:
        all_his: dict[str, list[str]] = main(run_dir, run_types, runs, frame_num, op_file)
    else:
        for frame_num in frame_nums:
            print(f"\n\n{frame_num}")
            all_his: dict[str, list[str]] = main(frame_num, run_types, runs, frame_num, op_file)
            for his, charges in all_his.items():
                print(f"{his}-{len(charges)}", end=" ")
