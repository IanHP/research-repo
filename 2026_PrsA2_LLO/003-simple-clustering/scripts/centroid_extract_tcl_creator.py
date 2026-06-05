from scipy.cluster.hierarchy import centroid
from pathlib import Path

SEED: int = 600411609

DIR_SCRIPT: Path = Path(__file__).parent.resolve()
DIR_STUDY: Path = Path(DIR_SCRIPT / ".." / ".." / "..").resolve()

def read_in_clusters(centroids_file: Path) -> tuple[list[int], list[list[str]]]:
    # get the different sections
    with open(centroids_file, "r") as f:
        all_text: str = f.read()
    all_text: list[str] = all_text.split("\n\n\n")
    centroids_str: str = all_text[1]
    cluster_types: str = all_text[2]

    # determine if a cluster is valid (only includes one run type)
    included_clusters: list[int] = []
    cluster_types_lst: list[str] = cluster_types.split("\n")[1:]
    for cluster_num, cluster_type in enumerate(cluster_types_lst):
        types: list[str] = cluster_type.split(",")[1:]
        type_set: set[str] = set()
        for type in types:
            type_set.add(type.split("-")[0])
        if(len(type_set) < 2):
            included_clusters.append(cluster_num)
    
    # for all included clusters, create representative list of centroids
    # indicies line up in included_clusters and centroids
    centroid_lst: list[str] = centroids_str.split("\n")[1:]
    centroids: list[list[str]] = []
    for cluster in included_clusters:
        centroid: str = centroid_lst[cluster].split(",")[1]
        centroid: list[str] = centroid.split("-")
        temp= centroid[1].split(" ")
        centroid.extend(temp)
        centroid.pop(1)
        centroids.append(centroid)
    
    return included_clusters, centroids


def split_clusters(centroids: list[list[str]]) -> dict[str, dict[str, list[int]]]:
    centroid_type_dict: dict[str, dict[str, list[int]]] = {}
    """Holds the centroids in every run_type, run_num as a double dictionary in a list"""
    for centroid in centroids:
        if(not centroid[0] in centroid_type_dict.keys()):
            centroid_type_dict[centroid[0]] = {}
        if(not centroid[1] in centroid_type_dict[centroid[0]].keys()):
            centroid_type_dict[centroid[0]][centroid[1]] = []

        centroid_type_dict[centroid[0]][centroid[1]].append(int(centroid[2]))

    return centroid_type_dict


def write_tcl(centroid_type_dict: dict[str, dict[str, list[int]]], traj_dir: Path,
              dcd_op: Path, op_tcl_dir: Path):
    # for each run_type create a different file
    for run_type, centroids in centroid_type_dict.items():
        to_print: list[str] = []
        all_pdb_names: list[str] = []
        run_type_dir: Path = traj_dir / run_type
        # for each run_num, import all data
        for run_num, frames in centroids.items():
            # import files
            to_print.append(f"### {run_type} - {run_num}")
            psf_dir: Path = run_type_dir / f"prsa2-llo-{run_type}-{run_num}_op.psf"
            to_print.append(f'mol new {{{str(psf_dir)}}}')
            dcd_dir: Path = run_type_dir / f"prsa2-llo-{run_type}-{run_num}_op.dcd"
            to_print.append(f'mol addfile {{{str(dcd_dir)}}} waitfor all\n')
            # create PDBs of the section
            for frame in frames:
                to_print.append(f'set c2 [atomselect top "all" frame {frame}]')
                to_print.append(f'$c2 writepdb {{{str(dcd_op / f"{run_type}_{run_num}_{frame}.pdb")}}}')
            to_print.append(f'mol delete all\n\n')
            all_pdb_names.extend([f"{run_type}_{run_num}_{item}" for item in frames])
        # create total DCD
        to_print.append(f'mol new {{{str(dcd_op / f"{all_pdb_names[0]}.pdb")}}}')
        for pdb_names in all_pdb_names[1:]:
            to_print.append(f'mol addfile {{{str(dcd_op / f"{pdb_names}.pdb")}}}')
        to_print.append(f'set c2 [atomselect top "all"]')
        to_print.append(f'animate write dcd {{{str(dcd_op / f"{run_type}_selected_frame.dcd")}}} \
sel $c2 skip 1 waitfor all')

        #write out to file
        op_tcl: Path = op_tcl_dir / f"create_dcd_{run_type}.tcl"
        with open(str(op_tcl),"w") as f:
            for line in to_print:
                f.write(f"{line}\n")










def main():
    """Main function. Will take in a cluster .csv, location of 
    main dcd / psf files of simulation, and create a .tcl script
    that will extract the centroid frames and create new dcd for each
    run type
    """

    # Files to read in / base vars
    traj_dir: Path = (DIR_STUDY / 'analysis' / '002-finding-descriptors' /
                      'data' / 'output_trajectories')
    """Trajectories directories"""
    centroids_file: Path = DIR_SCRIPT / ".." / "data" / "hdbscan_cluster_leaf_m25_eps.csv"
    """CSV with all centroids"""
    dcd_pdb_dir: Path = (DIR_STUDY / 'analysis' / '003-simple-clustering' / 
                         "data" / "centroid_data" / "hdb_leaf_m25_eps")
    """Where PDB/DCDs of centroids will be place"""
    op_tcl_dir: Path = DIR_SCRIPT
    """Directory where tcls will be placed"""
    run_types: list[str] = ["ph5","ph7"]
    """Run types. The different types of MD simulations created. Ex: ['ph5']"""
    run_nums: list[str] = ["r1","r2","r3"]
    """The number of runs done for each type. Ex: ['r1']"""

    # read in the centroids file. Get all valid clusters, and their centroids
    centroid_cluster, centroids = read_in_clusters(centroids_file)

    # create list of centroids for each pH
    centroid_type_dict: dict[str, dict[str, list[int]]] = split_clusters(centroids)
    
    # create file that creates pdbs of those specific frames. 1 file for ph5, 1 for ph7
    write_tcl(centroid_type_dict, traj_dir, dcd_pdb_dir, op_tcl_dir)







if __name__ == "__main__":
    main()









