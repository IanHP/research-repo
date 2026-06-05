from pathlib import Path

# CPPTRAJ TEMPLATE
"""# import parameter / topology of system
parm /ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/004-lm-prsa2-llo-ph5/simulations/01-prep/mol.prmtop

# import trajectory
trajin /ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/004-lm-prsa2-llo-ph5/simulations/04-prod/run-01/outputs/prod_npt_?.nc 1 last 20
trajin /ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/004-lm-prsa2-llo-ph5/simulations/04-prod/run-01/outputs/prod_npt_??.nc 1 last 20
trajin /ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/004-lm-prsa2-llo-ph5/simulations/04-prod/run-01/outputs/prod_npt_???.nc 1 last 20

# remove water / ions + output
strip :WAT,Na+,Cl- parmout prsa2-llo-ph5-r1_op.psf
trajout prsa2-llo-ph5-r1_op.dcd "dcd"

run"""

SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()


def main(dirs: list[Path], cpp_output_dir: Path, output_dir: Path, time_step: int):
    """Will take in cpptraj scripts that will take in each runs trajectories
    put them together with a specific time step

    Args:
        dirs (list[Path]): where all the trajectories are currently stored
        cpp_output_dir (Path): where the cpptraj files will be stored
        output_dir (Path): where the sliced / combined trajs will be stored
        time_step (int): how big a stride through the trajectory
    """    

    # put in the folder the holds output trajectories for each run
    # get list of trajectories in each folder. Index aligns to dir index above.
    all_files: list[list[Path]] = []
    for dir in dirs:
        all_files_temp: list[Path] = [f.absolute() for f in dir.iterdir() if f.suffix == ".nc"]
        all_files_temp.sort(key=lambda f: int(f.name.split("_")[2].split(".")[0]))
        all_files.append(all_files_temp)
    
    # write cpp output folder
    if not cpp_output_dir.is_dir():
        cpp_output_dir.mkdir(parents=True, exist_ok=True)

    # for each run, create a cpptraj that will create a valid trajectory
    for ind, dir in enumerate(dirs):
        base = """# import parameter / topology of system"""
        # get the parameter / topology of system
        prmtop_path = Path(dir / ".." / ".." / ".." / "01-prep" / "mol.prmtop").resolve()
        base = base + f"\nparm {str(prmtop_path)}"
        # import all the trajectory files
        base = base + "\n# import trajectory"
        for file in all_files[ind]:
            base = base + f"\ntrajin {file} 1 last {time_step}"
        # output file name
        ph_dir: Path = Path(dir / ".." / ".." / "..").resolve()
        ph: str = ph_dir.parent.name[-1]
        num: str = dir.parent.name[-1]
        file_name: str = f"prsa2-llo-ph{ph}-r{num}_op"
        psf_op_name: Path = Path(output_dir / f"ph{ph}" / f"{file_name}.psf").resolve()
        dcd_op_name: Path = Path(output_dir / f"ph{ph}" / f"{file_name}.dcd").resolve()
        base = (
            base
            + f"\n\n# remove water / ions + output"
            + f"\nstrip :WAT,Na+,Cl- parmout {str(psf_op_name)}"
            + f'\ntrajout {str(dcd_op_name)} "dcd"\n\nrun'
        )
        #print out cpptraj
        cpp_output_file: Path = Path(cpp_output_dir / f"export_{file_name}.cpp").resolve()
        with open(cpp_output_file, "w") as f:
            f.write(base)
    
        # create directory to hold output files
        if not dcd_op_name.parent.is_dir():
            dcd_op_name.parent.mkdir(parents=True, exist_ok=True)





if __name__ == "__main__":
    base_ph7_path: Path = Path("/ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/003-lm-prsa2-llo-ph7/simulations/04-prod").resolve()
    base_ph5_path: Path = Path("/ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/004-lm-prsa2-llo-ph5/simulations/04-prod").resolve()
    dirs: list[Path] = [
        Path(base_ph7_path / "run-01" / "outputs").resolve(),
        Path(base_ph7_path / "run-02" / "outputs").resolve(),
        Path(base_ph7_path / "run-03" / "outputs").resolve(),
        Path(base_ph5_path / "run-01" / "outputs").resolve(),
        Path(base_ph5_path / "run-02" / "outputs").resolve(),
        Path(base_ph5_path / "run-03" / "outputs").resolve(),
    ]
    cpp_output_dir: Path = Path(DIR_SCRIPT / "cpp").resolve()
    output_dir: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "output_trajectories").resolve()
    main(dirs, cpp_output_dir, output_dir, 20)