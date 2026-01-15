import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter
import os

def cleanTraj(u, currDirectory):
    try:
        selection_str = "not (resname TIP3 or resname HOH or resname SOD or resname CLA)"
        selected_atoms = u.select_atoms(selection_str)

        # Write the cleaned-up structure to new PSF and PDB files
        selected_atoms.write(os.path.join(currDirectory, 'step3_input.pdb'))

        # Write out the cleaned trajectory
        with mda.Writer(os.path.join(currDirectory, 'step5_1.dcd'), selected_atoms.n_atoms) as W:
            for ts in u.trajectory:
                W.write(selected_atoms)
    except Exception as e:
        print(f"Selection error: {e}")

def writeFrame(currMutation, currRun):
    currDirectory = os.path.join(currMutation, currRun)
    u = mda.Universe(os.path.join(currDirectory, 'init_coords/tog_step3_input.psf'),
                     os.path.join(currDirectory, 'init_coords/tog_step5.dcd'))

    output_dir = os.path.join(currMutation, f"{currRun}_frame")
    os.makedirs(output_dir, exist_ok=True)

    protein = u.select_atoms("protein and (not resname SOD and not resname CLA) or resname FAD or resname SAH")
    print(protein)

    fr_num = 1
    for frame in range(len(u.trajectory)):
        if(frame%5 == 0):
            u.trajectory[frame]  # Update the trajectory to the current frame
            filename = os.path.join(output_dir, f"frame_{fr_num:04d}.pdb")
            fr_num = fr_num + 1
            if os.path.isfile(filename):
                continue
            with PDBWriter(filename, multiframe=False) as pdb_writer:
                pdb_writer.write(protein)

    print("Finished extracting all frames.")


def main():
    
    inp = input("What is the run type name (ex: WT_v2)")
    Mutations = []
    while inp != "":
        Mutations.append(inp)
        inp = input("What is the run type name (ex: WT_v2)")
    print("The run types that will be opened are: ", Mutations)
    
    inp = int(input("How many runs of each type are there?"))
    Runs = []
    for i in range(1, inp+1):
        Runs.append(f"Run{i}")
    
    print("For each run type, the subfolders opened will be: ", Runs)
            
    for currMutation in Mutations:
        for currRun in Runs:
            writeFrame(currMutation, currRun)

if __name__ == "__main__":
    main()