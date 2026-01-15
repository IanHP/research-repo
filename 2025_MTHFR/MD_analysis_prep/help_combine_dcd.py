# THIS JUST PUTS DCD TO THE END OF THE PREVIOUS. NEED TO SPECIFY START AND END POINTS CAREFULLY
# First DCD should go from 0 --> start point of second DCD minus 1
# Second DCD should be entire thing

# doesnt actually work on windows, so just take number output. That is the first frame of the second dcd
# therefore DCD1 should go to one before that

from MDAnalysis import *
op = input("What is the output name: ")

if(op != "def"):
    psf = input("What is the psf name: ")
    inp1 = input("What is the first traj: ")
    tot_len = input("Total Length (2000): ")
    inp2 = input("What is the second traj: ")
else:
    op = "combo"
    psf = "step3_input.psf"
    inp1 = "step5_1.dcd"
    tot_len = "2000"
    inp2 = "step5_2.dcd"

# determine the first frame of dcd2
u1 = Universe(psf, inp1)
u2 = Universe(psf, inp2)

inp2_first = int(tot_len) - len(u2.trajectory)
print(inp2_first)

ag = u1.select_atoms("all")
with Writer("e_"+inp1, ag.n_atoms) as w:
    for ts in u1.trajectory[0:int(inp2_first)]:
        w.write(ag)
