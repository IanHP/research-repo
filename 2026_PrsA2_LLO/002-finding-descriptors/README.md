# 002 Finding Descriptors
General goal is to find a list of potentially descriptive HIS - charged residue charged atom distances across each frame
Initial idea was to use as features for a clustering algorithm, however also can be used for other things.

## Data
Output Trajectories/: the shortened / combined trajectories and there parameter/topology file. In NAMD format (.dcd, .psf)
- Created by 001
all_his_interactions.txt: Each line is a histidine, and the list of other charged residues it may interact with
- Created by 002
interaction_distances/: Holds each run's interaction distances for each frame. Stored in parquets
- Each file is all distances for each frame in 1 run for 1 histidine
- Created by 003
ranked_interactions.csv: holds each interaction's average in order from greatest difference from ph5 to ph7 to lowest
- Created by 004

## Scripts
### 001-export-trajectories
Goal is to take the fragmented .nc trajectories and turn into a single trajectory for each run with a certain time step
to make them smaller and easier to use 
- create_export_cpp.py: will create cpptraj files that will take in the fragmented trajectories and put them together with some
time step. Will need to put valid directory of fragmented trajectories and specified framestep (20 in this case)
- run_create_export_cpp.sh: will run the above script 
- cpp/: this will hold all the cpptraj input files that will actually create the single trajectories
- run_all_export.sh: will run all the cpptrajs

### 002-get-potential-interactions
Goal is find all HIS-charged residue (including other HIS) interactions that might exist. Done by going through all interactions
and a random set of frames, and checking if distances is <20. If yes, then store it.
- all_potential_interactions.py: script that takes in trajectories and finds potential interactions

### 003-distance-calc
Goal is to calculate all interactions found in previous section. Then any that do not go under 10a at least once across
either run will be thrown out. The rest will be written to parquet files, storing distance in each frame of each run
- each_frame_intraction_distance.py: script that does above

### 004-rank-interactions
Goal is to calculate average of each interaction, and difference between phs. If the average is below 15a for either it is kept.
Everything is output into a csv
- rank_interactions.py: will rank the interactions
- residue_conversion.py: small library that converts residue numbers between formats