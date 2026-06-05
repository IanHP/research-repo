# 004-distance-over-time
Will take in interactions and graph that change in distance (angstroms) over time (ns)

# data
marker_interactions.csv: interactions that may show general trends in protein b/w ph5 and ph7
real_interaction.csv: actual interactions that may be different b/w
ph5 and ph7

# figures
These are the figures created by the script. Y-axis is distance between charged atoms of residues, x-axis is time
- marker_interactions: all graphs for marker interactions
- real_interactions: all graphs for real interactions


# scripts
## 001-categorize
These categorize found interactions in 002
- categerize_interactions: goes through all ranked_interactions and determines if they are real / marker / none. Then places in those respecitive .csv in data

## 002-create-graph
These create the graphs
- distance_over_time: will take in interactions.txt and the parquets to create the graphs for specific interactions
- residue_conversion: library to help determine what protein a residue is in
