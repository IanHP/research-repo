from matplotlib.pylab import real
from pathlib import Path


## TEMPLATE
"""
set ambient, 0.25
set cartoon_transparency, 0.35
set field_of_view 20
set cartoon_discrete_colours, on
set antialias, 4


delete all
load C:\\Users\sodiu\OneDrive\Desktop\Lm-PrsA2-LLO-github\study\analysis\002-finding-descriptors\data\output_trajectories\ph5\prsa2-llo-ph5-r1_op.psf
load C:\\Users\sodiu\OneDrive\Desktop\Lm-PrsA2-LLO-github\study\analysis\002-finding-descriptors\data\output_trajectories\ph5\prsa2-llo-ph5-r1_op.dcd
dss

select prsa2_1, residue 1-273
select prsa2_2, residue 274-546
select llo, residue 547-1051

hide sticks
show cartoon
color yellow, prsa2_1
color orange, prsa2_2
color cyan, llo

select ha, residue 152+985
color purple, ha
show sticks, ha

center residue 1000, -1

set_view (\
    -0.536496997,   -0.835647225,    0.117732026,\
    -0.076919228,   -0.090499170,   -0.992914736,\
     0.840375423,   -0.541750610,   -0.015727460,\
    -0.000482872,    0.001356579, -447.645690918,\
    92.667434692,   42.113895416,   53.924823761,\
   307.347137451,  588.326293945,  -20.000000000 )

ray
png test.png

"""

SEED: int = 600411609 # test is 520586740
DIR_SCRIPT: Path = Path(__file__).parent.resolve()


def main():
    """Will create a pymol script that takes in each interaction in marker / 
    real interactions and creates a photo highlighting the residues
    """

    script_path: Path = DIR_SCRIPT / "pymol_script.pml"
    """Where the pymol script will be created"""
    photos_path: Path = Path(DIR_SCRIPT / ".." / ".." / "figures" / "photos").resolve()
    """Where the protein photos will be created"""
    marker_path: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "marker_interactions.csv").resolve()
    real_path: Path = Path(DIR_SCRIPT / ".." / ".." / "data" / "real_interactions.csv").resolve()
    """Where the marker / real interaction .csv are stored"""
    dcd_path: Path = Path(DIR_SCRIPT / ".." / ".." / ".." / "002-finding-descriptors" / "data" / "output_trajectories" / "ph5" / "prsa2-llo-ph5-r1_op.dcd").resolve()
    psf_path: Path = Path(DIR_SCRIPT / ".." / ".." / ".." / "002-finding-descriptors" / "data" / "output_trajectories" / "ph5" / "prsa2-llo-ph5-r1_op.psf").resolve()
    """Where the dcd / psf for ph5 run 1 are stored"""

    view_dict = {
        "prsa2_1_top": "set_view (\
    -0.536496997,   -0.835647225,    0.117732026,\
    -0.076919228,   -0.090499170,   -0.992914736,\
     0.840375423,   -0.541750610,   -0.015727460,\
    -0.000482872,    0.001356579, -447.645690918,\
    92.667434692,   42.113895416,   53.924823761,\
   307.347137451,  588.326293945,  -20.000000000 )",
        "prsa2_2_top":"set_view (\
     0.301102459,   -0.765755236,   -0.568288565,\
     0.251737416,   -0.510965407,    0.821905971,\
    -0.919748127,   -0.390538216,    0.038917623,\
    -0.000482872,    0.001356579, -447.645690918,\
    92.667434692,   42.113895416,   53.924823761,\
   307.347137451,  588.326293945,  -20.000000000 )",
        "prsa2_bot":"set_view (\
    -0.457111269,   -0.872142613,    0.174388438,\
     0.024727724,   -0.208451077,   -0.977712452,\
     0.889050603,   -0.442609996,    0.116848543,\
    -0.001306191,    0.001153737, -434.362762451,\
    87.091224670,   41.923095703,   56.279808044,\
   294.100921631,  575.080017090,  -20.000000000 )",
        "inner_prsa2_2":"set_view (\
    -0.536496997,   -0.835647225,    0.117732026,\
    -0.076919228,   -0.090499170,   -0.992914736,\
     0.840375423,   -0.541750610,   -0.015727460,\
    -0.000482872,    0.001356579, -447.645690918,\
    92.667434692,   42.113895416,   53.924823761,\
   307.347137451,  588.326293945,  -20.000000000 )",
        "inner_prsa2_1":"set_view (\
     0.301102459,   -0.765755236,   -0.568288565,\
     0.251737416,   -0.510965407,    0.821905971,\
    -0.919748127,   -0.390538216,    0.038917623,\
    -0.000482872,    0.001356579, -447.645690918,\
    92.667434692,   42.113895416,   53.924823761,\
   307.347137451,  588.326293945,  -20.000000000 )"
    }

    create_script(script_path, photos_path, marker_path, 
                  real_path, dcd_path, psf_path, view_dict)


def create_script(script_path: Path, photos_path: Path, marker_path: Path,
                  real_path: Path, dcd_path: Path, psf_path: Path, view_dict: dict[str, str]):
    """Creates pymol script that highlights all interactions

    Args:
        script_path (Path): where pymol script will be created
        photos_path (Path): Where protein photos will be created
        marker_path (Path): where the marker interaction .csv is stored
        real_path (Path): where real interaction .csv is stored
        dcd_path (Path): where the dcd of the simulation is stored
        psf_path (Path): where the psf of the simulation is stored
        view_dict (Dict): holds all view angles of protein

    """
    # header of script
    base: str = "set ambient, 0.25\nset cartoon_transparency, 0.35 \
                \nset field_of_view 20\nset cartoon_discrete_colours, on \
                \nset antialias, 4\n\ndelete all\n"

    base = base + f"load {str(psf_path)}\n"
    base = base + f"load {str(dcd_path)}\ndss\n\n"
    
    base = base + "select prsa2_1, residue 1-273\nselect prsa2_2, residue 274-546\n"
    base = base + "select llo, residue 547-1051\n\n\n"
    
    # each different interaction
    for (type, type_path) in [("marker", marker_path), ("real",real_path)]:

        with open(type_path, "r") as f:
            interactions: list[str] = [item.split(",")[0] for item in f.read().split("\n")]
            type_photos_path: Path = Path(photos_path / type).resolve()
        
        if not type_photos_path.is_dir():
            type_photos_path.mkdir(parents=True, exist_ok=True)

        for interaction in interactions:
            base = each_interaction_section(base, interaction, type_photos_path, view_dict)

    with open(script_path, "w") as f:
        f.write(base)


def each_interaction_section(base: str, interaction: str, photos_path: Path, view_dict: dict[str, str]) -> str:
    """Will take in an interaction and add its info to base

    Args:
        base (str): the current script
        interaction (str): the interaction being added
        photos_path (Path): where the photo will be placed
        view_dict (Dict): holds all view angles of protein

    Returns:
        str: updated script
    """
    # get interaction info
    try:
        residues: list[str] = [item.split("_")[1] for item in interaction.split(" - ")]
        # create path
        op_path = Path(photos_path / f"{"-".join(residues)}.png")
        # determine view angle
        if(residues[0] == "200" or residues[0] == "985" or residues[0] == "833"):
            viewpoint = view_dict["prsa2_1_top"]
        elif(residues[0] == "972" or residues[0] == "473"):
            viewpoint = view_dict["prsa2_2_top"]
        elif(residues[0] == "945"):
            viewpoint = view_dict["prsa2_bot"]
        elif(residues[0] == "395"):
            viewpoint = view_dict["inner_prsa2_2"]
        else:
            viewpoint = view_dict["prsa2_1_top"]
    except:
        return base

    
    base = base + "hide sticks\nshow cartoon\ncolor yellow, prsa2_1\n"
    base = base + "color orange, prsa2_2\ncolor cyan, llo\n\n"

    base = base + f"select ha, residue {residues[0]}+{residues[1]}\n"
    base = base + "color purple, ha\nshow sticks, ha\ncenter residue 1000, -1\n\n"
    
    base = base + f"{viewpoint}\n\n"
    base = base + f"ray\npng {op_path}\n\n\n"
    return base







if __name__ == "__main__":
    main()