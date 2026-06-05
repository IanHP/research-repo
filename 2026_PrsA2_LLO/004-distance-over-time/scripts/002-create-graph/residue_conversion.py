def prot_name_from_res_id(resnum: int | str):
    """Will take in a (ID from computer) residue number and return the protein it is apart of

    Args:
        resnum (int): residue number

    Returns:
        string: which protein it is apart of
    """

    # 1-273, 274-546, 547-1051
    if(isinstance(resnum, str)):
        resnum: int = int(resnum.split("_")[1])

    if resnum >= 1 and resnum <= 273:
        return "PrsA2_1"
    if resnum >= 274 and resnum <= 546:
        return "PrsA2_2"
    if resnum >= 547 and resnum <= 1051:
        return "LLO"

    raise Exception("Unknown res num")


def experiment_num_from_PC_num(resnum):
    """Will return experimental residue # using the
    pdb number

    Args:
        resnum (int): the residues pdb number

    Returns:
        int: the residue experimental number
    """
    prot_name = prot_name_from_res_id(resnum).split("_")

    if prot_name == "PrsA2_1":
        # experimental has length of 293
        # PC has length of 273
        return resnum + 20
    if prot_name == "PrsA2_2":
        return resnum + 20 - 273
    if prot_name == "LLO":
        # experiumental has length of 529
        # PC has length of 505
        return resnum + 24 - 546
