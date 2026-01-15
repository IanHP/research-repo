import matplotlib.pyplot as plt
import numpy as np

def create_ddg_func_graph(inp, inp2, ifRoll, x_range = ()):
    """Creates graph of functionality and ddg

    Args:
        inp (file path): location of the ddg data
        inp2 (file path): location of the functionality data
    """
    
    # read in the ddg scores csv
    aa_res_list, roll_aa_ddg_list, max_res = ddg_scores(inp, ifRoll)
    
    # read up functionality scores csv
    aa_res_list_f, roll_aa_func_list = func_scores(inp2, ifRoll)


    # create the graph
    fig, ax = plt.subplots(figsize=[8.4, 4.8])
    marker_dict = {"none":"none","point":".","circle":"o","triangle down":"v","octagon":"8","square":"s","star":"*"}
    ax.set_ybound((0, 660))
    ax.set_title("Average Mutation Impact on Residue Stability and Functionality")
    ax.set_xlabel("Residue #")

    # create the ddg graph
    ax.set_ylabel("ddG")
    marker_ = marker_dict["point"]
    ax.scatter([-100],[0],c="red",label="Fuctionality",s=36.0)
    ax.scatter(aa_res_list, y=roll_aa_ddg_list, c="#1F77B4", marker=marker_, s=36.0, label="ddG")
    ax.plot(aa_res_list, roll_aa_ddg_list, c="#1F77B4")

    #create the functionality graph
    ax2 = ax.twinx()
    ax2.set_ylabel("Fuctionality Scores")
    ax2.scatter(aa_res_list_f, y=roll_aa_func_list, c="red", marker=marker_, s=36.0, label="Functionality")
    ax2.plot(aa_res_list_f, roll_aa_func_list, c="red")
    ax2.set_ybound(0.0,1)

    # annotations
    ax.axvspan(1, 35, color='#FFFFC5', alpha=0.5)
    ax.annotate('Disordered\nRegion',
            xy=(2, -1.3),
            horizontalalignment='left', verticalalignment='top',
            fontsize=8)

    ax.axvspan(48, 296, color='#808080', alpha=0.5)
    ax.axvspan(296, 338, color='#808080', alpha=0.5)
    ax.annotate('Catalytic\nDomain',xy=(100, -1.3),horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('Hinge\nRegion',xy=(298, -1.2),horizontalalignment='left', verticalalignment='top',fontsize=8)

    
    ax.axvspan(338, 381, color='#5db371', alpha=0.5)
    ax.axvspan(381, 394, color='#5db371', alpha=0.5)
    ax.axvspan(394, 412, color='#5db371', alpha=0.5)
    ax.annotate('Linker\nRegion',xy=(340, -1.3),horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('L1',xy=(352, -1.2),horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('L2',xy=(382, -1.2),horizontalalignment='left', verticalalignment='top',fontsize=8)
    ax.annotate('L3',xy=(397, -1.2),horizontalalignment='left', verticalalignment='top',fontsize=8)

    ax.axvspan(412, 644, color='#808080', alpha=0.5)
    ax.annotate('Regulatory\nDomain',xy=(500, -1.3),horizontalalignment='left', verticalalignment='top',fontsize=8)
    
    ax.plot([-10, max_res+10], [-0.1,-0.1], marker = 'o')
    ax.set_ybound(-1.45,0.25)

    
    #finilize and show graph
    ax.legend()
    if len(x_range) != 2:
        ax.set_xbound(-5, max_res+5)
    else:
         ax.set_xbound(x_range[0],x_range[1])
    plt.show()





def ddg_scores(inp, ifRoll):
    # Read in data
    with open(inp, "r") as file:
        data = file.read()
    

    # split into lines, that split each line based on tab
    data = data.split("\n")
    split_data = []
    for line in data:
        if line.startswith("6fcx_fill_ligands.pdb"):
            split_data.append(line.split("\t"))
    


    # create dictionary of each residue --> list of all ddgun scores for each mut
    aa_ddg_dict = {}
    max_res = 0
    
    for line in split_data:  
        name = get_aa_name(line[2]) # extracts aa index from the text
        
        if name != -10: # -10 means there is no index found
            
            # if the index is not in keys, add it 
            if not name in aa_ddg_dict.keys():
                aa_ddg_dict[name] = []
                if(name > max_res): # check if it is highest
                    max_res = name
            
            aa_ddg_dict[name].append(float(line[3]))

    

    # calculate average for each amino acid
    aa_res_list = []
    aa_ddg_list = []
    
    # goes through each residue, calculates average and adds it to the list
    keys = list(aa_ddg_dict.keys())
    keys.sort()
    for aa in keys:
        ddg_list = aa_ddg_dict[aa]
        aa_res_list.append(aa)
        aa_ddg_list.append(sum(ddg_list) / len(ddg_list))


    
    # create rolling
        # This is "residue" rolling, meaning it only includes residues-5 and residues+5. If there is no
        # data for one of those residues, it is skipped
    roll_aa_ddg_list = []
    if ifRoll:
        for res_start in aa_res_list:
            temp = []

            for res in range(res_start-5, res_start+6):
                try:
                    ind = aa_res_list.index(res)
                    temp.append(aa_ddg_list[ind])
                except ValueError:
                    pass

            roll_aa_ddg_list.append(sum(temp) / len(temp))
    else:
        roll_aa_ddg_list = aa_ddg_list
    
    return aa_res_list, roll_aa_ddg_list, max_res





def func_scores(inp2, ifRoll):
    # get all text
    with open(inp2, "r") as file:
        text = file.read()
    
    # turn into list of list for every line / data point
    text = text.split("\n")
    split_text = []
    for line in text:
        line = line.replace('"',"")
        if line.startswith('p.'):
            split_text.append(line.split(","))
    


    # get dictionary. Key: residue #, value: all func scores
    aa_func_dict = {}
    for line in split_text:
        if not "=" in line[0] and not "Ter" in line[0]:
            aa_res = get_aa_name(line[0], 5)
            if not aa_res in aa_func_dict.keys():
                aa_func_dict[aa_res] = []
            aa_func_dict[aa_res].append(float(line[1]))
    


    # calculate average for each amino acid
    aa_res_list = []
    aa_func_list = []
    
    # goes through each residue, calculates average and adds it to the list
    keys = list(aa_func_dict.keys())
    keys.sort()
    for aa in keys:
        func_list = aa_func_dict[aa]
        aa_res_list.append(aa)
        aa_func_list.append(sum(func_list) / len(func_list))


    
    # create rolling
        # This is "residue" rolling, meaning it only includes residues-5 and residues+5. If there is no
        # data for one of those residues, it is skipped
    roll_aa_func_list = []
    if ifRoll:
        for res_start in aa_res_list:
            temp = []

            for res in range(res_start-5, res_start+6):
                try:
                    ind = aa_res_list.index(res)
                    temp.append(aa_func_list[ind])
                except ValueError:
                    pass

            roll_aa_func_list.append(sum(temp) / len(temp))
    else:
        roll_aa_func_list = aa_func_list

    
    return aa_res_list, roll_aa_func_list
    


def get_aa_name(unedited, start=1):
    op = ""
    for letter in unedited[start:]:
        if letter.isnumeric():
            op = op + letter
        else:
            break
    if op == "":
        return -10
    else:
        return int(op)




if __name__ == "__main__":
    inp = r"F:\Research Lab\ddGun - Functionality Analysis\ddgun output\ddg_op.txt"
    inp2 = r"F:\Research Lab\ddGun - Functionality Analysis\ddgun output\f25WT_t1_simple_aa_floored.csv"
    
    rolling = input("If rolling (t/f): ")
    if rolling == "t":
        roll = True
    else:
        roll = False
    
    start = input("start of range")
    if(start.isnumeric()):
        end = input("end of range")
    
    if not start.isnumeric() or not end.isnumeric():
        rg = ()
    else:
        rg = (float(start), float(end))
    
    create_csv(inp, inp2, roll, rg)